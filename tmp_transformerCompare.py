# import pandas as pd
# import numpy as np


# path_var_response = '../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv'
# path_bed_var = '../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6'



# bed_var = pd.read_csv(path_bed_var, sep = '\t', header=None)
# bed_var = bed_var.iloc[np.where((bed_var[0]).isin(['chr1', 'chr2', 'chr3']))]


# var_Y = pd.read_csv(path_var_response, sep = '\t', index_col = 'binID')
# var_Y = var_Y.loc[bed_var[3]]

# var_Y.to_csv('../external/tmp_transformer/chr1_3_varY.tsv', sep = '\t')



# Y_1k = pd.read_csv('../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/resTab_bed6Generated_1k_cnn_withInfo.tsv', 
#                    sep = '\t', index_col = 'binID')

# def read_bed(path_bed):
#     bed = pd.read_csv(path_bed, sep = '\t', header = None)
#     # excluded = bed.iloc[np.where((bed[0] == 'chrX') | (bed[0] == 'chrY') | (bed[0] == 'chrM'))]
    
#     # bed = bed.iloc[~bed.index.isin(excluded.index)]
#     bed['binID'] = bed[3]
#     bed = bed.set_index('binID')
    
#     return bed


# bed_1k_chr1_3 = read_bed('../external/tmp/bed_1kCNN_chr1_3.bed')
# Y_1k_chr1_3 = Y_1k.loc[bed_1k_chr1_3.index]

# Y_1k_chr1_3.to_csv('../external/tmp/Y_1k_chr1_3.tsv', sep = '\t')

# with h5py.File('../external/tmp_transformer/ftrMtrix_chr1_3.h5', 'r') as f:
#      all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
#      selected_cols = [col for col in all_features if col not in new_ftrs]
#      X = f['/X/block0_values'] 
#      print(X.shape)
#      X = X[:, np.where(np.isin(all_features,selected_cols))[0]]
#      print(X.shape)
#      Sc_data = standard_scale_data(X)
#      dump(Sc_data[1], open('../external/tmp_transformer/1kChr1_3_cnn_stdScaler_1373Ftrs.pkl', 'wb'))






import numpy as np
from tensorflow.keras.layers import Input,  Dropout, Dense
from tensorflow.keras.models import Model
import tensorflow as tf
import h5py
from tensorflow.keras.layers import LayerNormalization, MultiHeadAttention
from readFtrs_Rspns import set_gpu_memory_limit
from scipy.stats import spearmanr
import pandas as pd
from pickle import load
from sklearn.metrics import mean_squared_error
from AdjuscentBins.generators import data_generator, test_data_generator, prepare_test_dataY, get_validation_bins, create_info_test, create_info_train



gpu_fraction = 0.05
set_gpu_memory_limit(gpu_fraction)



# Positional Encoding
def get_positional_encoding(seq_length, d_model):
    pos_encoding = np.array([
        [pos / np.power(10000, 2 * (i // 2) / d_model) for i in range(d_model)]
        for pos in range(seq_length)
    ])
    pos_encoding[1:, 0::2] = np.sin(pos_encoding[1:, 0::2])
    pos_encoding[1:, 1::2] = np.cos(pos_encoding[1:, 1::2])
    pos_encoding = tf.cast(pos_encoding, dtype=tf.float32)
    pos_encoding = tf.expand_dims(pos_encoding, 0)  # Expand to 3D for batch compatibility
    return pos_encoding


# Transformer Encoder Layer
class TransformerEncoderLayer(tf.keras.layers.Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(TransformerEncoderLayer, self).__init__()
        self.mha = MultiHeadAttention(num_heads=num_heads, key_dim=d_model)
        self.ffn = tf.keras.Sequential([
            Dense(dff, activation='relu'),
            Dense(d_model)
        ])
        self.layernorm1 = LayerNormalization(epsilon=1e-6)
        self.layernorm2 = LayerNormalization(epsilon=1e-6)
        self.dropout1 = Dropout(rate)
        self.dropout2 = Dropout(rate)
        
    def call(self, x, training, mask=None):
        attn_output = self.mha(x, x, attention_mask=mask)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(x + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        out2 = self.layernorm2(out1 + ffn_output)
        return out2

def build_transformer_model(input_shape, num_layers=4, d_model=128, num_heads=4, dff=512, rate=0.1):
    seq_length = input_shape[0]
    num_features = input_shape[1]
    pos_encoding = get_positional_encoding(seq_length, d_model)
    
    inputs = Input(shape=input_shape)
    x = Dense(d_model)(inputs)
    x *= tf.sqrt(tf.cast(d_model, tf.float32))  # Scale embeddings
    x += pos_encoding  # Add positional encoding
    
    for _ in range(num_layers):
        x = TransformerEncoderLayer(d_model, num_heads, dff, rate)(x)
        
    outputs = Dense(1, activation='softplus')(x)
    model = Model(inputs=inputs, outputs=outputs)
    return model

def custom_poisson_loss(y_true, y_pred):
      mask = tf.cast(tf.not_equal(y_true, -1.0), tf.float32)  # Mask for available rates
      
      # Ensure y_pred is greater than 0 to avoid log(0)
      y_pred_safe = tf.maximum(y_pred, tf.keras.backend.epsilon())
      
      # Compute element-wise Poisson loss
      loss = y_pred_safe - y_true * tf.math.log(y_pred_safe)
      
      loss *= mask  # Element-wise multiplication, broadcasting loss to match mask shape
      
      # Return the loss without summing and averaging across the sequence
      return loss  # Shape [batch_size, sequence_length]

#############################################################################
#############################################################################
#############################################################################

nn_batch_size = 23
num_regions_per_sample = 100

# Model Configuration
input_shape = (num_regions_per_sample, 1372)  # 100 regions, each with 1500 features
transformer_model = build_transformer_model(input_shape)
transformer_model.compile(optimizer='adam', loss=custom_poisson_loss)





path_response_tr = '../external/tmp_transformer/Y_1k_chr1_3.tsv'

path_features = '../external/tmp_transformer/ftrMtrix_chr1_3.h5'

path_scaler = '../external/tmp_transformer/1kChr1_3_cnn_stdScaler_1372Ftrs.pkl'

path_validation_set_gbm = '../external/tmp_transformer/compareModels/GBM/rep_train_test/GBM_predTest1.tsv'

path_bed_validation = '../external/tmp_transformer/bed_var_chr1_3.bed'
path_bed_train = '../external/tmp_transformer/bed_1kCNN_chr1_3.bed'









validation_bins = get_validation_bins(path_validation_set_gbm, 
                                      path_bed_validation, path_bed_train) #[0]

train_info = create_info_train(path_response_tr, validation_bins)
total_n_samples = train_info.shape[0] # 2881044
print(f'Model training on {total_n_samples} bins')



           
from tensorflow.keras.callbacks import ModelCheckpoint

checkpoint_callback = ModelCheckpoint(filepath='../external/tmp_transformer/model_checkpoint_{epoch:02d}.h5',
                                      save_freq='epoch', period = 30)

transformer_model.fit(
         data_generator(path_features, train_info, path_scaler, nn_batch_size, num_regions_per_sample),
         steps_per_epoch=(total_n_samples // (num_regions_per_sample * nn_batch_size)),  # Total smaller batches per epoch
         epochs=900, callbacks=[checkpoint_callback])




print('*******************************')

#################################################################################

# # load the last model
from tensorflow.keras.models import load_model
with h5py.File('../external/tmp_transformer/model_checkpoint_420.h5', 'r') as f:
            # load the model
    transformer_model = load_model(f, custom_objects = {'TransformerEncoderLayer': TransformerEncoderLayer,
                                                        'custom_poisson_loss': custom_poisson_loss})


path_test_response = '../external/tmp_transformer/Y_1k_chr1_3.tsv'
path_bed_test = '../external/tmp_transformer/bed_1kCNN_chr1_3.bed'
path_test_features = '../external/tmp_transformer/ftrMtrix_chr1_3.h5'
test_on = 'validation_set'

info_test = create_info_test(path_test_response, path_bed_test, validation_bins)


middle_region_index = 50
y_pred = transformer_model.predict(test_data_generator(info_test, path_test_features,
                                                       path_test_response, 
                                                       path_bed_test, 
                                                       path_scaler,
                                                       nn_batch_size,
                                                       num_regions_per_sample,
                                                       middle_region_index, test_on))



Y_preds = y_pred[:, middle_region_index]
Y_obs, obs_df =prepare_test_dataY(info_test, nn_batch_size, 
                       num_regions_per_sample, middle_region_index, test_on)


mse = mean_squared_error(Y_obs, Y_preds)
print(f'Mean Squared Error for Middle Region: {mse}')
mae = np.mean(np.abs(Y_obs - Y_preds))
print(f'Mean Absolute Error for Middle Region: {mae}')

corr, p_value = spearmanr(Y_preds, Y_obs)
print(f'Spearman correlation for Middle Region: {corr}. p-value: {p_value}')


# spearmanr(Y_preds[np.where(obs_df['nMut'] != 0)], Y_obs[np.where(obs_df['nMut'] != 0)])


