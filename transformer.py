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





path_response_tr = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'
path_scaler = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler_1372Ftrs.pkl'
path_validation_set_gbm = '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/GBM_predTest5.tsv'
path_bed_validation = '../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6'
path_bed_train = '../external/database/bins/CNN/1k_window.bed'









validation_bins = get_validation_bins(path_validation_set_gbm, 
                                      path_bed_validation, path_bed_train)

train_info = create_info_train(path_response_tr, validation_bins)
total_n_samples = train_info.shape[0] # 2881044


transformer_model.fit(
    data_generator(path_features, train_info, path_scaler, nn_batch_size, num_regions_per_sample),
    steps_per_epoch=(total_n_samples // (num_regions_per_sample * nn_batch_size)),  # Total smaller batches per epoch
    epochs=400
)



#################################################################################
path_test_response = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_bed_test = '../external/database/bins/CNN/1k_window.bed'
path_test_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'




info_test = create_info_test(path_test_response, path_bed_test, validation_bins)
test_on = 'validation_set'

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
print(f'Spearman correlation for Middle Region: {mae}. p-value: {p_value}')


spearmanr(Y_preds[np.where(obs_df['nMut'] != 0)], Y_obs[np.where(obs_df['nMut'] != 0)])


