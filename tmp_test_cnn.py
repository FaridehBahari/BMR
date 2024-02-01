import pandas as pd
import numpy as np
import h5py
from readFtrs_Rspns import set_gpu_memory_limit
from sklearn.metrics import mean_squared_error
# from sklearn.preprocessing import StandardScaler
from pickle import load
import json
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, UpSampling1D, Dropout, Flatten, Dense
from tensorflow.keras.models import Model
import tensorflow as tf
from scipy.stats import spearmanr
from AdjuscentBins.generators import data_generator, test_data_generator, prepare_test_dataY



# df = pd.DataFrame({
#     'Column1': range(1, 5),
#     'Column2': range(5, 9),
#     'Column3': range(9, 13)
# })

# # Assuming 'df' is your DataFrame
# nrow, _ = df.shape

# # Define the other dimensions
# dim_2 = 2
# dim_1 = int(nrow/dim_2)  # variable size

# dim_3 = 3

# # Reshape the DataFrame to a 3D NumPy array
# array_3d = df.to_numpy().reshape((dim_1, dim_2, dim_3))



####################################################################


### from cnn.py


def build_model_from_config(config_path):
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)
        
    input_shape = tuple(config['input_shape'])
    inputs = Input(shape=input_shape)
    
    x = inputs
    for layer_conf in config['layers']:
        if layer_conf['type'] == 'Conv1D':
            # Select relevant parameters for Conv1D
            conv_params = {k: v for k, v in config['default_conv_params'].items()}
            conv_params.update({k: v for k, v in layer_conf.items() if k != 'type'})
            x = Conv1D(**conv_params)(x)
        elif layer_conf['type'] == 'MaxPooling1D':
            x = MaxPooling1D(**{k: v for k, v in layer_conf.items() if k != 'type'})(x)
        elif layer_conf['type'] == 'UpSampling1D':
            x = UpSampling1D(**{k: v for k, v in layer_conf.items() if k != 'type'})(x)
        elif layer_conf['type'] == 'Dropout':
            x = Dropout(**{k: v for k, v in layer_conf.items() if k != 'type'})(x)
            
    model = Model(inputs=inputs, outputs=x)
    return model



# def custom_poisson_loss(y_true, y_pred):
#     mask = tf.cast(tf.not_equal(y_true, -1.0), tf.float32)  # Mask for available rates
#     loss = tf.keras.losses.poisson(y_true, y_pred)  # Shape [batch_size]
    
#     # Expand the dimensions of loss to match the mask
#     loss = tf.expand_dims(loss, axis=-1)  # Now loss shape becomes [batch_size, 1]
    
#     loss *= mask  # Element-wise multiplication, broadcasting loss to match mask shape
#     # Return the loss without summing and averaging across the sequence
#     return loss  # Shape [batch_size, sequence_length]

def custom_poisson_loss(y_true, y_pred):
      mask = tf.cast(tf.not_equal(y_true, -1.0), tf.float32)  # Mask for available rates
      
      # Ensure y_pred is greater than 0 to avoid log(0)
      y_pred_safe = tf.maximum(y_pred, tf.keras.backend.epsilon())
      
      # Compute element-wise Poisson loss
      loss = y_pred_safe - y_true * tf.math.log(y_pred_safe)
      
      loss *= mask  # Element-wise multiplication, broadcasting loss to match mask shape
      
      # Return the loss without summing and averaging across the sequence
      return loss  # Shape [batch_size, sequence_length]


####### finish ###


########################################################################
gpu_fraction = 0.05
set_gpu_memory_limit(gpu_fraction)


# Build and Compile the Neural Network Model
model = build_model_from_config('model_config.json')
model.compile(optimizer='adam', loss=custom_poisson_loss)
print(model.summary())

# Train the Model Using the Generator
nn_batch_size = 23        # Neural network batch size
num_regions_per_sample = 100


path_response = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'
path_scaler = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler.pkl'

model.fit(
    data_generator(path_features, path_response, path_scaler, nn_batch_size, num_regions_per_sample),
    steps_per_epoch=(2881044 // (num_regions_per_sample * nn_batch_size)),  # Total smaller batches per epoch
    epochs=10
)


##################################################
# path_test_response = '../external/BMR/rawInput/responseTabs_bedtools/Pan_Cancer/tmp_test_100bp_cnn_withInfo.tsv'
# path_test_features = '../external/ftrMtrix/tmpTEST100bp_cnn_features.h5'

path_test_response = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_bed_test = '../external/database/bins/CNN/1k_window.bed'
path_test_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'





middle_region_index = 50
y_pred = model.predict(test_data_generator(path_test_features, path_test_response, 
                                           path_bed_test, path_scaler,
                                           nn_batch_size, num_regions_per_sample, 
                                           middle_region_index))

Y_preds = y_pred[:, middle_region_index]
Y_obs, obs_df =prepare_test_dataY(path_test_response,path_bed_test, nn_batch_size, 
                         num_regions_per_sample, middle_region_index)

mean_squared_error(Y_obs, Y_preds)
spearmanr(Y_preds[np.where(obs_df['nMut'] != 0)], Y_obs[np.where(obs_df['nMut'] != 0)])

###############################################################################



