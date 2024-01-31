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
def data_generator(path_features, path_response, path_scaler, nn_batch_size, num_regions_per_sample=100):
    info = pd.read_csv(path_response, sep='\t', index_col='binID')
    info['mutRate'] = np.where((info['PCAWG_test_genomic_elements'] == 0) |
                                (info['varSize_longer50'] == 0), # | 
                                # (info['chr'] == 'chrX') | 
                                # (info['chr'] == 'chrY') | 
                                # (info['chr'] == 'chrM'),
                                (info['nMut'] / (info['length'] * info['N'])),
                                -1)
    
    scaler = load(open(path_scaler, 'rb'))
    with h5py.File(path_features, 'r') as f:
        nrow = f['/X/block0_values'].shape[0]
        
            
        # print(start_vec[:3])
        # print(end_vec[:3])
        # print(start_vec[-3:])
        # print(end_vec[-3:])
        
        
        while True:
            
            indices = list(range(0, (nrow-num_regions_per_sample), 
                                 num_regions_per_sample))
            
            initial_start_vec = np.random.choice(indices, len(indices), replace=False)
            initial_end_vec = (initial_start_vec).copy() + num_regions_per_sample
            
            start_vec = initial_start_vec.copy()
            end_vec = initial_end_vec.copy()
            
            for i in range(0, len(end_vec)-nn_batch_size, nn_batch_size):
                
                # print(f'....  i:{i}')
                
                sts = start_vec[i : i + nn_batch_size]
                ends = end_vec[i: i + nn_batch_size]
                idx = [list(range(start, end)) for start, end in zip(sts, ends)]
                
                
                # Initialize a list to hold the subsets
                subsets = []
                tmp_binIDs_features = []
                tmp_binIDs_feature = []
                for sublist in idx:
                    # Make sure the sublist has precisely 100 indices
                    if len(sublist) == num_regions_per_sample:
                        subset = f['/X/block0_values'][sublist]
                        subsets.append(subset)
                        
                        tmp_binIDs_feature = f['/X/axis1'][sublist]
                        tmp_binIDs_features.append(tmp_binIDs_feature)
                    else:
                        print(f"One of the nn_batch_sizes does not contain {num_regions_per_sample} indices.")
                
                        
                # if i == 0:
                    # print(idx)
                    
                raw_data_batch_X = np.stack(subsets, axis=0)
                reshaped_raw_data_batch_X = raw_data_batch_X.reshape(-1, raw_data_batch_X.shape[2])
                
                scaled_data_batch_X = scaler.transform(reshaped_raw_data_batch_X)
                data_batch_X = scaled_data_batch_X.reshape(raw_data_batch_X.shape)
                
                binIDs_features = [item.decode('utf-8') for slist in tmp_binIDs_features for item in slist]
                info_subset = info.loc[binIDs_features]
                
                # Before proceeding, check the number of elements
                expected_num_elements = nn_batch_size * num_regions_per_sample
                actual_num_elements = len(binIDs_features)
                
                
                if actual_num_elements != expected_num_elements:
                    print(i)
                    print(idx)
                    # print(f'Expected number of elements: {expected_num_elements}')
                    # print(f'Actual number of elements: {actual_num_elements}')
                    # print(info_subset)
                    # print(data_batch_X.shape)
                    raise ValueError('number of elements not passed')
                   
                
                data_batch_Y = info_subset['mutRate'].values.reshape((nn_batch_size, num_regions_per_sample))
                
                yield data_batch_X, data_batch_Y

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

def test_data_generator(path_test_features, path_test_response, path_scaler, 
                        nn_batch_size, num_regions_per_sample, middle_region_index):
    info = pd.read_csv(path_test_response, sep='\t', index_col='binID')
    info['mutRate'] = (info.nMut / (info.length * info.N))
    bed = pd.read_csv(path_bed_test, sep='\t', index_col=3, header=None)
    info = info.loc[bed.index]
    
    scaler = load(open(path_scaler, 'rb'))
    
    with h5py.File(path_test_features, 'r') as f:
        pcawg_ov_binIdx = np.where(info['PCAWG_test_genomic_elements'] != 0)[0]
        n_testElems = pcawg_ov_binIdx.shape[0]
        
        print(f'Please wait! The mutRates of {n_testElems} elements are predicting')
        
        for i in range(0, n_testElems, nn_batch_size):
            chunk_indices = pcawg_ov_binIdx[i:i+nn_batch_size]
            
            subsets = []
            for idx in chunk_indices:
                if (idx - middle_region_index < 0) or (idx + middle_region_index > info.shape[0]):
                    print(f'there is not enough bins before/after element at index {idx}')
                    continue
                X_subset = f['/X/block0_values'][idx-middle_region_index:idx+middle_region_index]
                subsets.append(X_subset)
            if len(subsets) > 0:
                
                
                raw_data_batch_X = np.stack(subsets, axis=0)
                reshaped_raw_data_batch_X = raw_data_batch_X.reshape(-1, raw_data_batch_X.shape[2])
                
                scaled_data_batch_X = scaler.transform(reshaped_raw_data_batch_X)
                data_batch_X = scaled_data_batch_X.reshape(raw_data_batch_X.shape)
                
                
                # expected_shape = (nn_batch_size, num_regions_per_sample, f['/X/block0_values'].shape[1])
                # if data_X.shape != expected_shape:
                #     raise ValueError(f"Unexpected shape for data_X. Expected {expected_shape}, got {data_X.shape}.")
                
                
                yield data_batch_X


def prepare_test_dataY(path_test_response,path_bed_test, nn_batch_size, num_regions_per_sample, middle_region_index):
    info = pd.read_csv(path_test_response, sep = '\t', index_col = 'binID')
    info['mutRate'] = (info.nMut / (info.length * info.N))
    bed = pd.read_csv(path_bed_test, sep = '\t', index_col = 3, header = None)
    info = info.loc[bed.index]
    
    pcawg_ov_binIdx = np.where(info['PCAWG_test_genomic_elements'] != 0)[0]
    mask = info.iloc[pcawg_ov_binIdx]
    
    
    n_testElems = pcawg_ov_binIdx.shape[0]
    
    print(f'Please wait! The mutRates of  {n_testElems} elements are predicting')
    
    
    # Initialize an empty list to store the selected rows
    idx_test_windows = []
    selected_rows = []
    subsets = []
    # Iterate through every 100 rows
    for i in range(n_testElems):
        idx = pcawg_ov_binIdx[i]
        # Get the indices for the current chunk
        if (idx-middle_region_index < 0) | (idx + middle_region_index > info.shape[0]):
            row_removed = info.index[idx]
            print(f'there is int enough bins before/ after {row_removed}')
            continue
        chunk_indices = list(range(idx-middle_region_index, idx + middle_region_index))
        
        if len(chunk_indices) == num_regions_per_sample:
            'passed'
            
        else:
            raise ValueError(f"One of the nn_batch_sizes does not contain {num_regions_per_sample} indices.")
    
        # Append the selected rows to the list
        idx_test_windows.append(idx)
     
    info_subset = info.iloc[idx_test_windows]
    
    mut_rate_array = info_subset['mutRate'].values.reshape((len(info_subset), 1))
    
    return mut_rate_array, info_subset





middle_region_index = 50
y_pred = model.predict(test_data_generator(path_test_features, path_test_response,
                                           nn_batch_size, num_regions_per_sample, 
                                           middle_region_index))

Y_preds = y_pred[:, middle_region_index]
Y_obs, obs_df =prepare_test_dataY(path_test_response,path_bed_test, nn_batch_size, 
                         num_regions_per_sample, middle_region_index)

mean_squared_error(Y_obs, Y_preds)
spearmanr(Y_preds[np.where(obs_df['nMut'] != 0)], Y_obs[np.where(obs_df['nMut'] != 0)])

###############################################################################


