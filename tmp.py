

import pandas as pd
import numpy as np
import h5py
from readFtrs_Rspns import set_gpu_memory_limit

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)

path_response = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features.h5'



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
import json
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, UpSampling1D, Dropout, Flatten, Dense
from tensorflow.keras.models import Model
import tensorflow as tf

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



def custom_poisson_loss(y_true, y_pred):
    mask = tf.cast(tf.not_equal(y_true, -1.0), tf.float32)  # Mask for available rates
    loss = tf.keras.losses.poisson(y_true, y_pred)  # Shape [batch_size]
    
    # Expand the dimensions of loss to match the mask
    loss = tf.expand_dims(loss, axis=-1)  # Now loss shape becomes [batch_size, 1]
    
    loss *= mask  # Element-wise multiplication, broadcasting loss to match mask shape
    # Return the loss without summing and averaging across the sequence
    return loss  # Shape [batch_size, sequence_length]


####### finish ###
import pandas as pd
import numpy as np
import h5py

def data_generator(path_features, path_response, nn_batch_size, num_regions_per_sample=100):
    info = pd.read_csv(path_response, sep='\t', index_col='binID')
    info['mutRate'] = np.where((info['PCAWG_test_genomic_elements'] == 0) | (info['varSize_longer50'] == 0),
                                (info['nMut'] / (info['length'] * info['N'])),
                                -1)
    
    
    with h5py.File(path_features, 'r') as f:
        nrow = f['/X/block0_values'].shape[0]
        indices = list(range(0, nrow, num_regions_per_sample))
        
        initial_start_vec = np.random.choice(indices, len(indices), replace=False)
        initial_end_vec = (initial_start_vec).copy() + num_regions_per_sample
        
        # Handle case where the end index exceeds the number of rows
        if np.max(initial_end_vec) > nrow:
            valid_idx = np.where(initial_end_vec <= nrow)[0]
            start_vec = initial_start_vec[valid_idx]
            end_vec = initial_end_vec[valid_idx]
        else:
            start_vec = initial_start_vec.copy()
            end_vec = initial_end_vec.copy()
            
        print(start_vec[:3])
        print(end_vec[:3])
        print(start_vec[-3:])
        print(end_vec[-3:])
        
        
        while True:
            for i in range(0, len(end_vec)):
                print(f'!!!!!!!!!   ....  i:{i}')
                
                sts = start_vec[i : i + nn_batch_size]
                ends = end_vec[i: i + nn_batch_size]
                idx = [list(range(start, end)) for start, end in zip(sts, ends)]
                
                print(idx[:5])
                print('....')
                print(idx[-5:])
                # Flatten the list of lists
                data_batch_X = np.array([f['/X/block0_values'][s:e] for s, e in zip(sts, ends)])
                
                print('******2******')
                print(data_batch_X.shape)
                
                tmp_binIDs_features = [f['/X/axis1'][s:e] for s, e in zip(sts, ends)]
                binIDs_features = [item.decode('utf-8') for sublist in tmp_binIDs_features for item in sublist]
                info_subset = info.loc[binIDs_features]
                
                data_batch_Y = info_subset['mutRate'].values.reshape((nn_batch_size, num_regions_per_sample))
                print('******3******')
                print(data_batch_Y.shape)            
                yield data_batch_X, data_batch_Y  # Adjust for your model's needs





# Build and Compile the Neural Network Model
model = build_model_from_config('model_config.json')
model.compile(optimizer='adam', loss=custom_poisson_loss)
print(model.summary())

# Train the Model Using the Generator
nn_batch_size = 32        # Neural network batch size
num_regions_per_sample = 100


# model.fit(
#     data_generator(path_features, path_response, large_batch_size, nn_batch_size, num_regions_per_sample),
#     steps_per_epoch=(30000000 // (num_regions_per_sample * nn_batch_size)),  # Total smaller batches per epoch
#     epochs=1
# )


model.fit(
    data_generator(path_features, path_response, nn_batch_size, num_regions_per_sample),
    steps_per_epoch=2,  # Total smaller batches per epoch
    epochs=3
)



# import numpy as np
# import pandas as pd

# path_test_response = '../external/BMR/rawInput/responseTabs_bedtools/Pan_Cancer/tmp_test_100bp_cnn_withInfo.tsv'
# path_test_features = '../external/ftrMtrix/tmpTEST100bp_cnn_features.h5'


# def prepare_test_data(path_test_features, path_test_response, num_regions_per_sample):
#     info = pd.read_csv(path_test_response, sep = '\t', index_col = 'binID')
#     info['mutRate'] = (info.nMut / (info.length * info.N))

#     X = pd.read_hdf(path_test_features)
#     info_subset = info.loc[X.index]

#     mut_rate_array = info_subset['mutRate'].values.reshape((len(info_subset), 1))

#     n_samples = (int(X.shape[0]/num_regions_per_sample))


#     X = X.values
#     large_batch_dataX = X.reshape(n_samples, num_regions_per_sample, X.shape[1])

#     Y = mut_rate_array[:]
#     large_batch_dataY = Y.reshape(n_samples, num_regions_per_sample)
    
#     return large_batch_dataX, large_batch_dataY




# y_pred = model.predict(large_batch_dataX)


# reshaped_large_batch_dataY = large_batch_dataY.reshape((200, 100, 1))
# np.mean(np.abs(reshaped_large_batch_dataY - y_pred))




# middle_region_index = 50  # Adjust based on your data specifics
# y_test_middle_region = large_batch_dataY[:, middle_region_index]
# y_pred_middle_region = y_pred[:, middle_region_index]
# np.mean(np.abs(y_test_middle_region - y_pred_middle_region))


# from sklearn.metrics import mean_squared_error
# mean_squared_error(y_test_middle_region, y_pred_middle_region)



#############################################################
import pandas as pd
import numpy as np
import pybedtools

path_info = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_var_bed = '../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6'
path_1k_cnn_bed = '../external/database/bins/CNN/1k_window.bed'

var_bed = pd.read_csv(path_var_bed, sep = '\t', index_col=3, header=None)
var_bed = var_bed.iloc[np.where(var_bed[2] - var_bed[1] >=50)]

info = pd.read_csv(path_info, sep = '\t', index_col='binID')
info_var0 = info.iloc[np.where(info.varSize_longer50 == 0)]


info_bed = pd.read_csv(path_1k_cnn_bed, sep = '\t', index_col=3, header=None)
info_bed = info_bed.loc[info_var0.index]


var_bed_obj = pybedtools.BedTool.from_dataframe(var_bed)
info_bed_obj = pybedtools.BedTool.from_dataframe(info_bed)

var_bed_obj = var_bed_obj.subtract(info_bed_obj)

# Perform an intersection
var_bed_obj = var_bed_obj.intersect(info_bed_obj, wa=True, u=True)  # wa: Write the original entry in A for each overlap. u: Write original A entry once if any overlaps found.

# Create a set of names from df2 that intersect with df1
intersecting_names = set([feature.name for feature in intersection])

# Label df1: 1 if the row's name is in the intersection set; otherwise 0
df1['label'] = df1['name'].apply(lambda x: 1 if x in intersecting_names else 0)




# def data_generator(path_features, path_response, large_batch_size, nn_batch_size, num_regions_per_sample=100):
#     # Load response info
#     info = pd.read_csv(path_response, sep='\t', index_col='binID')
#     info['mutRate'] = np.where((info['PCAWG_test_genomic_elements'] == 0) | (info['varSize_longer50'] == 0),
#                                (info['nMut'] / (info['length'] * info['N'])),
#                                -1)
    
#     # Open HDF5 file for reading
#     with h5py.File(path_features, 'r') as f:
#         nrow = f['/X/block0_values'].shape[0]
#         indices = list(range(0, nrow, nn_batch_size))
        
#         # Choose random start positions for each region
#         initial_start_vec = np.random.choice(indices, len(indices), replace=False)
#         initial_end_vec = initial_start_vec + nn_batch_size
        
#         # Handle case where the end index exceeds the number of rows
#         if np.max(initial_end_vec) > nrow:
#             valid_idx = np.where(initial_end_vec < nrow)[0]
#             start_vec = initial_start_vec[valid_idx]
#             end_vec = initial_end_vec[valid_idx]
            
#         else:
#             start_vec = initial_start_vec
#             end_vec = initial_end_vec
        
#         print(start_vec[:3])
#         print(end_vec[:3])
#         print(start_vec[-3:])
#         print(end_vec[-3:])
        
#         while True:
#             for batch_offset in range(0, len(start_vec), large_batch_size):
#                 idx = []
                
#                 for batch_index in range(batch_offset, min(batch_offset + large_batch_size, len(start_vec))):
#                     start = start_vec[batch_index]
#                     end = min(end_vec[batch_index], start + num_regions_per_sample)
#                     idx.extend(range(start, end))
                
#                 if not idx:  # This should not happen but adds robustness
#                     break
                
                
#                 print('******1*****')
#                 print(idx[:5])
#                 print(idx[-5:])
                
#                 # Select subsets from the dataset for each start_vec
#                 subsets = [f['/X/block0_values'][s:e] for s, e in zip(start_vec[batch_offset:batch_offset + large_batch_size], end_vec[batch_offset:batch_offset + large_batch_size])]

#                 # Vertically stack the subsets until reaching large_batch_size
#                 large_batch_dataX = np.vstack(subsets)

#                 print('******2******')
#                 print(large_batch_dataX.shape)

#                 tmp_binIDs_features = [f['/X/axis1'][s:e] for s, e in zip(start_vec[batch_offset:batch_offset + large_batch_size], end_vec[batch_offset:batch_offset + large_batch_size])]
#                 binIDs_features = [item.decode('utf-8') for sublist in tmp_binIDs_features for item in sublist]
                
                
#                 info_subset = info.loc[binIDs_features]
                
#                 mut_rate_array = info_subset['mutRate'].values
#                 print('******3******')
#                 # Reshape data for the CNN
#                 reshaped_size = min(large_batch_size * num_regions_per_sample, len(idx)) // num_regions_per_sample
#                 large_batch_dataX = large_batch_dataX.reshape((reshaped_size, 
#                                                                num_regions_per_sample,
#                                                                f['/X/block0_values'].shape[1]))
#                 large_batch_dataY = mut_rate_array.reshape((reshaped_size, num_regions_per_sample))
                
#                 print('******4******')
#                 print(large_batch_dataX.shape)
                
#                 # Yield smaller batches
#                 for j in range(0, reshaped_size, nn_batch_size):
#                     print(f'====j:{j}====')
#                     nn_batch_dataX = large_batch_dataX[j:(j + nn_batch_size)]
#                     print(nn_batch_dataX.shape)
#                     nn_batch_dataY = large_batch_dataY[j:(j + nn_batch_size)]
#                     yield nn_batch_dataX, nn_batch_dataY




# model.fit(
#     data_generator(path_features, path_response, large_batch_size, nn_batch_size, num_regions_per_sample),
#     steps_per_epoch=6,  # Total smaller batches per epoch
#     epochs=10
# )

