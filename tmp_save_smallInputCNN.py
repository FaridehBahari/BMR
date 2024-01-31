import numpy as np
import matplotlib.pyplot as plt
import json
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, UpSampling1D, Dropout, Flatten, Dense
from tensorflow.keras.models import Model
import tensorflow as tf
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten
from readFtrs_Rspns import set_gpu_memory_limit
import pandas as pd
import pandas as pd
import numpy as np
import h5py
from readFtrs_Rspns import set_gpu_memory_limit
from sklearn.metrics import mean_squared_error
# from sklearn.preprocessing import StandardScaler
from pickle import load


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
        nrow = f['/X/block0_values'].shape[0]//2
        
            
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
                
                return data_batch_X, data_batch_Y
            
            


nn_batch_size = 2000        # Neural network batch size
num_regions_per_sample = 100


path_response = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'
path_scaler = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler.pkl'


X_train, Y_train = data_generator(path_features, path_response, path_scaler, nn_batch_size, num_regions_per_sample)

X_train.shape
(2000, 100, 1500)
Y_train.shape
(2000, 100)
np.save('../../../../Projects/bahari_work/tmp_X_train.npy', X_train)
np.save('../../../../Projects/bahari_work/tmp_Y_train.npy', Y_train)

##################################################################################################################


def data_generator_test(path_features, path_response, path_scaler, nn_batch_size, num_regions_per_sample=100):
    info = pd.read_csv(path_response, sep='\t', index_col='binID')
    info['mutRate'] = info['nMut'] / (info['length'] * info['N'])
    
    scaler = load(open(path_scaler, 'rb'))
    with h5py.File(path_features, 'r') as f:
        nrow = f['/X/block0_values'].shape[0]
        
            
        # print(start_vec[:3])
        # print(end_vec[:3])
        # print(start_vec[-3:])
        # print(end_vec[-3:])
        
        
        while True:
            
            indices = list(range(nrow//2, (nrow-num_regions_per_sample), 
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
                
                return data_batch_X, data_batch_Y
            
            


X_test, Y_test = data_generator_test(path_features, path_response, path_scaler, nn_batch_size, num_regions_per_sample)


X_test.shape
Y_test.shape

np.save('../../../../Projects/bahari_work/tmp_X_test.npy', X_test)
np.save('../../../../Projects/bahari_work/tmp_Y_test.npy', Y_test)

#################################################################
