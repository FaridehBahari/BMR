import numpy as np
import pandas as pd
import h5py
from tensorflow.keras.utils import Sequence
from tensorflow.keras.models import Sequential
from readFtrs_Rspns import set_gpu_memory_limit

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)


class HDF5BatchGenerator(Sequence):
    
    
    def __init__(self, dataset_ftrs, dataset_bins, path_response, file_path, batch_size, num_regions_per_sample):
        self.dataset_ftrs = dataset_ftrs
        self.dataset_bins = dataset_bins
        self.path_response = path_response
        self.file_path = file_path
        self.batch_size = batch_size
        self.num_regions_per_sample = num_regions_per_sample
        self.info = pd.read_csv(self.path_response, sep='\t', index_col='binID')
        self.info['mutRate'] = np.where((self.info['PCAWG_test_genomic_elements'] == 0) |
                                    (self.info['varSize_longer50'] == 0) | 
                                    (self.info['chr'] == 'chrX') | 
                                    (self.info['chr'] == 'chrY') | 
                                    (self.info['chr'] == 'chrM'),
                                    (self.info['nMut'] / (self.info['length'] * self.info['N'])),
                                    -1)
        
        with h5py.File(self.file_path, 'r') as f:
            self.n_samples = f[self.dataset_ftrs].shape[0]
            
            
    def __len__(self):
        return int(np.floor(self.n_samples / (self.batch_size * self.num_regions_per_sample))) 
    
    def __getitem__(self, idx):
        start = idx * self.batch_size * self.num_regions_per_sample
        end = (idx + 1) * self.batch_size * self.num_regions_per_sample
        
        # start = np.random.randint(0, self.n_samples, size=batch_size * num_regions_per_sample)
        # end = start + batch_size * num_regions_per_sample 
        
        with h5py.File(self.file_path, 'r') as f:
            batch_x = f[self.dataset_ftrs][start:end] 
            tmp_binIDs_feature = f[self.dataset_bins][start:end]
        # print('......1.......')    
        # Process your data (e.g., normalization) here, and reshape if necessary for your CNN
        # batch_x = preprocess(batch_x)
        batch_x = batch_x.reshape(batch_size, num_regions_per_sample, batch_x.shape[1])
        # Get the corresponding batch of labels (assumes labels are stored in a similar fashion)
        
        # print('......2.......')
        binIDs_features = np.array([val.decode('utf-8') for val in tmp_binIDs_feature])
        info_subset = self.info.loc[binIDs_features]
        
        batch_y = info_subset['mutRate'].values.reshape((batch_size, num_regions_per_sample)) 
        
        return batch_x, batch_y

# Usage
batch_size = 23
num_regions_per_sample = 100

path_response = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'




# Build and Compile the Neural Network Model
model = build_model_from_config('model_config.json')
model.compile(optimizer='adam', loss=custom_poisson_loss)
print(model.summary())

# Train the Model Using the Generator
nn_batch_size = 23        # Neural network batch size
num_regions_per_sample = 100
dataset_ftrs = '/X/block0_values'
dataset_bins = '/X/axis1'


# Instantiate generator
training_generator = HDF5BatchGenerator(dataset_ftrs, dataset_bins, path_response, path_features, 
                                        batch_size, num_regions_per_sample)



# Train model on dataset
model.fit(training_generator, epochs=100)
