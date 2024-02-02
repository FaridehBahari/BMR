import pandas as pd
import numpy as np
import h5py
# from sklearn.preprocessing import StandardScaler
from pickle import load
from AdjuscentBins.generators import read_bed
import pybedtools




def get_validation_bins(path_validation_set_gbm, path_bed_validation, path_bed_train):
    
    var_validation = pd.read_csv(path_validation_set_gbm, sep = '\t', index_col = 'binID')
    var_bed = read_bed(path_bed_validation)

    var_bed = var_bed.iloc[var_bed.index.isin(var_validation.index)]
    var_bed = pybedtools.BedTool.from_dataframe(var_bed)

    train_bed = pybedtools.BedTool(path_bed_train)

    # Perform an intersection
    intersection_tr= train_bed.intersect(var_bed, wa=True, u=True).to_dataframe() # wa: Write the original entry in A for each overlap. u: Write original A entry once if any overlaps found.

    non_train_set_binIDs = np.unique(intersection_tr.name)
    
    return non_train_set_binIDs
    

def create_info_train(path_response, validation_bins = None):
    
    info = pd.read_csv(path_response, sep='\t', index_col='binID')
    info['mutRate'] = np.where((info['PCAWG_test_genomic_elements'] == 0)|
                               (~info['chr'].isin(['chrX', 'chrY', 'chrM'])),
                               
                                (info['nMut'] / (info['length'] * info['N'])),
                                -1)
    
    if validation_bins is not None:
        info['mutRate'] = np.where(info.index.isin(validation_bins),
                                    -1, info['mutRate'])
    
    return info



def create_info_test(path_test_response, path_bed_test, validation_bins = None): 
    
    info = pd.read_csv(path_test_response, sep='\t', index_col='binID')
    info['mutRate'] = (info.nMut / (info.length * info.N))
    bed = read_bed(path_bed_test)
    info = info.loc[bed.index]
    
    if validation_bins is not None:
        info['test_on'] = np.where(info.index.isin(validation_bins), 1, 0)
    
    info = info[~info['chr'].isin(['chrX', 'chrY', 'chrM'])]
    
    return info
    
    

def data_generator(path_features, info, path_scaler, nn_batch_size, num_regions_per_sample=100):
    
    new_ftrs = ['APOBEC3A', 'E001-DNAMethylSBS', 'E002-DNAMethylSBS', 'E003-DNAMethylSBS', 
     'E004-DNAMethylSBS', 'E005-DNAMethylSBS', 'E006-DNAMethylSBS', 
     'E007-DNAMethylSBS', 'E008-DNAMethylSBS', 'E009-DNAMethylSBS', 
     'E010-DNAMethylSBS', 'E011-DNAMethylSBS', 'E012-DNAMethylSBS', 
     'E013-DNAMethylSBS', 'E014-DNAMethylSBS', 'E015-DNAMethylSBS',
     'E016-DNAMethylSBS', 'E017-DNAMethylSBS', 'E018-DNAMethylSBS', 
     'E019-DNAMethylSBS', 'E020-DNAMethylSBS', 'E021-DNAMethylSBS',
     'E022-DNAMethylSBS', 'E023-DNAMethylSBS', 'E024-DNAMethylSBS', 
     'E025-DNAMethylSBS', 'E026-DNAMethylSBS', 'E027-DNAMethylSBS',
     'E028-DNAMethylSBS', 'E029-DNAMethylSBS', 'E030-DNAMethylSBS', 
     'E031-DNAMethylSBS', 'E032-DNAMethylSBS', 'E033-DNAMethylSBS',
     'E034-DNAMethylSBS', 'E035-DNAMethylSBS', 'E036-DNAMethylSBS', 
     'E037-DNAMethylSBS', 'E038-DNAMethylSBS', 'E039-DNAMethylSBS', 
     'E040-DNAMethylSBS', 'E041-DNAMethylSBS', 'E042-DNAMethylSBS', 
     'E043-DNAMethylSBS', 'E044-DNAMethylSBS', 'E045-DNAMethylSBS', 
     'E046-DNAMethylSBS', 'E047-DNAMethylSBS', 'E048-DNAMethylSBS', 
     'E049-DNAMethylSBS', 'E050-DNAMethylSBS', 'E051-DNAMethylSBS', 
     'E052-DNAMethylSBS', 'E053-DNAMethylSBS', 'E054-DNAMethylSBS',
     'E055-DNAMethylSBS', 'E056-DNAMethylSBS', 'E057-DNAMethylSBS',
     'E058-DNAMethylSBS', 'E059-DNAMethylSBS', 'E061-DNAMethylSBS',
     'E062-DNAMethylSBS', 'E063-DNAMethylSBS', 'E065-DNAMethylSBS',
     'E066-DNAMethylSBS', 'E067-DNAMethylSBS', 'E068-DNAMethylSBS',
     'E069-DNAMethylSBS', 'E070-DNAMethylSBS', 'E071-DNAMethylSBS', 
     'E072-DNAMethylSBS', 'E073-DNAMethylSBS', 'E074-DNAMethylSBS', 
     'E075-DNAMethylSBS', 'E076-DNAMethylSBS', 'E077-DNAMethylSBS', 
     'E078-DNAMethylSBS', 'E079-DNAMethylSBS', 'E080-DNAMethylSBS', 
     'E081-DNAMethylSBS', 'E082-DNAMethylSBS', 'E083-DNAMethylSBS', 
     'E084-DNAMethylSBS', 'E085-DNAMethylSBS', 'E086-DNAMethylSBS', 
     'E087-DNAMethylSBS', 'E088-DNAMethylSBS', 'E089-DNAMethylSBS', 
     'E090-DNAMethylSBS', 'E091-DNAMethylSBS', 'E092-DNAMethylSBS', 
     'E093-DNAMethylSBS', 'E094-DNAMethylSBS', 'E095-DNAMethylSBS', 
     'E096-DNAMethylSBS', 'E097-DNAMethylSBS', 'E098-DNAMethylSBS', 
     'E099-DNAMethylSBS', 'E100-DNAMethylSBS', 'E101-DNAMethylSBS', 
     'E102-DNAMethylSBS', 'E103-DNAMethylSBS', 'E104-DNAMethylSBS', 
     'E105-DNAMethylSBS', 'E106-DNAMethylSBS', 'E107-DNAMethylSBS', 
     'E108-DNAMethylSBS', 'E109-DNAMethylSBS', 'E110-DNAMethylSBS', 
     'E111-DNAMethylSBS', 'E112-DNAMethylSBS', 'E113-DNAMethylSBS', 
     'E114-DNAMethylSBS', 'E115-DNAMethylSBS', 'E116-DNAMethylSBS', 
     'E117-DNAMethylSBS', 'E118-DNAMethylSBS', 'E119-DNAMethylSBS', 
     'E120-DNAMethylSBS', 'E121-DNAMethylSBS', 'E122-DNAMethylSBS', 
     'E123-DNAMethylSBS', 'E124-DNAMethylSBS', 'E125-DNAMethylSBS', 
     'E126-DNAMethylSBS', 'E127-DNAMethylSBS', 'E128-DNAMethylSBS', 
     'E129-DNAMethylSBS'
     # , 'primates_phastCons46way', 
     # 'primates_phyloP46way', 'vertebrate_phastCons46way'
     ]
    
    
    scaler = load(open(path_scaler, 'rb'))
    with h5py.File(path_features, 'r') as f:
        nrow = f['/X/block0_values'].shape[0]
        
        # just include DP features
        all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
        selected_cols = [col for col in all_features if col not in new_ftrs]
        
        
        
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
                        subset = subset[:, np.where(np.isin(all_features,selected_cols))[0]]
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
                

def test_data_generator(info, path_test_features, path_test_response, path_bed_test, path_scaler, 
                        nn_batch_size, num_regions_per_sample, middle_region_index, test_on): # test_on can be: 'PCAWG_test_genomic_elements' or 'validation_set'
    
    new_ftrs = ['APOBEC3A', 'E001-DNAMethylSBS', 'E002-DNAMethylSBS', 'E003-DNAMethylSBS', 
     'E004-DNAMethylSBS', 'E005-DNAMethylSBS', 'E006-DNAMethylSBS', 
     'E007-DNAMethylSBS', 'E008-DNAMethylSBS', 'E009-DNAMethylSBS', 
     'E010-DNAMethylSBS', 'E011-DNAMethylSBS', 'E012-DNAMethylSBS', 
     'E013-DNAMethylSBS', 'E014-DNAMethylSBS', 'E015-DNAMethylSBS',
     'E016-DNAMethylSBS', 'E017-DNAMethylSBS', 'E018-DNAMethylSBS', 
     'E019-DNAMethylSBS', 'E020-DNAMethylSBS', 'E021-DNAMethylSBS',
     'E022-DNAMethylSBS', 'E023-DNAMethylSBS', 'E024-DNAMethylSBS', 
     'E025-DNAMethylSBS', 'E026-DNAMethylSBS', 'E027-DNAMethylSBS',
     'E028-DNAMethylSBS', 'E029-DNAMethylSBS', 'E030-DNAMethylSBS', 
     'E031-DNAMethylSBS', 'E032-DNAMethylSBS', 'E033-DNAMethylSBS',
     'E034-DNAMethylSBS', 'E035-DNAMethylSBS', 'E036-DNAMethylSBS', 
     'E037-DNAMethylSBS', 'E038-DNAMethylSBS', 'E039-DNAMethylSBS', 
     'E040-DNAMethylSBS', 'E041-DNAMethylSBS', 'E042-DNAMethylSBS', 
     'E043-DNAMethylSBS', 'E044-DNAMethylSBS', 'E045-DNAMethylSBS', 
     'E046-DNAMethylSBS', 'E047-DNAMethylSBS', 'E048-DNAMethylSBS', 
     'E049-DNAMethylSBS', 'E050-DNAMethylSBS', 'E051-DNAMethylSBS', 
     'E052-DNAMethylSBS', 'E053-DNAMethylSBS', 'E054-DNAMethylSBS',
     'E055-DNAMethylSBS', 'E056-DNAMethylSBS', 'E057-DNAMethylSBS',
     'E058-DNAMethylSBS', 'E059-DNAMethylSBS', 'E061-DNAMethylSBS',
     'E062-DNAMethylSBS', 'E063-DNAMethylSBS', 'E065-DNAMethylSBS',
     'E066-DNAMethylSBS', 'E067-DNAMethylSBS', 'E068-DNAMethylSBS',
     'E069-DNAMethylSBS', 'E070-DNAMethylSBS', 'E071-DNAMethylSBS', 
     'E072-DNAMethylSBS', 'E073-DNAMethylSBS', 'E074-DNAMethylSBS', 
     'E075-DNAMethylSBS', 'E076-DNAMethylSBS', 'E077-DNAMethylSBS', 
     'E078-DNAMethylSBS', 'E079-DNAMethylSBS', 'E080-DNAMethylSBS', 
     'E081-DNAMethylSBS', 'E082-DNAMethylSBS', 'E083-DNAMethylSBS', 
     'E084-DNAMethylSBS', 'E085-DNAMethylSBS', 'E086-DNAMethylSBS', 
     'E087-DNAMethylSBS', 'E088-DNAMethylSBS', 'E089-DNAMethylSBS', 
     'E090-DNAMethylSBS', 'E091-DNAMethylSBS', 'E092-DNAMethylSBS', 
     'E093-DNAMethylSBS', 'E094-DNAMethylSBS', 'E095-DNAMethylSBS', 
     'E096-DNAMethylSBS', 'E097-DNAMethylSBS', 'E098-DNAMethylSBS', 
     'E099-DNAMethylSBS', 'E100-DNAMethylSBS', 'E101-DNAMethylSBS', 
     'E102-DNAMethylSBS', 'E103-DNAMethylSBS', 'E104-DNAMethylSBS', 
     'E105-DNAMethylSBS', 'E106-DNAMethylSBS', 'E107-DNAMethylSBS', 
     'E108-DNAMethylSBS', 'E109-DNAMethylSBS', 'E110-DNAMethylSBS', 
     'E111-DNAMethylSBS', 'E112-DNAMethylSBS', 'E113-DNAMethylSBS', 
     'E114-DNAMethylSBS', 'E115-DNAMethylSBS', 'E116-DNAMethylSBS', 
     'E117-DNAMethylSBS', 'E118-DNAMethylSBS', 'E119-DNAMethylSBS', 
     'E120-DNAMethylSBS', 'E121-DNAMethylSBS', 'E122-DNAMethylSBS', 
     'E123-DNAMethylSBS', 'E124-DNAMethylSBS', 'E125-DNAMethylSBS', 
     'E126-DNAMethylSBS', 'E127-DNAMethylSBS', 'E128-DNAMethylSBS', 
     'E129-DNAMethylSBS'
     # , 'primates_phastCons46way', 
     # 'primates_phyloP46way', 'vertebrate_phastCons46way'
     ]
    
    
    
    
    
    scaler = load(open(path_scaler, 'rb'))
    
    with h5py.File(path_test_features, 'r') as f:
        
        # test_on can be: 'PCAWG_test_genomic_elements' or 'validation_set
        if test_on == 'PCAWG_test_genomic_elements':
            pcawg_ov_binIdx = np.where(info['PCAWG_test_genomic_elements'] != 0)[0]
            
        elif test_on == 'validation_set':
            pcawg_ov_binIdx = np.where(info['test_on'] == 1)[0]
            
        n_testElems = pcawg_ov_binIdx.shape[0]
        
        all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
        selected_cols = [col for col in all_features if col not in new_ftrs]
        Ftrs = np.where(np.isin(all_features,selected_cols))[0]
        
        print(f'Please wait! The mutRates of {n_testElems} elements are predicting')
        
        for i in range(0, n_testElems, nn_batch_size):
            chunk_indices = pcawg_ov_binIdx[i:i+nn_batch_size]
            
            subsets = []
            for idx in chunk_indices:
                if (idx - middle_region_index < 0) or (idx + middle_region_index > info.shape[0]):
                    print(f'there is not enough bins before/after element at index {idx}')
                    continue
                X_subset = f['/X/block0_values'][idx-middle_region_index:idx+middle_region_index]
                X_subset = X_subset[:, Ftrs]
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


def read_bed(path_bed):
    bed = pd.read_csv(path_bed, sep = '\t', header = None)
    # excluded = bed.iloc[np.where((bed[0] == 'chrX') | (bed[0] == 'chrY') | (bed[0] == 'chrM'))]
    
    # bed = bed.iloc[~bed.index.isin(excluded.index)]
    bed['binID'] = bed[3]
    bed = bed.set_index('binID')
    
    return bed

def prepare_test_dataY(info_test, nn_batch_size, 
                       num_regions_per_sample, middle_region_index):
    
    info = info_test
    
    pcawg_ov_binIdx = np.where(info['PCAWG_test_genomic_elements'] != 0)[0]
    mask = info.iloc[pcawg_ov_binIdx]
    
    
    n_testElems = pcawg_ov_binIdx.shape[0]
    
    print(f'Please wait! The mutRates of  {n_testElems} elements are predicting')
    
    
    # Initialize an empty list to store the selected rows
    idx_test_windows = []
    
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

