import pandas as pd
import numpy as np
import h5py
# from sklearn.preprocessing import StandardScaler
from pickle import load
import pybedtools



def read_bed(path_bed):
    bed = pd.read_csv(path_bed, sep = '\t', header = None)
    # excluded = bed.iloc[np.where((bed[0] == 'chrX') | (bed[0] == 'chrY') | (bed[0] == 'chrM'))]
    
    # bed = bed.iloc[~bed.index.isin(excluded.index)]
    bed['binID'] = bed[3]
    bed = bed.set_index('binID')
    
    return bed


def get_validation_bins(path_validation_set_gbm, path_bed_validation, path_bed_train):
    
    var_validation = pd.read_csv(path_validation_set_gbm, sep = '\t', index_col = 'binID')
    var_binIDs_validation = var_validation.index
    var_bed = read_bed(path_bed_validation)
    
    var_bed = var_bed.iloc[var_bed.index.isin(var_binIDs_validation)]
    var_bed = pybedtools.BedTool.from_dataframe(var_bed)
    
    train_bed = pybedtools.BedTool(path_bed_train)
    
    # Perform an intersection
    intersection_tr= train_bed.intersect(var_bed, wa=True, u=True).to_dataframe() # wa: Write the original entry in A for each overlap. u: Write original A entry once if any overlaps found.
    non_train_set_binIDs = np.unique(intersection_tr.name)
    
    return non_train_set_binIDs, var_binIDs_validation


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




path_validation_set_gbm = '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/GBM_predTest5.tsv'
path_bed_validation = '../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6'
path_bed_train = '../external/database/bins/CNN/1k_window.bed'


fixed_bins_ov_var, bins_var = get_validation_bins(path_validation_set_gbm, path_bed_validation, path_bed_train)


path_test_response = '../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv'
path_bed_test = '../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6'
path_fixed_bins_response = '../external/BMR/rawInput/responseTabs_cnn/Pan_Cancer/1k_cnn_withInfo.tsv'
path_bed_fixed_bins = '../external/database/bins/CNN/1k_window.bed'



def create_info_test2(path_test_response, path_bed_test, 
                     path_fixed_bins_response, path_bed_fixed_bins,
                     fixed_bins_ov_var, bins_var): 
    
    info_test = pd.read_csv(path_test_response, sep='\t', index_col='binID')
    info_test['mutRate'] = (info_test.nMut / (info_test.length * info_test.N))
    
    bed_test= read_bed(path_bed_test)
    bed_test = bed_test[~bed_test[0].isin(['chrX', 'chrY', 'chrM'])]
    bed_test = bed_test.iloc[np.where(bed_test.index.isin(bins_var))]
    info_test = info_test.loc[bed_test.index]
    
    info_fixed_bins = pd.read_csv(path_fixed_bins_response, sep='\t', index_col='binID')
    bed_train = read_bed(path_bed_fixed_bins)
    info_fixed_bins['mutRate'] = (info_fixed_bins.nMut / (info_fixed_bins.length * info_fixed_bins.N))
    bed_train = bed_train.iloc[~(bed_train.index.isin(fixed_bins_ov_var))]
    bed_train = bed_train[~bed_train[0].isin(['chrX', 'chrY', 'chrM'])]
    
    info_fixed_bins = info_fixed_bins.loc[bed_train.index]
    info_fixed_bins['mutRate'] = (info_fixed_bins.nMut / (info_fixed_bins.length * info_fixed_bins.N))
    
    bed = pybedtools.BedTool.from_dataframe(pd.concat([bed_train , bed_test], axis=0)).sort().to_dataframe()
    info = pd.concat([info_fixed_bins, info_test], axis=0)
    info = info.loc[bed.name]
    
    return info


info = create_info_test2(path_test_response, path_bed_test, 
                     path_fixed_bins_response, path_bed_fixed_bins,
                     fixed_bins_ov_var, bins_var)

path_fixed_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'
path_test_features = '../external/ftrMtrix/var_features.h5'
path_scaler = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler_1372Ftrs.pkl'
nn_batch_size = 23
num_regions_per_sample = 100
middle_region_index = 50

# def test_data_generator2(info, bins_var, path_test_features, path_fixed_features, path_scaler, 
#                         nn_batch_size, num_regions_per_sample, middle_region_index):    
                
#     new_ftrs = ['APOBEC3A', 'E001-DNAMethylSBS', 'E002-DNAMethylSBS', 
#                 'E003-DNAMethylSBS', 
#                  'E004-DNAMethylSBS', 'E005-DNAMethylSBS', 'E006-DNAMethylSBS', 
#                  'E007-DNAMethylSBS', 'E008-DNAMethylSBS', 'E009-DNAMethylSBS', 
#                  'E010-DNAMethylSBS', 'E011-DNAMethylSBS', 'E012-DNAMethylSBS', 
#                  'E013-DNAMethylSBS', 'E014-DNAMethylSBS', 'E015-DNAMethylSBS',
#                  'E016-DNAMethylSBS', 'E017-DNAMethylSBS', 'E018-DNAMethylSBS', 
#                  'E019-DNAMethylSBS', 'E020-DNAMethylSBS', 'E021-DNAMethylSBS',
#                  'E022-DNAMethylSBS', 'E023-DNAMethylSBS', 'E024-DNAMethylSBS', 
#                  'E025-DNAMethylSBS', 'E026-DNAMethylSBS', 'E027-DNAMethylSBS',
#                  'E028-DNAMethylSBS', 'E029-DNAMethylSBS', 'E030-DNAMethylSBS', 
#                  'E031-DNAMethylSBS', 'E032-DNAMethylSBS', 'E033-DNAMethylSBS',
#                  'E034-DNAMethylSBS', 'E035-DNAMethylSBS', 'E036-DNAMethylSBS', 
#                  'E037-DNAMethylSBS', 'E038-DNAMethylSBS', 'E039-DNAMethylSBS', 
#                  'E040-DNAMethylSBS', 'E041-DNAMethylSBS', 'E042-DNAMethylSBS', 
#                  'E043-DNAMethylSBS', 'E044-DNAMethylSBS', 'E045-DNAMethylSBS', 
#                  'E046-DNAMethylSBS', 'E047-DNAMethylSBS', 'E048-DNAMethylSBS', 
#                  'E049-DNAMethylSBS', 'E050-DNAMethylSBS', 'E051-DNAMethylSBS', 
#                  'E052-DNAMethylSBS', 'E053-DNAMethylSBS', 'E054-DNAMethylSBS',
#                  'E055-DNAMethylSBS', 'E056-DNAMethylSBS', 'E057-DNAMethylSBS',
#                  'E058-DNAMethylSBS', 'E059-DNAMethylSBS', 'E061-DNAMethylSBS',
#                  'E062-DNAMethylSBS', 'E063-DNAMethylSBS', 'E065-DNAMethylSBS',
#                  'E066-DNAMethylSBS', 'E067-DNAMethylSBS', 'E068-DNAMethylSBS',
#                  'E069-DNAMethylSBS', 'E070-DNAMethylSBS', 'E071-DNAMethylSBS', 
#                  'E072-DNAMethylSBS', 'E073-DNAMethylSBS', 'E074-DNAMethylSBS', 
#                  'E075-DNAMethylSBS', 'E076-DNAMethylSBS', 'E077-DNAMethylSBS', 
#                  'E078-DNAMethylSBS', 'E079-DNAMethylSBS', 'E080-DNAMethylSBS', 
#                  'E081-DNAMethylSBS', 'E082-DNAMethylSBS', 'E083-DNAMethylSBS', 
#                  'E084-DNAMethylSBS', 'E085-DNAMethylSBS', 'E086-DNAMethylSBS', 
#                  'E087-DNAMethylSBS', 'E088-DNAMethylSBS', 'E089-DNAMethylSBS', 
#                  'E090-DNAMethylSBS', 'E091-DNAMethylSBS', 'E092-DNAMethylSBS', 
#                  'E093-DNAMethylSBS', 'E094-DNAMethylSBS', 'E095-DNAMethylSBS', 
#                  'E096-DNAMethylSBS', 'E097-DNAMethylSBS', 'E098-DNAMethylSBS', 
#                  'E099-DNAMethylSBS', 'E100-DNAMethylSBS', 'E101-DNAMethylSBS', 
#                  'E102-DNAMethylSBS', 'E103-DNAMethylSBS', 'E104-DNAMethylSBS', 
#                  'E105-DNAMethylSBS', 'E106-DNAMethylSBS', 'E107-DNAMethylSBS', 
#                  'E108-DNAMethylSBS', 'E109-DNAMethylSBS', 'E110-DNAMethylSBS', 
#                  'E111-DNAMethylSBS', 'E112-DNAMethylSBS', 'E113-DNAMethylSBS', 
#                  'E114-DNAMethylSBS', 'E115-DNAMethylSBS', 'E116-DNAMethylSBS', 
#                  'E117-DNAMethylSBS', 'E118-DNAMethylSBS', 'E119-DNAMethylSBS', 
#                  'E120-DNAMethylSBS', 'E121-DNAMethylSBS', 'E122-DNAMethylSBS', 
#                  'E123-DNAMethylSBS', 'E124-DNAMethylSBS', 'E125-DNAMethylSBS', 
#                  'E126-DNAMethylSBS', 'E127-DNAMethylSBS', 'E128-DNAMethylSBS', 
#                  'E129-DNAMethylSBS'
#                  # , 'primates_phastCons46way', 
#                  # 'primates_phyloP46way', 'vertebrate_phastCons46way'
#      ]
    
    
#     n_testElems = len(bins_var)
    
    
#     scaler = load(open(path_scaler, 'rb'))
    
#     with h5py.File(path_test_features, 'r') as f:
        
#         all_test_binIDs = np.array([val.decode('utf-8') for val in f['/X/axis1'][:]])
        
                    
        
        
#         with h5py.File(path_fixed_features, 'r') as fixed_f: 
            
#             all_train_binIDs = np.array([val.decode('utf-8') for val in fixed_f['/X/axis1'][:]])
        
        

        
#         all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
#         selected_cols = [col for col in all_features if col not in new_ftrs]
#         Ftrs = np.where(np.isin(all_features,selected_cols))[0]
        
#         print(f'Please wait! The mutRates of {n_testElems} elements are predicting')
        
#         for i in range(0, n_testElems, nn_batch_size):
            
            
#             middle_idx = np.where(info.index.isin(bins_var))[0]
                      
#             chunk_indices = middle_idx[i:i+nn_batch_size]
            
#             subsets = []
            
#             for idx in chunk_indices:
#                 if (idx - middle_region_index < 0) or (idx + middle_region_index > info.shape[0]):
#                     print(f'there is not enough bins before/after element at index {idx}: {info.loc[idx]}')
#                     continue
                
#                 test_binID = np.where(np.isin(all_test_binIDs, info.index[idx]))[0]
#                 X_middle = f['/X/block0_values'][test_binID]
#                 X_middle = X_middle[:, Ftrs]
                
#                 up_bins_fixed = idx - (middle_region_index)+1
#                 up_bins_fixed = np.where(np.isin(all_train_binIDs, info.index[range(up_bins_fixed,up_bins_fixed+middle_region_index-1)]))[0]
                
#                 down_bins_fixed = idx + middle_region_index
#                 down_bins_fixed = np.where(np.isin(all_train_binIDs, info.index[range(idx+1, down_bins_fixed)]))[0]
                
                
                
#                 X_train_up = fixed_f['/X/block0_values'][up_bins_fixed]
#                 X_train_up = X_train_up[:, Ftrs]
                
#                 X_train_down = fixed_f['/X/block0_values'][down_bins_fixed]
#                 X_train_down = X_train_down[:, Ftrs]
                
#                 X_subset = np.vstack((X_train_up, X_middle, X_train_down))
#                 subsets.append(X_subset)
                
#             if len(subsets) > 0:
                
                
#                 raw_data_batch_X = np.stack(subsets, axis=0)
#                 reshaped_raw_data_batch_X = raw_data_batch_X.reshape(-1, raw_data_batch_X.shape[2])
                
#                 scaled_data_batch_X = scaler.transform(reshaped_raw_data_batch_X)
#                 data_batch_X = scaled_data_batch_X.reshape(raw_data_batch_X.shape)
                
                
#                 expected_shape = (nn_batch_size, num_regions_per_sample, len(Ftrs))
#                 if data_batch_X.shape != expected_shape:
#                     raise ValueError(f"Unexpected shape for data_X. Expected {expected_shape}, got {data_batch_X.shape}.")
                
                
#                 yield data_batch_X





def test_data_generator2(info, bins_var, path_test_fixed_features, path_scaler, 
                        nn_batch_size, num_regions_per_sample, middle_region_index):    
                
    new_ftrs = ['APOBEC3A', 'E001-DNAMethylSBS', 'E002-DNAMethylSBS', 
                'E003-DNAMethylSBS', 
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
    
    
    n_testElems = len(bins_var)
    
    
    scaler = load(open(path_scaler, 'rb'))
    
    with h5py.File(path_test_fixed_features, 'r') as f:
        
        all_test_binIDs = np.array([val.decode('utf-8') for val in f['/X/axis1'][:]]) 
        
        all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']]) 
        
        selected_cols = [col for col in all_features if col not in new_ftrs]
        Ftrs = np.where(np.isin(all_features,selected_cols))[0]
        
        
        print(f'Please wait! The mutRates of {n_testElems} elements are predicting')
        
        middle_idxs = np.where(info.index.isin(bins_var))[0]
        
        for i in range(0, n_testElems, nn_batch_size):
            
            
                      
            chunk_indices = middle_idxs[i:i+nn_batch_size]
            
            subsets = []
            
            for idx in chunk_indices:
                if (idx - middle_region_index < 0) or (idx + middle_region_index > info.shape[0]):
                    print(f'there is not enough bins before/after element at index {idx}: {info.loc[idx]}')
                    continue
                
                 
                test_binID = np.where(np.isin(all_test_binIDs, info.index[idx]))[0]
                X_subset = f['/X/block0_values'][(test_binID-middle_region_index)[0]:(test_binID+middle_region_index)[0]]
                X_subset = X_subset[:, Ftrs] 
                
                subsets.append(X_subset)
                
            if len(subsets) > 0:
                
                
                raw_data_batch_X = np.stack(subsets, axis=0)
                reshaped_raw_data_batch_X = raw_data_batch_X.reshape(-1, raw_data_batch_X.shape[2])
                
                scaled_data_batch_X = scaler.transform(reshaped_raw_data_batch_X)
                data_batch_X = scaled_data_batch_X.reshape(raw_data_batch_X.shape)
                
                
                expected_shape = (nn_batch_size, num_regions_per_sample, len(Ftrs))
                if data_batch_X.shape != expected_shape:
                    # if i != len(range(0, n_testElems, nn_batch_size)) -1 :
                        raise ValueError(f"Unexpected shape for data_X. Expected {expected_shape}, got {data_batch_X.shape}.")
                
                
            yield data_batch_X


path_test_fixed_features = '../../../../Projects/bahari_work/ftrMtrix/cnn/tmp_1k_varTest5.h5'

y_pred = transformer_model.predict(test_data_generator2(info, bins_var, 
                                                        path_test_fixed_features, 
                                                        path_scaler, 
                       1, num_regions_per_sample,
                        middle_region_index))






Y_preds = y_pred[:, middle_region_index]



obs = info.loc[bins_var]
Y_obs = (obs.mutRate).values.reshape(-1, 1)   #(111171,)


mse = mean_squared_error(Y_obs, Y_preds)
print(f'Mean Squared Error for Middle Region: {mse}')
mae = np.mean(np.abs(Y_obs - Y_preds))
print(f'Mean Absolute Error for Middle Region: {mae}')

corr, p_value = spearmanr(Y_preds, Y_obs)
print(f'Spearman correlation for Middle Region: {corr}. p-value: {p_value}')


spearmanr(Y_preds[np.where(obs['nMut'] != 0)], Y_obs[np.where(obs['nMut'] != 0)])


