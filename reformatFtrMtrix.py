import pandas as pd
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import StandardScaler
import pandas as pd
import h5py
from pickle import dump
import numpy as np

def read_bed(path_bed):
    bed = pd.read_csv(path_bed, sep = '\t', header = None)
    # excluded = bed.iloc[np.where((bed[0] == 'chrX') | (bed[0] == 'chrY') | (bed[0] == 'chrM'))]
    
    # bed = bed.iloc[~bed.index.isin(excluded.index)]
    bed['binID'] = bed[3]
    bed = bed.set_index('binID')
    
    return bed



def convert_tsv_to_h5(input_path, output_path):
    # Read TSV file
    X = pd.read_csv(input_path, sep='\t', index_col='binID')
    
    # Handle missing values
    X = X.fillna(0)
    
    # Save as H5 file
    X.to_hdf(output_path, key='/X', mode='w', data_columns=True)





def scale_data(X, scaler=None):
    """ Scale X with robust scaling.
    
    Args:
        X (np.array): feature matrix indexed by binID.
        scaler (RobustScaler): pre-trained scaler. Default is None
        
    Returns:
        np.array: normalized feature matrix.
        RobustScaler: robust scaler fitted with training data,
            only returned when there is no pre-trained scaler.
    
    """
    if scaler is not None:
        return scaler.transform(X)
    else:
        scaler = RobustScaler(copy=False)
        scaler.fit(X)
        return scaler.transform(X), scaler





def standard_scale_data(X, scaler=None):
    """ Scale X with robust scaling.
    
    Args:
        X (np.array): feature matrix indexed by binID.
        scaler (RobustScaler): pre-trained scaler. Default is None
        
    Returns:
        np.array: normalized feature matrix.
        RobustScaler: robust scaler fitted with training data,
            only returned when there is no pre-trained scaler.
    
    """
    if scaler is not None:
        return scaler.transform(X)
    else:
        scaler = StandardScaler(copy=False)
        scaler.fit(X)
        return scaler.transform(X), scaler



# reformat FtrMtrixes
# convert_tsv_to_h5('../external/ftrMtrix/1M_features.tsv', '../external/ftrMtrix/1M_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/pcawg_features.tsv', '../external/ftrMtrix/pcawg_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/var_features.tsv', '../external/ftrMtrix/var_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/100k_features.tsv', '../external/ftrMtrix/100k_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/50k_features.tsv', '../external/ftrMtrix/50k_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/10k_features.tsv', '../external/ftrMtrix/10k_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/1k_features.tsv', '../external/ftrMtrix/1k_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/var_dpIntergenic_features.tsv',
#                   '../external/ftrMtrix/var_dpIntergenic_features.h5')
# convert_tsv_to_h5('../external/ftrMtrix/bedtools/1M_features.tsv',
#                   '../external/ftrMtrix/bedtools/1M_features.h5')

# convert_tsv_to_h5('../external/ftrMtrix/bedtools/FullSet_features.tsv',
#                   '../external/ftrMtrix/bedtools/FullSet_features.h5')

# convert_tsv_to_h5('../external/ftrMtrix/bedtools/100k_features.tsv',
#                   '../external/ftrMtrix/bedtools/100k_features.h5')

# convert_tsv_to_h5('../external/ftrMtrix/bedtools/50k_features.tsv',
#                   '../external/ftrMtrix/bedtools/50k_features.h5')

#############################################################################







# with h5py.File('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5', 'r') as f:
#      X = f['/X/block0_values'] 
#      print(X.shape)
#      Sc_data = scale_data(X)
#      dump(Sc_data[1], open('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler.pkl', 'wb'))


# save the scaler for 1372Ftrs
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


with h5py.File('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5', 'r') as f:
     all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
     selected_cols = [col for col in all_features if col not in new_ftrs]
     X = f['/X/block0_values'] 
     print(X.shape)
     X = X[:, np.where(np.isin(all_features,selected_cols))[0]]
     print(X.shape)
     Sc_data = scale_data(X)
     dump(Sc_data[1], open('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler_1372Ftrs.pkl', 'wb'))





with h5py.File('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5', 'r') as f:
     all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
     selected_cols = [col for col in all_features if col not in new_ftrs]
     X = f['/X/block0_values'] 
     print(X.shape)
     X = X[:, np.where(np.isin(all_features,selected_cols))[0]]
     print(X.shape)
     Sc_data = standard_scale_data(X)
     dump(Sc_data[1], open('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_stdScaler_1372Ftrs.pkl', 'wb'))




###########################################################################
path_bed_1k_cnn = '../external/database/bins/CNN/1k_window.bed'

path_ftrs = '../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5'
path_callable_percent = '../external/output/1k_cnn/pcawg_callable_var.info'

bed = read_bed(path_bed_1k_cnn)
X = pd.read_hdf(path_ftrs, index_col = 'binID')
X =X.loc[bed.index]

callable_percent = pd.read_csv(path_callable_percent, sep = '\t', index_col='binID')
callable_percent = callable_percent.loc[bed.index]
callable_df = callable_percent['callable']

new_X = pd.concat([X,callable_df], axis = 1)
new_X = new_X.loc[bed.index]
new_X.to_hdf('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder_callFtr.h5', key='/X', mode='w'
, data_columns=True)


with h5py.File('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder_callFtr.h5', 'r') as f:
     all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
     selected_cols = [col for col in all_features if col not in new_ftrs]
     X = f['/X/block0_values'] 
     print(X.shape)
     X = X[:, np.where(np.isin(all_features,selected_cols))[0]]
     print(X.shape)
     Sc_data = standard_scale_data(X)
     dump(Sc_data[1], open('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_stdScaler_1373FtrsCALL.pkl', 'wb'))


with h5py.File('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder_callFtr.h5', 'r') as f:
     all_features = np.array([val.decode('utf-8') for val in f['/X/axis0']])
     selected_cols = [col for col in all_features if col not in new_ftrs]
     X = f['/X/block0_values'] 
     print(X.shape)
     X = X[:, np.where(np.isin(all_features,selected_cols))[0]]
     print(X.shape)
     Sc_data = scale_data(X)
     dump(Sc_data[1], open('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler_1373FtrsCALL.pkl', 'wb'))


