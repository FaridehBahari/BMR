import pandas as pd


def convert_tsv_to_h5(input_path, output_path):
    # Read TSV file
    X = pd.read_csv(input_path, sep='\t', index_col='binID')
    
    # Handle missing values
    X = X.fillna(0)
    
    # Save as H5 file
    X.to_hdf(output_path, key='/X', mode='w', data_columns=True)

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
convert_tsv_to_h5('../external/ftrMtrix/bedtools/1M_features.tsv',
                  '../external/ftrMtrix/bedtools/1M_features.h5')

convert_tsv_to_h5('../external/ftrMtrix/bedtools/FullSet_features.tsv',
                  '../external/ftrMtrix/bedtools/FullSet_features.h5')

convert_tsv_to_h5('../external/ftrMtrix/bedtools/100k_features.tsv',
                  '../external/ftrMtrix/bedtools/100k_features.h5')

convert_tsv_to_h5('../external/ftrMtrix/bedtools/50k_features.tsv',
                  '../external/ftrMtrix/bedtools/50k_features.h5')

#############################################################################

from sklearn.preprocessing import RobustScaler

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


import pandas as pd
import h5py
from pickle import dump


with h5py.File('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5', 'r') as f:
     X = f['/X/block0_values'] 
     print(X.shape)
     Sc_data = scale_data(X)
     dump(Sc_data[1], open('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_RobustScaler.pkl', 'wb'))


# save the scaler
