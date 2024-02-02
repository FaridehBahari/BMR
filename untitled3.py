import pandas as pd


fixed = pd.read_hdf('../../../../Projects/bahari_work/ftrMtrix/cnn/1k_cnn_features_bedOrder.h5')

val = pd.read_hdf('../../../../Projects/bahari_work/ftrMtrix/var_features.h5')

all_df = pd.concat([fixed, val], axis=0)

all_df = all_df.loc[info.index]


all_df.to_hdf('../../../../Projects/bahari_work/ftrMtrix/cnn/tmp_1k_varTest5.h5', key='/X', mode='w', data_columns=True))