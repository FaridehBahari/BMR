import pandas as pd
import numpy as np
import h5py

path = '../../BMR_proj/external/rawInput/train_feature.hdf5'

f = h5py.File(path, 'r')
X = f['/X/block0_values'] # feature values
ftr_names = f['/X/block0_items'][:] # feature names
ftr_names = np.array([val.decode('utf-8') for val in ftr_names])

bins = f['/X/axis1'][:]
binID = np.array([val.decode('utf-8') for val in bins])

X_dp = pd.DataFrame(X, index=binID, columns=ftr_names)   
f.close()



# Remove the 'GERP' column
X_dp = X_dp.drop('GERP', axis=1, errors='ignore')

column_mapping = {
    
    'pphastCons46way': 'primates_phastCons46way',
    'pphyloP46way': 'primates_phyloP46way',
    'vphastCons46way': 'vertebrate_phastCons46way',
    'A549_small': 'ENCFF460FMH',
    'Caki2_small': 'ENCFF648KNC',
      'G401_small':'ENCFF820URQ',
    'LNCaPcloneFGC_small': 'ENCFF713UVF',
    'NCI-H460_small': 'ENCFF807HEX',
      'Panc1_small': 'ENCFF783PDP',
     'RPMI-7951_small': 'ENCFF760OPW',
      'SJCRH30_small': 'ENCFF025LMM',
    'SK-MEL-5_small': 'ENCFF364HBJ' ,
    'SK-N-DZ_small': 'ENCFF110KIR',
    'SK-N-MC_small': 'ENCFF978RIC',
     'T47D_small': 'ENCFF411JKH'
}

# Rename columns based on the mapping
X_dp.rename(columns=column_mapping, inplace=True)

path2 = '../external/ftrMtrix/var_features.h5'
f = h5py.File(path2, 'r')
X = f['/X/block0_values'] # feature values
ftr_names = f['/X/block0_items'][:] # feature names
ftr_names = np.array([val.decode('utf-8') for val in ftr_names])

bins = f['/X/axis1'][:]
binID = np.array([val.decode('utf-8') for val in bins])

x_var = pd.DataFrame(X, index=binID, columns=ftr_names)   
f.close()

 


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
 'E129-DNAMethylSBS' ]
x_var = x_var.drop(new_ftrs, axis=1, errors='ignore')



#######################################3
y_var = pd.read_csv('../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv',
                   sep = '\t',
                   index_col = 'binID')

y_var = y_var.iloc[np.where(y_var.nMut !=0)]
x_var=x_var.loc[y_var.index]

y_dp = pd.read_csv('../../BMR_proj/external/rawInput/Pan_Cancer_train_y.tsv',
                   sep = '\t',
                   index_col = 'binID')


# just include mutated bins
y_dp = y_dp.iloc[np.where(y_dp.nMut !=0)]
X_dp =X_dp.loc[y_dp.index]


Y_all = pd.concat([y_dp, y_var], axis = 0)

x_all = pd.concat([x_var, X_dp], axis=0)
x_all = x_all.fillna(0)
X_all = x_all.loc[Y_all.index]


################ save feature and response tables #############
Y_all.to_csv('../external/BMR/rawInput/responseTabs/Pan_Cancer/var_dpIntergenic_.tsv',
             sep = '\t')
# Save the DataFrame to an HDF5 file with the specified structure
file_path = '../external/ftrMtrix/var_dpIntergenic_features.h5'
X_all.to_hdf(file_path, key='/X', mode='w',
          # , format='table',
          data_columns=True)
print(x_all)
# # Confirm the file has been saved correctly
# loaded_df = pd.read_hdf(file_path, key='/X')

