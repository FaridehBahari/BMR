import numpy as np
import h5py
import pandas as pd


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


path = '../external/ftrMtrix/var_features.h5'


X = pd.read_hdf(path)
X.shape


columns_to_exclude = [col for col in new_ftrs if col in X.columns]
X_selected = X.drop(columns=columns_to_exclude, errors='ignore') 


X_selected.shape
# Save as H5 file
X_selected.to_hdf('../external/ftrMtrix/var_1372features.h5', key='/X', mode='w', data_columns=True)

