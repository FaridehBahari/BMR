# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 00:06:09 2023

@author: Farideh
"""
# import numpy as np
# from readFtrs_Rspns import read_response

# path_Y_all_intergenic = '../external/rawInput/Pan_Cancer_train_y.tsv'
# S = 100000
# all_Y_intergenic = read_response(path_Y_all_intergenic)
# bins = np.random.choice(all_Y_intergenic.index.unique(), size=S, replace=False)
# Y_validate = all_Y_intergenic.loc[bins]
# Y_validate.to_csv('../external/rawInput/Pan_Cancer_validate_y.tsv', sep = '\t')

import numpy as np
from readFtrs_Rspns import read_response
import os

def create_validation_sets(data_path, num_folds, output_folder):
    os.makedirs(output_folder, exist_ok= True)
    all_Y_intergenic = read_response(data_path)
    indices = all_Y_intergenic.index.unique()
    num_indices = len(indices)
    fold_size = num_indices // num_folds
    
    # Create a set of all indices to sample from
    unsampled_indices = set(indices)
    
    for fold in range(num_folds):
        fold_indices = np.random.choice(list(unsampled_indices), size=fold_size, replace=False)
        unsampled_indices.difference_update(fold_indices)
        
        Y_validate = all_Y_intergenic.loc[fold_indices]
        output_path = f'{output_folder}/Pan_Cancer_validate_y_fold_{fold + 1}.tsv'
        Y_validate.to_csv(output_path, sep='\t')

data_path = '../external/rawInput/Pan_Cancer_train_y.tsv'
num_folds = 10  # You can adjust the number of folds as needed
output_folder = '../external/rawInput/validation_sets_10folds'
create_validation_sets(data_path, num_folds, output_folder)
