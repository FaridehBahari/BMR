# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 00:06:09 2023

@author: Farideh
"""


import numpy as np
from readFtrs_Rspns import read_response
import os

def create_validation_sets(data_path, num_folds, output_folder):
    os.makedirs(output_folder, exist_ok= True)
    all_Y_intergenic = read_response(data_path)
    all_Y_intergenic = all_Y_intergenic.iloc[np.where(all_Y_intergenic.length >= 20)]
    indices = all_Y_intergenic.index.unique()
    num_indices = len(indices)
    fold_size = num_indices // num_folds
    
    # Create a set of all indices to sample from
    unsampled_indices = set(indices)
    
    for fold in range(num_folds):
        np.random.seed(0)
        fold_indices = np.random.choice(list(unsampled_indices), 
                                        size=fold_size, replace=False)
        unsampled_indices.difference_update(fold_indices)
        
        Y_validate = all_Y_intergenic.loc[fold_indices]
        output_path = f'{output_folder}/Pan_Cancer_validate_y_fold_{fold + 1}.tsv'
        Y_validate.to_csv(output_path, sep='\t')

data_path = '../external/BMR/rawInput/responseTabs/Pan_Cancer/responseTab_intgnic_intervals.tsv'
num_folds = 10  
output_folder = '../external/BMR/procInput/val_sets_10folds'
create_validation_sets(data_path, num_folds, output_folder)
