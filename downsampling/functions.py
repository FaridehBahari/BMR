import numpy as np
from readFtrs_Rspns import read_response
import os

def downSample_bins(Y_train, Y_val, n_sample, output_folder):
    os.makedirs(output_folder, exist_ok= True)
    # Y_train = read_response(path_Y_train) #(867266, 6)
    # Y_val = read_response(path_Y_val) # (86726, 6)
    val_bins = Y_val.index
    Y_train = Y_train.drop(val_bins) #(780540, 6)
     
    Y_train = Y_train[Y_train['nMut'] != 0] # (773937, 6)
    
   
       
    tr_indices = np.random.choice(list(Y_train.index), size=n_sample, replace=False)
    
    Y_train = Y_train.loc[tr_indices]
    output_path = f'{output_folder}/Pan_Cancer_train_y_Dsamples_{n_sample}.tsv'
    print(Y_train.shape)
    Y_train.to_csv(output_path, sep='\t')
    

n_sample = 50
output_folder = '../external/rawInput/downSampling'
path_Y_train = '../external/rawInput/Pan_Cancer_train_y.tsv'
path_Y_val = '../external/rawInput/validation_sets_10folds/Pan_Cancer_validate_y_fold_1.tsv'

downSample_bins(path_Y_train, path_Y_val, n_sample, output_folder)
