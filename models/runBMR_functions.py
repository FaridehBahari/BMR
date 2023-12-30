#################################
#         TODO: preds are count not rate

from collections import namedtuple
from readFtrs_Rspns import create_TestTrain_TwoSources, get_latest_commit_hash, read_response
import os
import sys 
import ast
import pickle
import pandas as pd
import numpy as np
from models.GLM_functions import GLM_model_info
from models.NN_functions import nn_model_info
from models.GBM_functions import gbm_model_info
from sklearn.model_selection import train_test_split
from readFtrs_Rspns import scale_train, scale_test, load_data, read_feature

###########

def generate_train_valvar_sets(path_var_intervals_Y, path_Y_train, 
                               path_train_info, val_size):
    
    # read all variable-size bins
    var_interval_response = read_response(path_var_intervals_Y)
    var_interval_response = var_interval_response.iloc[np.where(var_interval_response.length >= 20)]
    
    # sample from variable-size bins to have validation set
    val_indices = np.random.choice(var_interval_response.index, 
                                    size=val_size, replace=False)
    Y_val = var_interval_response.loc[val_indices]
    
    # read all fixed-size bins
    Y_train = read_response(path_Y_train)
    train_info = pd.read_csv(path_train_info, sep = '\t', index_col='binID')
    train_Y_annotated = pd.concat([Y_train, train_info], axis=1)

    # remove validation bins from train set
    filtered_train_Y = train_Y_annotated[~train_Y_annotated['orig_name'].str.contains('|'.join(Y_val.index))]
    
    return filtered_train_Y, Y_val



def load_data_sim(sim_setting):
    path_X_test = sim_setting['path_X_test']
    path_X_train = sim_setting['path_X_train']
    path_X_val = sim_setting['path_X_validate']
    path_train_info = sim_setting['path_train_info']
    path_Y_test = sim_setting['path_Y_test']
    path_Y_train = sim_setting['path_Y_train']
    path_Y_val = sim_setting['path_Y_validate']
    val_size = 800
    scale = ast.literal_eval(sim_setting['scale'])
    split_intergenic = ast.literal_eval(sim_setting['split_intergenic'])
    DSmpl = ast.literal_eval(sim_setting['DSmpl'])
    n_sample = sim_setting['n_sample']
    remove_unMutated = ast.literal_eval(sim_setting['remove_unMutated'])
    
    
    X_test, Y_test = load_data(path_X_test, path_Y_test)
    
    Y_train, Y_val = generate_train_valvar_sets(path_Y_val, path_Y_train, 
                                   path_train_info, val_size)
    
    X_train = read_feature(path_X_train, Y_train.index)
    X_val = read_feature(path_X_val, Y_val.index)
    
    # use_features = np.nan
    use_features = np.load('dp_ftrs.npy', allow_pickle=True)
    if use_features is not None:
        X_train = X_train[use_features]
        X_test = X_test[use_features]
        X_val = X_val[use_features]
    
    # reorder test columns based on train columns
    X_test = X_test[X_train.columns]
    X_val = X_val[X_train.columns]
    
    
    
    
    if scale:
        X_train, meanSc, sdSc = scale_train(X_train)
        X_test = scale_test(X_test, meanSc, sdSc)
        X_val = scale_test(X_val, meanSc, sdSc)
    
    
    
    if remove_unMutated:
        Y_train = Y_train[Y_train['nMut'] != 0]
        X_train = X_train.loc[Y_train.index]
        
        Y_test = Y_test[Y_test['nMut'] != 0]
        X_test = X_test.loc[Y_test.index]

            
    if split_intergenic:
        Y_val = read_response(path_Y_val)
        if remove_unMutated:
            Y_val = Y_val[Y_val.nMut != 0]
        X_val = X_train.loc[Y_val.index]
        if (X_val.index != Y_val.index).all():
            raise ValueError('X_val and Y_val indexes are not the same')
        val_bins = Y_val.index
        X_train = X_train.drop(val_bins)
        Y_train = Y_train.drop(val_bins)
        
        X_test = pd.concat([X_test, X_val], axis=0)
        Y_test = pd.concat([Y_test, Y_val], axis=0)
    
    if DSmpl:
        np.random.seed(0)
        tr_indices = np.random.choice(list(Y_train.index), size=n_sample, replace=False)
        Y_train = Y_train.loc[tr_indices]
        print(f'Down sampling was performed... number of training bins: {Y_train.shape[0]}')
        X_train = X_train.loc[Y_train.index]
    
    X_test = pd.concat([X_test, X_val], axis=0)
    Y_test = pd.concat([Y_test, Y_val], axis=0)
    
    if (Y_test.index != X_test.index).all():
        raise ValueError('X_test and Y_test indexes are not the same')
    if (Y_train.index != X_train.index).all():
        raise ValueError('X_train and Y_train indexes are not the same')
    return X_train, Y_train, X_test, Y_test



def fit_model(X_train, Y_train, X_test, Y_test, run_func, predict_func,
              make_pred = True, *args):
    
    model = run_func(X_train, Y_train, args[0])
    length_test_elems = Y_test['length']
    length_train_elems = Y_train['length']
    
    if make_pred:
        predRates_test = predict_func(model, X_test, length_test_elems)
        predRates_train = predict_func(model, X_train, length_train_elems)
    else:
        predRates_train = predRates_test = pd.DataFrame({'predRate': None}, 
                                                        index=['None'])
    
    
    # create a named tuple
    DataTuple = namedtuple('DataTuple', ['model', 'predRates_train', 
                                         'predRates_test', 'run_func',
                                         'predict_func'])
    
    fitted_Model = DataTuple(model=model, predRates_train=predRates_train,
                     predRates_test = predRates_test, run_func = run_func,
                     predict_func = predict_func)
    
    return  fitted_Model


        
def write_readme_file(information, filename):
    # Call the function to get the latest commit hash
    latest_commit_hash = get_latest_commit_hash()
    information['git_commit']= latest_commit_hash
    with open(filename, 'w') as f:
        # Write the dictionary to the file
        f.write('information\n\n')
        for key, value in information.items():
            f.write(f'- {key}: {value}\n')
            


def RUN_BMR(sim_setting,  X_train, Y_train, X_test, Y_test, make_pred = True,
            overwrite = True):
    
    models = sim_setting['models']
    base_dir = sim_setting['base_dir']
        
    for key in models:
       
        m = models[key]
        name = m['save_name']
        
        os.makedirs(f'{base_dir}/{name}/', exist_ok= True)
        readme_file_name = f'{base_dir}/{name}/README.md'
        print(f'@@@@  model: {name}  @@@@')
        params = m['Args']
        save_path_model = f'{base_dir}/{name}/'
        params['path_save'] = f'{save_path_model}models_interval/'
        # check_file_func = m['check_file_func']
        # file_check = check_file_func(base_dir, name)
        # if not os.path.exists(file_check) or sim_setting['overwrite']:
        if not os.path.exists(readme_file_name) or overwrite:
            write_readme_file(m, readme_file_name)
            fitted_Model = fit_model(X_train, Y_train, X_test, Y_test,
                                     m['run_func'], m['predict_func'], make_pred, m['Args'])
            save_func = m['save_func']
            save_func(fitted_Model, base_dir, name, save_model = True)
        
        
        print("=============================")


def RUN_pairRank(sim_setting,  X_train, Y_train, X_test, Y_test, overwrite = True):
    
    models = sim_setting['models']
    base_dir = sim_setting['base_dir']
      
    for key in models:
       
        m = models[key]
        save_name = m['save_name']
        save_path_model = f'{base_dir}/{save_name}/'
        os.makedirs(f'{base_dir}/{save_name}/', exist_ok= True)
        readme_file_name = f'{base_dir}/{save_name}/README.md'
        print(f'@@@@  model: {save_name}  @@@@')
        if not os.path.exists(readme_file_name) or overwrite:
            write_readme_file(m, readme_file_name)
            
            run_func = m['run_func']
            params = m['Args']
            params['path_save'] = f'{save_path_model}models_interval/'
            model_data = run_func(X_train, Y_train, params)
            model = model_data['model']
            params = model_data['NN_hyperparams']
            
            model.save(f'{save_path_model+save_name}_model.h5')
            
            # Save the dictionary to a file using pickle
            with open(f'{base_dir}/{save_name}/{save_name}_params.pkl', 'wb') as f: 
                pickle.dump(params, f)
        
        print("=============================")
        
        

import shutil
import configparser
from simulation_settings import load_sim_settings, config_get

def config_save(sim_file):
    sim_setting = load_sim_settings(sim_file)
    base_dir = sim_setting['base_dir']
    
    sim_config = configparser.ConfigParser()
    sim_config.read(sim_file)
    for model_name in sim_config['models']:
        print(model_name)
        config_file = sim_config['models'][model_name]
        config_model = configparser.ConfigParser()
        config_model.read(config_file)
        save_name = config_get(config_model, 'main', 'method',config_file)
        os.makedirs(f'{base_dir}/{save_name}/', exist_ok= True)
        shutil.copy(config_file, f'{base_dir+ save_name }/{os.path.basename(config_file)}')
        shutil.copy(sim_file, f'{base_dir+ save_name }/{os.path.basename(sim_file)}')
        


def save_train_ids(sim_file, Y_train, model_number):
    sim_setting = load_sim_settings(sim_file)
    base_dir = sim_setting['base_dir']
    model_name = list((sim_setting['models']).keys())[0]
    os.makedirs(f'{base_dir}/{model_name}_{model_number}/', exist_ok= True)
    # Save the DataFrame index to a .npy file
    index_array = Y_train.index.to_numpy()
    np.save(f'{base_dir}/{model_name}_{model_number}/{model_name}_{model_number}_trainBins.npy',
            index_array, allow_pickle=True)
        