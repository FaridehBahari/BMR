#################################
#         TODO: preds are count not rate

from collections import namedtuple
from readFtrs_Rspns import create_TestTrain_TwoSources, get_latest_commit_hash, read_response
import os
import ast
import pickle
import pandas as pd
import numpy as np
# from models.GLM_functions import GLM_model_info
# from models.NN_functions import nn_model_info
# from models.GBM_functions import gbm_model_info
from readFtrs_Rspns import scale_train, scale_test, load_data
from performance.assessModels import assess_model
from models.repeated_train_test import save_metrics_summary

import platform
if platform.system() == 'Linux':
    from pybedtools import BedTool
    import pybedtools
    
###########
def load_data_sim_2(sim_setting):
    
    path_X_train = sim_setting['path_X_train']
    path_X_val = sim_setting['path_X_validate']
    
    
    path_Y_train = sim_setting['path_Y_train']
    path_Y_val = sim_setting['path_Y_validate']
    
    scale = ast.literal_eval(sim_setting['scale'])
    
    DSmpl = ast.literal_eval(sim_setting['DSmpl'])
    n_sample = sim_setting['n_sample']
    remove_unMutated = ast.literal_eval(sim_setting['remove_unMutated'])
    
    # load all train 
    X_tr_cmplt, Y_tr_cmplt = load_data(path_X_train, path_Y_train)
    
    # load all val 
    X_val_cmplt, Y_val_cmplt = load_data(path_X_val, path_Y_val)
    
    # restrict features to DP features
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
    
    columns_to_exclude = [col for col in new_ftrs if col in X_tr_cmplt.columns]
    X_tr_cmplt = X_tr_cmplt.drop(columns=columns_to_exclude, errors='ignore') #
    X_val_cmplt = X_val_cmplt.drop(columns=columns_to_exclude, errors='ignore')
    
    # reorder val columns based on train columns
    X_val_cmplt = X_val_cmplt[X_tr_cmplt.columns]
    
    if scale:
        X_tr_cmplt, meanSc, sdSc = scale_train(X_tr_cmplt)
        X_val_cmplt = scale_test(X_val_cmplt, meanSc, sdSc)
    
    
    
    if remove_unMutated:
        Y_tr_cmplt = Y_tr_cmplt[Y_tr_cmplt['nMut'] != 0]
        X_tr_cmplt = X_tr_cmplt.loc[Y_tr_cmplt.index]
        
        Y_val_cmplt = Y_val_cmplt[Y_val_cmplt['nMut'] != 0]
        X_val_cmplt = X_val_cmplt.loc[Y_val_cmplt.index]
    
    
    if DSmpl:
        np.random.seed(0)
        tr_indices = np.random.choice(list(Y_tr_cmplt.index), size=n_sample, replace=False)
        Y_tr_cmplt = Y_tr_cmplt.loc[tr_indices]
        print(f'Down sampling was performed... number of training bins: {Y_tr_cmplt.shape[0]}')
        X_tr_cmplt = X_tr_cmplt.loc[Y_tr_cmplt.index]
    
    
    if (Y_val_cmplt.index != X_val_cmplt.index).all():
        raise ValueError('X_val and Y_val indexes are not the same')
    if (X_tr_cmplt.index != Y_tr_cmplt.index).all():
        raise ValueError('X_train and Y_train indexes are not the same')
        
    return X_tr_cmplt, Y_tr_cmplt, X_val_cmplt, Y_val_cmplt



# def generate_train_valvar_sets(path_var_intervals_Y, path_Y_train, 
#                                path_train_info, val_size, seed_value):
    
#     # read all variable-size bins
#     var_interval_response = read_response(path_var_intervals_Y)
#     var_interval_response = var_interval_response.iloc[np.where(var_interval_response.length >= 20)]
    
#     # sample from variable-size bins to have validation set
#     np.random.seed(seed_value)
#     val_indices = np.random.choice(var_interval_response.index, 
#                                     size=val_size, replace=False)
#     Y_val = var_interval_response.loc[val_indices]
    
#     # read all fixed-size bins
#     Y_train = read_response(path_Y_train)
#     train_info = pd.read_csv(path_train_info, sep = '\t', index_col='binID')
#     train_Y_annotated = pd.concat([Y_train, train_info], axis=1)

#     # remove validation bins from train set
#     filtered_train_Y = train_Y_annotated[~train_Y_annotated['orig_name'].str.contains('|'.join(Y_val.index))]
    
#     return filtered_train_Y, Y_val

def load_data_sim(sim_setting):
    
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

    
    
    path_X_test = sim_setting['path_X_test']
    path_X_train = sim_setting['path_X_train']
    path_Y_test = sim_setting['path_Y_test']
    path_Y_train = sim_setting['path_Y_train']
    scale = ast.literal_eval(sim_setting['scale'])
    DSmpl = ast.literal_eval(sim_setting['DSmpl'])
    n_sample = sim_setting['n_sample']
    remove_unMutated = ast.literal_eval(sim_setting['remove_unMutated'])
    
    
    X_train, Y_train, X_test, Y_test = create_TestTrain_TwoSources(path_X_train, 
                                                               path_Y_train, 
                                                               path_X_test, 
                                                               path_Y_test,
                                                               scale)
    
    columns_to_exclude = [col for col in new_ftrs if col in X_train.columns]
    X_train = X_train.drop(columns=columns_to_exclude, errors='ignore') #

    X_test = X_test.drop(columns=columns_to_exclude, errors='ignore')
    
    
    if remove_unMutated:
        Y_train = Y_train[Y_train['nMut'] != 0]
        X_train = X_train.loc[Y_train.index]
        
        Y_test = Y_test[Y_test['nMut'] != 0]
        X_test = X_test.loc[Y_test.index]
    
    if DSmpl:
        
        np.random.seed(40)
        tr_indices = np.random.choice(list(Y_train.index), size=n_sample, replace=False)
        Y_train = Y_train.loc[tr_indices]
        print(f'Down sampling was performed... number of training bins: {Y_train.shape[0]}')
        X_train = X_train.loc[Y_train.index]
    
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


########################################################################################


def sample_train_valvar(var_interval_response, Y_train, 
                               train_info, val_size, seed_value):
    
    # restrict to > 20nt length variable-size bins
    
    var_interval_response = var_interval_response.iloc[np.where(var_interval_response.length >= 20)]
    
    # sample from variable-size bins to have validation set
    np.random.seed(seed_value)
    val_indices = np.random.choice(var_interval_response.index, 
                                    size=val_size, replace=False)
    Y_val = var_interval_response.loc[val_indices]
    
    # annotate fixed-size bins
    train_Y_annotated = pd.concat([Y_train, train_info], axis=1)

    # remove validation bins from train set
    filtered_train_Y = train_Y_annotated[~train_Y_annotated['orig_name'].str.contains('|'.join(Y_val.index))]
    
    return filtered_train_Y, Y_val


def sample_validations(Y_train, val_size, seed_value):
    
    Y_train_L = Y_train.iloc[np.where(Y_train.length >= 20)]
    
    # sample from variable-size bins to have validation set
    np.random.seed(seed_value)
    val_indices = np.random.choice(Y_train_L.index, 
                                    size=val_size, replace=False)
    Y_val = Y_train_L.loc[val_indices]
    
    # remove validation bins from train set
    filtered_train_Y = Y_train.loc[~Y_train.index.isin(val_indices)]
    
    return filtered_train_Y, Y_val

def val_IDs_fixedElems(bed_tr, bed_val, seed_value, val_size):
    
    bed_val['binID'] = bed_val[3]
    bed_val = bed_val.set_index('binID')
    
    bed_val = bed_val.iloc[np.where(bed_val[2] - bed_val[1] >= 20)]

    np.random.seed(seed_value)
    val_indices = np.random.choice(bed_val.index, 
                                    size=val_size, replace=False)

    bed_val = bed_val.loc[val_indices]
    bedObj_val = BedTool.from_dataframe(bed_val)
    
    bed_tr['binID'] = bed_tr[3]
    bed_tr = bed_tr.set_index('binID')
    bedObj_tr = BedTool.from_dataframe(bed_tr)
    intersection_val = bedObj_tr.intersect(bedObj_val).to_dataframe()
    val_set_binIDs = np.unique(intersection_val.name)
    
    intersection_tr = bedObj_val.intersect(bedObj_tr).to_dataframe()
    train_set_binIDs = np.unique(intersection_tr.name)
    
    return val_set_binIDs, train_set_binIDs


def sample_train_val_fixedSize(Y_train, Y_val, bed_tr, bed_var, seed_value, val_size):
    
    train_set_binIDs, val_set_binIDs = val_IDs_fixedElems(bed_tr, bed_var, seed_value, val_size)
    
    Y_train = Y_train.loc[train_set_binIDs ]
    Y_val = Y_val.loc[val_set_binIDs]
    
    return Y_train, Y_val
    


def repeated_train_test(sim_setting,  X_tr_cmplt, Y_tr_cmplt, X_val_cmplt, Y_val_cmplt,
            make_pred = True, overwrite = True):
    
    
    fixed_size_train = ast.literal_eval(sim_setting['fixed_size_train'])
    path_train_info = sim_setting['path_train_info']
    path_bed_tr = sim_setting['path_bed_tr']
    path_bed_var = sim_setting['path_bed_var']
    
    models = sim_setting['models']
    base_dir = sim_setting['base_dir']
    
    val_size = X_tr_cmplt.shape[0]//5 
    Nr_pair_acc = sim_setting['Nr_pair_acc']

    if path_train_info != '':
        train_info = pd.read_csv(path_train_info, sep = '\t', index_col='binID')
        train_info = train_info.loc[Y_tr_cmplt.index]
    elif fixed_size_train:
        bed_tr = pd.read_csv(path_bed_tr, sep = '\t', header = None)
        bed_val = pd.read_csv(path_bed_var, sep = '\t', header = None)
        
    seed_values = [1, 5, 14, 10, 20, 30, 40, 50, 60, 70, 80, 90, 77, 100, 110]
    
    
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
            
            for i in range(10):
                print(f'.......... repeat number {i+1} of train-test for evaluation of the {name} ......')
                seed_value = np.random.seed(seed_values[i])

                if os.path.exists(f'{save_path_model}/rep_train_test/{name}_M{i+1}_assessment.tsv'):
                    print(f"Skipping iteration {i+1} as the file already exists.")
                    continue

                if path_train_info != '':
                    Y_train, Y_test = sample_train_valvar(Y_val_cmplt, Y_tr_cmplt, 
                                               train_info, val_size, seed_value)
                elif fixed_size_train: 
                    Y_train, Y_test = sample_train_val_fixedSize(Y_tr_cmplt, bed_tr,
                                                                 bed_val, seed_value, val_size)
                    
                else: 
                    Y_train, Y_test = sample_validations(Y_tr_cmplt, val_size, seed_value)

                X_train = X_tr_cmplt.loc[Y_train.index]
                X_test = X_val_cmplt.loc[Y_test.index]
                print(X_train.shape)
                print(X_test.shape)

                common_indices = X_test.index.intersection(X_train.index)

                if not common_indices.empty:
                    raise ValueError(f"Common indices found between X_test and X_train:{common_indices}")
                else:
                    print("No common indices found between X_test and X_train.")

                fitted_Model = fit_model(X_train, Y_train, X_test, Y_test,
                                         m['run_func'], m['predict_func'], make_pred, m['Args'])
                save_func = m['save_func']
                itr = i+1
                save_func(fitted_Model, base_dir, name, iteration=itr, save_model=True)

                Y_pred = fitted_Model.predRates_test
                Y_obs = Y_test.nMut/(Y_test.N * Y_test.length)
                assessments = assess_model(Y_pred, Y_obs, Nr_pair_acc, name, per_element=False)

                path_assessments = f'{save_path_model}/rep_train_test/{name}_M{i+1}_assessment.tsv'
                assessments.to_csv(path_assessments, sep='\t')
        
        print("=============================")
        dir_path = f'{save_path_model}/rep_train_test/'
        save_metrics_summary(dir_path)







# seed_value = 7
# val_size = 5000


# path_bed_val = '../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6'
# bed_tr = pd.read_csv(path_bed_tr, sep = '\t', header = None)


# path_bed_tr = '../external/database/bins/proccessed_bedtools/callable_50k_intergenic_bins.bed12'
# bed_val = pd.read_csv(path_bed_val, sep = '\t', index_col = 3, header = None)



