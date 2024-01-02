# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 16:13:00 2023

@author: Farideh
"""
import pandas as pd
import sys
import os
import ast
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
sys.path.append('models/')
from readFtrs_Rspns import read_response
from scipy.stats import spearmanr
import numpy as np
from sklearn.metrics import mean_squared_error
from simulation_settings import load_sim_settings_perBatchPerf
from models.runBMR_functions import load_data
import h5py

from tensorflow.keras.models import load_model

import warnings
def read_pred(path_pred):
    Y_pred = pd.read_csv(path_pred, sep = "\t", header=0, index_col='binID',
                         usecols=['binID', 'predRate'])
    return Y_pred

def read_obs(path_Y):
    Y = read_response(path_Y)
    Y = Y[Y['nMut'] != 0]
    Y_obs = Y['obsRates']
    Y_obs = Y_obs.to_frame(name = 'obsRates')
    return Y_obs

def split_by_element(df, element_type):
    df_element = df.loc[df.index.str.contains(element_type)]
    return df_element


def calc_corr(pred, obs):
    corr, p_value = spearmanr(pred, obs)
    
    return corr


def generate_pair_indices(N_pairs, N_tot):
    
    # retain uniqe rows
    pairs = np.unique(np.column_stack((np.random.choice(N_tot, N_pairs, replace=True),
                                       np.random.choice(N_tot, N_pairs, replace=True))), axis=0)
    
    # just retain unique pairs per row (pair i,i is not alloweded)
    non_duplicate_rows = pairs[~np.all(pairs[:, 1:] == pairs[:, :-1], axis=1)]
    
    return non_duplicate_rows 

def calc_Pairs_acc(Nr_pair_acc, obs, pred):
    
    pairs = generate_pair_indices(Nr_pair_acc, obs.shape[0])
    pred_res = pred.iloc[pairs[:, 0]].values > pred.iloc[pairs[:, 1]].values
    obs_res = obs.iloc[pairs[:, 0]].values > obs.iloc[pairs[:, 1]].values
    
    return np.mean(pred_res == obs_res)


def assess_model_element_type(Y_pred, Y_obs, Nr_pair_acc, model_name, elem):
    
    if elem == "intergenic":
        obs_elem = split_by_element(Y_obs,  r'[v]\d')
        pred_elems = split_by_element(Y_pred,  r'[v]\d')
    else:
        obs_elem = split_by_element(Y_obs, elem)
        pred_elems = split_by_element(Y_pred, elem)
    
    pred_elem = pred_elems[~np.isinf(pred_elems)]
    pred_elem = pred_elem[~np.isnan(pred_elem)]
    if pred_elem.shape[0] != pred_elems.shape[0]:
        warnings.warn(f"{elem} prediction contains NaN or infinite values. These elements will be discarded.", stacklevel=1)
    obs_elem = obs_elem.loc[pred_elem.index]
    
    corr_elem =  calc_corr(pred_elem, obs_elem)
    acc_elem = calc_Pairs_acc(Nr_pair_acc, obs_elem, pred_elem)
    mse_elem = mean_squared_error(obs_elem, pred_elem)
    
    return corr_elem, acc_elem, mse_elem
    

def assess_model(Y_pred, Y_obs, Nr_pair_acc, model_name, per_element=True):
    
    acc_name = f"acc_{model_name}"
    mse_name = f"mse_{model_name}"
    corr_name  = f"corr_{model_name}"
    
    if per_element:
        elems = ["gc19_pc.cds", "enhancers", "gc19_pc.3utr", "gc19_pc.5utr",
                 "gc19_pc.promCore", "gc19_pc.ss", "lncrna.ncrna", 
                 "lncrna.promCore"]
        if sum(Y_pred.index.str.contains(r'[v]\d')) != 0:
            elems.append("intergenic")
        
        results = []
        for elem in elems:
            corr_elem, acc_elem, mse_elem = assess_model_element_type(Y_pred, 
                                                                      Y_obs, Nr_pair_acc,
                                                                      model_name,
                                                                      elem)
            
            results.append({'Element': elem, acc_name : acc_elem,
                            corr_name : corr_elem, 
                            mse_name : mse_elem})
        performances = pd.DataFrame(results).set_index('Element').pivot_table(index=None, columns='Element')
        
    else:
        Y_obs = Y_obs.loc[Y_pred.index]
        corr =  calc_corr(Y_pred, Y_obs)
        acc = calc_Pairs_acc(Nr_pair_acc, Y_obs, Y_pred)
        mse = mean_squared_error(Y_obs, Y_pred)
        results = {'Element': 'train', acc_name: [acc],
                   corr_name: [corr], mse_name: [mse]}
        performances = pd.DataFrame(results).set_index('Element').pivot_table(index=None, columns='Element')
    
    return performances


def assess_models(sim_setting):
    models = sim_setting['models']
    path_Y_train = sim_setting['path_Y_train']
    path_Y_test = sim_setting['path_Y_test']
    
    base_dir = sim_setting['base_dir']
    Nr_pair_acc = sim_setting['Nr_pair_acc']
    
    
    Y_obs_all_intergenic = read_obs(path_Y_train)
    Y_obs_all_elems = read_obs(path_Y_test)
    
    acc_all = []
    corr_all = []
    mse_all = []
    
    for key in models:
        m = models[key]
        # load train
        save_name = m['save_name']
        print(save_name)
        
        Y_obs_unseen = Y_obs_all_elems.copy()
        Y_obs_seen = Y_obs_all_intergenic.copy()
            
        path_pred_unseen = f'{base_dir}/{save_name}/{save_name}_predTest.tsv'
        Y_pred_unseen = read_pred(path_pred_unseen)
        Y_pred_unseen = Y_pred_unseen.loc[Y_obs_unseen.index]

        if (Y_pred_unseen.index != Y_obs_unseen.index).all():
            raise ValueError('index mismatch')
        assessments_test = assess_model(Y_pred_unseen, Y_obs_unseen, 
                                        Nr_pair_acc, save_name, per_element=True)
        
        
        
        path_pred_seen = f'{base_dir}/{save_name}/{save_name}_predTrain.tsv'
        Y_pred_seen = read_pred(path_pred_seen)
        Y_pred_seen = Y_pred_seen.loc[Y_obs_seen.index]
        
        
        if (Y_pred_seen.index != Y_obs_seen.index).all():
            ValueError('index mismatch')
         
         
        assessments_train = assess_model(Y_pred_seen, Y_obs_seen, Nr_pair_acc,
                                         save_name, per_element=False)
        
        assessments = pd.concat([assessments_test, assessments_train], axis=1)
        
        assessments.to_csv(f'{base_dir}/{save_name}/{save_name}_assessments.tsv', sep='\t')
        
        acc_all.append(assessments.loc['acc_'+save_name])
        corr_all.append(assessments.loc['corr_'+save_name])
        mse_all.append(assessments.loc['mse_'+save_name])
        
        print("=========================")
    
    acc_all_df = pd.concat(acc_all, axis=1)
    corr_all_df = pd.concat(corr_all, axis=1)
    mse_all_df = pd.concat(mse_all, axis=1)
    
    acc_all_df.to_csv(f'{base_dir}/acc_all.tsv', sep='\t')
    corr_all_df.to_csv(f'{base_dir}/corr_all.tsv', sep='\t')
    mse_all_df.to_csv(f'{base_dir}/mse_all.tsv', sep='\t')
    

def assess_perBatch(dir_path):
    setting_config = 'sim_setting.ini'
    param_config = 'dev.ini'
        
    sim_setting = load_sim_settings_perBatchPerf(dir_path, setting_config,
                                                 param_config)
    save_name = list(sim_setting['models'].keys())[0]
    sim_params = sim_setting['models'][save_name]
    predict_func = sim_params['predict_func']
    base_dir = sim_setting['base_dir']
    directory_path = f'{base_dir + save_name}/models_interval/'
    split_intergenic = ast.literal_eval(sim_setting['split_intergenic'])
    # List all files in the directory
    model_names = [os.path.join(directory_path, file) for file in os.listdir(directory_path)]

    X_train, Y_train, X_test, Y_test = load_data(sim_setting)
    Y_obs_all_intergenic = read_obs(sim_setting['path_Y_train'])
    Y_obs_all_elems = read_obs(sim_setting['path_Y_test'])
    Y_obs_val = read_obs(sim_setting['path_Y_validate'])
    bins_val = Y_obs_val.index
    if split_intergenic:
        Y_obs_unseen = pd.concat([Y_obs_all_elems, Y_obs_val], axis= 0)
        Y_obs_seen = Y_obs_all_intergenic.drop(bins_val)
    else:
        Y_obs_unseen = Y_obs_all_elems
        Y_obs_seen = Y_obs_all_intergenic
    
    N = np.unique(Y_train.N)[0]
    batch_indexes = []
    for batch in range(sim_params['Args']['epochs'], 0, -1):
        if batch % (sim_params['Args']['save_interval']) == 0:
            btch_idx = f'batch_{batch}'
            batch_indexes.append(btch_idx)
    
    batch_indexes = [batch_idx for batch_idx in batch_indexes if any(batch_idx in model_name for model_name in model_names)]
    acc_file = f'{base_dir}/{save_name}/perBatch_assessments.tsv'
    # acc_file = f'{base_dir}/{save_name}/perElement_accuracies.tsv'
    # corr_file = f'{base_dir}/{save_name}/perElement_correlations.tsv'
    # MSE_file = f'{base_dir}/{save_name}/perElement_MSEs.tsv'
    if os.path.exists(acc_file):
        accuracies = pd.read_csv(acc_file, sep='\t', index_col=0)
        # Filter out batches that are already in the accuracies DataFrame
        existing_batches = set(accuracies.index.str.replace('acc_', ''))
        batch_indexes = [batch for batch in batch_indexes if batch not in existing_batches]
    else:
        accuracies = pd.DataFrame()
    
    # If no new batches to assess, return immediately
    if not batch_indexes:
        print("All batches already assessed.")
        return
    
    for batch_id in batch_indexes:
        
        model_name = os.path.join(directory_path, f'{batch_id}_model.h5')
        # Check if the model file exists before proceeding
        if not os.path.exists(model_name):
            print(f"Model file for {batch_id} not found. Skipping.")
            continue

        print(batch_id)
        # print(model_name)

        # load the model
        with h5py.File(model_name, 'r') as f:
            # load the model
            model = load_model(f)

        model_data = {'model': model, 
                      'N': N,
                      'NN_hyperparams': sim_params['Args']}
        
        # if classic_NN:
            # model_data['response_type'] = sim_params['Args']['response']
        model_data['response_type'] = 'rate'
        
        predRates_test = predict_func(model_data, X_test, "")
        predRates_train = predict_func(model_data, X_train, "")
        
        
        assessDF_test = assess_model(predRates_test, Y_obs_unseen,
                                     sim_setting['Nr_pair_acc'], 
                                     batch_id, per_element=True)
        
        
        assessDF_train = assess_model(predRates_train, Y_obs_seen,
                                     sim_setting['Nr_pair_acc'], 
                                     batch_id, per_element=False)
        
        
        

        assessDF = pd.concat([assessDF_test, assessDF_train], axis=1)
        accuracies = pd.concat([assessDF, accuracies], axis=0)
        
        accuracies.to_csv(f'{base_dir}/{save_name}/perBatch_assessments.tsv',
                          sep='\t')
        print('****************')
        
        
###################################################################3
def load_all_obsRates():
    path_Y_test = '../external/rawInput/Pan_Cancer_test_y.tsv'
    path_Y_train = '../external/rawInput/Pan_Cancer_train_y.tsv'
    
    Y_test = read_response(path_Y_test) 
    Y_train = read_response(path_Y_train) 
    Y_all = pd.concat([Y_test, Y_train], axis=0)
    return Y_all


def extract_Y_obs(path_preds):
    Yobs_all = load_all_obsRates()
    Y_pred = read_pred(path_preds)
    Y_obs = Yobs_all.loc[Y_pred.index]
    Y_obs = Y_obs[Y_obs.nMut != 0]
    Y_obs = pd.DataFrame(Y_obs.obsRates)
    return Y_obs
 

    


def assess_models_new(sim_setting):
    models = sim_setting['models']
    base_dir = sim_setting['base_dir']
    Nr_pair_acc = sim_setting['Nr_pair_acc']
    # model_name = list(models.keys())[0]
    # save_name = sim_setting['models'][model_name]['save_name']
    
    for key in models:
       
        m = models[key]
        save_name = m['save_name']   
        
        print(f'assessment for {save_name}')
        path_predTest = f'{base_dir}/{save_name}/{save_name}_predTest.tsv'

        Y_pred_test = read_pred(path_predTest)
        Y_obs_test = extract_Y_obs(path_predTest)

        path_predTrain = f'{base_dir}/{save_name}/{save_name}_predTrain.tsv'
        Y_pred_train = read_pred(path_predTrain)
        Y_obs_train = extract_Y_obs(path_predTrain)


        test_ensemble = assess_model(Y_pred_test, Y_obs_test, Nr_pair_acc = Nr_pair_acc, 
                     model_name = save_name, per_element = True)

        train_ensemble = assess_model(Y_pred_train, Y_obs_train, Nr_pair_acc = Nr_pair_acc, 
                     model_name = save_name, per_element = False)

        assessments = pd.concat([test_ensemble, train_ensemble], axis=1)

        assessments.to_csv(f'{base_dir}/{save_name}/{save_name}_assessments.tsv', sep='\t')
    
    


def assess_model_number_n(sim_setting, base_dir, save_name):
    Nr_pair_acc = sim_setting['Nr_pair_acc']
    
    path_predTest = f'{base_dir}/{save_name}/{save_name}_predTest.tsv'

    Y_pred_test = read_pred(path_predTest)
    Y_obs_test = extract_Y_obs(path_predTest)

    path_predTrain = f'{base_dir}/{save_name}/{save_name}_predTrain.tsv'
    Y_pred_train = read_pred(path_predTrain)
    Y_obs_train = extract_Y_obs(path_predTrain)


    test_ensemble = assess_model(Y_pred_test, Y_obs_test, Nr_pair_acc = Nr_pair_acc, 
                 model_name = save_name, per_element = True)

    train_ensemble = assess_model(Y_pred_train, Y_obs_train, Nr_pair_acc = Nr_pair_acc, 
                 model_name = save_name, per_element = False)

    assessments = pd.concat([test_ensemble, train_ensemble], axis=1)

    assessments.to_csv(f'{base_dir}/{save_name}/{save_name}_assessments.tsv', sep='\t')
