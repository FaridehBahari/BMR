# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 14:02:50 2023

@author: Farideh
"""
import configparser
import pandas as pd
import os
import numpy as np
import pickle
from readFtrs_Rspns import read_fi
from sklearn.ensemble import RandomForestRegressor
from readFtrs_Rspns import save_preds_tsv

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def build_RF_params(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    
    max_depth_str = config['train']['max_depth'] 
    max_depth = int(max_depth_str) if max_depth_str != 'None' else None
    # Load the hyperparameters from the config file
    params = {
        'n_estimators': config.getint('train', 'n_estimators'),
        'max_depth': max_depth,
        'n_threads': config.getint('train', 'n_threads')
    }
    return params


def run_rf(X_train, Y_train, RF_params):
    
    # use_features = read_fi('../external/rawInput/DP_results/Rndm_Frst.feature_importance.tsv',
    #                        cutoff=0.5)
    # X_train = X_train[use_features]

    rf = RandomForestRegressor(n_estimators=RF_params['n_estimators'], 
                               max_depth=RF_params['max_depth'],
                               verbose= 2,
                               n_jobs= RF_params['n_threads']) # n_estimators: number of trees
    rf.fit(X_train, Y_train.obsRates)
    feature_names = X_train.columns
    
    model_data = {'model' : rf,
                  'cols' : feature_names}
    
    return model_data

def predict_rf(model, X_test, length_elems):
    
    X_test = X_test[model['cols']]
    binID = X_test.index
    pred_test = model['model'].predict(X_test) 
    prediction_df = pd.DataFrame({'predRate': pred_test.ravel()}, 
                                 index=binID) 
    return prediction_df


def RF_model_info(save_name, *args):
    params = build_RF_params(args[0])
    model_dict = {"save_name" : save_name,
                  "Args" : params,
                  "run_func": run_rf,
                  "predict_func": predict_rf,
                  "save_func": save_rf,
                  "check_file_func": check_file_rf
                  }
    
    return model_dict


def save_rf(fitted_Model, path_save, save_name, save_model = True): 
    
    save_preds_tsv(fitted_Model, path_save, save_name)
    
    if save_model:
        M = fitted_Model.model['model']
        # Save the model using pickle
        save_path_model = f'{path_save}/{save_name}/{save_name}_model.pkl'
        with open(save_path_model, 'wb') as f: 
            pickle.dump(M, f)
    

def check_file_rf(path_save, save_name):
    file_name = f'{path_save}/{save_name}/{save_name}_model.pkl'
    return file_name