# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:33:07 2023

@author: Farideh
"""
import pandas as pd
from readFtrs_Rspns import save_preds_tsv
import configparser

def build_DP_params(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    # Load the hyperparameters from the config file
    params = {
        'path_DP_fitted': config.get('train', 'path_DP_fitted'),
        'path_DP_predicted': config.get('train', 'path_DP_predicted')
    }
    return params


def run_DP(X_train, Y_train, params):
    res_test = pd.read_csv(params['path_DP_predicted'], sep = '\t', header=0, 
                           index_col='binID',
                    usecols=['binID', 'length', 'nMut', 'nSample', 'N', 'nPred'])
    res_train = pd.read_csv(params['path_DP_fitted'], sep = '\t', header=0, 
                           index_col='binID',
                    usecols=['binID', 'length', 'nMut', 'nSample', 'N', 'nPred'])
    results = pd.concat([res_test, res_train], axis= 0)
    return results


def predict_DP(model, X, length_elems):
    
    binID = X.index
    pred_DP = model['nPred']/(model['N'] * model['length'])
    pred_DP = pred_DP.loc[binID]
    prediction_df = pd.DataFrame({'predRate': pred_DP.ravel()}, 
                                 index=binID) 
    return prediction_df

def DP_model_info(save_name, *args):
    
    params = build_DP_params(args[0])
    model_dict = {"save_name" : save_name,
                  "Args" : params,
                  "run_func": run_DP,
                  "predict_func": predict_DP,
                  "save_func": save_DP,
                  "check_file_func": check_file_DP
                  }
    
    return model_dict


def save_DP(fitted_Model, path_save, save_name, save_model = False): 
    
    save_preds_tsv(fitted_Model, path_save, save_name)
    
    if save_model:
        print("DriverPower was used to train and predict BMR...")


def check_file_DP(path_save, save_name):
    file_name = f'{path_save}/{save_name}/{save_name}_predTrain.tsv'
    return file_name