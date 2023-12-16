import pandas as pd
import os

from readFtrs_Rspns import read_response
from performance.assessModels import assess_model


def load_all_obsRates():
    path_Y_test = '../external/rawInput/Pan_Cancer_test_y.tsv'
    
    path_Y_train = '../external/rawInput/Pan_Cancer_train_y.tsv'
    
    Y_test = read_response(path_Y_test) 
    Y_train = read_response(path_Y_train) 
    Y_all = pd.concat([Y_test, Y_train], axis=0)
    return Y_all

def read_pred(path_pred):
    Y_pred = pd.read_csv(path_pred, sep = "\t", header=0, index_col='binID',
                         usecols=['binID', 'predRate'])
    return Y_pred

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
    save_name = list(models.keys())[0]

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

