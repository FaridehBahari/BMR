import pandas as pd
import os
import glob
import h5py
import re

from readFtrs_Rspns import read_response
from performance.assessModels import assess_model
from tensorflow.keras.models import  load_model
from models.runBMR_functions import load_data_new
from simulation_settings import load_sim_settings
from models.siamese_new import predict_rankNN, save_siamese
from collections import namedtuple


def get_files_not_predicted(first_directory, second_directory):
    # first_dir is the model dir and the second one is the pred dir.
    # List the files in the first and second directories
    
    if os.path.exists(second_directory):
        first_files = os.listdir(first_directory)
        second_files = os.listdir(second_directory)
        
        # Extract the base file names without the suffix from the second directory
        second_base_names = set(file.split('_modelNumber')[0] for file in second_files)
        
        # Find the file names in the first directory that don't have a corresponding base name in the second directory
        files_not_in_second_directory = [file for file in first_files if file.split('_modelNumber')[0] not in second_base_names]
        
        # Create the full file paths for the files not in the second directory
        full_paths = [os.path.join(first_directory, file_name) for file_name in files_not_in_second_directory]
        
    else:
        full_paths = [os.path.join(first_directory, d) for d in os.listdir(first_directory) if os.path.isdir(os.path.join(first_directory, d))] 
    
    print(f'There are {len(full_paths)} files for prediction...')
    
    return full_paths


def save_preds_batch_n(dir_path, n_batches, sim_file, path_save):
    
    saved_model = f'batch_{n_batches}_model.h5'
    
    model_path = os.path.join(dir_path, 'models_interval', saved_model)
    # Load the model 
    if os.path.isfile(model_path):
        with h5py.File(model_path, 'r') as f:
            model = load_model(f)
            
        # load test and train data
        sim_setting = load_sim_settings(sim_file)
        model_name = list((sim_setting['models']).keys())[0]
        sim_setting['models'][model_name]['save_name'] = dir_path.split('/')[-1]
        n_sample = None
        model_number = int((dir_path.split('_')[-1]).split('/')[0]) 
        X_train, Y_train, X_test, Y_test = load_data_new(sim_setting, n_sample, model_number)
        
        model_data = {'model': model,
                      'N': Y_train.N[0],
                      'NN_hyperparams': sim_setting['models'][model_name]['Args']
            }
        
        length_test_elems = Y_test['length']
        length_train_elems = Y_train['length']
        predRates_test = predict_rankNN(model_data, X_test, length_test_elems)
        predRates_train = predict_rankNN(model_data, X_train, length_train_elems)
        
        DataTuple = namedtuple('DataTuple', ['model', 'predRates_train', 
                                             'predRates_test'])
        
        fitted_Model = DataTuple(model=model, predRates_train=predRates_train,
                         predRates_test = predRates_test)
        
        save_name = sim_setting['models'][model_name]['save_name']
        print(f'@@@@@@@@**** {save_name} ****@@@@@@@@@')
        save_name = f'{save_name}_modelNumber{n_batches}'
        save_siamese(fitted_Model, path_save, save_name, save_model = False)


def load_all_obsRates():
    path_Y_test = '../external/rawInput/Pan_Cancer_test_y.tsv'
    
    path_Y_train = '../external/rawInput/Pan_Cancer_train_y.tsv'
    
    Y_test = read_response(path_Y_test) 
    Y_train = read_response(path_Y_train) 
    Y_all = pd.concat([Y_test, Y_train], axis=0)
    return Y_all


def compute_average_prediction(file_paths, path_save, save_name, pred_on):
    os.makedirs(f'{path_save+save_name}/', exist_ok= True)
    all_data = pd.DataFrame()
    
    for file_path in file_paths:
        df = pd.read_csv(file_path, delimiter='\t', index_col='binID')
        
        all_data = pd.concat([all_data, df], axis = 1)
    
    # Compute the average of all prediction rates
    all_data = pd.DataFrame(all_data.mean(axis=1))
    all_data.columns = ['predRate']
    # Save the result to a new file
    if pred_on == 'test':
        all_data.to_csv(f'{path_save+save_name}/{save_name}_ensemble_predTest.tsv', sep = '\t')
    elif pred_on == 'train':
        all_data.to_csv(f'{path_save+save_name}/{save_name}_ensemble_predTrain.tsv', sep = '\t')
    else:
        ValueError('pred_on should be test or train')
        

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


def get_all_Preds_dir(base_dir, file_pattern):
    # Initialize an empty list to store the file paths
    file_paths = []
    # Walk through the root directory and its subdirectories
    for root, dirs, files in os.walk(base_dir):
        for file in glob.glob(os.path.join(root, file_pattern)):
            file_paths.append(file)
            
    return file_paths
      

def save_ensemble_assessments(base_dir, path_save, save_name, n_models):
    print(f'====={save_name}=====')
    file_paths_test = get_all_Preds_dir(base_dir, '*_predTest.tsv')
    file_paths_train =  get_all_Preds_dir(base_dir, '*_predTrain.tsv')
    
    if len(file_paths_test) != len(file_paths_train):
        raise ValueError('the length of test preds is not compatible with train preds')
    
    print(f'there are {len(file_paths_train)} prediction files')
    
    if n_models > len(file_paths_train):
        raise ValueError('the number of provided models exceed the existing prediction files')
    
    else:
        file_paths_test = file_paths_test[0:n_models]        
        file_paths_train = file_paths_train[0:n_models]
        
        compute_average_prediction(file_paths_test, path_save, save_name, pred_on = 'test')
        compute_average_prediction(file_paths_train, path_save, save_name, pred_on = 'train')
        path_test_ensemble = f'{path_save+save_name}/{save_name}_ensemble_predTest.tsv'
        path_train_ensemble = f'{path_save+save_name}/{save_name}_ensemble_predTrain.tsv'
        Y_pred_test = read_pred(path_test_ensemble)
        Y_obs_test = extract_Y_obs(path_test_ensemble)
        Y_pred_train = read_pred(path_train_ensemble)
        Y_obs_train = extract_Y_obs(path_train_ensemble)
        test_ensemble = assess_model(Y_pred_test, Y_obs_test, Nr_pair_acc = 1000000, 
                     model_name = 'ensemble', per_element = True)
        
        train_ensemble = assess_model(Y_pred_train, Y_obs_train, Nr_pair_acc = 1000000, 
                     model_name = 'ensemble', per_element = False)
        
        assessments = pd.concat([test_ensemble, train_ensemble], axis=1)
        
        assessments.to_csv(f'{path_save}/{save_name}/{save_name}_ensemble_assessments.tsv', sep='\t')
        print('=====================================================')


