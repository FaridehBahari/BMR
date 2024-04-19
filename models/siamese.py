# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 20:12:20 2023

@author: Farideh
"""
import time
import configparser
import os
import pandas as pd
import numpy as np
import pickle
import h5py
from readFtrs_Rspns import split_by_element, save_preds_tsv, set_gpu_memory_limit, load_regulatory_elems
from tensorflow.keras.layers import Dense, BatchNormalization, Concatenate, LeakyReLU, PReLU
from tensorflow.keras import layers, models
from tensorflow.keras.models import Sequential, load_model
# from simulation_settings import load_sim_settings
# from tens import keras

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def build_rank_params(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    
    # Load the hyperparameters from the config file
    param = {
        'activation1': config.get('architecture', 'activation1'),
        'activation_out': config.get('architecture', 'activation_out'),
        'dropout': config.getfloat('architecture', 'dropout'),
        'dropout_cd': config.getfloat('architecture', 'dropout_cd'),
        'batch_norm_cd': config.getboolean('architecture', 'batch_norm_cd'),
        'batch_norm': config.getboolean('architecture', 'batch_norm'),
        'P_unique': config.getfloat('training', 'P_unique'),
        'n_inputs': config.getint('main', 'n_inputs'),
        'concat_dense': [int(x.strip()) for x in config.get('architecture','concat_dense').split(',')],
        'save_interval': config.getint('training', 'save_interval'),
        'loss': config.get('training', 'loss'),
        'optimizer': config.get('training', 'optimizer'),
        'metrics': [x.strip() for x in config.get('training', 'metrics').split(',')],
        'epochs': config.getint('training', 'epochs'),
        'batch_size': config.getint('training', 'batch_size'),
        'pred_type': config.get('prediction', 'pred_type'),
        'architecture': [int(x.strip()) for x in config.get('architecture', 'architecture').split(',')]
    }
    return param


# def build_model_architecture_siamese2(hyperparams, n_ftrs):
#     n = hyperparams['n_inputs']
#     n_ftrs = (n_ftrs,)
    
#     if hyperparams['loss'] == 'mse':
#         nr_outputs = 2**(n-1)

#     elif hyperparams['loss'] == 'poisson':
#         nr_outputs = 2**(n)
    
    
#     shared_network = models.Sequential()
#     for i, units in enumerate(hyperparams['architecture']):
#         shared_network.add(layers.Dense(units, input_shape=n_ftrs if i == 0 else ()))
        
#         if hyperparams['batch_norm']:
#             shared_network.add(BatchNormalization())
        
#         if hyperparams['activation1'] == 'leaky_relu':
#             shared_network.add(LeakyReLU(alpha=0.2))
#         elif hyperparams['activation1'] == 'PReLU':
#             shared_network.add(PReLU())
#         else:
#             shared_network.add(layers.Activation(hyperparams['activation1']))
        
#         if hyperparams['dropout'] != 0:
#             shared_network.add(layers.Dropout(hyperparams['dropout']))
    
#     input_layer = layers.Input(shape=(n*n_ftrs[0],))
    
#     # Define Lambda layers to split the inputs
#     inputs = []
#     for i in range(n):
#         inputs.append(layers.Lambda(lambda x: x[:, i*n_ftrs[0]:(i+1)*n_ftrs[0]])(input_layer))
    
#     # Connect the inputs to the shared sub-network (twin)
#     outputs = []
#     for i in range(n):
#         outputs.append(shared_network(inputs[i]))
        
#     if n != 1:
#         concat = layers.Concatenate()(outputs)
#     else:
#         concat = outputs
    
#     if hyperparams['concat_dense']== [0]:
#         output = Dense(nr_outputs, activation=hyperparams['activation_out'])(concat)
#     else:
#         for i, units in enumerate(hyperparams['concat_dense']):
#             if i == 0:
#                 if units != [0]:
                    
#                     dense = Dense(units)(concat)
#                     if hyperparams['batch_norm_cd']:
#                           dense = BatchNormalization()(dense)
                          
#                     if hyperparams['activation1'] == 'leaky_relu':
#                         dense = LeakyReLU(alpha=0.2)(dense)
                    
#                     elif hyperparams['activation1'] == 'PReLU':
#                         dense = PReLU()(dense)
#                     else:
#                         dense = layers.Activation(hyperparams['activation1'])(dense) 
                    
#             else:
#                 if hyperparams['batch_norm']:
#                     dense = BatchNormalization()(dense)
#                 dense = Dense(units)(dense)
                
                
#                 if hyperparams['activation1'] == 'leaky_relu' :
#                     dense = LeakyReLU(alpha=0.2)(dense)
#                 elif hyperparams['activation1'] == 'PReLU' :
#                     dense = PReLU()(dense)
#                 else:
#                    dense = layers.Activation(hyperparams['activation1'])(dense) 
#             if hyperparams['dropout_cd'] != 0:
#                 dense = layers.Dropout(hyperparams['dropout_cd'])(dense) 
            
#         output = Dense(nr_outputs, activation=hyperparams['activation_out'])(dense)
    
#     model = models.Model(inputs=input_layer, outputs=output)
    
#     model.compile(loss=hyperparams['loss'], optimizer=hyperparams['optimizer'], 
#                   metrics=hyperparams['metrics'])
#     return model


def build_model_architecture_siamese2(hyperparams, n_ftrs):
    n = hyperparams['n_inputs']
    n_ftrs = (n_ftrs,)
    
    if hyperparams['loss'] == 'mse':
        nr_outputs = 2**(n-1)
        
    elif hyperparams['loss'] == 'poisson':
        nr_outputs = 2**(n)
    
    # Define the shared sub-network (twin)
    shared_network = models.Sequential()
    for i, units in enumerate(hyperparams['architecture']):
        if i == 0:
            shared_network.add(layers.Dense(units, input_shape = n_ftrs))
        else:
            shared_network.add(layers.Dense(units))
            
        if hyperparams['batch_norm']:
            shared_network.add(BatchNormalization())
            
        if hyperparams['activation1'] == 'leaky_relu' :
            shared_network.add(LeakyReLU(alpha=0.2))
            
        elif hyperparams['activation1'] == 'PReLU' :
            shared_network.add(PReLU())
            
        else:
            shared_network.add(layers.Activation(hyperparams['activation1']))    
        
        if hyperparams['dropout'] != 0:
            shared_network.add(layers.Dropout(hyperparams['dropout']))
    
    input_layer = layers.Input(shape=(n*n_ftrs[0],))
    
    # Define Lambda layers to split the inputs
    inputs = []
    for i in range(n):
        inputs.append(layers.Lambda(lambda x: x[:, i*n_ftrs[0]:(i+1)*n_ftrs[0]])(input_layer))
    
    # Connect the inputs to the shared sub-network (twin)
    outputs = []
    for i in range(n):
        outputs.append(shared_network(inputs[i]))
        
    if n != 1:
        concat = layers.Concatenate()(outputs)
    else:
        concat = outputs
    
    for i, units in enumerate(hyperparams['concat_dense']):
            if i == 0:
                if units != [0]:
                    
                    dense = Dense(units)(concat)
                    if hyperparams['batch_norm_cd']:
                          dense = BatchNormalization()(dense)
                          
                    if hyperparams['activation1'] == 'leaky_relu':
                        dense = LeakyReLU(alpha=0.2)(dense)
                    
                    elif hyperparams['activation1'] == 'PReLU':
                        dense = PReLU()(dense)
                    else:
                        dense = layers.Activation(hyperparams['activation1'])(dense) 
                    
            else:
                if hyperparams['batch_norm']:
                    dense = BatchNormalization()(dense)
                dense = Dense(units)(dense)
                
                
                if hyperparams['activation1'] == 'leaky_relu' :
                    dense = LeakyReLU(alpha=0.2)(dense)
                elif hyperparams['activation1'] == 'PReLU' :
                    dense = PReLU()(dense)
                else:
                    dense = layers.Activation(hyperparams['activation1'])(dense) 
            if hyperparams['dropout_cd'] != 0:
                dense = layers.Dropout(hyperparams['dropout_cd'])(dense) 
            
    output = Dense(nr_outputs, activation=hyperparams['activation_out'])(dense)
    
    model = models.Model(inputs=input_layer, outputs=output)
    
    model.compile(loss=hyperparams['loss'], optimizer=hyperparams['optimizer'], 
                  metrics=hyperparams['metrics'])
    return model


def get_distinct_indices(n, batch_size, n_input):
    indices = np.random.choice(n, size=(n_input, batch_size), replace=True)
    
    return indices.T

def generate_outputs(vector):
    n = len(vector)  # Number of columns
    rows = 2**(n-1)  # Number of rows
    sign_matrix = np.zeros((rows, n), dtype=int)
    
    for i in range(rows):
        binary = bin(i)[2:].zfill( n )
        for j in range ( n ):
            sign_matrix[i, j] = int(binary[j]) * 2 - 1
    
    sign_matrix = -sign_matrix
    result = np.dot(sign_matrix, vector)
    
    return result

def generate_outputs_pois(vector):
    n = len(vector)  # Number of columns
    rows = 2**(n)  # Number of rows
    ZO_matrix = np.zeros((rows, n), dtype=int)
    
    for i in range(rows):
        binary = bin(i)[2:].zfill( n )
        for j in range ( n ):
            ZO_matrix[i, j] = int(binary[j])
    result = np.dot(ZO_matrix, vector)
    return result
            

def generate_siam_batch_samples(X_train, Y_train, batch_size, n_input, P_unique,
                                loss_type): # loss_type can be mse or poisson 
    
    indices = get_distinct_indices(X_train.shape[0], batch_size, n_input)
        
    if P_unique != 0.0:
       # Get unique values from the initial indices
       unique_values = np.unique(indices) 
       
       # select P_unique% of the unique_values to concat to the initial indices
       unique_values_subset = np.random.choice(unique_values,
                        int(P_unique * batch_size),
                        replace = False)
       # Create new rows filled with each unique value
       new_rows = np.array([[value] * indices.shape[1] for value in unique_values_subset]) 
       # Concatenate the new rows to the original array
       indices = np.concatenate((indices, new_rows))
    
    X_batch = X_train.iloc[indices.flatten()].values.reshape(-1, n_input * X_train.shape[1])
    
    if loss_type == 'mse':
        generate_out_func = generate_outputs
        Y_train['obsRate'] = np.log((Y_train['nMut']) / (Y_train['length']))
        Y = Y_train['obsRate'].iloc[indices.flatten()].values.reshape(-1, n_input)
    elif loss_type == 'poisson':
        generate_out_func = generate_outputs_pois
        Y_tr = (Y_train['nMut'])/(Y_train['length'])
        Y = Y_tr.iloc[indices.flatten()].values.reshape(-1, n_input)
    
    Y_batch = []    
    for Y_row in Y:
        outputs_row = generate_out_func(Y_row)
        Y_batch.append(outputs_row)
    
    Y_batch = np.array(Y_batch)
    
    return X_batch, Y_batch


def iterate_training(init_batch, model,
                     X_train, Y_train, NN_hyperparams):
    save_path_interval = NN_hyperparams['path_save']
    for batch in range(init_batch, NN_hyperparams['epochs']+1):
        
        
        # prepare training data
        X_batch, Y_batch = generate_siam_batch_samples(X_train, Y_train,
                                                       NN_hyperparams['batch_size'],
                                                       NN_hyperparams['n_inputs'],
                                                       NN_hyperparams['P_unique'],
                                                       NN_hyperparams['loss'])
        
        # Train on batch
        metrics = model.train_on_batch(X_batch, Y_batch)
        
        if (batch) % NN_hyperparams['save_interval'] == 0:
            print(f'====== batch: {batch} =======')
            print(f'batch_size = {X_batch.shape[0]}')
            lastBatch_npy = f'{save_path_interval}lastBatch.npy'
            np.save(lastBatch_npy, batch)
            model_save_name = f'{save_path_interval}batch_{batch}_model.h5'
            model.save(model_save_name)
    return metrics
            

def calculate_r1_pois(outputs):
    rows = outputs.shape[0]
    B = outputs[range(int(rows/2))]
    A = outputs[range(int(rows/2), rows)]
    r1 = A - B
    
    r1 = np.mean(r1)
    return r1


def run_rankNN_iteration2(X_train, Y_train, NN_hyperparams):
    
    n_ftrs = int(X_train.shape[1])
    save_path_interval = NN_hyperparams['path_save']
    os.makedirs(save_path_interval, exist_ok= True)
    
    
    
    #print("***0*******")
    lastBatch_npy = f'{save_path_interval}lastBatch.npy'
    
    if not os.path.exists(lastBatch_npy):
        # Bulid initial model
        model = build_model_architecture_siamese2(NN_hyperparams, n_ftrs)
        print(model.summary())
        init_batch = 1
    else:
        last_saved_batch = np.load(lastBatch_npy)
        init_batch = last_saved_batch + 1
        last_saved_model = f'{save_path_interval}batch_{last_saved_batch}_model.h5'
        n_batch = NN_hyperparams['epochs'] 
        desired_model = f'{save_path_interval}batch_{n_batch}_model.h5'
        if n_batch < int(last_saved_model):
            
            # load desired model
            with h5py.File(desired_model, 'r') as f:
                # load the model
                model = load_model(f)
        else:
            # load the last model
            with h5py.File(last_saved_model, 'r') as f:
            # load the model
                 model = load_model(f)
            
    if init_batch < NN_hyperparams['epochs']+1 or init_batch == NN_hyperparams['epochs']:
        
        metrics = iterate_training(init_batch, model,
                             X_train, Y_train, NN_hyperparams)
        metrics_used  = str(NN_hyperparams['metrics'])
        print(f'training loss and training {metrics_used}: {metrics}')
    else:
        print('all batches were completed...')   
    
    model_data = {'model': model,
                  'NN_hyperparams': NN_hyperparams,
                  'N': Y_train.N[0]
                  }
    # use the model for feature extraction
    # sim_file = 'configs/rate_based/fb_sim_setting_elemSp.ini' 
    
    # sim_setting = load_sim_settings(sim_file)
    # X_regLmnt, Y_regLmnt = load_regulatory_elems(sim_setting)
    
    # pretrained_model = model
    # layer_name = 'sequential'
    # encoder_model = keras.Model(inputs=pretrained_model.get_layer(layer_name).input,
    #                             outputs=pretrained_model.get_layer(layer_name).output)
    
    # # Use the encoder model to obtain the encoded data
    # encoded_data = encoder_model.predict(X_regLmnt)
    
    # encoded_data = pd.DataFrame(encoded_data)
    # encoded_data.index = X_regLmnt.index
    
    # encoded_data.to_csv(f'{save_path_interval}X_regLmnt_extracted.tsv')
    # print(encoded_data.shape)
    
    
    return model_data

# from sklearn.model_selection import train_test_split
# def run_rankNN_iteration2(X_train, Y_train, NN_hyperparams):
    
#     n_ftrs = int(X_train.shape[1])
#     save_path_interval = NN_hyperparams['path_save']
#     os.makedirs(save_path_interval, exist_ok= True)
#     # Bulid initial model
#     model = build_model_architecture_siamese2(NN_hyperparams, n_ftrs)
#     print(model.summary())
#     #print("***0*******")
    
#     # Split the training data into a training and validation set
#     X_train, X_val, Y_train, Y_val = train_test_split(X_train, Y_train,
#                                                       test_size=0.1)
    
#     best_loss = float('inf')  # Start with a high error
#     patience_count = 0  # Initial patience count
#     patience = 5  # Patience limit
    
#     # prepare validation data
#     X_val, Y_val = generate_siam_batch_samples(X_val, Y_val,
#                                                    100000,
#                                                    NN_hyperparams['n_inputs'],
#                                                    NN_hyperparams['P_unique'])
    
#     for batch in range(1, NN_hyperparams['epochs']+1):
#         # prepare training data
#         X_batch, Y_batch = generate_siam_batch_samples(X_train, Y_train,
#                                                               NN_hyperparams['batch_size'],
#                                                               NN_hyperparams['n_inputs'],
#                                                               NN_hyperparams['P_unique'])
        
#         # Train on batch
#         metrics = model.train_on_batch(X_batch, Y_batch)
        
#         # Calculate loss on validation set
#         validation_loss, _ = model.evaluate(X_val, Y_val, verbose=0)  # assuming you have a single loss function
#         if validation_loss > best_loss:
#            patience_count += 1
#            print(f'val_loss has not improved from {best_loss:.5f}. Patience counter: {patience_count}/{patience}')
           
#            if patience_count >= patience:
#                print('Early stopping due to validation loss has not improved.')
#                break
               
#         else:
#            best_loss = validation_loss
#            patience_count = 0
#            print(f'val_loss improved to {best_loss:.5f}. Resetting patience counter.')
        
#         if (batch) % NN_hyperparams['save_interval'] == 0:
#             print(f'====== batch: {batch} =======')
#             print(f'batch_size = {X_batch.shape[0]}')
            
#             model_save_name = f'{save_path_interval}batch_{batch}_model.h5'
#             model.save(model_save_name)
            
        
#     metrics_used  = str(NN_hyperparams['metrics'])
#     print(f'training loss and training {metrics_used}: {metrics}')
    
#     model_data = {'model': model,
#                   'NN_hyperparams': NN_hyperparams,
#                   'N': Y_train.N[0]
#                   }
#     return model_data

def calculate_r1(outputs):
    r1 = np.mean(outputs)
    return r1

def generate_siam_testDat(X_test, n_input, pred_type):
    
    if pred_type == 'random_pairs':
        # Create a column to be added at the beginning of indices
        fixed_column = np.arange(X_test.shape[0])
        
        # generate indices array for next n_input-1 elements
        indices = get_distinct_indices(X_test.shape[0], 
                                       X_test.shape[0], n_input-1)
        
        # Stack the fixed first column with the indices array
        indices = np.hstack((fixed_column.reshape(-1, 1), indices))
        test_data = X_test.iloc[indices.flatten()].values.reshape(-1, n_input * X_test.shape[1])
        
    elif pred_type == 'given_elements':
        all_indices = np.arange(X_test.shape[0])
        indices = np.tile(np.array([all_indices]).transpose(), (n_input))
        test_data = X_test.iloc[indices.flatten()].values.reshape(-1, n_input * X_test.shape[1])
        
    return test_data

# Function to check if elements exist in the feature matrix index
def check_pattern_in_index(df, patterns):
    for pattern in patterns:
        if any(pattern in index for index in df.index):
            return True
    return False


def predict_rankNN(model, feature_matrix, elem_length):
    N = model['N']
    NN_hyperparams = model['NN_hyperparams']
    pred_type = NN_hyperparams['pred_type']
    
    if NN_hyperparams['loss'] == 'mse':
        calc_r_func = calculate_r1
    elif NN_hyperparams['loss'] == 'poisson':
        calc_r_func = calculate_r1_pois    
    
    M = model['model'] 
    elems = ["gc19_pc.cds", "enhancers", "gc19_pc.3utr", "gc19_pc.5utr",
              "gc19_pc.promCore", "gc19_pc.ss", "lncrna.ncrna", 
              "lncrna.promCore"]
    if sum(feature_matrix.index.str.contains(r'[v]\d')) != 0:
        elems.append("intergenic")
    
    # Check if it is X_train or X_test
    # use_X_test = check_pattern_in_index(feature_matrix, elems)
    use_X_test = (sum(feature_matrix.index.str.contains(r'[v]\d')) != feature_matrix.shape[0])
    if use_X_test:
        predictions = pd.DataFrame()
        for elem in elems:
            print(elem)
            # split the X_test to element-types
            if elem == "intergenic":
                test_elem = split_by_element(feature_matrix, r'[v]\d')
            else:
                test_elem = split_by_element(feature_matrix, elem)
            
            
            # generate test_data 
            test_data = generate_siam_testDat(test_elem,
                                              NN_hyperparams['n_inputs'],
                                              pred_type)
            
            # make prediction on the pairs of an element-type
            prediction = M.predict(test_data, verbose = 0)
            # report r1
            pred_elem = []    
            for row in prediction:
                pred_row = calc_r_func(row)
                pred_elem.append(pred_row)
            
            if NN_hyperparams['loss'] == 'mse':
                pred_elem = np.exp(pred_elem)/N 
                
            elif NN_hyperparams['loss'] == 'poisson':
                pred_elem = pred_elem/N
            
            predElem_df = pd.DataFrame({'predRate': pred_elem}, 
                                         index=test_elem.index)
            predictions = pd.concat([predictions, predElem_df],
                                              axis = 0)
    else:
       st_time = time.time()
       train_data = generate_siam_testDat(feature_matrix, 
                                          NN_hyperparams['n_inputs'],
                                          pred_type)
       end_t = time.time()
       print(f'generatate_siam_testDat took {end_t - st_time} seconds')
       
       print("============")
       st_time3 = time.time()
       prediction = M.predict(train_data, verbose = 0)
       end_t3 = time.time()
       print(f'M.predict took {end_t3 - st_time3} seconds')
       
       
       st_time0 = time.time()
       
       # report r1
       pred = []    
       for row in prediction:
           pred_row = calc_r_func(row)
           pred.append(pred_row)
           
       if NN_hyperparams['loss'] == 'mse':
           pred = np.exp(pred)/N 
           
       elif NN_hyperparams['loss'] == 'poisson':
           pred = pred/N
           
       predictions = pd.DataFrame({'predRate': pred}, 
                                    index=feature_matrix.index) 
       end_t0 = time.time()
       print(f'report r1 took {end_t0 - st_time0} seconds')
       
    return predictions

def save_nn(fitted_Model, path_save, save_name, iteration = None, save_model = True): 
    
    save_preds_tsv(fitted_Model, path_save, save_name, iteration)
    
    if save_model:
        M = fitted_Model.model['model']
        # Save the model 
        if iteration is not None:
            save_path_model = f'{path_save}/{save_name}/rep_train_test/{save_name}_model_{iteration}.h5'
        else:
            save_path_model = f'{path_save}/{save_name}/{save_name}_model.h5'
            
        M.save(save_path_model)
        
def save_siamese(fitted_Model, path_save, save_name, iteration = None, save_model = True): 
    
    save_preds_tsv(fitted_Model, path_save, save_name, iteration)
    
    if save_model:
        model = fitted_Model.model['model']
        params = fitted_Model.model['NN_hyperparams']
        
        # Save the model 
        if iteration is not None:
            save_path_model = f'{path_save}/{save_name}/rep_train_test/{save_name}_model_{iteration}.h5'
            save_path_param = f'{path_save}/{save_name}/rep_train_test/{save_name}_params_{iteration}.pkl'
        else:
            save_path_model = f'{path_save}/{save_name}/{save_name}_model.h5'
            save_path_param = f'{path_save}/{save_name}/{save_name}_params.pkl'
            
        model.save(save_path_model)
            
    # Save the dictionary to a file using pickle
        with open(save_path_param, 'wb') as f: 
            pickle.dump(params, f)
        
    
    

def pair_rank_info_siamese2(save_name, *args):
    params = build_rank_params(args[0])
    model_dict = {"save_name" : save_name,
                  "Args" : params,
                  "run_func": run_rankNN_iteration2,
                  "predict_func": predict_rankNN,
                  "save_func": save_siamese
                  # "check_file_func": check_file_nn
                  }
    return model_dict

