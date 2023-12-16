import numpy as np
import configparser
import pandas as pd
import os
import h5py
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Activation, Dense, Dropout, LeakyReLU, BatchNormalization, PReLU, Input, Add
from tensorflow.keras.models import Sequential, load_model, Model
from tensorflow.keras.losses import Poisson
from keras import regularizers
from tensorflow.keras.models import  load_model
from models.runBMR_functions import  load_data, config_save
from simulation_settings import load_sim_settings
from readFtrs_Rspns import split_by_element
from scipy.stats import spearmanr
from performance.assessModels import assess_model


# Generate bootstrap samples
def bootstrap_samples(X_elem, Y_elem):
    n_samples = Y_elem.shape[0]
    indices = np.random.choice(X_elem.index.unique(), size=n_samples, replace=True)
    X_samples = X_elem.loc[indices]
    y_samples = Y_elem.loc[indices]
    return X_samples, y_samples, np.unique(indices)


def generate_test_train_bootstrapSamples(X_elem, Y_elem):
    
    X_samples_train, y_samples_train, seen_bins = bootstrap_samples(X_elem, Y_elem)

    unseen_bins = Y_elem.index[np.where(~(Y_elem.index).isin(np.unique(seen_bins)))]
    y_samples_test = Y_elem.loc[np.unique(unseen_bins)]
    X_samples_test = X_elem.loc[np.unique(unseen_bins)]


    # check the number of bins in test and train data
    if len(np.unique(y_samples_train.index)) + len(np.unique(y_samples_test.index)) != Y_elem.shape[0]:
        raise ValueError('number of sample test_train indices is incompatible with the input data indices ')
        
    return X_samples_train, y_samples_train, X_samples_test, y_samples_test



def build_NN_params(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    
    # Load the hyperparameters from the config file
    param = {
        'activation1': config.get('architecture', 'activation1'),
        'activation_out': config.get('architecture', 'activation_out'),
        'res_connection': config.getboolean('architecture', 'res_connection'),
        'dropout': config.getfloat('architecture', 'dropout'),
        'learning_rate': config.getfloat('training', 'learning_rate'),
        'loss': config.get('training', 'loss'),
        'optimizer': config.get('training', 'optimizer'),
        'metrics': [x.strip() for x in config.get('training', 'metrics').split(',')],
        'epochs': config.getint('training', 'epochs'),
        'batch_size': config.getint('training', 'batch_size'),
        'save_interval': config.getint('training', 'save_interval'),
        'response': config.get('main', 'response'),
        'architecture': [int(x.strip()) for x in config.get('architecture', 'architecture').split(',')]
    }
    return param



def build_model_architecture_bootstrap(hyperparams, n_ftrs):
    
    optimizer = hyperparams['optimizer']
    learning_rate = hyperparams['learning_rate']
    if optimizer == 'adam':
        opt = tf.keras.optimizers.Adam(learning_rate=learning_rate)
    elif optimizer == 'sgd':
        opt = tf.keras.optimizers.SGD(learning_rate=learning_rate)

    inputs = Input(shape=(n_ftrs,))
    x = inputs  # Initialize x with inputs

    if hyperparams['architecture'] == [0]:
        x = Dense(1, activation=hyperparams['activation_out'])(x)
    else:
        for i, units in enumerate(hyperparams['architecture']):
            print(i)
            if i == 0:
                x = Dense(units)(x)
                x = BatchNormalization()(x)
                
                if hyperparams['activation1'] == 'leaky_relu':
                    x = LeakyReLU(alpha=0.2)(x)
                elif hyperparams['activation1'] == 'PReLU':
                    x = PReLU()(x)
                else:
                    
                    x = Activation(hyperparams['activation1'])(x)
            else:
                residual_connection = x  # Store the current output for the residual connection
                x = BatchNormalization()(x)
                x = Dense(units)(x)

                if hyperparams['res_connection']:
                    # Implement the residual connection by adding the residual_connection
                    x = Add()([x, residual_connection])

                if hyperparams['activation1'] == 'leaky_relu':
                    x = LeakyReLU(alpha=0.2)(x)
                elif hyperparams['activation1'] == 'PReLU':
                    x = PReLU()(x)
                else:
                    x = Activation(hyperparams['activation1'])(x)


            if hyperparams['dropout'] != 0:
                x = Dropout(hyperparams['dropout'])(x)

        x = Dense(1, activation=hyperparams['activation_out'])(x)

    model = Model(inputs=inputs, outputs=x)
    model.compile(loss=hyperparams['loss'], optimizer=opt)

    return model

def run_nn_bootstrap(X_samples_train, Y_samples_train, NN_hyperparams):
        
    if NN_hyperparams['response'] == "rate":
        Y_tr = (Y_samples_train['nMut']) / (Y_samples_train['length'])
    elif NN_hyperparams['response'] == "count":
        Y_tr = Y_samples_train['nMut']
    else:
        raise ValueError("error")  # Use "raise" to raise an exception
    
    
    n_ftrs = int(X_samples_train.shape[1])
    
    # Build initial model
    model = build_model_architecture_bootstrap(NN_hyperparams, n_ftrs)
    print(model.summary())
    
    model.fit(X_samples_train, Y_tr, epochs=NN_hyperparams['epochs'],
              batch_size=NN_hyperparams['batch_size'], verbose=1) 
   
    model_data = {'model': model,
                  'N': Y_train.N[0],
                  'response_type': NN_hyperparams['response']}

    return model_data




def predict_nn_bootstrap(model_data, X_samples_test, length_elems):
    
    response = model_data['response_type']
    prediction_df = None
    M = model_data['model']
    prediction = M.predict(X_samples_test, verbose = 1)
    N = model_data['N']
    
    if response == "rate":
        prediction = prediction/ N  
    elif response == "count":
        log_offset = N * length_elems
        prediction = prediction[:,0] /np.exp(log_offset)
    else:
        ValueError("error")
    
    prediction_df = pd.DataFrame({'predRate': prediction.ravel()}, 
                                 index=X_samples_test.index)
    return prediction_df

def nn_model_info(save_name, *args):
    NN_hyperparams = build_NN_params(args[0])
    model_dict = {"save_name" : save_name,
                  "Args" : NN_hyperparams,
                  "run_func": run_nn_bootstrap,
                  "predict_func": predict_nn_bootstrap
                  # "save_func": save_nn_bootstrap,
                  # "check_file_func": check_file_nn
                  }
    
    return model_dict


def fit_per_element_bootstrap_nn(X_test, Y_test, elem, NN_hyperparams, n_bootstrap):
    
    pred_samples = pd.DataFrame()
    for n in range(n_bootstrap):
        print(f'bootstrap sample number: {n}')
        X_elem = split_by_element(X_test, elem)
        Y_elem = Y_test.loc[X_elem.index]
        
        # X_elem = X_elem.iloc[0:10]
        # Y_elem = Y_elem.loc[X_elem.index]

        # X_val = split_by_element(X_test, r'[b]\d')
        # Y_val = Y_test.loc[Y_test.index]


        X_samples_train, Y_samples_train, X_samples_test, Y_samples_test = generate_test_train_bootstrapSamples(X_elem, Y_elem)
        length_elems = Y_samples_test.length
        
        # train model and make prediction on a sample:
        model_data = run_nn_bootstrap(X_samples_train, Y_samples_train, 
                                      NN_hyperparams)
        pred_sample = predict_nn_bootstrap(model_data, X_samples_test, length_elems)
        obs_sample = Y_samples_test.nMut/(Y_samples_test.N * Y_samples_test.length)
        corr, p_value = spearmanr(pred_sample, obs_sample)
        print('=============================')
        print(f'obs-pred spearman corr: {corr}') 
        
        pred_samples = pd.concat([pred_samples, pred_sample], axis = 1)
        
        print('=============================')
        ensemble_preds = pd.DataFrame(pred_samples.mean(axis=1))
        tmp_ensemble_obs = Y_test.loc[ensemble_preds.index]
        ensemble_obs = tmp_ensemble_obs.nMut/(tmp_ensemble_obs.N * tmp_ensemble_obs.length)
        
        corr_ensemble, p_value = spearmanr(ensemble_preds, ensemble_obs)
        print('=============================')
        print(f'ensemble obs-pred spearman corr: {corr_ensemble}') 
        print('############ Job Done ##############')
    return ensemble_preds
    
    
    
   

#############################################

sim_file = 'configs/rate_based/sim_setting.ini'

sim_setting = load_sim_settings(sim_file)
# config_save(sim_file)
X_train, Y_train, X_test, Y_test = load_data(sim_setting)

elems = ["gc19_pc.cds", "enhancers", "gc19_pc.3utr", "gc19_pc.5utr",
          "gc19_pc.promCore", "gc19_pc.ss", "lncrna.ncrna", 
          "lncrna.promCore"]

elem = elems[0]


##############################
base_dir = sim_setting['base_dir']
models = sim_setting['models']
model_name = list(models.keys())[0]
m = models[model_name]
name = m['save_name']
NN_hyperparams = m['Args']
save_path_model = f'{base_dir}/{model_name}/'
NN_hyperparams['path_save'] = f'{save_path_model}models_interval/'
Nr_pair_acc = sim_setting['Nr_pair_acc']

############3

n_bootstrap = 2
pred_ensemble_bootstraps = fit_per_element_bootstrap_nn(X_test, Y_test,
                                                        elem, NN_hyperparams, n_bootstrap) 
# pred_ensemble_bootstraps.to_csv(f'{base_dir}/{model_name}/{model_name}_{elem}_ensemble_bootstraps.tsv')

obs_ensemble_bootstraps = X_test.loc[pred_ensemble_bootstraps.index]

assess_model(pred_ensemble_bootstraps, obs_ensemble_bootstraps, Nr_pair_acc,
             model_name, per_element=True)