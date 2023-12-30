import pandas as pd
import os
import numpy as np
import statsmodels.api as sm
import pickle
from readFtrs_Rspns import read_fi
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

# os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
from readFtrs_Rspns import save_preds_tsv
import configparser
def build_GLM_params(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    # Load the hyperparameters from the config file
    params = {
        
    }
    return params


def run_glm(X_train, Y_train, param):
    
    use_features = read_fi('../external/rawInput/DP_results/Binom_GLM.feature_importance.tsv',
                           cutoff=0.5)
    X_train = X_train[use_features]
    # Add const manually. sm.add_constant cannot add 1 for shape (1, n)
    X_train = np.c_[X_train, np.ones(X_train.shape[0])]
    # make two columns response (# success, # failure)
    y_binom = np.zeros((Y_train.shape[0], 2), dtype=np.int_)
    y_binom[:, 0] = Y_train.nMut
    y_binom[:, 1] = Y_train.length * Y_train.N - Y_train.nMut
    glm = sm.GLM(y_binom, X_train, family=sm.families.Binomial())
    model = glm.fit()
    
    model_data = {"model": model,
                  "use_features": use_features
                  }
    
    return model_data

def predict_glm(model, X_test, length_elems):
    use_features = model['use_features']
    binID = X_test.index
    X_test = X_test[use_features]
    X_test = np.c_[X_test, np.ones(X_test.shape[0])]
    M = model['model']
    pred_test = M.predict(X_test) 
    prediction_df = pd.DataFrame({'predRate': pred_test.ravel()}, 
                                 index=binID) 
    return prediction_df


def GLM_model_info(save_name, *args):
    params = build_GLM_params(args[0])
    model_dict = {"save_name" : save_name,
                  "Args" : params,
                  "run_func": run_glm,
                  "predict_func": predict_glm,
                  "save_func": save_glm,
                  # "check_file_func": check_file_glm
                  }
    
    return model_dict


def save_glm(fitted_Model, path_save, save_name, save_model = True): 
    
    save_preds_tsv(fitted_Model, path_save, save_name)
    
    if save_model:
        M = fitted_Model.model['model']
        # Save the model using pickle
        save_path_model = f'{path_save}/{save_name}/{save_name}_model.pkl'
        with open(save_path_model, 'wb') as f: 
            pickle.dump(M, f)
    

def check_file_glm(path_save, save_name):
    file_name = f'{path_save}/{save_name}/{save_name}_model.pkl'
    return file_name