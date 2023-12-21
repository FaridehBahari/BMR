import pandas as pd
import os
import numpy as np
import statsmodels.api as sm
import pickle
from scipy.special import logit
from sklearn.linear_model import LassoCV, RandomizedLasso
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

# os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
from readFtrs_Rspns import save_preds_tsv
import configparser

def build_GLM_params(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    # Load the hyperparameters from the config file
    params = {
        'fi_cut': config.getfloat('main', 'fi_cut')
    }
    return params


def read_fi(path, cutoff=0.8):
    """Read feature importance table in TSV format.

    Feature importance table must contain two columns: name and importance

    Args:
        path (str): path to the file.
        cutoff (float): cutoff of feature selection.

    Returns:
        list: useful features. Return None if path is None.

    """
    fi = pd.read_csv(path, sep='\t', header=0, index_col='name',
                     usecols=['name', 'importance'])
    assert len(fi.index.values) == len(fi.index.unique()), \
        "Feature name in feature importance table is not unique."
    keep = (fi.importance >= cutoff).values
    use_features = fi.loc[keep]
    use_features = use_features.index.values
    
    return use_features



def run_lasso(X, y, max_iter=3000, cv=5, n_threads=1):
    """ Implement LassoCV in sklearn
    
    Args:
        X (np.array): scaled X.
        y (pd.df): four columns response table. 
        max_iter (int): max iteration. 
        cv (int): CV fold.
        n_threads (int): Number of threads to use for parallel computing.

    Returns:
        float: trained alpha value.

    """
    # generate logit response
    y_logit = logit((y.nMut + 0.5) / (y.length * y.N))
    # sub-sampling X and y (300,000)
    use_ix = np.random.choice(y_logit.shape[0], 300000, replace=False)
    Xsub = X[use_ix, :]
    ysub = y_logit[use_ix]
    reg = LassoCV(max_iter=max_iter, cv=cv, copy_X=False, n_jobs=n_threads)
    lassocv = reg.fit(Xsub, ysub)
    return lassocv.alpha_

def run_rndlasso(X, y, alpha,
    n_resampling=500, sample_fraction=0.1, n_threads=1):
    """  Implement Randomized Lasso in sklearn

    Args:
        X (np.array): scaled X. 
        y (pd.df): four columns response table. 
        alpha (float): parameter trained from lassoCV 
        n_resampling (int): number of times for resampling 
        sample_fraction (float): fraction of data to use at each resampling

    Returns:
        np.array: feature importance scores

    """
    
    # generate logit response
    y_logit = logit((y.nMut + 0.5) / (y.length * y.N))
    reg = RandomizedLasso(alpha=alpha,
                          n_resampling=n_resampling,
                          sample_fraction=sample_fraction,
                          selection_threshold=1e-3,
                          max_iter=3000,
                          normalize=False,
                          n_jobs=n_threads)
    rndlasso = reg.fit(X, y_logit)
    fi_scores = rndlasso.scores_
    return fi_scores

def save_fi(fi_scores, feature_names, project_name, out_dir):
    """ Save feature importance results to disk
    
    Args:
        fi_scores (np.array): array of feature importance scores.
        feature_names (np.array): array of feature names.
        project_name (str): name of the project, prefix of the output file.
        out_dir (str): output directory. 

    Returns:
        pd.DF: feature importance table
        
    """
    res = pd.DataFrame({'name':feature_names, 'importance':fi_scores}, columns=['name', 'importance'])
    path = os.path.join(out_dir, project_name+'.feature_importance.tsv')
    res.to_csv(path, sep='\t', index=False)
    return res


def run_glm(X_train, Y_train, param):
    fi_cut = param['fi_cut']
    path = param['path_save']
    parts = path.split('/')
    model_name = parts[-3]
    
    feature_names = X_train.columns.values
    
    # Run lasso to get alpha
    alpha = run_lasso(X_train, Y_train)
    # Run rnd lasso to get feature importance
    fi_scores = run_rndlasso(X_train, Y_train, alpha)
    fi = save_fi(fi_scores, feature_names, model_name, path)
    # Remove unimportant features
    keep = (fi.importance >= fi_cut).values
    use_features = fi.name.values[keep]
    X_train = X_train[:, np.isin(feature_names, use_features)]
    
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