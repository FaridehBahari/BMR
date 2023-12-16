# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 12:08:51 2023

@author: Farideh
"""

# Import matplotlib library
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from performance.assessModels import read_obs, read_pred, split_by_element
from scipy.stats import spearmanr

def create_unseen_obs_pred_df(path_Y_test, path_Y_val, path_pred_unseen):
    Y_obs_all_elems = read_obs(path_Y_test)
    Y_obs_val = read_obs(path_Y_val)
    
    Y_obs_unseen = pd.concat([Y_obs_all_elems, Y_obs_val], axis= 0)
    Y_pred_unseen = read_pred(path_pred_unseen)
    
    obs_elem = split_by_element(Y_obs_unseen,  r'[b]\d')
    pred_elem = split_by_element(Y_pred_unseen,  r'[b]\d')
    
    pred_elem = pred_elem.loc[obs_elem.index]
    
    return obs_elem, pred_elem


def plot_scatterPlot(obs, pred, model_name, color):
    obs = obs["obsRates"]
    pred = pred["predRate"]
    
    min_value = np.min([pred, obs]) 
    max_value = np.max([pred, obs]) 
    
    print(f'spearman corr: {spearmanr(pred, obs)}')
    label = model_name
    fig, ax = plt.subplots()
    ax.scatter(pred, obs, label=label, color=color)
    
    ax.plot(np.unique(pred), 
              np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)),
              color = 'black')
    
    # Set both axes to log scale
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    # Set the axis limits to ensure the same intervals
    ax.set_xlim([min_value, max_value])
    ax.set_ylim([min_value, max_value])

    # Add labels and legend
    ax.set_xlabel("Prediction")
    ax.set_ylabel("Observed")
       
    ax.legend()
    

# path_pred_unseen = '../external/output/N2_60k_800_do6_cd5025/N2_60k_800_do6_cd5025_predTest.tsv'
# path_pred_unseen = '../external/output/nn_pois/nn_pois_predTest.tsv'
# path_pred_unseen = '../external/output/N2_60k_800_do6_cd5025/N2_60k_800_do6_cd5025_predTest.tsv'
path_pred_unseen = '../external/output/GBM/GBM_predTest.tsv'
# path_pred_unseen = '../external/output/nn_classic_small/nn_classic_small_predTest.tsv'
path_pred_unseen = '../external/output/220_epochs50/220_epochs50_predTest.tsv'
path_Y_val = '../external/rawInput/validation_sets_10folds/Pan_Cancer_validate_y_fold_1.tsv'
path_Y_test = '../external/rawInput/Pan_Cancer_test_y.tsv'

obs, pred = create_unseen_obs_pred_df(path_Y_test, path_Y_val, 
                                      path_pred_unseen)

model_name = os.path.basename(path_pred_unseen)[:-13]
color = 'orange'
plot_scatterPlot(obs, pred, model_name, color)






# import seaborn as sn
# sn.regplot(x = pred["predRate"], y = obs["obsRates"])
