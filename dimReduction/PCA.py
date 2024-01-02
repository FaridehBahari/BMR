# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 16:38:05 2023

@author: Farideh
"""
import os
import pandas as pd 
from sklearn.decomposition import PCA 
# import seaborn as sns # to plot the heatmap
from readFtrs_Rspns import paste0

# #Check the Co-relation between features without PCA
# sns.heatmap(scaled_data.corr())

# def train_PCA(n_dim, X_train):
    
#     pca_model = PCA(n_components = n_dim)
#     pca_model.fit(X_train)
    
#     return pca_model

# def dim_reduction_PCA(pca_model, subject_X): 
#     data_pca = pca_model.transform(subject_X)
#     data_pca = pd.DataFrame(data_pca, index = subject_X.index, 
#                           columns = paste0('F', range(data_pca.shape[1])))


def train_PCA(X_train, variance_threshold):
    pca_model = PCA(n_components=None)  # Use None to get all components initially
    pca_model.fit(X_train)

    # Calculate the cumulative explained variance
    explained_variance_ratio_cumsum = pca_model.explained_variance_ratio_.cumsum()

    # Determine the number of components that explain the specified threshold of variance
    n_components = (explained_variance_ratio_cumsum >= variance_threshold).argmax() + 1

    # Now, fit the PCA model with the selected number of components
    pca_model = PCA(n_components=n_components)
    pca_model.fit(X_train)

    return pca_model

def dim_reduction_PCA(pca_model, subject_X): 
    data_pca = pca_model.transform(subject_X)
    data_pca = pd.DataFrame(data_pca, index=subject_X.index, 
                            columns=paste0('F', range(data_pca.shape[1])))
    
    return data_pca



def save_PCA_reduced(X_train, X_test, path_save, save_name, variance_threshold = 0.8):
    os.makedirs(f'{path_save}/{save_name}/', exist_ok= True)
    pca_model = train_PCA(X_train, variance_threshold)
    pca_reduced_train = dim_reduction_PCA(pca_model, X_train)
    pca_reduced_test = dim_reduction_PCA(pca_model, X_test)

    pca_reduced_train.to_csv(f'{path_save}/{save_name}/train_PCAreduced.tsv',
                             sep = '\t')
    pca_reduced_test.to_csv(f'{path_save}/{save_name}/test_PCAreduced.tsv', 
                            sep = '\t')