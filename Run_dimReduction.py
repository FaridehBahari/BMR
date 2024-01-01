import os
import pandas as pd
from readFtrs_Rspns import create_TestTrain_TwoSources, set_gpu_memory_limit
from dimReduction.autoencoder import train_autoencoder_model, dim_reduction_AE
from dimReduction.PCA import train_PCA, dim_reduction_PCA
import matplotlib.pyplot as plt
import pickle

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)


os.makedirs('../external/procInput/', exist_ok= True)
########### load data ################   
path_X_test = '../external/rawInput/test_feature.hdf5'
path_X_train = '../external/rawInput/train_feature.hdf5'
path_Y_train = '../external/rawInput/Pan_Cancer_train_y.tsv'
path_Y_test = '../external/rawInput/Pan_Cancer_test_y.tsv'

X_train, Y_train, X_test, Y_test = create_TestTrain_TwoSources(path_X_train, 
                                                           path_Y_train, 
                                                           path_X_test, 
                                                           path_Y_test)

############### autoencoder for dim reduction ##############
# n_dim = 130
# params = {'activation_encoded': 'relu',
#           'activation_decoded': 'relu',
#           'activation_bottleneck': 'relu',
#           'loss': 'mse',
#           'optimizer': 'adam',
#           'epochs': 2000,
#           'batch_size': 64}

# # Train the autoencoder and get the history
# encoder, history, ae_model = train_autoencoder_model(n_dim, X_train, params)

# # save the plot and reduced test and train Xs
# path_save = '../external/procInput/AEreduced_bs64_2000Ep_dim130/'
# os.makedirs(path_save, exist_ok= True)

# ae_model.save(f'{path_save}_model.h5')
        
# # save the model and its params
# with open(f'{path_save}_params.pkl', 'wb') as f: 
#     pickle.dump(params, f)


# # Plot the loss curve
# plt.plot(history.history['loss'])
# plt.title('Model Loss')
# plt.xlabel('Epoch')
# plt.ylabel('Loss')
# # plt.show()
# # Save the plot to a file
# plt.savefig(f'{path_save}loss_curve.png')

# # Apply the trained encoder to reduce dimensionality
# reduced_test = dim_reduction_AE(encoder, X_test)
# reduced_train = dim_reduction_AE(encoder, X_train)
# # save
# reduced_test.to_csv(f'{path_save}test_AEreduced.tsv', sep = '\t')
# reduced_train.to_csv(f'{path_save}train_AEreduced.tsv', sep = '\t')
# print('files saved successfuly...')

# ###############
# import pandas as pd
# path_save = '../external/procInput/AEreduced_bs64_2000Ep_dim130/'
# reduced_test = pd.read_csv(f'{path_save}test_AEreduced.tsv', sep = '\t', 
#                            index_col='binID')
# non_zero_columns_test = reduced_test.loc[:, (reduced_test != 0).any()]




# reduced_train = pd.read_csv(f'{path_save}train_AEreduced.tsv', sep = '\t', 
#                            index_col='binID')
# non_zero_columns_train = reduced_train.loc[:, (reduced_train != 0).any()]



# columns_not_in_test = non_zero_columns_train.columns.difference(non_zero_columns_test.columns)

# # Columns in non_zero_columns_test but not in non_zero_columns_train
# columns_not_in_train = non_zero_columns_test.columns.difference(non_zero_columns_train.columns)

# # Display the results
# print("Columns not in non_zero_columns_test:", columns_not_in_test)
# print("Columns not in non_zero_columns_train:", columns_not_in_train)



# new_reduced_test = reduced_test.loc[:, (non_zero_columns_train.columns)]
# new_reduced_train = non_zero_columns_train.copy()

# new_reduced_test.to_csv(f'{path_save}test_AEreduced_75F.tsv', sep = '\t')
# new_reduced_train.to_csv(f'{path_save}train_AEreduced_75F.tsv', sep = '\t')
############### PCA for dim reduction ##############
variance_threshold = 0.8  # .9...370pc-------------------.85...221pc----------------.8...130pc
pca_model = train_PCA(X_train, variance_threshold)
pca_reduced_train = dim_reduction_PCA(pca_model, X_train)
pca_reduced_test = dim_reduction_PCA(pca_model, X_test)

pca_reduced_train.to_csv('../external/procInput/train_PCAreduced.tsv', sep = '\t')
pca_reduced_test.to_csv('../external/procInput/test_PCAreduced.tsv', sep = '\t')