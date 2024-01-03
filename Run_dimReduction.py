import time
import os
import platform
import argparse
from readFtrs_Rspns import set_gpu_memory_limit
from models.runBMR_functions import  load_data_sim, config_save
from simulation_settings import load_sim_settings
from dimReduction.PCA import save_PCA_reduced

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)


if platform.system() == 'Windows':
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    sim_file = 'configs/quick/dimReduction.ini'
    # sim_file = 'configs/rate_based/dimReduction.ini'
else:
    parser = argparse.ArgumentParser(description='Train different models for rate prediction on intergenic bins')
    parser.add_argument('sim_file', type=str, help='the path to the simulation setting config')
    args = parser.parse_args()
    sim_file = args.sim_file

config_save(sim_file)
########### load data ################      
st_time = time.time()
sim_setting = load_sim_settings(sim_file)

X_train, Y_train, X_test, Y_test = load_data_sim(sim_setting)
X_train = X_train.fillna(0)


############### PCA for dim reduction ##############
save_name = list(sim_setting['models'].keys())[0]
save_PCA_reduced(X_train, X_test, sim_setting['base_dir'], save_name, 0.8)

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
