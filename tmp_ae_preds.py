import pickle
from tensorflow.keras.models import load_model
from readFtrs_Rspns import read_response
from models.runBMR_functions import load_data
import os
import argparse
import platform
from simulation_settings import load_sim_settings

import h5py


sim_file = 'configs/rate_based/sim_settingDSmpl.ini'   
sim_setting = load_sim_settings(sim_file)
# remve_nM=0 False, split_intergenic=False

X_train, Y_train, X_test, Y_test = load_data(sim_setting)


X_train.shape
X_test.shape
###### load ae model
path_model = '../external/procInput/AEreduced_bs64_2000Ep_dim130/_model.h5'
path_params = '../external/procInput/AEreduced_bs64_2000Ep_dim130/_params.pkl'

with h5py.File(path_model, 'r') as f:
    # load the model
    model = load_model(f)


# Load the dictionary from the .pkl file
with open(path_params, 'rb') as file:
    params = pickle.load(file)

preds = model.predict(X_train)



import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# Convert preds to a NumPy array
preds_array = np.array(preds)

# Calculate and print metrics
mse = mean_squared_error(X_train.values, preds_array)
mae = mean_absolute_error(X_train.values, preds_array)
r2 = r2_score(X_train.values, preds_array)


# Calculate the correlation matrix
preds_array = np.array(preds)
correlation_matrix = np.corrcoef(X_train.values.flatten(), preds_array.flatten())

# Extract the correlation coefficient between observed and predicted values
correlation_coefficient = correlation_matrix[0, 1]

# for X_train:
print(f'Mean Squared Error: {mse}') # Mean Squared Error: 0.3745056484469502
print(f'Mean Absolute Error: {mae}') # Mean Absolute Error: 0.367274813585284
print(f'R-squared Score: {r2}') # R-squared Score: 0.6254939197293006
print(f'Correlation Coefficient: {correlation_coefficient}') # Correlation Coefficient: 0.8101391181247778



# for X_test
preds = model.predict(X_test)
# Convert preds to a NumPy array
preds_array = np.array(preds)

# Calculate and print metrics
mse = mean_squared_error(X_test.values, preds_array)
mae = mean_absolute_error(X_test.values, preds_array)
r2 = r2_score(X_test.values, preds_array)


# Calculate the correlation matrix
preds_array = np.array(preds)
correlation_matrix = np.corrcoef(X_test.values.flatten(), preds_array.flatten())

# Extract the correlation coefficient between observed and predicted values
correlation_coefficient = correlation_matrix[0, 1]

print(f'Mean Squared Error: {mse}') # Mean Squared Error: 24.759909770296893

print(f'Mean Absolute Error: {mae}') # Mean Absolute Error: 1.2746613971348946

print(f'R-squared Score: {r2}') # R-squared Score: 0.46866942753767504

print(f'Correlation Coefficient: {correlation_coefficient}') # Correlation Coefficient: 0.8878968523662699
