import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten
from keras.losses import Poisson
from scipy.stats import spearmanr
from readFtrs_Rspns import set_gpu_memory_limit

gpu_fraction = 0.05
set_gpu_memory_limit(gpu_fraction)

# Define the model architecture
def build_basic_nn_model(input_shape):
    model = Sequential()
    model.add(Flatten(input_shape=input_shape))  # Flatten the input
    model.add(Dense(64, activation='relu'))      # First hidden layer with 64 neurons
    model.add(Dense(64, activation='relu'))      # Second hidden layer with 64 neurons
    model.add(Dense(1, activation='softplus'))     # Output layer for regression
    
    return model


def prepare_training_data(X_simulated, y_simulated):
    """
    Prepares training data from simulated data by filtering out instances with missing mutation rates.
    
    Parameters:
    X_simulated (numpy.ndarray): Simulated feature data.
    y_simulated (numpy.ndarray): Simulated target data (mutation rates).
    
    Returns:
    X_train, y_train: Filtered training data.
    """
    
    # Filter out instances where mutation rates are missing
    # Assuming missing values are represented as None, NaN, or a specific marker (e.g., -1)
    valid_indices = (y_simulated != -1)
    X_filtered = X_simulated[valid_indices]
    y_filtered = y_simulated[valid_indices]
    
    return X_filtered, y_filtered



def prepare_test_data(X_simulated, y_simulated, middle_region_index=50):
    """
    Prepares test data from simulated data by filtering out instances where mutation rates are missing
    (indicated by a value of -1) and extracts the middle region for evaluation.
    
    Parameters:
    X_simulated (numpy.ndarray): Simulated feature data.
    y_simulated (numpy.ndarray): Simulated target data (mutation rates).
    middle_region_index (int): Index of the middle region to focus on in the evaluation.
    
    Returns:
    X_test, y_test_middle_region: Filtered test data and corresponding middle region mutation rates.
    """
    
    # Filter out instances where mutation rates are missing (indicated by -1)
    valid_indices = (y_simulated != -1)
    X_filtered = X_simulated[valid_indices]
    y_filtered = y_simulated[valid_indices]
    
    # Extract the middle region values from y_filtered
    y_test_middle_region = y_filtered[:, middle_region_index]
    
    return X_filtered, y_test_middle_region


#################################################################################
#################################################################################
X_simulated = np.load('../../../../Projects/bahari_work/tmp_X_train.npy')
y_simulated = np.load('../../../../Projects/bahari_work/tmp_Y_train.npy')
 
X_simulated_test = np.load('../../../../Projects/bahari_work/tmp_X_test.npy')
y_simulated_test = np.load('../../../../Projects/bahari_work/tmp_Y_test.npy')

middle_region_index = 50


X_train, y_train = prepare_training_data(X_simulated, y_simulated)



# Assuming X_train.shape[1:] gives the shape of a single input sample
model_simple = build_basic_nn_model(X_train.shape[1:])

#model.compile(optimizer='adam',
 #             loss='mean_squared_error',  # Mean Squared Error for regression
  #            metrics=['mean_absolute_error'])  # Including MAE for interpretability



model_simple.compile(optimizer='adam',
              loss=Poisson(),  # Poisson loss for regression
              metrics=['mean_absolute_error'])

history = model_simple.fit(X_train, y_train,
                    epochs=50,             # Number of epochs
                    batch_size=32,          # Batch size
                    validation_split=0.2)   # Split part of the training data for validation


# Extract features of the middle region (index 50) from each sequence
X_test_middle_region = X_simulated_test[:, middle_region_index, :]

# Now, X_test_middle_region has the shape (num_samples, features_of_single_bin)
# Predict using the model
y_pred = model_simple.predict(X_test_middle_region)

# Since y_pred will have a single column, there's no need to extract a specific column for the middle region
y_pred_middle_region = y_pred.flatten()
y_test_middle_region = y_simulated_test[:, middle_region_index]

# Evaluate the performance on the middle region
# Calculating mean absolute error
mae = np.mean(np.abs(y_test_middle_region - y_pred_middle_region))
print(f'Mean Absolute Error for Middle Region: {mae}')

from scipy.stats import spearmanr
corr, p_value = spearmanr(y_pred_middle_region, y_test_middle_region)