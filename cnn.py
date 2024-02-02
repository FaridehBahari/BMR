
import numpy as np
import matplotlib.pyplot as plt
import json
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, UpSampling1D, Dropout, Flatten, Dense
from tensorflow.keras.models import Model
import tensorflow as tf


def simulate_data(num_samples=1000, input_shape=(100, 1000), influential_feature_count=10, unavailable_rate_probability=0.1):
    X_simulated = np.random.normal(size=(num_samples, *input_shape))

    # Selecting random indices for influential features
    influential_features_indices = np.random.choice(input_shape[1], size=influential_feature_count, replace=False)
    X_influence = np.sum(X_simulated[:, :, influential_features_indices], axis=2)

    # Base mutation rates influenced by X
    y_base = np.random.poisson(lam=np.clip(X_influence / 10, 0.1, 10), size=(num_samples, input_shape[0]))

    # Weak influence from neighboring regions
    neighborhood_influence = (np.roll(y_base, 1, axis=1) + np.roll(y_base, -1, axis=1)) / 20
    y_simulated = y_base + neighborhood_influence

    # Introduce unavailable mutation rates
    mask = np.random.rand(num_samples, input_shape[0]) > unavailable_rate_probability
    y_simulated = np.where(mask, y_simulated, -1.0)  # Set unavailable rates to -1

    return X_simulated, y_simulated



def build_model_from_config(config_path):
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)

    input_shape = tuple(config['input_shape'])
    inputs = Input(shape=input_shape)

    x = inputs
    for layer_conf in config['layers']:
        if layer_conf['type'] == 'Conv1D':
            # Select relevant parameters for Conv1D
            conv_params = {k: v for k, v in config['default_conv_params'].items()}
            conv_params.update({k: v for k, v in layer_conf.items() if k != 'type'})
            x = Conv1D(**conv_params)(x)
        elif layer_conf['type'] == 'MaxPooling1D':
            x = MaxPooling1D(**{k: v for k, v in layer_conf.items() if k != 'type'})(x)
        elif layer_conf['type'] == 'UpSampling1D':
            x = UpSampling1D(**{k: v for k, v in layer_conf.items() if k != 'type'})(x)
        elif layer_conf['type'] == 'Dropout':
            x = Dropout(**{k: v for k, v in layer_conf.items() if k != 'type'})(x)

    model = Model(inputs=inputs, outputs=x)
    return model

def custom_poisson_loss(y_true, y_pred):
    mask = tf.cast(tf.not_equal(y_true, -1.0), tf.float32)  # Mask for available rates
    loss = tf.keras.losses.poisson(y_true, y_pred)  # Shape [batch_size]

    # Expand the dimensions of loss to match the mask
    loss = tf.expand_dims(loss, axis=-1)  # Now loss shape becomes [batch_size, 1]

    loss *= mask  # Element-wise multiplication, broadcasting loss to match mask shape

    # Return the loss without summing and averaging across the sequence
    return loss  # Shape [batch_size, sequence_length]


# Simulate test data
X_simulated_test, y_simulated_test = simulate_data(num_samples=500, unavailable_rate_probability=0)
print(X_simulated_test)
print(y_simulated_test)


# Simulate data
X_simulated, y_simulated = simulate_data(num_samples=1000)

print(y_simulated.shape)
y_simulated[1,0:99]
# Function to build model from config
model = build_model_from_config('../model_config.json')

model.compile(optimizer='adam', loss=custom_poisson_loss)
# Train the model on the simulated data
#model.summary()
model.fit(X_simulated, y_simulated, epochs=200, batch_size=32)



# Extract the middle region values from y_test
# Assuming each sequence in y_test is 100 units long
middle_region_index = 50  # Adjust based on your data specifics
y_test_middle_region = y_simulated_test[:, middle_region_index]

# Predict using the model
y_pred = model.predict(X_simulated_test)
y_pred_middle_region = y_pred[:, middle_region_index]

# Evaluate the performance on the middle region
# You might use a custom evaluation metric as appropriate
# For instance, calculating mean absolute error
mae = np.mean(np.abs(y_test_middle_region - y_pred_middle_region))
print(f'Mean Absolute Error for Middle Region: {mae}')

# Compare to simple NN
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

import numpy as np

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


# Example usage
X_train, y_train = prepare_training_data(X_simulated, y_simulated)

import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten

# Define the model architecture
def build_basic_nn_model(input_shape):
    model = Sequential()
    model.add(Flatten(input_shape=input_shape))  # Flatten the input
    model.add(Dense(64, activation='relu'))      # First hidden layer with 64 neurons
    model.add(Dense(64, activation='relu'))      # Second hidden layer with 64 neurons
    model.add(Dense(1, activation='relu'))     # Output layer for regression

    return model


# Assuming X_train.shape[1:] gives the shape of a single input sample
model = build_basic_nn_model(X_train.shape[1:])

#model.compile(optimizer='adam',
 #             loss='mean_squared_error',  # Mean Squared Error for regression
  #            metrics=['mean_absolute_error'])  # Including MAE for interpretability

from keras.losses import Poisson

model.compile(optimizer='adam',
              loss=Poisson(),  # Poisson loss for regression
              metrics=['mean_absolute_error'])

history = model.fit(X_train, y_train,
                    epochs=50,             # Number of epochs
                    batch_size=32,          # Batch size
                    validation_split=0.2)   # Split part of the training data for validation


import numpy as np


# Extract features of the middle region (index 50) from each sequence
X_test_middle_region = X_simulated_test[:, middle_region_index, :]

# Now, X_test_middle_region has the shape (num_samples, features_of_single_bin)
# Predict using the model
y_pred = model.predict(X_test_middle_region)

# Since y_pred will have a single column, there's no need to extract a specific column for the middle region
y_pred_middle_region = y_pred.flatten()

# Evaluate the performance on the middle region
# Calculating mean absolute error
mae = np.mean(np.abs(y_test_middle_region - y_pred_middle_region))
print(f'Mean Absolute Error for Middle Region: {mae}')

