
import numpy as np
import matplotlib.pyplot as plt
import json
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, UpSampling1D, Dropout, Flatten, Dense
from tensorflow.keras.models import Model
import tensorflow as tf
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten
from readFtrs_Rspns import set_gpu_memory_limit
from scipy.stats import spearmanr


gpu_fraction = 0.05
set_gpu_memory_limit(gpu_fraction)



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
      
      # Ensure y_pred is greater than 0 to avoid log(0)
      y_pred_safe = tf.maximum(y_pred, tf.keras.backend.epsilon())
      
      # Compute element-wise Poisson loss
      loss = y_pred_safe - y_true * tf.math.log(y_pred_safe)
      
      loss *= mask  # Element-wise multiplication, broadcasting loss to match mask shape
      
      # Return the loss without summing and averaging across the sequence
      return loss  # Shape [batch_size, sequence_length]
  
    
#############################################################
############################################################# 

X_simulated = np.load('../../../../Projects/bahari_work/tmp_X_train.npy')
y_simulated = np.load('../../../../Projects/bahari_work/tmp_Y_train.npy')
 
X_simulated_test = np.load('../../../../Projects/bahari_work/tmp_X_test.npy')
y_simulated_test = np.load('../../../../Projects/bahari_work/tmp_Y_test.npy')


print(X_simulated.shape)
print(X_simulated)
print(y_simulated.shape)
print(y_simulated)

# Function to build model from config
model = build_model_from_config('model_config.json')

model.compile(optimizer='adam', loss=custom_poisson_loss)
# Train the model 
model.summary()
model.fit(X_simulated, y_simulated, epochs=400, batch_size=32)

# Extract the middle region values from y_test
# Assuming each sequence in y_test is 100 units long
middle_region_index = 50  # Adjust based on your data specifics
y_test_middle_region = y_simulated_test[:, middle_region_index]
#print(y_test_middle_region)
# Predict using the model
y_pred = model.predict(X_simulated_test)
y_pred_middle_region = y_pred[:, middle_region_index]
#print(y_pred_middle_region)


# Evaluate the performance on the middle region
mae = np.mean(np.abs(y_test_middle_region - y_pred_middle_region))
print(f'Mean Absolute Error for Middle Region: {mae}')

corr, p_value = spearmanr(y_pred_middle_region, y_test_middle_region)
