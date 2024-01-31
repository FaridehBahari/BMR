import numpy as np
import matplotlib.pyplot as plt
import json
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, UpSampling1D, Dropout, Flatten, Dense
from tensorflow.keras.models import Model
import tensorflow as tf
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten
from tensorflow.keras.layers import Input, Dense, LayerNormalization, Dropout, MultiHeadAttention
from readFtrs_Rspns import set_gpu_memory_limit
from scipy.stats import spearmanr

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)



# Positional Encoding
def get_positional_encoding(seq_length, d_model):
    pos_encoding = np.array([
        [pos / np.power(10000, 2 * (i // 2) / d_model) for i in range(d_model)]
        for pos in range(seq_length)
    ])
    pos_encoding[1:, 0::2] = np.sin(pos_encoding[1:, 0::2])
    pos_encoding[1:, 1::2] = np.cos(pos_encoding[1:, 1::2])
    pos_encoding = tf.cast(pos_encoding, dtype=tf.float32)
    pos_encoding = tf.expand_dims(pos_encoding, 0)  # Expand to 3D for batch compatibility
    return pos_encoding


# Transformer Encoder Layer
class TransformerEncoderLayer(tf.keras.layers.Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(TransformerEncoderLayer, self).__init__()
        self.mha = MultiHeadAttention(num_heads=num_heads, key_dim=d_model)
        self.ffn = tf.keras.Sequential([
            Dense(dff, activation='relu'),
            Dense(d_model)
        ])
        self.layernorm1 = LayerNormalization(epsilon=1e-6)
        self.layernorm2 = LayerNormalization(epsilon=1e-6)
        self.dropout1 = Dropout(rate)
        self.dropout2 = Dropout(rate)
        
    def call(self, x, training, mask=None):
        attn_output = self.mha(x, x, attention_mask=mask)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(x + attn_output)
        ffn_output = self.ffn(out1)
        ffn_output = self.dropout2(ffn_output, training=training)
        out2 = self.layernorm2(out1 + ffn_output)
        return out2

def build_transformer_model(input_shape, num_layers=4, d_model=128, num_heads=4, dff=512, rate=0.1):
    seq_length = input_shape[0]
    num_features = input_shape[1]
    pos_encoding = get_positional_encoding(seq_length, d_model)
    
    inputs = Input(shape=input_shape)
    x = Dense(d_model)(inputs)
    x *= tf.sqrt(tf.cast(d_model, tf.float32))  # Scale embeddings
    x += pos_encoding  # Add positional encoding
    
    for _ in range(num_layers):
        x = TransformerEncoderLayer(d_model, num_heads, dff, rate)(x)
        
    outputs = Dense(1, activation='softplus')(x)
    model = Model(inputs=inputs, outputs=outputs)
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


#############################################################################
#############################################################################
X_simulated = np.load('../../../../Projects/bahari_work/tmp_X_train.npy')
y_simulated = np.load('../../../../Projects/bahari_work/tmp_Y_train.npy')
 
X_simulated_test = np.load('../../../../Projects/bahari_work/tmp_X_test.npy')
y_simulated_test = np.load('../../../../Projects/bahari_work/tmp_Y_test.npy')



# Model Configuration
input_shape = (100, 1500)  # 100 regions, each with 1500 features
transformer_model = build_transformer_model(input_shape)
transformer_model.compile(optimizer='adam', loss=custom_poisson_loss)


transformer_model.fit(X_simulated, y_simulated, epochs=400, batch_size=32)

# Evaluate the performance on the middle region
middle_region_index = 50
y_test_middle_region = y_simulated_test[:, middle_region_index]
#print(y_test_middle_region)
y_pred = transformer_model.predict(X_simulated_test)
y_pred_middle_region = y_pred[:, middle_region_index, :]
#print(y_pred_middle_region)

mae = np.mean(np.abs(y_test_middle_region - y_pred_middle_region))
print(f'Mean Absolute Error for Middle Region: {mae}')


corr, p_value = spearmanr(y_pred_middle_region, y_test_middle_region)