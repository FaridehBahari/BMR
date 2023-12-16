from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
import pandas as pd
from readFtrs_Rspns import paste0



def train_autoencoder_model(encoding_dim, train_X, params):
    # Define the shape of the input data
    input_shape = train_X.shape[1]

    # Define the input layer
    input_data = Input(shape=(input_shape,))

    # Define the encoded layer
    encoded = Dense(250, activation=params['activation_encoded'])(input_data)
    encoded = Dense(200, activation=params['activation_encoded'])(encoded)
    encoded = Dense(150, activation=params['activation_encoded'])(encoded)

    # bottleneck
    latent_space = Dense(encoding_dim, activation=params['activation_bottleneck'])(encoded)

    # Define the decoded layer
    decoded = Dense(150, activation=params['activation_decoded'])(latent_space)
    decoded = Dense(200, activation=params['activation_decoded'])(decoded)
    decoded = Dense(250, activation=params['activation_decoded'])(decoded)
    decoded = Dense(input_shape, activation=params['activation_decoded'])(decoded)

    # Create the autoencoder model
    autoencoder = Model(input_data, decoded)

    # Compile the model
    autoencoder.compile(params['optimizer'], params['loss'])

    # Train the autoencoder using the training data
    history = autoencoder.fit(train_X, train_X, epochs=params['epochs'], 
                              batch_size=params['batch_size'],
                              shuffle=True)
    
    encoder = Model(input_data, latent_space)
    return encoder, history, autoencoder


def dim_reduction_AE(encoder, subject_X):
    # Encode the training and test data using the encoded layer
    encoded_subject_X = encoder.predict(subject_X)

    # convert reduced feature matrices to dataframe and save them
    reduced_X = pd.DataFrame(encoded_subject_X, index=subject_X.index, 
                             columns=paste0('F', range(encoded_subject_X.shape[1])))
    return reduced_X
