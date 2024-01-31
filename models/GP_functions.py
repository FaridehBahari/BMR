import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

# Let's assume 'df' is your pandas DataFrame.
# Replace this line with your actual DataFrame loading step if needed.
df_raw = pd.read_csv('../external/BMR/output/bin_size_effect/var_size_nonMutInclude_longer50/GBM/chr22_dat.csv')

# remove NAs
df = df_raw.dropna(subset=['predRate'])

# We will use the middle point between 'start' and 'end' for our 'X' axis in the Gaussian process.
X = np.array((df['start'] + df['end']) / 2).reshape(-1, 1)
y = np.array(df['predRate'])

# Define Gaussian Process Model
# We combine the RBF kernel with a constant kernel.
kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)

# Fit to data using Maximum Likelihood Estimation of the parameters
gp.fit(X, y)

# Make the prediction on the meshed x-axis (ask for MSE as well)
y_pred, sigma = gp.predict(X, return_std=True)

# Add the predictions back to the dataframe
df['smoothed_predRate'] = y_pred

# If you want to examine the confidence intervals you could use sigma to create them
df['confidence_interval_lower'] = y_pred - 2 * sigma
df['confidence_interval_upper'] = y_pred + 2 * sigma

#######################################################################

import statsmodels.api as sm

# Use the LOWESS (Locally Weighted Scatterplot Smoothing) from statsmodels
lowess = sm.nonparametric.lowess
smoothed = lowess(y, X, frac=0.025)  # You can adjust the frac value to change the level of smoothing

# Add the smoothed values back to the dataframe
df['smoothed_predRate'] = smoothed[:, 1]

# Now 'df' has an additional column 'smoothed_predRate' with the LOESS smoothed rates.

