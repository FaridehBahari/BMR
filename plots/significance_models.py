from models.repeated_train_test import compare_two_models
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# List of directory paths to compare
dir_paths = [
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/RF/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_poisLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_mseLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/RF/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_poisLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_mseLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/RF/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_poisLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_mseLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/RF/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_poisLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_mseLoss/rep_train_test/'
    
]

model_names = ['Variabl size XGBoost', 'Variabl size RF','Variabl size NN (Poisson loss)','Variabl size NN (MSE loss)',
               '10k XGBoost', '10k RF','10k NN (Poisson loss)','10k NN (MSE loss)',
               '50k XGBoost', '50k RF','50k NN (Poisson loss)','50k NN (MSE loss)',
               '100k XGBoost', '100k RF','100k NN (Poisson loss)','100k NN (MSE loss)',
               '1M XGBoost', '1M RF','1M NN (Poisson loss)','1M NN (MSE loss)'
               ]

# Metric to compare
performance_metric = 'corr'

# Initialize a DataFrame to store p-values with model names as both index and columns
p_value_matrix = pd.DataFrame(index=model_names, columns=model_names)

# Calculate the p-value for each pair of directory paths
for i, dir_path1 in enumerate(dir_paths):
    for j, dir_path2 in enumerate(dir_paths):
        if i < j:  # Only compute for pairs (avoid redundant comparisons)
            p_value = compare_two_models(dir_path1, dir_path2, performance_metric)
            p_value_matrix.loc[model_names[i], model_names[j]] = p_value
            p_value_matrix.loc[model_names[j], model_names[i]] = p_value  # Symmetric matrix

# Display or save the p-value matrix
print(p_value_matrix)
# Optionally save to CSV
p_value_matrix.to_csv("../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_model_binGeneration.csv")

###############################################################################
# List of directory paths to compare
dir_paths = [
    '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/GBM/rep_train_test/',
  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/GBM/rep_train_test/',
  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/GBM/rep_train_test/',
  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/GBM/rep_train_test/',
  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/GBM/rep_train_test/',
  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/GBM/rep_train_test/',
  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/GBM/rep_train_test/',
  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/RF/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/RF/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/RF/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/RF/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/RF/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/RF/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/RF/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/nn_poisLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/nn_poisLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/nn_poisLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/nn_poisLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/nn_poisLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/nn_poisLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/nn_poisLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/nn_mseLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/nn_mseLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/nn_mseLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/nn_mseLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/nn_mseLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/nn_mseLoss/rep_train_test/',
'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/nn_mseLoss/rep_train_test/'
    
]

model_names = ['XGBoost (Full set) ', 'XGBoost (DS1M) ','XGBoost (DS800k) ','XGBoost (DS600k) ','XGBoost (DS300k) ','XGBoost (DS100k) ','XGBoost (DS50k) ',
               'RF (Full set) ', 'RF (DS1M) ','RF (DS800k) ','RF (DS600k) ','RF (DS300k) ','RF (DS100k) ','RF (DS50k) ',
               'NN (Poisson loss, Full set) ', 'NN (Poisson loss, DS1M) ','NN (Poisson loss, DS800k) ','NN (Poisson loss, DS600k) ','NN (Poisson loss, DS300k) ','NN (Poisson loss, DS100k) ','NN (Poisson loss, DS50k) ',
               'NN (MSE loss, Full set) ', 'NN (MSE loss, DS1M) ','NN (MSE loss, DS800k) ','NN (MSE loss, DS600k) ','NN (MSE loss, DS300k) ','NN (MSE loss, DS100k) ','NN (MSE loss, DS50k) '
               ]

# Metric to compare
performance_metric = 'corr'

# Initialize a DataFrame to store p-values with model names as both index and columns
p_value_matrix = pd.DataFrame(index=model_names, columns=model_names)

# Calculate the p-value for each pair of directory paths
for i, dir_path1 in enumerate(dir_paths):
    for j, dir_path2 in enumerate(dir_paths):
        if i < j:  # Only compute for pairs (avoid redundant comparisons)
            p_value = compare_two_models(dir_path1, dir_path2, performance_metric)
            p_value_matrix.loc[model_names[i], model_names[j]] = p_value
            p_value_matrix.loc[model_names[j], model_names[i]] = p_value  # Symmetric matrix

# Display or save the p-value matrix
print(p_value_matrix)
# Optionally save to CSV
p_value_matrix.to_csv("../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_model_DS.csv")


###############################################################################
# List of directory paths to compare
dir_paths = [
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/rep_train_test/',
    
    "../external/BMR/output/dimReduction_effect/PCA/GBM/rep_train_test/",
    "../external/BMR/output/dimReduction_effect/PCA/RF/rep_train_test/",
    "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/rep_train_test/",
    "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/rep_train_test/", 
   
   "../external/BMR/output/dimReduction_effect/AE/GBM/rep_train_test/",
   "../external/BMR/output/dimReduction_effect/AE/RF/rep_train_test/",
   "../external/BMR/output/dimReduction_effect/AE/nn_poisLoss/rep_train_test/",
   "../external/BMR/output/dimReduction_effect/AE/nn_mseLoss/rep_train_test/"
    
]

model_names = ['All features XGBoost', 'All features RF','All features NN (Poisson loss)','All features NN (MSE loss)',
               'PCA XGBoost', 'PCA RF','PCA NN (Poisson loss)','PCA NN (MSE loss)',
               'AE XGBoost', 'AE RF','AE NN (Poisson loss)','AE NN (MSE loss)'
               ]
# Metric to compare
performance_metric = 'corr'

# Initialize a DataFrame to store p-values with model names as both index and columns
p_value_matrix = pd.DataFrame(index=model_names, columns=model_names)

# Calculate the p-value for each pair of directory paths
for i, dir_path1 in enumerate(dir_paths):
    for j, dir_path2 in enumerate(dir_paths):
        if i < j:  # Only compute for pairs (avoid redundant comparisons)
            p_value = compare_two_models(dir_path1, dir_path2, performance_metric)
            p_value_matrix.loc[model_names[i], model_names[j]] = p_value
            p_value_matrix.loc[model_names[j], model_names[i]] = p_value  # Symmetric matrix

# Display or save the p-value matrix
print(p_value_matrix)
# Optionally save to CSV
p_value_matrix.to_csv("../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_model_dimRed.csv")

###############################################################################

dir_paths = ['../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_rep/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_rep/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100_rep/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100_rep/rep_train_test/']


model_names = ['#mutations >= 1',
               'All',
               '#mutations >= 1 & Length > 100',
               'Length > 100'
               ]
# Metric to compare
performance_metric = 'corr'

# Initialize a DataFrame to store p-values with model names as both index and columns
p_value_matrix = pd.DataFrame(index=model_names, columns=model_names)

# Calculate the p-value for each pair of directory paths
for i, dir_path1 in enumerate(dir_paths):
    for j, dir_path2 in enumerate(dir_paths):
        if i < j:  # Only compute for pairs (avoid redundant comparisons)
            p_value = compare_two_models(dir_path1, dir_path2, performance_metric)
            p_value_matrix.loc[model_names[i], model_names[j]] = p_value
            p_value_matrix.loc[model_names[j], model_names[i]] = p_value  # Symmetric matrix

# Display or save the p-value matrix
print(p_value_matrix)
# Optionally save to CSV
p_value_matrix.to_csv("../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_mutStatus.csv")


###############################################################################
import pandas as pd
from statsmodels.stats.multitest import multipletests


def convert_pvals_to_qvals(path_pval_matrix):
    # Load the CSV file into a DataFrame
    pval_df = pd.read_csv(path_pval_matrix, index_col=0)
    
    # Convert all values to numeric, coercing errors to NaN
    pval_df = pval_df.apply(pd.to_numeric, errors='coerce')
    
    # Flatten the DataFrame values to apply FDR adjustment
    pval_flat = pval_df.values.flatten()
    non_nan_mask = ~np.isnan(pval_flat)  # Mask of non-NaN values
    
    # Apply FDR correction only to non-NaN values
    _, qval_flat_non_nan, _, _ = multipletests(pval_flat[non_nan_mask], alpha=0.05, method='fdr_bh')
    
    # Create a new array to store q-values, inserting NaNs where necessary
    qval_flat = np.full(pval_flat.shape, np.nan)
    qval_flat[non_nan_mask] = qval_flat_non_nan
    
    # Reshape back to the original DataFrame format
    qval_df = pd.DataFrame(qval_flat.reshape(pval_df.shape), columns=pval_df.columns, index=pval_df.index)
    
    # Define the output path for q-values file
    output_path = path_pval_matrix.replace("p_value", "q_value")
    
    # Save the q-value DataFrame to a new CSV file
    qval_df.to_csv(output_path)
    
    return qval_df

# Example usage
path_pval_matrixs = ["../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_model_binGeneration.csv", 
                     '../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_model_DS.csv',
                     "../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_model_dimRed.csv",
                     '../external/BMR/output/Res_reviewerComments/significance_of_models/p_value_matrix_mutStatus.csv']
for path_pval_matrix in path_pval_matrixs:
    qval_df = convert_pvals_to_qvals(path_pval_matrix)



##############################################################################
##############################################################################
##############################################################################
def plot_p_value_heatmap(path_matrix, path_save):
    # Load the p-value matrix
    p_value_matrix = pd.read_csv(path_matrix, index_col=0)
    
    # Convert p-value entries to floats if necessary
    p_value_matrix = p_value_matrix.apply(pd.to_numeric, errors='coerce')
    
    
    # Create a mask for the upper triangle
    mask = np.triu(np.ones_like(p_value_matrix, dtype=bool))
    
    # Set up the matplotlib figure and axis explicitly
    fig, ax = plt.subplots(figsize=(14, 10))  # Adjust size as needed
    
    # Create the heatmap with the mask
    heatmap = sns.heatmap(
        p_value_matrix,
        mask=mask,  # Apply the mask to avoid redundant comparisons
        annot=True,  # Display p-values in each cell
        fmt=".2f",   # Show values with 3 decimal places
        cmap="coolwarm",  # Color map for visual effect
        cbar_kws={'label': 'q-value'},  # Label for color bar
        square=False,  # Rectangular shape
        xticklabels=True,
        yticklabels=True,
        annot_kws={"size": 10},  # Font size for cell annotations
        vmin=0, vmax=1,
        ax=ax  # Pass the axis to sns.heatmap
    )
    
    # Increase color bar label font size
    colorbar = heatmap.collections[0].colorbar
    colorbar.ax.yaxis.label.set_size(18)  # Set the color bar label font size
    
    # Enhance the font size of model names (tick labels)
    plt.xticks(rotation=45, ha='right', fontsize=17)
    plt.yticks(rotation=0, fontsize=18)
    
    # Remove x and y axis labels
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    # Title
    plt.title("")
    
    # Adjust layout for readability and save
    plt.tight_layout()
    plt.savefig(path_save, dpi=300)  # Save with high resolution
    plt.close()  # Close the plot to free memory




################################################################################


# Example usage
plot_p_value_heatmap(
    path_matrix="../external/BMR/output/Res_reviewerComments/significance_of_models/q_value_matrix_model_binGeneration.csv",
    path_save="../external/BMR/output/Res_reviewerComments/significance_of_models/q_value_matrix_model_binGeneration.png"
)



plot_p_value_heatmap(
    path_matrix="../external/BMR/output/Res_reviewerComments/significance_of_models/q_value_matrix_model_DS.csv",
    path_save="../external/BMR/output/Res_reviewerComments/significance_of_models/q_value_matrix_model_DS.png"
)

plot_p_value_heatmap(
    path_matrix="../external/BMR/output/Res_reviewerComments/significance_of_models/q_value_matrix_model_dimRed.csv",
    path_save="../external/BMR/output/Res_reviewerComments/significance_of_models/q_value_matrix_model_dimRed.png"
)
