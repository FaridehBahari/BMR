import pandas as pd
import os
import re
import numpy as np
from performance.assessModels import assess_model

M_name = 'GBM0_rep'

# Directory containing the files
directory = f'../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/{M_name}/rep_train_test/'
y = pd.read_csv('../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv', sep='\t', header=0, index_col='binID',
                usecols=['binID', 'length', 'nMut', 'nSample', 'N'])

y = y[y['nMut'] != 0]
# List all files in the directory
files = os.listdir(directory)

# Filter files matching the pattern pred_test{some number}.tsv
pattern = re.compile(fr'{M_name}_predTest\d+\.tsv')
matching_files = [f for f in files if pattern.match(f)]

print("Matching files:", matching_files)



# Function to read a .tsv file
def read_tsv(file_path):
    return pd.read_csv(file_path, sep='\t', index_col='binID')

# Read all matching files into a list of DataFrames
dataframes = [read_tsv(os.path.join(directory, f)) for f in matching_files]


# Merge each df with y to add the length column
def add_length_column(df, y):
    return df.merge(y[['length']], left_index=True, right_index=True, how='left')

# Split each DataFrame based on the length column into specified bins
def split_by_length(df):
    bins = [0, 20, 50, 100 , 200, 400, 1500, float('inf')]
    labels = ['0-20', '20-50', '50-100', '100-200', '200-400', '400-1500', '1500<']
    df['length_bin'] = pd.cut(df['length'], bins=bins, labels=labels, right=False)
    return {label: df[df['length_bin'] == label] for label in labels}

# Function to run assess_model on each subset
def run_assess_model(df, y, model_name):
    results = {}
    subsets = split_by_length(df)
    for label, subset in subsets.items():
        if subset.empty:
            results[label] = np.nan  # Assign "NA" when the subset is empty
        else:
            Y_pred = subset['predRate']
            Y_obs = y.loc[Y_pred.index] 
            Y_obs = Y_obs['nMut'] / (Y_obs['length'] * Y_obs['N'])
            Nr_pair_acc = 100000
            assessment_result = assess_model(Y_pred, Y_obs, Nr_pair_acc, model_name, per_element=False)
            # Fill all values in the DataFrame with "NA" if the DataFrame is empty
            if assessment_result.empty:
                assessment_result = assessment_result.fillna(np.nan)
            results[label] = assessment_result
    return results






# Apply the merging and assessment to each DataFrame
merged_dataframes = [add_length_column(df, y) for df in dataframes]
all_results = []

for i, df in enumerate(merged_dataframes):
    model_name = f"{M_name}_{i+1}"
    results = run_assess_model(df, y, model_name)
    all_results.append(results)

# Print the results
for i, result in enumerate(all_results):
    print(f"Results for {M_name}_{i+1}:")
    for length_bin, metrics in result.items():
        print(f"  Length bin {length_bin}: {metrics}")
        

# Extract the corr values and calculate the mean corr for each length interval
length_bins =  ['0-20', '20-50', '50-100', '100-200', '200-400', '400-1500', '1500<']
mean_corrs = {bin: [] for bin in length_bins}

for i, result in enumerate(all_results):
    n = str(i+1)
    for length_bin, metrics in result.items():
        # Check if metrics is a float (indicating "NA" or missing data)
        if isinstance(metrics, float):
            # Append "NA" or missing data to the list
            mean_corrs[length_bin].append(np.nan)
        else:
            # Append the value to the list after accessing it with .loc
            mean_corrs[length_bin].append(metrics.loc['corr_' + M_name + '_' + n])


# Calculate the mean corr for each length bin
mean_corrs = {bin: sum(corrs) / len(corrs) for bin, corrs in mean_corrs.items()}

mean_corrs_extracted = {}
for bin, corr in mean_corrs.items():
    if isinstance(corr, float):
        mean_corrs_extracted[bin] = corr
    else:
        mean_corrs_extracted[bin] = corr.values[0]

print(mean_corrs_extracted)


# Convert the extracted dictionary to a DataFrame with length bins as the index
mean_corrs_df = pd.DataFrame(mean_corrs_extracted, index=['corr']).transpose()
np.mean(mean_corrs_df['corr'])


mean_corrs_df.to_csv(f'{directory}_length_intervals_assessment.tsv')


pd.read_csv('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_rep/rep_train_test/_length_intervals_assessment.tsv')
pd.read_csv('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100_rep/rep_train_test/_length_intervals_assessment.tsv')
pd.read_csv('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_rep/rep_train_test/_length_intervals_assessment.tsv')
pd.read_csv('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100_rep/rep_train_test/_length_intervals_assessment.tsv')