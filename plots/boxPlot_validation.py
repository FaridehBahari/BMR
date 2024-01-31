import os
import pandas as pd
import matplotlib.pyplot as plt

def extract_bin_size_from_path(directory_path):
    # Extract bin size from directory path
    parts = directory_path.split('/')
    bin_size_index = parts.index('dimReduction_effect') + 2
    return parts[bin_size_index]

def plot_metric_boxplot(directory_paths, metric):
    all_data = []

    for directory_path in directory_paths:
        data_frames = []

        # Read data from TSV files in the current directory
        files = os.listdir(directory_path)
        for file in files:
            if file.endswith("_assessment.tsv"):
                metric_data = pd.read_csv(os.path.join(directory_path, file),
                                          sep='\t', index_col=0, header=None)
                metric_data = metric_data.T

                metric_name = [s for s in metric_data.columns if metric in str(s)]
                mDat = pd.DataFrame(metric_data[metric_name].values.astype(float))
                data_frames.append(mDat)

        if data_frames:
            # Combine data for the current directory
            combined_data = pd.concat(data_frames, axis=0)

            # Reset the index before concatenating
            combined_data.reset_index(drop=True, inplace=True)

            # Extract bin size from the current directory path
            bin_size = extract_bin_size_from_path(directory_path)

            # Add combined data and bin size to the list
            all_data.append((combined_data, bin_size))

    if all_data:
        # Combine data for plotting
        combined_data = pd.concat([df[0] for df in all_data], axis=1)
        combined_data.columns = [df[1] for df in all_data]

        # Plotting the box plot for all values
        plt.figure(figsize=(12, 8))
        combined_data.boxplot(grid=False)
        plt.ylabel(metric)
        plt.title(f'{metric} Distribution Comparison between Models')
        plt.show()
    else:
        print(f"No data found for metric '{metric}' in the specified directories.")

# Example usage with multiple directory paths
directory_paths_PCA = [
    "../external/BMR/output/dimReduction_effect/PCA/GBM/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/RF/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/NN_PoisLoss/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/NN_MSEloss/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/GLM/rep_train_test/"
    
]
metric = 'corr'

plot_metric_boxplot(directory_paths_PCA, metric)



directory_paths_AE = ["../external/BMR/output/dimReduction_effect/AE/GBM/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/AE/RF/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/AE/NN_PoisLoss/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/AE/NN_MSEloss/rep_train_test/",
                         "../external/BMR/output/dimReduction_effect/AE/GLM/rep_train_test/"]





####################################################
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# First, let's enhance extract_bin_size_from_path to extract both bin size and method
def extract_info_from_path(directory_path):
    parts = directory_path.split('/')
    method = parts[parts.index('dimReduction_effect') + 1]
    model = parts[parts.index('dimReduction_effect') + 2]
    return method, model

# Next, we modify plot_metric_boxplot to handle grouping
def plot_metric_boxplot(directory_paths, metric):
    all_data = []

    for directory_path in directory_paths:
        files = os.listdir(directory_path)
        for file in files:
            if file.endswith("_assessment.tsv"):
                metric_data = pd.read_csv(os.path.join(directory_path, file),
                                          sep='\t', index_col=0, header=None).T
                
                method, model = extract_info_from_path(directory_path)
                
                metric_name = [s for s in metric_data.columns if metric in str(s)]
                metric_data['Method'] = method  # Add a column for the method (PCA or AE)
                metric_data['Model'] = model  # Add a column for the model
                metric_data['Value'] = metric_data[metric_name].values.astype(float) # Consolidate metric under 'Value'
            
                all_data.append(metric_data[['Model', 'Method', 'Value']])
    
    if all_data:
        # Concatenate all data frames
        combined_data = pd.concat(all_data, ignore_index=True)
        
        # Use seaborn to plot the grouped box plot
        plt.figure(figsize=(12, 8))
        sns.boxplot(data=combined_data, x='Model', y='Value', hue='Method')
        plt.ylabel(metric)
        plt.title(f'{metric} Distribution Comparison between Models and Methods')
        plt.show()
    else:
        print(f"No data found for metric '{metric}' in the specified directories.")



directory_paths_PCA = [
    "../external/BMR/output/dimReduction_effect/PCA/GBM/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/RF/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/NN_PoisLoss/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/NN_MSEloss/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/PCA/GLM/rep_train_test/"
    
]
metric = 'corr'




directory_paths_AE = ["../external/BMR/output/dimReduction_effect/AE/GBM/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/AE/RF/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/AE/NN_PoisLoss/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/AE/NN_MSEloss/rep_train_test/",
                         "../external/BMR/output/dimReduction_effect/AE/GLM/rep_train_test/"]


directory_paths_orig = ["../external/BMR/output/dimReduction_effect/original/GBM/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/original/RF/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/original/NN_PoisLoss/rep_train_test/",
                          "../external/BMR/output/dimReduction_effect/original/NN_MSEloss/rep_train_test/"]


plot_metric_boxplot(directory_paths_PCA + directory_paths_AE + directory_paths_orig, metric)
