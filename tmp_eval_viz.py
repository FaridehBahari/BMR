import os
import pandas as pd
from scipy.stats import ttest_rel

def extract_bin_size_from_path(directory_path):
    # Extract bin size from directory path
    parts = directory_path.split('/')
    bin_size_index = parts.index('bin_size_effect') + 1
    return parts[bin_size_index]

def perform_paired_t_test(directory_path1, directory_path2, metric):
    data1 = pd.DataFrame()
    data2 = pd.DataFrame()
    
    # Read data from TSV files in directory 1
    files1 = os.listdir(directory_path1)
    for file in files1:
        if file.endswith("_assessment.tsv"):
            metric_data = pd.read_csv(os.path.join(directory_path1, file),
                                      sep='\t', index_col=0, header=None)
            metric_data = metric_data.T
            
            metric_name = [s for s in metric_data.columns if metric in str(s)]
            mDat = pd.DataFrame(metric_data[metric_name].values.astype(float))
            data1 = pd.concat([data1, mDat], axis=0)
            
    # Read data from TSV files in directory 2
    files2 = os.listdir(directory_path2)
    for file in files2:
        if file.endswith("_assessment.tsv"):
            metric_data = pd.read_csv(os.path.join(directory_path2, file), sep='\t', 
                                      index_col=0, header=None)
            metric_data = metric_data.T
            metric_name = [s for s in metric_data.columns if metric in str(s)]
            mDat = pd.DataFrame(metric_data[metric_name].values.astype(float))
            data2 = pd.concat([data2, mDat], axis=0)
            
    # Perform paired t-test
    t_statistic, p_value = ttest_rel(data1[0], data2[0])

    # Extract bin size from directory paths
    bin_size_folder1 = extract_bin_size_from_path(directory_path1)
    bin_size_folder2 = extract_bin_size_from_path(directory_path2)
    
    mean1 = data1.mean().values[0]
    mean2 = data2.mean().values[0]
    
    if metric != 'mse':
        better_folder = bin_size_folder1 if mean1 > mean2 else bin_size_folder2
    else:
        better_folder = bin_size_folder1 if mean1 < mean2 else bin_size_folder2

    # Additional information
    result = {
        't_statistic': t_statistic,
        'p_value': p_value,
        f'mean_{bin_size_folder1}': mean1,
        f'mean_{bin_size_folder2}': mean2,
        f'var_{bin_size_folder1}': data1.var().values[0],
        f'var_{bin_size_folder2}': data2.var().values[0],
        'better': better_folder
    }
    
    return result

# Example usage:
directory_path1 = '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/'
directory_path2 = '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/'
metric = 'mse'

result = perform_paired_t_test(directory_path1, directory_path2, metric)
print(result)
########################################################################
import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_metric_boxplot(directory_path1, directory_path2, metric):
    data1 = pd.DataFrame()
    data2 = pd.DataFrame()

    # Read data from TSV files in directory 1
    files1 = os.listdir(directory_path1)
    for file in files1:
        if file.endswith("_assessment.tsv"):
            metric_data = pd.read_csv(os.path.join(directory_path1, file),
                                      sep='\t', index_col=0, header=None)
            metric_data = metric_data.T

            metric_name = [s for s in metric_data.columns if metric in str(s)]
            mDat = pd.DataFrame(metric_data[metric_name].values.astype(float))
            data1 = pd.concat([data1, mDat], axis=0)

    # Read data from TSV files in directory 2
    files2 = os.listdir(directory_path2)
    for file in files2:
        if file.endswith("_assessment.tsv"):
            metric_data = pd.read_csv(os.path.join(directory_path2, file), sep='\t',
                                      index_col=0, header=None)
            metric_data = metric_data.T
            metric_name = [s for s in metric_data.columns if metric in str(s)]
            mDat = pd.DataFrame(metric_data[metric_name].values.astype(float))
            data2 = pd.concat([data2, mDat], axis=0)

    bin_size_folder1 = extract_bin_size_from_path(directory_path1)
    bin_size_folder2 = extract_bin_size_from_path(directory_path2)
     
    # Combine data for plotting
    combined_data = pd.concat([data1, data2], axis=1)
    combined_data.columns = [bin_size_folder1, bin_size_folder2]

    # Plotting the box plot for all values
    plt.figure(figsize=(10, 6))
    combined_data.boxplot(grid=False)
    plt.ylabel(metric)
    plt.title(f'{metric} Distribution Comparison between models')
    plt.show()

# Example usage:
directory_path1 = '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/'
directory_path2 = '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/'
metric = 'corr'

plot_metric_boxplot(directory_path1, directory_path2, metric)

###################################################################
import os
import pandas as pd
import matplotlib.pyplot as plt

def extract_bin_size_from_path(directory_path):
    # Extract bin size from directory path
    parts = directory_path.split('/')
    bin_size_index = parts.index('bin_size_effect') + 1
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
directory_paths = [
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/',
    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/'
    
]
metric = 'corr'

plot_metric_boxplot(directory_paths, metric)
