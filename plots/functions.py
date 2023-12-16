import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns


def plot_perBatch_assessments(file_paths, path_to_save, measurement, horizontal_line_GBM):
    os.makedirs(path_to_save, exist_ok= True)
    all_data = pd.DataFrame()
    directory_names = [os.path.basename(os.path.dirname(file_path)) for file_path in file_paths]
    
    for file_path in file_paths:
        df = pd.read_csv(file_path, delimiter='\t', index_col=0)
        intergenic_df = df['intergenic']
        intergenic_df = intergenic_df.loc[intergenic_df.index.str.contains(f'{measurement}')]
        all_data = pd.concat([all_data, intergenic_df], axis=1)
        
    numeric_index = all_data.index.str.extract(r'(\d+)', expand=False)
    all_data.columns = directory_names
    
    fig, ax = plt.subplots()
    
    for column in all_data.columns:
        ax.plot(numeric_index, all_data[column], label=column)
        
    # Get the path of the GBM file
    tsv_path = horizontal_line_GBM
    tsv_directory = os.path.dirname(tsv_path)
    
    # Get the directory name for the legend label
    legend_label = os.path.basename(tsv_directory)
    
    # Read the TSV file to extract the value for the horizontal line
    tsv_df = pd.read_csv(tsv_path, delimiter='\t', index_col=0)
    horizontal_line_value = tsv_df.loc[tsv_df.index.str.contains(f'{measurement}'), 'intergenic'].values[0]
    
    # Add the horizontal line with the specified value and legend label
    ax.axhline(y=horizontal_line_value, color='red', linestyle='--', label=legend_label)
    
    ax.set_xlabel('Batch Number')
    ax.set_ylabel(f'{measurement} Value')
    if measurement == 'mse':
        ax.set_yscale('log')
    ax.legend()
    plt.xticks(rotation='vertical', ha='right')
    
    fig_size = plt.gcf().get_size_inches() 
    sizefactor = 2.8 # Set a zoom factor
    # Modify the current size by the factor
    plt.gcf().set_size_inches(sizefactor * fig_size) 
    plt.tight_layout()
    plt.savefig(f'{path_to_save}/{legend_label}_{measurement}')
    plt.close()
# file_paths = ['../external/output/cd8_N2/perBatch_assessments.tsv',
#               '../external/output/N2_60k/perBatch_assessments.tsv', 
#               '../external/output/N3_60k/perBatch_assessments.tsv',
#               '../external/output/N4_60k/perBatch_assessments.tsv',
#               '../external/output/N5_60k/perBatch_assessments.tsv'] 

file_paths = ['../external/output/check_classic_NNs/savingIntervals/nn_classic_small_100Epochs_correctRates_learningRate0.00001/perBatch_assessments.tsv']
horizontal_line_GBM = '../external/output/GBM/GBM_assessments.tsv'


path_to_save = '../external/plots/perbatch/'
measurement = 'mse'
plot_perBatch_assessments(file_paths, path_to_save, measurement, horizontal_line_GBM)
