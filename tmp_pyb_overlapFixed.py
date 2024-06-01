import pandas as pd
import numpy as np
import pickle
from pybedtools import BedTool
import pybedtools
from readFtrs_Rspns import  load_data, read_response

def val_IDs_fixedElems(bed_tr, bed_val, seed_value, val_size):
    
    
    
    bed_val = bed_val.iloc[np.where(bed_val[2] - bed_val[1] >= 20)]
    
    np.random.seed(seed_value)
    tr_indices = np.random.choice(bed_tr.index, size=val_size, replace=False)
    bed_tr = bed_tr.loc[tr_indices]
    
    bedObj_val = BedTool.from_dataframe(bed_val)
    
    
    bedObj_tr = BedTool.from_dataframe(bed_tr)
    
    intersection_tr= bedObj_tr.intersect(bedObj_val).to_dataframe()
    non_train_set_binIDs = np.unique(intersection_tr.name)
    
    intersection_val = bedObj_val.intersect(bedObj_tr).to_dataframe()
    val_set_binIDs = np.unique(intersection_val.name)
    
    
    return val_set_binIDs, non_train_set_binIDs


# Function to process all bed_tr files and save results
def process_and_save_results(bed_val, bed_tr, seed_values, val_size, output_file,
                             bin_size,  n_repeat = 10):
    
    results = {}
    for i in range(n_repeat):
        print(i)
        val_set_binIDs, non_train_set_binIDs = val_IDs_fixedElems(bed_tr, bed_val, seed_values[i], val_size)
        results[f'{bin_size}_{i}'] = {
            'val_set_binIDs': val_set_binIDs,
            'non_train_set_binIDs': non_train_set_binIDs
        }
    
    with open(output_file, 'wb') as f:
        pickle.dump(results, f)


def select_groups_from_dict(dictionary, keys_to_include):
    
    # Create an empty list to store the values
    included_values = []
    
    # Iterate through the original dictionary
    for key, value in dictionary.items():
        # Check if the key should be included
        if key in keys_to_include:
            # Extend the list with the values
            included_values.extend(value)
            
    return included_values


def get_features_category(category, path_featureURLs = '../external/database/all_feature_URLs.xlsx'):
    
    # # Load feature groups from Excel file
    # feature_groups_df = pd.read_excel(path_featureURLs)  
    # feature_groups = feature_groups_df.groupby('Group Name')['Feature Name'].apply(list).to_dict()
    # nucleotide_content = ['ACA', 'ACC', 'ACG', 'ACT', 'ATA', 'ATC', 'ATG', 'ATT',
    #                       'CCA', 'CCC', 'CCG', 'CCT', 'CTA', 'CTC', 'CTG', 'CTT', 
    #                       'GCA', 'GCC', 'GCG', 'GCT', 'GTA', 'GTC', 'GTG', 'GTT', 
    #                       'TCA', 'TCC', 'TCG', 'TCT', 'TTA', 'TTC', 'TTG', 'TTT', 
    #                       'TA5p', 'TC5p', 'TG5p', 'TT5p', 'CA5p', 'CC5p', 'CG5p', 
    #                       'CT5p', 'AT3p', 'CT3p', 'GT3p', 'AC3p', 'GC3p', 'TC3p']
    
    
    # # Adding 'nucleotide content'and 'APOBEC' key to the dictionary
    # feature_groups['nucleotide content'] = nucleotide_content
    # feature_groups['APOBEC'] = ['APOBEC3A']
    
    
    
    # # # File path from where to load the dictionary
    # # file_path = '../external/procInput/ftrs_dict.pickle'
    # # # Save the dictionary to disk
    # # with open(file_path, 'wb') as file:
    # #     pickle.dump(feature_groups, file)
    
    
    # Load the dictionary from disk
    with open('../external/procInput/ftrs_dict.pickle', 'rb') as file:
        feature_groups = pickle.load(file)
    
    features = select_groups_from_dict(feature_groups, category)
    
    return features




#####################################################################################

bin_sizes = ['1M', '100k', '50k', '10k']
for bin_size in bin_sizes:
    
    # Parameters
    seed_values = [1, 5, 14, 10, 20, 30, 40, 50, 60, 70, 80, 90, 77, 100, 110]
    
    path_X_train = f'../external/ftrMtrix/{bin_size}_features.h5'
    path_Y_train = f'../external/BMR/rawInput/responseTabs/Pan_Cancer/{bin_size}_bins.tsv'
    
    remove_unMutated = True
    
    ftrs = get_features_category(['HiC'])
    X_tr_cmplt, Y_tr_cmplt = load_data(path_X_train, path_Y_train,
                                           use_features=ftrs)
    
    if remove_unMutated:
        Y_tr_cmplt = Y_tr_cmplt[Y_tr_cmplt['nMut'] != 0]
        X_tr_cmplt = X_tr_cmplt.loc[Y_tr_cmplt.index]
        
    # Load your data (replace with your actual data loading method)
    bed_val = pd.read_csv('../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6', sep = '\t', header = None)
    bed_tr = pd.read_csv(f'../external/database/bins/proccessed/intergenic_fixed{bin_size}.bed6', sep = '\t', header = None)
    
    bed_tr['binID'] = bed_tr[3]
    bed_tr = bed_tr.set_index('binID')
    bed_tr = bed_tr.loc[Y_tr_cmplt.index]
    
    bed_val['binID'] = bed_val[3]
    bed_val = bed_val.set_index('binID')
    Y_val_cmplt = read_response('../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv')
    Y_val_cmplt = Y_val_cmplt.iloc[np.where(Y_val_cmplt.nMut !=0)]
    bed_val = bed_val.loc[Y_val_cmplt.index]
    
    print(bed_val.shape)
    print(bed_tr.shape)
    
    val_size = X_tr_cmplt.shape[0]//5 
    print(val_size)
    base_dir = '../external/BMR/procInput/fixedSize_trainValIDs/'
    
    output_file = f'{base_dir}fixed{bin_size}_IDs.pkl'
    
    # Process and save results
    process_and_save_results(bed_val, bed_tr, seed_values, val_size, output_file, bin_size)



#########################################################################

def extract_bin_size(path_bed_tr):
    
    if path_bed_tr == '../external/database/bins/proccessed/intergenic_fixed1M.bed6':
        bin_size = '1M'
    elif path_bed_tr == '../external/database/bins/proccessed/intergenic_fixed100k.bed6':
        bin_size = '100k'
    elif path_bed_tr == '../external/database/bins/proccessed/intergenic_fixed50k.bed6':
        bin_size = '50k'
    elif path_bed_tr == '../external/database/bins/proccessed/intergenic_fixed10k.bed6':
        bin_size = '10k'
        
    return bin_size


def sample_train_val_fixedSize2(Y_train, Y_val, bed_tr, bed_var, 
                               path_bed_tr, iteration):
    
    bin_size = extract_bin_size(path_bed_tr)
    saved_file = '../external/BMR/procInput/fixedSize_trainValIDs/fixed{bin_size}_IDs.pkl'
    pickle_file_path = saved_file
    with open(pickle_file_path, 'rb') as f:
         results = pickle.load(f)
    
    val_set_binIDs = results[f'{bin_size}_0']['val_set_binIDs']
    non_train_set_binIDs = results[f'{bin_size}_0']['non_train_set_binIDs']
    
    Y_train = Y_train[~Y_train.index.isin(non_train_set_binIDs) ]
    Y_val = Y_val.loc[val_set_binIDs]
    
    return Y_train, Y_val



def repeated_train_test(sim_setting,  X_tr_cmplt, Y_tr_cmplt, X_val_cmplt, Y_val_cmplt,
            make_pred = True, overwrite = True, n_repeat = 10):
    
    
    path_bed_tr = sim_setting['path_bed_tr']
    path_bed_var = sim_setting['path_bed_var']
    
    models = sim_setting['models']
    base_dir = sim_setting['base_dir']
    
    Nr_pair_acc = sim_setting['Nr_pair_acc']
    
    bed_tr = pd.read_csv(path_bed_tr, sep = '\t', header = None)
    bed_tr['binID'] = bed_tr[3]
    bed_tr = bed_tr.set_index('binID')
    bed_tr = bed_tr.loc[Y_tr_cmplt.index]
    
    bed_val = pd.read_csv(path_bed_var, sep = '\t', header = None)
    bed_val['binID'] = bed_val[3]
    bed_val = bed_val.set_index('binID')
    bed_val = bed_val.loc[Y_val_cmplt.index]
    
    for key in models:
       
        m = models[key]
        name = m['save_name']
        
        os.makedirs(f'{base_dir}/{name}/', exist_ok= True)
        readme_file_name = f'{base_dir}/{name}/README.md'
        print(f'@@@@  model: {name}  @@@@')
        params = m['Args']
        save_path_model = f'{base_dir}/{name}/'
        
        # check_file_func = m['check_file_func']
        # file_check = check_file_func(base_dir, name)
        # if not os.path.exists(file_check) or sim_setting['overwrite']:
        if not os.path.exists(readme_file_name) or overwrite:
            write_readme_file(m, readme_file_name)
            
            for i in range(n_repeat):
                params['path_save'] = f'{save_path_model}rep_train_test/models{i+1}/'
                print(f'.......... repeat number {i+1} of train-test for evaluation of the {name} ......')
                                
                if os.path.exists(f'{save_path_model}/rep_train_test/{name}_M{i+1}_assessment.tsv'):
                    print(f"Skipping iteration {i+1} as the file already exists.")
                    continue
                
                Y_train, Y_test = sample_train_val_fixedSize2(Y_tr_cmplt, Y_val_cmplt,
                                                              bed_tr, bed_val,
                                                              path_bed_tr, i)
                
                X_train = X_tr_cmplt.loc[Y_train.index]
                X_test = X_val_cmplt.loc[Y_test.index]
                
                print(X_train.shape)
                print(X_test.shape)
                
                common_indices = X_test.index.intersection(X_train.index)
                
                if not common_indices.empty:
                    raise ValueError(f"Common indices found between X_test and X_train:{common_indices}")
                else:
                    print("No common indices found between X_test and X_train.")
                    
                fitted_Model = fit_model(X_train, Y_train, X_test, Y_test,
                                         m['run_func'], m['predict_func'], make_pred, m['Args'])
                save_func = m['save_func']
                itr = i+1
                                
                save_func(fitted_Model, base_dir, name, iteration=itr, save_model=True)
                
                Y_pred = fitted_Model.predRates_test
                Y_obs = Y_test.nMut/(Y_test.N * Y_test.length)
                assessments = assess_model(Y_pred, Y_obs, Nr_pair_acc, name, per_element=False)
                
                path_assessments = f'{save_path_model}/rep_train_test/{name}_M{i+1}_assessment.tsv'
                assessments.to_csv(path_assessments, sep='\t')
        
        print("=============================")
        dir_path = f'{save_path_model}/rep_train_test/'
        save_metrics_summary(dir_path)

