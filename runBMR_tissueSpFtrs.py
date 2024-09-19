import ast
import numpy as np
import pandas as pd
from simulation_settings import load_sim_settings
from readFtrs_Rspns import create_TestTrain_TwoSources
from models.runBMR_functions import  config_save, RUN_BMR
from performance.assessModels import assess_models

# # Load the Excel files
# urls_df = pd.read_excel('../external/tmp/all_feature_URLs.xlsx')
# tissues_df = pd.read_excel('../external/tmp/ftrs_histology.xlsx')

# # Initialize a column for the matched pattern/cell
# urls_df['matched_pattern'] = ''

# # Perform a substring search for each pattern in 'pattern/cell' within 'Feature Name'
# for index, row in tissues_df.iterrows():
#     pattern = row['pattern/cell']
#     # Use str.contains to find rows where the pattern is in 'Feature Name'
#     urls_df.loc[urls_df['Feature Name'].str.contains(pattern, case=False, na=False), 'matched_pattern'] = pattern

# # Perform the left join using the 'matched_pattern' column
# merged_df = pd.merge(urls_df, tissues_df, left_on='matched_pattern', right_on='pattern/cell', how='left')

# merged_df = merged_df[['Feature Name', 'URL', 'Group Name', 'md5sum', 
#        'pattern/cell', 'origin', 'description ']]

# # Save or display the merged result
# merged_df.to_csv('../external/tmp/merged_output.csv', index=False)
# print(merged_df)


def load_data_tissueSpFtrs(sim_setting, cohort, 
                           path_tissueSpFtrs = '../external/tmp/merged_output.csv'):
    
    all_ftrs_df = pd.read_csv(path_tissueSpFtrs)
    
    tissue_specific_features = {
        'Skin-Melanoma': ['skin', 'melanoma'], #101
        # 'Lymph-BNHL': ['blood'],
        # 'Lymph-CLL': ['blood'],
        # 'Lymph-NOS': ['blood'],
        'Panc-Endocrine': ['pancreas'], #24
        'Panc-AdenoCA': ['pancreas'], #24
        # 'Myeloid-AML': ['blood'],
        # 'Myeloid-MPN': ['blood'],
        # 'Myeloid-MDS': ['blood'],
        'Liver-HCC': ['liver'], #33
        'Prost-AdenoCA': ['prostate'], # 4
        'Bone-Epith': ['bone'], #13
        'Bone-Osteosarc': ['bone'], #13
        'Bone-Leiomyo': ['bone'], #13
        'Bone-Cart': ['bone'], #13
        'Lung-SCC': ['lung', 'lung carcinoma'], #83
        'Lung-AdenoCA': ['lung', 'lung carcinoma'], #83
        'Breast-AdenoCa': ['breast'], #37
        'Breast-LobularCa': ['breast'],
        'Breast-DCIS': ['breast'],
        'Stomach-AdenoCA': ['GI_STOMACH'], #35
        'ColoRect-AdenoCA': ['GI_COLON', 'GI_RECTUM'], #53
        # 'Bladder-TCC': [],
        'Head-SCC': ['brain', 'cerebellum', 'frontal cortex'], #127
        'CNS-Medullo': ['brain', 'cerebellum', 'frontal cortex'], #127
        'CNS-Oligo': ['brain', 'cerebellum', 'frontal cortex'], #127
        'CNS-PiloAstro':  ['brain', 'cerebellum', 'frontal cortex'], #127
        'CNS-GBM': ['brain', 'cerebellum', 'frontal cortex'], #127
        'Kidney-RCC': ['kidney'], #13
        'Kidney-ChRCC': ['kidney'], #13
        'Cervix-AdenoCA': ['uterus', 'cervix'], #20
        'Cervix-SCC': ['uterus', 'cervix'], #20
        'Uterus-AdenoCA': ['uterus', 'cervix'], #20
        # 'Biliary-AdenoCA': [],
        # 'Eso-AdenoCa': [],
        'Ovary-AdenoCA': ['ovary'],
        'Thy-AdenoCA': ['thymus'], #17
        
    }
    
    x = tissue_specific_features[cohort]
    
    filtered_df = all_ftrs_df[all_ftrs_df['origin'].isin(x)]
    tissue_sp_ftrs = (filtered_df[['Feature Name']])
    
    genomeFtrs = ['ACA', 'ACC', 'ACG', 'ACT', 'ATA', 'ATC', 'ATG', 'ATT',
                  'CCA', 'CCC', 'CCG', 'CCT', 'CTA', 'CTC', 'CTG', 'CTT', 
                  'GCA', 'GCC', 'GCG', 'GCT', 'GTA', 'GTC', 'GTG', 'GTT', 
                   'TCA', 'TCC', 'TCG', 'TCT', 'TTA', 'TTC', 'TTG', 'TTT', 
                   'TA5p', 'TC5p', 'TG5p', 'TT5p', 'CA5p', 'CC5p', 'CG5p', 
                    'CT5p', 'AT3p', 'CT3p', 'GT3p', 'AC3p', 'GC3p', 'TC3p',
                    'vertebrate_phastCons46way',
                     'primates_phastCons46way',
                     'primates_phyloP46way']
    
    # Convert the 'Feature Name' column to a list
    tissue_sp_ftrs_list = tissue_sp_ftrs['Feature Name'].tolist()
    
    # Combine the list from the DataFrame with nucleotide_content
    combined_list = tissue_sp_ftrs_list + genomeFtrs
    
    tissue_sp_ftrs = combined_list
    
    path_X_test = sim_setting['path_X_test']
    path_X_train = sim_setting['path_X_train']
    path_Y_test = sim_setting['path_Y_test']
    path_Y_train = sim_setting['path_Y_train']
    scale = ast.literal_eval(sim_setting['scale'])
    DSmpl = ast.literal_eval(sim_setting['DSmpl'])
    n_sample = sim_setting['n_sample']
    remove_unMutated = ast.literal_eval(sim_setting['remove_unMutated'])
    
    X_train, Y_train, X_test, Y_test = create_TestTrain_TwoSources(path_X_train, 
                                                                   path_Y_train, 
                                                                   path_X_test, 
                                                                   path_Y_test,
                                                                   scale,
                                                                   use_features = tissue_sp_ftrs)
    X_train = X_train.loc[:, tissue_sp_ftrs]
    X_test = X_test.loc[:, tissue_sp_ftrs]
    
    if remove_unMutated:
        Y_train = Y_train[Y_train['nMut'] != 0]
        X_train = X_train.loc[Y_train.index]
        
        Y_test = Y_test[Y_test['nMut'] != 0]
        X_test = X_test.loc[Y_test.index]
    
    if DSmpl:
        
        np.random.seed(40)
        tr_indices = np.random.choice(list(Y_train.index), size=n_sample, replace=False)
        Y_train = Y_train.loc[tr_indices]
        print(f'Down sampling was performed... number of training bins: {Y_train.shape[0]}')
        X_train = X_train.loc[Y_train.index]
    
    if (Y_test.index != X_test.index).all():
        raise ValueError('X_test and Y_test indexes are not the same')
    if (Y_train.index != X_train.index).all():
        raise ValueError('X_train and Y_train indexes are not the same')
        
    return X_train, Y_train, X_test, Y_test



# List of cohorts
# 'Pancan-no-skin-melanoma-lymph', 
cohorts = [ "Liver-HCC",  "ColoRect-AdenoCA", "Uterus-AdenoCA", "Kidney-RCC",
    # "Lymph-BNHL","Lymph-CLL", "Myeloid-MDS", "Lymph-NOS","Myeloid-AML", "Myeloid-MPN",
    "Lung-SCC", "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
    "Breast-AdenoCa", "Breast-DCIS", "Breast-LobularCa", 
    # "Biliary-AdenoCA", "Eso-AdenoCa", "Bladder-TCC",
    "CNS-GBM", "CNS-Medullo", "CNS-Oligo", "CNS-PiloAstro",   
    "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
    "Thy-AdenoCA",  "Bone-Leiomyo",    
     "Cervix-SCC", "Cervix-AdenoCA", 
     "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc", "Bone-Cart"
]

# Path to the template ini file
template_ini_path = 'configs/rate_based/sim_setting_tissueSpFtrs.ini'

# Read the template ini file
with open(template_ini_path, 'r') as file:
    ini_content = file.read()

# Loop over the cohorts and run the algorithm
for cohort in cohorts:
    print(cohort)
    # Replace the cohort name in the ini content
    modified_ini_content = ini_content.replace('Skin-Melanoma', cohort)
    
    # Save the modified ini content to a new file
    modified_ini_path = f'configs/rate_based/sim_setting_tissueSpFtrs_{cohort}.ini'
    with open(modified_ini_path, 'w') as file:
        file.write(modified_ini_content)
    
    # Run the algorithm with the modified ini file
    sim_file = modified_ini_path
    sim_setting = load_sim_settings(sim_file)
    config_save(sim_file)
    
    X_train, Y_train, X_test, Y_test = load_data_tissueSpFtrs(sim_setting, cohort)
    print(X_train.shape)
    print(X_test.shape)
    RUN_BMR(sim_setting, X_train, Y_train, X_test, Y_test, make_pred=True) 
    print('************') 
    assess_models(sim_setting)
    

