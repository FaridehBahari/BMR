import pickle
import xgboost as xgb
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests

path_tissueSpFtrs = '../external/tmp/merged_output.csv'

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
cohorts = ['Skin-Melanoma', 'Panc-Endocrine', 'Panc-AdenoCA', 'Liver-HCC', 
           'Prost-AdenoCA', 'Bone-Osteosarc', 'Bone-Leiomyo',
           'Lung-SCC', 'Lung-AdenoCA', 'Breast-AdenoCa',
           'Stomach-AdenoCA', 'ColoRect-AdenoCA', 'Head-SCC',
           'CNS-Medullo','CNS-GBM', 'Kidney-RCC',
           'Uterus-AdenoCA',
           'Ovary-AdenoCA'] # 'Cervix-AdenoCA', 'Cervix-SCC','CNS-Oligo','CNS-PiloAstro','Kidney-ChRCC','Thy-AdenoCA'

df = pd.DataFrame()
for cohort in cohorts:
    # Load the GBM model
    path = f'../external/BMR/output/reviewerComments/{cohort}/GBM/GBM_model.pkl'
    with open(path, 'rb') as file:
        model = pickle.load(file)
        
    # After training the model (assuming your model is named 'model')
    importance = model.get_score(importance_type='weight')  # Can also use 'gain' or 'cover'
    
    # Sort the importance dictionary by value (importance score) in descending order
    sorted_importance = {k: v for k, v in sorted(importance.items(), key=lambda item: item[1], reverse=True)}
    
    # Print the top 50 features
    top_50 = dict(list(sorted_importance.items())[:50])
            
    x = tissue_specific_features[cohort]
    
    filtered_df = all_ftrs_df[all_ftrs_df['origin'].isin(x)]
    tissue_sp_ftrs = (filtered_df[['Feature Name']])
    
    # Load tissue-specific features into a list
    tissue_sp_ftrs_list = tissue_sp_ftrs['Feature Name'].tolist()
    
    # Initialize the list to store the data for each cohort
    data = []
    
    # Collect all top 50 features as a comma-separated string
    top_50_str = ', '.join(top_50.keys())
    
    # Collect tissue-specific features along with their rank
    tissue_sp_ranks = []
    for i, feature in enumerate(top_50.keys(), 1):  # Enumerate to get the rank
        if feature in tissue_sp_ftrs_list:
            tissue_sp_ranks.append(f"{feature}: {i}")
            
    # Combine tissue-specific features and their ranks as a comma-separated string
    tissue_sp_ranks_str = ', '.join(tissue_sp_ranks)
    
    k = len(tissue_sp_ranks) # Number of successes observed in top 50
    M = 1372 # total number of features
    n = len(tissue_sp_ftrs) # number of tissue-specific features in allof the feature set
    N = 50 # top 50 important features
    
    # Append the data for the current cohort
    data.append({
        'cohort': cohort,
        'Top 50 important features': top_50_str,
        'rank of the tissue_sp_ftr in the top 50 important feature': tissue_sp_ranks_str,
        'total number of tissue-specific features': len(tissue_sp_ftrs),
        'enrichment score': hypergeom.sf(k-1, M, n, N),
        'fdr': multipletests(hypergeom.sf(k-1, M, n, N), method='fdr_bh')[1][0]
    })
    
    df = pd.concat([df, pd.DataFrame(data)], axis=0)

# Display the DataFrame
# print(df)
df.shape

import os
os.makedirs('../external/BMR/output/Res_reviewerComments/tissueSp_FtrImpo/', exist_ok= True)
df.to_csv('../external/BMR/output/Res_reviewerComments/tissueSp_FtrImpo/tissueSp_FtrImpo.tsv', sep = '\t')
