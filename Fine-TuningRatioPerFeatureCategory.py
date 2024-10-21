import os
import pandas as pd
import numpy as np

def save_importance_ratio_finetuned(oneGroup_eMET_path, oneGroup_intergenic_path, save_dir):
    
    categories = [ 'nucleotide content', 'HiC', 
                  'DNA_accessibility', 'Epigenetic_mark', 'RNA_expression', 'Replication_timing',
                  'conservation'] 
    os.makedirs(save_dir,  exist_ok= True)
    all_ratio_dfs = pd.DataFrame(columns=["Element", "Ratio", "feature_category"])
    
    for category in categories:
        print(category) 
        
        if category == 'nucleotide content':
            category_GBM = 'nucleotideContext'
        else:
           category_GBM = category
        
        path_one_group_ass_eMET = oneGroup_eMET_path.replace('conservation', category)
        path_one_group_ass_int = oneGroup_intergenic_path.replace('conservation', category_GBM)
        
        one_group_df_eMET = pd.read_csv(path_one_group_ass_eMET, sep="\t", index_col=0)
        one_group_df_int = pd.read_csv(path_one_group_ass_int, sep="\t", index_col=0)
        
        # Extract correlation values from the dataframes
        eMET_corr = one_group_df_eMET.loc['corr_eMET']
        intergenic_model_corr = one_group_df_int.loc[f'corr_GBM_{category_GBM}']
        
        # Create a dataframe to store the ratios
        ratios_df = pd.DataFrame(columns=["Element", "Ratio"])
        
        # Calculate and store the ratios for each element type
        for element in eMET_corr.index:
            if element in intergenic_model_corr.index:
                ratio = np.abs(eMET_corr[element] - intergenic_model_corr[element]) / np.max([eMET_corr[element], intergenic_model_corr[element]])
                ratios_df.loc[len(ratios_df)] = [element, ratio]
               
        # # Save the ratios to a new TSV file
        # ratios_df.to_csv(f'{save_dir}/importanceRatios_{category}.tsv', sep="\t", index=False)
        # ratio_df = pd.read_csv(f'{save_dir}/importanceRatios_{category}.tsv', sep="\t")
        ratios_df["feature_category"] = category
        all_ratio_dfs = pd.concat([all_ratio_dfs, ratios_df], axis = 0)
        
    all_ratio_dfs.to_csv(f"{save_dir}/importanceRatios.csv", sep=",")   




###############################################################################
oneGroup_intergenic_path = '../external/BMR/output/featureImportance/GBM_conservation/GBM_conservation_assessments.tsv'
oneGroup_eMET_path = '../external/BMR/output/eMET_GroupImportance/conservation/GBM/eMET_assessment.tsv'

save_dir = '../external/BMR/output/Res_reviewerComments/finetunedFtrImpo/'

save_importance_ratio_finetuned(oneGroup_eMET_path, oneGroup_intergenic_path, save_dir)