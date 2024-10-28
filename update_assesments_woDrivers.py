import time
import os
import platform
import argparse
from readFtrs_Rspns import set_gpu_memory_limit
from models.runBMR_functions import  RUN_BMR, load_data_sim, config_save, repeated_train_test
from performance.assessModels import assess_models
from simulation_settings import load_sim_settings
import pandas as pd
import shutil

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)

sim_files = ['../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/RF/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_poisLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_mseLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/RF/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_poisLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_mseLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/RF/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_poisLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_mseLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/RF/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_poisLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_mseLoss/sim_setting.ini',
                    '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/GBM/sim_setting.ini',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/GBM/sim_setting.ini',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/GBM/sim_setting.ini',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/GBM/sim_setting.ini',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/GBM/sim_setting.ini',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/GBM/sim_setting.ini',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/GBM/sim_setting.ini',
                  "../external/BMR/output/dimReduction_effect/PCA/GBM/sim_setting.ini",
                  "../external/BMR/output/dimReduction_effect/PCA/RF/sim_setting.ini",
                  "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/sim_setting.ini",
                  "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/sim_setting.ini", 
                 
                 "../external/BMR/output/dimReduction_effect/AE/GBM/sim_setting.ini",
                 "../external/BMR/output/dimReduction_effect/AE/RF/sim_setting.ini",
                 "../external/BMR/output/dimReduction_effect/AE/nn_poisLoss/sim_setting.ini",
                 "../external/BMR/output/dimReduction_effect/AE/nn_mseLoss/sim_setting.ini"
                 ]

for sim_file_path in sim_files:
    print(sim_file_path)
    st_time = time.time()
    sim_setting = load_sim_settings(sim_file_path)
    parts = sim_file_path.split('/')
    extracted_path_base = '/'.join(parts[:-1]) + '/' + parts[-2]
    extracted_path = extracted_path_base + '_assessmentsWithDrivers.tsv'
    
    if os.path.exists(extracted_path_base+'_assessments.tsv'):
        
        ass = pd.read_csv(extracted_path_base+'_assessments.tsv', sep = '\t')
        ass.to_csv( extracted_path, sep = '\t')
    
    assess_models(sim_setting)
    end_t = time.time()
    
    print('************')
    print(f'total time = {end_t - st_time} seconds')

    

############################## cancer-specific intergenic #############################
import time
import os
import platform
import argparse
from readFtrs_Rspns import set_gpu_memory_limit
from models.runBMR_functions import  RUN_BMR, load_data_sim, config_save, repeated_train_test
from performance.assessModels import assess_models
from simulation_settings import load_sim_settings
import pandas as pd
import shutil

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)

cohorts = ['Pancan-no-skin-melanoma-lymph', "Liver-HCC", "Bladder-TCC" ,
           "ColoRect-AdenoCA" , "Lymph-BNHL",
    "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
    "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
    "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",         
    "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
    "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",    
    "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
    "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc"
]

for cohort in cohorts:
    print(cohort)
    if cohort == 'Pancan-no-skin-melanoma-lymph':
        sim_file_path = f'../external/BMR/output/reviewerComments/{cohort}/GBM/sim_setting_iDriver.ini'
    else:
        sim_file_path = f'../external/BMR/output/reviewerComments/{cohort}/GBM/sim_setting_iDriver_{cohort}.ini'
        
    print(sim_file_path)
    st_time = time.time()
    sim_setting = load_sim_settings(sim_file_path)
    parts = sim_file_path.split('/')
    extracted_path_base = '/'.join(parts[:-1]) + '/' + parts[-2]
    extracted_path = extracted_path_base + '_assessmentsWithDrivers.tsv'
    
    if os.path.exists(extracted_path_base+'_assessments.tsv'):
        
        ass = pd.read_csv(extracted_path_base+'_assessments.tsv', sep = '\t')
        ass.to_csv( extracted_path, sep = '\t')
    
    assess_models(sim_setting)
    end_t = time.time()
    
    print('************')
    print(f'total time = {end_t - st_time} seconds')

       
        
    

    
############################## cancer-specific eMETs ##################################
import pandas as pd
from performance.assessModels import assess_model
# Load the CSV file into a DataFrame
df = pd.read_csv('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', sep=',')
filtered_df = df[(df['in_CGC'] | df['in_CGC_literature'] | df['in_CGC_new'] | df['in_oncoKB'] | df['in_pcawg'])]
drivers = filtered_df['PCAWG_IDs']


cohorts = ['Pancan-no-skin-melanoma-lymph', "Liver-HCC", "Bladder-TCC" ,
           "ColoRect-AdenoCA" , "Lymph-BNHL",
    "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
    "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
    "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",         
    "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
    "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",    
    "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
    "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc"
]

for cohort in cohorts:
    print(cohort)
    obs_pred_rates_path = f'../external/BMR/output/reviewerComments/{cohort}/eMET/eMET_100_predTest.tsv'
    
    Y_pred = pd.read_csv(obs_pred_rates_path, sep = "\t", header=0, index_col='binID')
    
    obs_pred_rates = Y_pred.loc[~(Y_pred.index).isin(drivers)]
    
    model_name = 'eMET'
    Nr_pair_acc = 100000
    assessment = assess_model(obs_pred_rates.predRate, obs_pred_rates.obs_rates, 
                  Nr_pair_acc, model_name, per_element=True)
    
    assessment.to_csv(f'../external/BMR/output/reviewerComments/{cohort}/eMET/eMET_assessment.tsv', sep = '\t')

#################################################################################
import pandas as pd
from performance.assessModels import assess_model
# Load the CSV file into a DataFrame
df = pd.read_csv('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', sep=',')
filtered_df = df[(df['in_CGC'] | df['in_CGC_literature'] | df['in_CGC_new'] | df['in_oncoKB'] | df['in_pcawg'])]
drivers = filtered_df['PCAWG_IDs']


obs_pred_rates_path = '../external/BMR/output/dimReduction_effect/AE/eMET/eMET_100_predTest.tsv'

Y_pred = pd.read_csv(obs_pred_rates_path, sep = "\t", header=0, index_col='binID')

obs_pred_rates = Y_pred.loc[~(Y_pred.index).isin(drivers)]

model_name = 'eMET'
Nr_pair_acc = 100000
assessment = assess_model(obs_pred_rates.pred_rates, obs_pred_rates.obs_rates,Nr_pair_acc, model_name, per_element=True)

assessment.to_csv('../external/BMR/output/dimReduction_effect/AE/eMET/eMET_assessment.tsv', sep = '\t')

############################## cancer-specific intergenic SNVs #############################
import time
import os
import platform
import argparse
from readFtrs_Rspns import set_gpu_memory_limit
from models.runBMR_functions import  RUN_BMR, load_data_sim, config_save, repeated_train_test
from performance.assessModels import assess_models
from simulation_settings import load_sim_settings
import pandas as pd
import shutil

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)

cohorts = ['Pancan-no-skin-melanoma-lymph',
            "Liver-HCC", "Bladder-TCC" ,
            "ColoRect-AdenoCA" , "Lymph-BNHL",
    "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
    "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
    "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",         
    "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
    "Bone-Leiomyo", "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
    "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc"
]

for cohort in cohorts:
    print(cohort)
    if cohort == 'Pancan-no-skin-melanoma-lymph':
        sim_file_path = f'../external/BMR/output/Res_reviewerComments/SNV_nonSNV/{cohort}/nonSNV/GBM/sim_setting_iDriver.ini'
    else:
        sim_file_path = f'../external/BMR/output/Res_reviewerComments/SNV_nonSNV/{cohort}/nonSNV/GBM/sim_setting_iDriver_{cohort}.ini'
        
    print(sim_file_path)
    st_time = time.time()
    sim_setting = load_sim_settings(sim_file_path)
    parts = sim_file_path.split('/')
    extracted_path_base = '/'.join(parts[:-1]) + '/' + parts[-2]
    extracted_path = extracted_path_base + '_assessmentsWithDrivers.tsv'
    
    if os.path.exists(extracted_path_base+'_assessments.tsv'):
        
        ass = pd.read_csv(extracted_path_base+'_assessments.tsv', sep = '\t')
        ass.to_csv( extracted_path, sep = '\t')
    
    assess_models(sim_setting)
    end_t = time.time()
    
    print('************')
    print(f'total time = {end_t - st_time} seconds')



############################## cancer-specific eMET SNV ##################################
import pandas as pd
from performance.assessModels import assess_model
# Load the CSV file into a DataFrame
df = pd.read_csv('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', sep=',')
filtered_df = df[(df['in_CGC'] | df['in_CGC_literature'] | df['in_CGC_new'] | df['in_oncoKB'] | df['in_pcawg'])]
drivers = filtered_df['PCAWG_IDs']


cohorts = ['Pancan-no-skin-melanoma-lymph',
            "Liver-HCC", "Bladder-TCC" ,
            "ColoRect-AdenoCA" , "Lymph-BNHL",
    "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
    "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
    "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",         
    "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
    "Bone-Leiomyo", "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
    "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc"
]

for cohort in cohorts:
    print(cohort)
    obs_pred_rates_path = f'../external/BMR/output/Res_reviewerComments/SNV_nonSNV/{cohort}/nonSNV/eMET/eMET_100_predTest.tsv'
    
    Y_pred = pd.read_csv(obs_pred_rates_path, sep = "\t", header=0, index_col='binID')
    
    obs_pred_rates = Y_pred.loc[~(Y_pred.index).isin(drivers)]
    
    model_name = 'eMET'
    Nr_pair_acc = 100000
    assessment = assess_model(obs_pred_rates.predRate, obs_pred_rates.obs_rates, 
                  Nr_pair_acc, model_name, per_element=True)
    
    assessment.to_csv(f'../external/BMR/output/Res_reviewerComments/SNV_nonSNV/{cohort}/nonSNV/eMET/eMET_assessment.tsv', sep = '\t')

############################### tissue-specific vs all featres (GBM) #########################################

import time
import os
import platform
import argparse
from readFtrs_Rspns import set_gpu_memory_limit
from models.runBMR_functions import  RUN_BMR, load_data_sim, config_save, repeated_train_test
from performance.assessModels import assess_models
from simulation_settings import load_sim_settings
import pandas as pd
import shutil

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)

cohorts = ["Liver-HCC", "ColoRect-AdenoCA" ,
                      "Uterus-AdenoCA" , "Kidney-RCC", "Lung-SCC",
                      "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
                      "Breast-AdenoCa",
                      "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
                      "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc"]

for cohort in cohorts:
    print(cohort)
    sim_file_path = f'../external/BMR/output/Res_reviewerComments/tissueSpFtrs/{cohort}/GBM/sim_setting_tissueSpFtrs_{cohort}.ini'
    
    print(sim_file_path)
    st_time = time.time()
    sim_setting = load_sim_settings(sim_file_path)
    parts = sim_file_path.split('/')
    extracted_path_base = '/'.join(parts[:-1]) + '/' + parts[-2]
    extracted_path = extracted_path_base + '_assessmentsWithDrivers.tsv'
    
    if os.path.exists(extracted_path_base+'_assessments.tsv'):
        
        ass = pd.read_csv(extracted_path_base+'_assessments.tsv', sep = '\t')
        ass.to_csv( extracted_path, sep = '\t')
    
    assess_models(sim_setting)
    end_t = time.time()
    
    print('************')
    print(f'total time = {end_t - st_time} seconds')


######################## GroupFtr_importance ###################################
#### Intergenic >>> # addedd save_name argument to the assess_models function and commented the save_name = m['save_name'] in the function


save_names = ['GBM_conservation', 'GBM_Epigenetic_mark', 'GBM_APOBEC', 
              'GBM_RNA_expression', 'GBM_Replication_timing', 'GBM_DNA_methylation',
              'GBM_DNA_accessibility', 'GBM_nucleotideContext', 'GBM_HiC']

for save_name in save_names:
    sim_file_path = f'../external/BMR/output/featureImportance/{save_name}/sim_setting.ini'
    print(sim_file_path)
    st_time = time.time()
    sim_setting = load_sim_settings(sim_file_path)
    parts = sim_file_path.split('/')
    extracted_path_base = '/'.join(parts[:-1]) + '/' + parts[-2]
    extracted_path = extracted_path_base + '_assessmentsWithDrivers.tsv'
    
    if os.path.exists(extracted_path_base+'_assessments.tsv'):
        
        ass = pd.read_csv(extracted_path_base+'_assessments.tsv', sep = '\t')
        ass.to_csv( extracted_path, sep = '\t')
    
    assess_models(sim_setting, save_name)
    end_t = time.time()
    
    print('************')
    print(f'total time = {end_t - st_time} seconds')

    




#### eMET
import pandas as pd
from performance.assessModels import assess_model
# Load the CSV file into a DataFrame
df = pd.read_csv('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', sep=',')
filtered_df = df[(df['in_CGC'] | df['in_CGC_literature'] | df['in_CGC_new'] | df['in_oncoKB'] | df['in_pcawg'])]
drivers = filtered_df['PCAWG_IDs']


ftrs = ['conservation', 'Epigenetic_mark', 'RNA_expression', 
        'Replication_timing', 'DNA_accessibility', 'nucleotide content', 'HiC']

for ftr in ftrs:
    print(ftr)
    obs_pred_rates_path = f'../external/BMR/output/eMET_GroupImportance/{ftr}/GBM/GBM_100_predTest.tsv'
    
    Y_pred = pd.read_csv(obs_pred_rates_path, sep = "\t", header=0, index_col='binID')
    
    obs_pred_rates = Y_pred.loc[~(Y_pred.index).isin(drivers)]
    
    model_name = 'eMET'
    Nr_pair_acc = 100000
    assessment = assess_model(obs_pred_rates.pred_rates, obs_pred_rates.obs_rates, 
                  Nr_pair_acc, model_name, per_element=True)
    
    assessment.to_csv(f'../external/BMR/output/eMET_GroupImportance/{ftr}/GBM/eMET_assessment.tsv', 
                      sep = '\t')

#################### importance ratio
import pandas as pd
oneGroup_eMET_path = '../external/BMR/output/eMET_GroupImportance/conservation/GBM/eMET_assessment.tsv'
path_full_model_eMET = "../external/BMR/output/TL/eMET/eMET_ensemble_bootstraps100_assessment.tsv"

categories = [ 'nucleotide content', 'HiC', 
              'DNA_accessibility', 'Epigenetic_mark', 'RNA_expression', 'Replication_timing',
              'conservation'] 

    
all_ratio_dfs = pd.DataFrame(columns=["Element", "Ratio", "feature_category"])
full_model_df = pd.read_csv(path_full_model_eMET, 
                            sep=",", index_col=0, skipinitialspace=True)
    
for category in categories:
    
    path_one_group_ass = oneGroup_eMET_path.replace('conservation', category)
    
    one_group_df = pd.read_csv(path_one_group_ass, sep="\t", index_col=0)
        
    # Extract correlation values from the dataframes
    one_group_corr = one_group_df.loc['corr_eMET']
    full_model_corr = full_model_df.loc['corr_eMET']
    
    # Create a dataframe to store the ratios
    ratios_df = pd.DataFrame(columns=["Element", "Ratio"])
    
    # Calculate and store the ratios for each element type
    for element in one_group_corr.index:
        if element in full_model_corr.index:
            ratio = one_group_corr[element] / full_model_corr[element]
            ratios_df.loc[len(ratios_df)] = [element, ratio]
           
    ratios_df["feature_category"] = category
    
    all_ratio_dfs = pd.concat([all_ratio_dfs, ratios_df], axis = 0)


all_ratio_dfs.to_csv('../external/BMR/output/eMET_GroupImportance/importanceRatios.csv', sep=",") 

################################# repeated train-test DS ######################################
# Update to include MAE

import pandas as pd
from performance.assessModels import assess_model, read_obs, read_pred
# import re

Y_obs_all_intergenic = read_obs('../external/BMR/rawInput/responseTabs_bedtools/Pan_Cancer/var_bins.tsv', True)
# DS = 'DS50k'
# model = 'nn_poisLoss'
n_repeat = 10

SampleSizes = ['FullSet', 'DS1M', 'DS800k', 'DS600k', 'DS300k', 'DS100k', 'DS50k']
models = ['GBM', 'RF', 'nn_poisLoss', 'nn_mseLoss']
for DS in SampleSizes:
    for model in models:
        pred_paths = []
        for i in range(n_repeat):
            rep_number = i+1
            pred_path = f'../external/BMR/output/with_RepliSeq_HiC/DownSampling/{DS}/{model}/rep_train_test/{model}_predTest{rep_number}.tsv'
            
            print(pred_path)
            
            Y_pred = read_pred(pred_path)
            Y_obs_unseen = Y_obs_all_intergenic.loc[Y_pred.index]
            
            parts = pred_path.split('/')
            
            model_name = parts[-3]
            Nr_pair_acc = 100000
            assessment = assess_model(Y_pred, Y_obs_unseen, 
                          Nr_pair_acc, model_name, per_element=False)
            
            # rep_number = re.findall(r'\d+', pred_path)[0]
            save_ass = '/'.join(parts[:-1]) + '/'  + model_name + '_M' + str(rep_number) + '_assessment.tsv'
            assessment.to_csv(save_ass, sep = '\t')




###################### update summary repeated train-test DS###################
from models.repeated_train_test import save_metrics_summary

SampleSizes = ['FullSet', 'DS1M', 'DS800k', 'DS600k', 'DS300k', 'DS100k', 'DS50k']
models = ['GBM', 'RF', 'nn_poisLoss', 'nn_mseLoss']
for DS in SampleSizes:
    for model in models:
        dir_path = f'../external/BMR/output/with_RepliSeq_HiC/DownSampling/{DS}/{model}/rep_train_test/'
        # save the model metrics summary:
        save_metrics_summary(dir_path)

################################# update train-test assessmnents without drivers ######################################
# Update to include MAE

import pandas as pd
from performance.assessModels import assess_model, read_obs, read_pred
# import re

Y_obs_all_intergenic = read_obs('../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv', True)

n_repeat = 10

binSizes = [ '1M', '100k', '50k', '10k'] #'var_size',
models = ['GBM', 'RF', 'nn_poisLoss', 'nn_mseLoss']
for BIN in binSizes:
    for model in models:
        pred_paths = []
        for i in range(n_repeat):
            rep_number = i+1
            pred_path = f'../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/{BIN}/{model}/rep_train_test/{model}_predTest{rep_number}.tsv'
            
            print(pred_path)
            
            Y_pred = read_pred(pred_path)
            Y_obs_unseen = Y_obs_all_intergenic.loc[Y_pred.index]
            
            parts = pred_path.split('/')
            
            model_name = parts[-3]
            Nr_pair_acc = 100000
            assessment = assess_model(Y_pred, Y_obs_unseen, 
                          Nr_pair_acc, model_name, per_element=False)
            
            # rep_number = re.findall(r'\d+', pred_path)[0]
            save_ass = '/'.join(parts[:-1]) + '/'  + model_name + '_M' + str(rep_number) + '_assessment.tsv'
            assessment.to_csv(save_ass, sep = '\t')



###################### update summary repeated train-test bin-model###################
from models.repeated_train_test import save_metrics_summary

bins = ['var_size', '1M', '100k', '50k', '10k']
models = ['GBM', 'RF', 'nn_poisLoss', 'nn_mseLoss']
for Bin in bins:
    for model in models:
        dir_path = f'../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/{Bin}/{model}/rep_train_test/'
        # save the model metrics summary:
        save_metrics_summary(dir_path)
