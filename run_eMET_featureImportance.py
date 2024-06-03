from simulation_settings import load_sim_settings
from models.eMET_functions import one_group_importance_eMET, save_importance_ratio_dfs


path_full_model_eMET = "../external/BMR/output/TL/eMET/eMET_ensemble_bootstraps100_assessment.tsv"
path_ann_pcawg_IDs = '../external/BMR/procInput/ann_PCAWG_ID_complement.csv'
sim_file = 'configs/rate_based/sim_setting.ini'
n_bootstrap = 100
path_oneGroup_intergenic = '../external/BMR/output/featureImportance/GBM_feature_group/GBM_feature_group_model.pkl'


sim_setting = load_sim_settings(sim_file)

one_group_importance_eMET(sim_file, path_ann_pcawg_IDs, path_oneGroup_intergenic, n_bootstrap)

save_importance_ratio_dfs(sim_setting, path_full_model_eMET)

