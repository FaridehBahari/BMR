import time
import os
import platform
import sys 
import argparse
from readFtrs_Rspns import set_gpu_memory_limit
from models.runBMR_functions import  RUN_BMR, load_data, config_save
from performance.assessModels import assess_models, assess_models_new
# from simulation_settings import create_quickSimSet, create_fullSimSet
from models.NN_functions import nn_model_info
from models.RF_functions import RF_model_info
from models.GBM_functions import  gbm_model_info
from models.DP_functions import DP_model_info
from models.siamese_new import pair_rank_info_siamese2
from simulation_settings import load_sim_settings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

gpu_fraction = 0.7
set_gpu_memory_limit(gpu_fraction)


if platform.system() == 'Windows':
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    # sim_file = 'configs/siam_O1O3/sim_setting.ini'
    # sim_file = 'configs/rate_based/sim_setting.ini'
    sim_file = 'configs/rate_based/sim_settingDSmpl.ini'
else:
    parser = argparse.ArgumentParser(description='Train different models for rate prediction on intergenic bins')
    parser.add_argument('sim_file', type=str, help='the path to the simulation setting config')
    args = parser.parse_args()
    sim_file = args.sim_file
    
st_time = time.time()
sim_setting = load_sim_settings(sim_file)
config_save(sim_file)
X_train, Y_train, X_test, Y_test = load_data(sim_setting)


RUN_BMR(sim_setting, X_train, Y_train, X_test, Y_test, make_pred=True)

print('************')


# assess_models(sim_setting) 
assess_models_new(sim_setting)
end_t = time.time()
print(f'total time = {end_t - st_time} seconds')