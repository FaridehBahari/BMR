import time
import os
import platform
import argparse
from readFtrs_Rspns import set_gpu_memory_limit
from models.runBMR_functions import  config_save, load_data_sim_2, repeated_train_test
from simulation_settings import load_sim_settings


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

gpu_fraction = 0.2
set_gpu_memory_limit(gpu_fraction)


if platform.system() == 'Windows':
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    sim_file = 'configs/quick/sim_setting.ini'
    # sim_file = 'configs/rate_based/sim_setting.ini'
else:
    parser = argparse.ArgumentParser(description='Train different models for rate prediction on intergenic bins')
    parser.add_argument('sim_file', type=str, help='the path to the simulation setting config')
    args = parser.parse_args()
    sim_file = args.sim_file
    
st_time = time.time()
sim_setting = load_sim_settings(sim_file)
config_save(sim_file)

##########################
print('repeated train and test for model evaluation is starting ...')
end_t = time.time()
X_tr_cmplt, Y_tr_cmplt, X_val_cmplt, Y_val_cmplt = load_data_sim_2(sim_setting)
print(X_tr_cmplt.shape)
print(X_val_cmplt.shape)
repeated_train_test(sim_setting,  X_tr_cmplt, Y_tr_cmplt, X_val_cmplt, Y_val_cmplt)

end_t2 = time.time()

print('************ Job Done *************')
print(f'total time = {end_t2 - end_t} seconds')
