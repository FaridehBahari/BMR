import os
import argparse
from inference.infer import perform_burdenTest

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


parser = argparse.ArgumentParser(description='Perform burden test for a specific BMR setting based on the path to the directory containing predictions')
parser.add_argument('dir_path', type=str, help='the path to the parent directory containing sim_setting and predictions')
args = parser.parse_args()
dir_path = args.dir_path
    
perform_burdenTest(dir_path)