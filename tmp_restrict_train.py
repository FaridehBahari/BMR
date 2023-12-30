# from readFtrs_Rspns import read_response
# import pandas as pd 

# path_y_validate = '../external/BMR/procInput/val_sets_10folds/Pan_Cancer_validate_y_fold_1.tsv'
# path_Y_train = '../external/BMR/rawInput/responseTabs/Pan_Cancer/1M_bins.tsv'
# path_train_info = '../external/database/bins/proccessed/intergenic_fixed1M_info.tsv'

# Y_val = read_response(path_y_validate)
# Y_train = read_response(path_Y_train)

# train_info = pd.read_csv(path_train_info, sep = '\t', index_col='binID')


# train_Y_annotated = pd.concat([Y_train, train_info], axis=1)

# filtered_train_Y = train_Y_annotated[~train_Y_annotated['orig_name'].str.contains('|'.join(Y_val.index))]

# ###########################################################
# import pandas as pd

# # Example DataFrames
# data_Y_val = {'nMut': [1, 0, 1, 3, 6],
#               'nSample': [1, 0, 1, 3, 6],
#               'length': [497, 168, 20, 2253, 429],
#               'N': [2253, 2253, 2253, 2253, 2253],
#               'obsRates': [8.930637e-07, 0.0, 2.219263e-05, 5.910155e-07, 6.207729e-06],
#               'offset': [1119741, 378504, 45060, 5076009, 966537]}

# data_train_Y = {'nMut': [3528, 4619, 12359, 14426, 7969],
#                 'nSample': [10, 15, 20, 25, 30],
#                 'length': [100, 200, 300, 400, 500],
#                 'N': [1000, 1500, 2000, 2500, 3000],
#                 'obsRates': [1e-05, 2e-05, 3e-05, 4e-05, 5e-05],
#                 'offset': [5000, 10000, 15000, 20000, 25000],
#                 'orig_name': ['b1__v0,b1__v1,b1__v2,b1__v3,b1__v4',
#                               'b10__v6511,b10__v6512,b10__v6513,b10__v6514',
#                               'b100__v41606,b100__v41607,b100__v41608,b100__v41609',
#                               'b1000__v414201,b1000__v414202,b1000__v414203,b1000__v414204',
#                               'b1001__v414344,b1001__v414345,b1001__v414346,b1001__v414347']}

# Y_val = pd.DataFrame(data_Y_val, index=['v0', 'v305760', 'v725230', 'v414202', 'v464826'])
# train_Y_annotated = pd.DataFrame(data_train_Y, index=['b1', 'b10', 'b100', 'b1000', 'b1001'])

# # Filter based on the condition
# filtered_train_Y = train_Y_annotated[~train_Y_annotated['orig_name'].str.contains('|'.join(Y_val.index))]

# # Display the filtered DataFrame
# print(filtered_train_Y)
##############################################################################
import numpy as np
from readFtrs_Rspns import read_response
import pandas as pd 

val_size = 800
path_var_intervals_Y = '../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv'

# np.random.seed(13)



path_Y_train = '../external/BMR/rawInput/responseTabs/Pan_Cancer/1M_bins.tsv'
path_train_info = '../external/database/bins/proccessed/intergenic_fixed1M_info.tsv'





def generate_train_valvar_sets(path_var_intervals_Y, path_Y_train, 
                               path_train_info, val_size):
    
    # read all variable-size bins
    var_interval_response = read_response(path_var_intervals_Y)
    var_interval_response = var_interval_response.iloc[np.where(var_interval_response.length >= 20)]
    
    # sample from variable-size bins to have validation set
    val_indices = np.random.choice(var_interval_response.index, 
                                    size=val_size, replace=False)
    Y_val = var_interval_response.loc[val_indices]
    
    # read all fixed-size bins
    Y_train = read_response(path_Y_train)
    train_info = pd.read_csv(path_train_info, sep = '\t', index_col='binID')
    train_Y_annotated = pd.concat([Y_train, train_info], axis=1)

    # remove validation bins from train set
    filtered_train_Y = train_Y_annotated[~train_Y_annotated['orig_name'].str.contains('|'.join(Y_val.index))]
    
    return filtered_train_Y, Y_val


Y_train, Y_val = generate_train_valvar_sets(path_var_intervals_Y, path_Y_train, 
                               path_train_info, val_size)


Y_val['length'].sum() #2723723
Y_train['length'].sum() #1802264277.0
