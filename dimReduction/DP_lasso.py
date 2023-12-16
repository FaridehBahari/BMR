import pandas as pd 
from readFtrs_Rspns import read_fi, create_TestTrain_TwoSources

def dim_reduction_DP(subject_X, cutoff):
    use_features = read_fi('../external/rawInput/DP_results/Binom_GLM.feature_importance.tsv',
                           cutoff)
    subject_X = subject_X[use_features]
    
    return subject_X

#########################
path_X_test = '../external/rawInput/test_feature.hdf5'
path_X_train = '../external/rawInput/train_feature.hdf5'
path_Y_train = '../external/rawInput/Pan_Cancer_train_y.tsv'
path_Y_test = '../external/rawInput/Pan_Cancer_test_y.tsv'

X_train, Y_train, X_test, Y_test  = create_TestTrain_TwoSources(path_X_train, 
                                                                path_Y_train, 
                                                                path_X_test,
                                                                path_Y_test, 
                                                                scale = True)

train_reduced_5 = dim_reduction_DP(X_train, .5)
train_reduced_8 = dim_reduction_DP(X_train, .8)

test_reduced_5 = dim_reduction_DP(X_test, .5)
test_reduced_8 = dim_reduction_DP(X_test, .8)

train_reduced_5.to_csv('../external/procInput/train_5_lasso.tsv', sep = '\t')
train_reduced_8.to_csv('../external/procInput/train_8_lasso.tsv', sep = '\t')
test_reduced_5.to_csv('../external/procInput/test_5_lasso.tsv', sep = '\t')
test_reduced_8.to_csv('../external/procInput/test_8_lasso.tsv', sep = '\t')