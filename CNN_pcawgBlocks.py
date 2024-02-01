# bedtools subtract -a ../external/database/bins/CNN/1k_window.bed -b ../external/database/bins/CNN/pcawg_blocks_all.bed -A | cat - ../external/database/bins/CNN/pcawg_blocks_all.bed | bedtools sort -i - > ../external/database/bins/CNN/sorted_1kwindow_no_pcawg_Acombined.bed


import pandas as pd
import numpy as np

# path_bed_blocks = '../external/database/bins/CNN/pcawg_blocks_all.bed'
path_bed_windows = '../external/database/bins/CNN/1k_window.bed'
path_test_1k = '../external/database/bins/CNN/sorted_1kwindow_no_pcawg_Acombined.bed'


def read_bed(path_bed):
    bed = pd.read_csv(path_bed, sep = '\t', header = None)
    excluded = bed.iloc[np.where((bed[0] == 'chrX') | (bed[0] == 'chrY') | (bed[0] == 'chrM'))]
    
    bed = bed.iloc[~bed.index.isin(excluded.index)]
    bed['binID'] = bed[3]
    bed = bed.set_index('binID')
    
    return bed


# blocks = read_bed(path_bed_blocks)
trainWindow = read_bed(path_bed_windows)
testWindow = read_bed(path_test_1k)

# (testWindow[2] - testWindow[1]).describe()
# Out[43]: 
# count    3.399278e+06
# mean     7.927370e+02
# std      4.199073e+02
# min      1.000000e+00
# 25%      1.000000e+03
# 50%      1.000000e+03
# 75%      1.000000e+03
# max      5.361300e+04
# dtype: float64