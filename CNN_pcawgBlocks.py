# bedtools subtract -a ../external/database/bins/CNN/1k_window.bed -b ../external/database/bins/CNN/pcawg_blocks_all.bed -A | cat - ../external/database/bins/CNN/pcawg_blocks_all.bed | bedtools sort -i - > ../external/database/bins/CNN/sorted_1kwindow_no_pcawg_Acombined.bed

import pandas as pd

path_bed_blocks = '../external/database/bins/CNN/pcawg_blocks_all.bed'
path_bed_windows = '../external/database/bins/CNN/1k_window.bed'
path_test_1k = '../external/database/bins/CNN/sorted_1kwindow_no_pcawg_Acombined.bed'


def read_bed(path_bed):
    bed = pd.read_csv(path_bed, sep = '\t', header = None)
    bed['binID'] = bed[3]
    bed = bed.set_index('binID')
    
    return bed


blocks = read_bed(path_bed_blocks)
trainWindow = read_bed(path_bed_windows)
testWindow = read_bed(path_test_1k)


