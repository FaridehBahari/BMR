#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 17:34:15 2023

@author: bahari
"""
import pandas as pd
import os


# df = pd.read_csv('../../iDriver/extdata/procInput/ann_PCAWG_ID_complement.csv', sep=',',
#                  header=0, index_col='PCAWG_IDs',
#                 usecols=['PCAWG_IDs', "length", "in_CGC","in_CGC_literature",
#                          "in_CGC_new","in_oncoKB","in_pcawg","tissues", 
#                          "type_of_element"])
# df = df.rename_axis("binID")
# os.mkdir('../external/PCAWG_annotations/', exist_ok= True)
# df.to_csv('../external/PCAWG_annotations/ann_PCAWG_IDs.tsv',
#                   sep='\t')

def annotate_results(path_result):
    ann_df = pd.read_csv('../external/PCAWG_annotations/ann_PCAWG_IDs.tsv', 
                         sep = '\t', header=0, index_col='binID',
                        usecols=['binID', "length", "in_CGC","in_CGC_literature",
                                 "in_CGC_new","in_oncoKB","in_pcawg","tissues", 
                                 "type_of_element"])
    
    result = pd.read_csv(path_result, sep = '\t', header=0, index_col='binID')
    
    ann_result = pd.merge(result, ann_df, left_index=True,
                           right_index=True, how='inner')
    
    directory_name = os.path.dirname(path_result)
   
    ann_result.to_csv(f'{directory_name}/ann_result.tsv', sep = '\t')
    
 
    
    


path_results = [#'../external/output/GBM/inference/GBM_inference.tsv',
                #'../external/output/GLM/inference/GLM_inference.tsv',
                #'../external/output/RF/inference/RF_inference.tsv',
                #'../external/output/check_classic_NNs/savingIntervals/leakyRelu/check_DO/five100_leakyRelu_do2/',
                # '../../iDriver/extdata/output/tmp_PanCancer_burden.result.tsv',
                '../external/output/small_BatchSize/500_50fiveLayers_do3_bs64/inference/500_50fiveLayers_do3_bs64_inference.tsv']

for path_result in path_results:
    annotate_results(path_result)


