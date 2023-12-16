rm(list = ls())

library(data.table)
library(openxlsx)
library(dplyr)
library(pROC)
library(PRROC)

source("benchmark/functions_benchmark.R")
source("benchmark/functions_annotate_elements.R")

path_procPCAWG_res <- "../extdata/procInput/PCAWG_results/"
path_proccessedGE <- "../extdata/procInput/"
tissue <-  "Pancan-no-skin-melanoma-lymph"
path_save_benchRes <- "../../BMR_proj/external/output/GBM_transferLearning/GBM/benchmark/binom/"

########### 1) load the pRes object and annotated PCAWG_IDs ###########
ann_PCAWG_ID_complement <- fread(paste0(path_proccessedGE, "ann_PCAWG_ID_complement.csv"))

load(paste0(path_procPCAWG_res, tissue, ".RData"))
########### 3) add new result(s) to pRes for computing measure stats and saving the results ##########
PATHs_newResults <- c("../../BMR_proj/external/output/GBM/inference/GBM_inference.tsv",
                      "../../BMR_proj/external/output/GBM_bootstrap/GBM/inference.tsv",
                      "../../BMR_proj/external/output/GBM_transferLearning/GBM/inference.tsv")

newRESULTS <- c("GBM", "elemSp_GBM",
                "transferLearning_GBM")

gbm <- fread(PATHs_newResults[1])
gbm$p_value = gbm$raw_p_binom
fwrite(gbm, '../../BMR_proj/external/output/GBM/inference/GBM_inference.tsv', sep = '\t')

gbm_elemSp <- fread(PATHs_newResults[2])
gbm_elemSp$p_value = gbm_elemSp$raw_p_binom
fwrite(gbm_elemSp, '../../BMR_proj/external/output/GBM_bootstrap/GBM/inference.tsv', sep = '\t')


gbm_TL <- fread(PATHs_newResults[3])
gbm_TL$p_value = gbm_TL$raw_p_binom
fwrite(gbm_TL, '../../BMR_proj/external/output/GBM_transferLearning/GBM/inference.tsv', sep = '\t')


sig_definition_methods <- c("fdr", "fixedNumberOfElems")
sig_definition_thresholds <- c(.1, 200)
########## use PCAWG drivers/ CGC as gold standard for defining true positives
#based_on =  "in_CGC_new"
based_on =  "all_Databases"
save_Measures_inSingletable(path_save_benchRes, PATHs_newResults,
                            pRes, sig_definition_methods,
                            sig_definition_thresholds, ann_PCAWG_ID_complement,
                            tissue, based_on,
                            compareRank_new = list(TRUE, "transferLearning_GBM"))
