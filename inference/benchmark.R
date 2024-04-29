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
path_save_benchRes <- "../../make_features/external/BMR/benchmark/"

########### 1) load the pRes object and annotated PCAWG_IDs ###########
ann_PCAWG_ID_complement <- fread(paste0(path_proccessedGE, "ann_PCAWG_ID_complement.csv"))
load(paste0(path_procPCAWG_res, tissue, ".RData"))


########### 2) add new result(s) to pRes for computing measure stats and saving the measurment results ##########
newRESULTS <- c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL")
PATHs_newResults <- c("../../make_features/external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/inference/GBM_inference.tsv",
                      "../../make_features/external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/inference/RF_inference.tsv",
                      "../../make_features/external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/inference/nn_mseLoss_inference.tsv",
                      "../../make_features/external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/inference/nn_poisLoss_inference.tsv",
                      "../../make_features/external/BMR/output/TL/GBM/inference/GBM_inference.tsv")

sig_definition_methods <- c("fdr", "fixedNumberOfElems")
sig_definition_thresholds <- c(.1, 100)

add_newMethod_to_pRes <- function(pRes, path_newRes){
  
  newMethod <- data.frame(fread(path_newRes))
  idx_pVal <- which(colnames(newMethod) %in% c("pValues", "p_value",
                                               "raw_p", "pvals", 
                                               "raw_p_nBinom")) #  raw_p_binom
  idx_ID <- which(colnames(newMethod) %in% c("binID", "PCAWG_ID"))
  
  newMethod <- data.frame("PCAWG_IDs" = newMethod[,idx_ID], 
                          "p_value" =  newMethod[,idx_pVal],
                          "q_value" = rep(NA, nrow(newMethod)),
                          "fdr" =  p.adjust(newMethod[,idx_pVal], "fdr"))
  
  newMethod <- newMethod[which(!is.na(newMethod$p_value)),]
  
  pRes$newMethod <- newMethod
  class(pRes) = "pRes"
  pRes
}

## use PCAWG drivers/ CGC as gold standard for defining true positives

driver_based_on <- c("in_CGC_new", "all", "in_pcawg", "in_oncoKB", "any")

for (based_on in driver_based_on) {
  save_Measures_inSingletable(path_save_benchRes, newRESULTS, PATHs_newResults, pRes, 
                              sig_definition_methods,
                              sig_definition_thresholds,
                              ann_PCAWG_ID_complement, tissue, based_on #,
                              # compareRank_new = list(TRUE, "GBM")
  )
  
  
}




