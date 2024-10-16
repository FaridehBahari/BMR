# Some driver mutations may only affect gene expression in early cancer stages and be undetectable by 
# the expression analysis. On the other hand, passenger mutations could potentially affect expression 
# without affecting cell survival. However, the much higher expression correlation signal among double-hit
# than single-hit mutations in CTCF binding sites is compatible with a selective enrichment for functional 
# impact and hence presence of driver mutations. However, mutational mechanisms may also correlate 
# with expression in some cases (see below).22,23


# rm(list = ls())
# 
# library(dplyr)
# library(data.table)
# 
# dp <- fread('../external/BMR/output/Res_reviewerComments/driverDiscovery/DP_allCohorts.tsv')
# eMET <- fread('../external/BMR/output/Res_reviewerComments/driverDiscovery/eMET_allCohorts.tsv')
# 
# df <- rbind(dp, eMET)
ann <- fread('../external/BMR/procInput/ann_PCAWG_ID_complement.csv')
ann <- ann[, c("PCAWG_IDs", "in_CGC_literature", "in_CGC_new",
               "in_oncoKB", "in_pcawg", "tissues")]
# df <- left_join(df, ann, by =c('PCAWG_ID' = 'PCAWG_IDs'))
# 
# # Function to get shared and unique PCAWG_IDs for each cohort
# get_shared_unique_ids <- function(data) {
#   as.list(data %>%
#             group_by(cohort) %>%
#             summarize(
#               shared = list(intersect(PCAWG_ID[tools == "DriverPower"], PCAWG_ID[tools == "eMET"])),
#               unique_DriverPower = list(setdiff(PCAWG_ID[tools == "DriverPower"], PCAWG_ID[tools == "eMET"])),
#               unique_eMET = list(setdiff(PCAWG_ID[tools == "eMET"], PCAWG_ID[tools == "DriverPower"]))
#             ))
# }
# 
# extract_novel_candidates <- function(updated_list){
#   extract_unReported <- function(df){
#     df <- as.data.frame(df)
#     
#     df <- (df[df$in_CGC_literature == F &
#                 df$in_CGC_new == F &
#                 df$in_oncoKB == F
#               & df$in_pcawg == F,])
#     df <- df[order(df$p_value, decreasing = F),]
#     df[, c('PCAWG_ID', 'cohort')]
#   }
#   
#   x = lapply(updated_list, function(s){
#     if(nrow(s) != 0){
#       s = extract_unReported(s)
#     } else {
#       s = c()
#     }
#     s
#   }
#   )
#   
#   novel_candidates <- do.call(rbind, x)
#   rownames(novel_candidates) <- c()
#   novel_candidates[!duplicated(novel_candidates),]
#   
# }
# 
# 
# 
# 
# # Apply the function
# result <- get_shared_unique_ids(df)
# 
# cohorts <- result$cohort
# 
# shared_drivers <- result$shared
# names(shared_drivers) <- cohorts
# 
# 
# DP_unique <- result$unique_DriverPower
# names(DP_unique) <- cohorts
# 
# eMET_unique <- result$unique_eMET
# names(eMET_unique) <- cohorts
# 
# 
# 
# # Assuming df contains a column 'PCAWG_IDs' and the list eMET_unique has vectors for each cohort
# df <- as.data.frame(df)
# update_tools_list <- function(T_unique, df) {
#   # Iterate over each element in the list
#   lapply(names(T_unique), function(cohort) {
#     # Extract the vector for each cohort
#     vec <- T_unique[[cohort]]
#     df_cohort <- df[df$cohort == cohort,]
#     
#     # If the vector is non-empty, filter df based on matching PCAWG_IDs
#     if (length(vec) > 0) {
#       filtered_df <- df_cohort[df_cohort$PCAWG_ID %in% vec, ]
#       return(filtered_df)
#     } else {
#       return(data.frame())  # Return an empty dataframe if the vector is empty
#     }
#   })
# }
# 
# # Apply the function to eMET_unique and store the result
# updated_eMET_list <- update_tools_list(eMET_unique, df)
# names(updated_eMET_list) <- cohorts
# 
# 
# updated_dp_list <- update_tools_list(DP_unique, df)
# names(updated_dp_list) <- cohorts
# 
# 
# updated_shared <- update_tools_list(shared_drivers, df)
# names(updated_shared) <- cohorts
# 
# 
# 
# 
# eMET_novel <- extract_novel_candidates(updated_eMET_list)
# dp_novel <- extract_novel_candidates(updated_dp_list)
# 
# shared_novel <- extract_novel_candidates(updated_shared)
# 
# 
# 
# path_genhancer_PCAWG <- '../../iDriver/extdata/rawInput/map.enhancer.gene.txt.gz'
# gene_enhancer <- data.frame(fread(path_genhancer_PCAWG, 
#                                   header = F, col.names = c("enhancer_pos", "ENSG_id")))
# ann_PCAWG_ID <- fread('../../iDriver/extdata/procInput/ann_PCAWG_ID_complement.csv')
# extract_ensembleIDs <- function(PCAWG_IDs){
#   
#   x <- unlist(lapply(strsplit(PCAWG_IDs, "::"), function(s){
#     unlist(lapply(s[4], function(z){
#       unlist(lapply(strsplit(z, "\\."), function(b){
#         b[1]
#       }))
#     }))
#   }))
#   
#   return(x)
# }
# 
# enhancerID <- 'enhancers::chr14:106445600-106455000::NA::NA'
# enhancer_pos <- unlist(lapply(strsplit(enhancerID, "::"), 
#                               function(s){
#                                 s[2]
#                               }))
# 
# gene_enhancer <- gene_enhancer[which(gene_enhancer$enhancer_pos == enhancer_pos),]
# ensemblIDs_enhancers <- strsplit(gene_enhancer$ENSG_id, ";")
# names(ensemblIDs_enhancers) <- gene_enhancer$enhancer_pos
# 
# gene_enhancer_long <- data.frame(enhancer_pos = rep(gene_enhancer$enhancer_pos, 
#                                                     unlist(lapply(ensemblIDs_enhancers, 
#                                                                   length))),
#                                  ENSG_ids = unlist(ensemblIDs_enhancers))
# 
# CDSs <- ann_PCAWG_ID[which(ann_PCAWG_ID$type_of_element == "gc19_pc.cds"),]
# CDSs$ENSG_ids <- extract_ensembleIDs(CDSs$PCAWG_IDs)
# CDSs <- CDSs[,c('PCAWG_IDs', 'ENSG_ids')]
# 
# CDS_enhancer <- left_join(gene_enhancer_long, CDSs, by = "ENSG_ids")
# unique(CDS_enhancer$PCAWG_IDs)

# # all ENSG_ids for the enhancers::chr14:106445600-106455000::NA::NA are for immunoglobulin heavy chain genes
# novel_enhancer <- enhancerMAp[grepl('chr14:106445600-106455000', enhancerMAp$`chr3:69993400-69994200`),]$`ENSG00000187098;ENSG00000187098`
##################################################################################################
rm(list = ls())
library(data.table)
library(dplyr)
library(readxl)

############################## Functions ##############################
select_cohort <- function(path_donorInfo, 
                          cohort, exclude_lymph_melanoma = TRUE,
                          exclude_hyper_mutated = TRUE){
  
  donorInfo <- fread(path_donorInfo)
  
  if (exclude_lymph_melanoma) {
    exceptions <- c("Skin-Melanoma", "SKCM-US",
                    "Lymph-NOS", 
                    "Lymph-CLL", "CLLE-ES",
                    "Lymph-BNHL", "MALY-DE", "DLBC-US")
    donorInfo <- donorInfo[-which(donorInfo$cohort1 %in% exceptions),] # 2279 donors
  }
  
  if (exclude_hyper_mutated) {
    donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),] 
  }
  
  if (!cohort %in% c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph')) {
    donorInfo <- donorInfo[which(donorInfo$cohort1 == cohort),]
  } 
  
  donorInfo <- donorInfo[,c("D_id","freq" )]
  colnames(donorInfo) <- c('donor_id', 'totalMutNrs')
  donorInfo
}


generate_samplesInfo <- function(path_sample_sheet, path_IDs_maf){
  # load samples info and include white list samples
  pcawg_supp <- fread(path_sample_sheet)
  pcawg_supp <- pcawg_supp[which(grepl('RNA-Seq', pcawg_supp$library_strategy) & pcawg_supp$donor_wgs_exclusion_white_gray == 'Whitelist'),] #1359
  
  tumor_barcodes <- fread(path_IDs_maf)
  aliquot_ids <- pcawg_supp[,c('icgc_donor_id', 'aliquot_id')]
  
  sampleInfo <- left_join(aliquot_ids, tumor_barcodes,  by = c('icgc_donor_id' = 'donor_id'))
  sampleInfo
}


load_cancer_specific_data <- function(path_donorInfo, cohort, mutData, sampleInfo){
  
  exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
  
  donorInfo <- select_cohort(path_donorInfo, cohort, exclude_lymph_melanoma,
                             exclude_hyper_mutated = TRUE)
  
  mutData <- mutData[which(mutData$donor_id %in% donorInfo$donor_id),]
  sampleInfo <- sampleInfo[which(sampleInfo$icgc_donor_id %in% donorInfo$donor_id),]
  list(mutData, sampleInfo)
}

load_cohort_specific_CNV_RNA <- function(count_matrix, cn, sampleInfo){
  
  # Load count matrix
  n_RNAseq = sum(colnames(count_matrix) %in% (sampleInfo$aliquot_id))
  print(paste0('#samples with RNAseq data: ', n_RNAseq))
  
  # Load gene-level CNVs
  n_CN = sum(colnames(cn) %in% (sampleInfo$Tumor_Sample_Barcode))
  print(paste0('#samples with CN data: ', n_CN))  
  
  # restrict to samples with RNAseq data
  sp_withRNAseq <- which(colnames(count_matrix) %in% (sampleInfo$aliquot_id))
  restricted_count_matrix <- count_matrix[, sp_withRNAseq]
  rownames(restricted_count_matrix) <- count_matrix$feature
  
  cns_withRNAseqData <- which(colnames(cn) %in% (sampleInfo$Tumor_Sample_Barcode))
  cns_withRNAseqData <- as.data.frame(cn)[,cns_withRNAseqData]
  rownames(cns_withRNAseqData) <- cn$`Gene Symbol`
  
  
  # Create a mapping of Tumor_Sample_Barcode and aliquot_id to icgc_donor_id
  barcode_to_donor <- setNames(sampleInfo$icgc_donor_id, sampleInfo$Tumor_Sample_Barcode)
  aliquot_to_donor <- setNames(sampleInfo$icgc_donor_id, sampleInfo$aliquot_id)
  
  # Update column names of cns_withRNAseqData using Tumor_Sample_Barcode
  new_colnames_cns <- barcode_to_donor[colnames(cns_withRNAseqData)]
  colnames(cns_withRNAseqData) <- ifelse(!is.na(new_colnames_cns), new_colnames_cns, colnames(cns_withRNAseqData))
  
  # Update column names of restricted_count_matrix using aliquot_id
  new_colnames_restricted <- aliquot_to_donor[colnames(restricted_count_matrix)]
  colnames(restricted_count_matrix) <- ifelse(!is.na(new_colnames_restricted), new_colnames_restricted, colnames(restricted_count_matrix))
  
  list('restricted_count_matrix' = restricted_count_matrix, 
       'cns_withRNAseqData' = cns_withRNAseqData,
       'n_CN' = n_CN,
       'n_RNAseq' = n_RNAseq)
}


restrict_CN_RNA_forCandidate <- function(cohort_specific_CN, cohort_specific_expression, candidate){
  pattern <- strsplit(candidate, '::')[[1]][3]
  
  
  candidate_cn <- cohort_specific_CN[which(rownames(cohort_specific_CN) == pattern),]
  candidate_cn <- candidate_cn[,which(!is.na(candidate_cn[1,]))]
  candidate_rna <- cohort_specific_expression[grep(paste0(':', pattern, ':'), rownames(cohort_specific_expression)),]
  candidate_rna <- candidate_rna[,which(!is.na(candidate_rna[1,]))]
  
  # Get the common column names
  common_cols <- intersect(colnames(candidate_cn), colnames(candidate_rna))
  N_with_both_RNA_CN = length(common_cols)
  print(paste0('#samples with both CN and RNA data for "', candidate, '" is ', N_with_both_RNA_CN))
  # Subset and sort columns for both dataframes
  candidate_cn_subset <- candidate_cn[, common_cols[order(common_cols)]]
  dim(candidate_cn_subset)
  candidate_rna_subset <- candidate_rna[, common_cols[order(common_cols)]]
  dim(candidate_rna_subset)
  
  list('candidate_rna_subset' = candidate_rna_subset, 
       'candidate_cn_subset' = candidate_cn_subset,
       'N_with_both_RNA_CN' = N_with_both_RNA_CN)
}

generate_data_glm <- function(candidate, mutData, donorInfo, cohort, candidate_rna_subset, candidate_cn_subset){
  mutated_IDs <- mutData[which(mutData$PCAWG_ID == candidate),]
  Donors <- colnames(candidate_cn_subset)
  MUT <- ifelse(Donors %in% mutated_IDs$donor_id, 'mutated', 'unmutated')
  
  totSamples <- nrow(donorInfo[which(donorInfo$cohort1 == cohort ),])
  
  print(paste0('total number of samples in the cohort: ', totSamples))
  nMutated = length(unique(mutated_IDs$donor_id))
  print(paste0('number of mutated samples in the cohort: ', nMutated))
  
  nMutated_withRNAseq_CN = sum(MUT == 'mutated')
  print(paste0(nMutated_withRNAseq_CN, ' form ', length(MUT), ' samples have RNAseq and CN data'))
  print('-------------------------------------------------')
  final_df <- as.data.frame(t(rbind(candidate_rna_subset, candidate_cn_subset, MUT)))
  colnames(final_df) <- c('FPKM_UQ', 'CN', 'MUT')
  
  final_df$FPKM_UQ <- as.numeric(final_df$FPKM_UQ)
  final_df$CN <- as.numeric(final_df$CN)
  
  tissue <- donorInfo[which(donorInfo$D_id %in% rownames(final_df)), c('D_id', 'cohort1')]
  
  
  final_df$D_id <- rownames(final_df)
  final_df <- left_join(final_df, tissue, by = 'D_id')
  colnames(final_df) <- c("FPKM_UQ", "CN", "MUT", "donor_id", "cohort")
  
  print(head(final_df))
  
  final_df$MUT <- ifelse(final_df$MUT == "mutated", 1, 0)
  final_df$CN <- as.numeric(final_df$CN)  
  final_df$cohort <- as.factor(final_df$cohort)  
  
  # Assign SCNA values based on CN
  final_df$SCNA <- ifelse(final_df$CN == 2, 0, ifelse(final_df$CN > 2, 1, -1))
  list('final_df' = final_df, 
       'totSamples' = totSamples,
       'nMutated' = nMutated, 
       'nMutated_withRNAseq_CN' = nMutated_withRNAseq_CN)
}


compute_pVal_quassiPois_CN <- function(final_df, cohort, Test = 'LRT'){ #Test = 'LRT' or 'F'
  if (cohort != 'Pancan-no-skin-melanoma-lymph') {
    # Fit the quasi-Poisson GLM
    glm_model <- glm(FPKM_UQ ~ MUT + CN , 
                     data = final_df, 
                     family = quasipoisson())
    
    # Perform likelihood ratio test to assess the significance of MUT
    # Reduced model without MUT
    glm_model_reduced <- glm(FPKM_UQ ~ CN , 
                             data = final_df, 
                             family = quasipoisson())
  } else {
    # Fit the quasi-Poisson GLM
    glm_model <- glm(FPKM_UQ ~ MUT + CN + cohort , 
                     data = final_df, 
                     family = quasipoisson())
    
    # Perform likelihood ratio test to assess the significance of MUT
    # Reduced model without MUT
    glm_model_reduced <- glm(FPKM_UQ ~ CN + cohort , 
                             data = final_df, 
                             family = quasipoisson())
  }
  
  # Likelihood ratio test (using anova)
  lrt <- anova(glm_model_reduced, glm_model, test = Test)
  p_value <- lrt$`Pr(>Chi)`[2]
  p_value
  # x <- summary(glm_model)
  # x[["coefficients"]]
  
}


compute_pVal_quassiPois_SCNA <- function(final_df, cohort, Test = 'LRT'){ #Test = 'LRT' or 'F'
  glm_model <- glm(FPKM_UQ ~ MUT + SCNA, 
                   data = final_df, 
                   family = quasipoisson())
  
  # Perform likelihood ratio test to assess the significance of MUT
  # Reduced model without MUT
  glm_model_reduced <- glm(FPKM_UQ ~ SCNA , 
                           data = final_df, 
                           family = quasipoisson())
  
  # Likelihood ratio test (using anova)
  lrt <- anova(glm_model_reduced, glm_model, test = Test)
  p_value <- lrt$`Pr(>Chi)`[2]
  p_value
  # x <- summary(glm_model)
  # x[["coefficients"]]
}

################################## Inputs ########################################
path_sample_sheet <- '../../pcawg_sample_sheet.tsv'
path_IDs_maf <- '../../iDriver/extdata/sampleID_donorID.tsv'
path_donorInfo <- '../../iDriver/extdata/procInput/iDriverInputs/donorInfo.tsv'
path_CN <- '../../gene_level_calls/all_samples.consensus_CN.by_gene.170214.txt.gz'
path_FPKM_UQ <- '../../expression/tophat_star_fpkm_uq.v2_aliquot_gl_ncg.tsv.gz'

all_mutData <- fread('../../iDriver/extdata/procInput/iDriverInputs/mutData.tsv')
all_sampleInfo <- generate_samplesInfo(path_sample_sheet, path_IDs_maf)
donorInfo <- fread(path_donorInfo)

complete_count_matrix <- as.data.frame(fread(path_FPKM_UQ))
complete_cn <- fread(path_CN) # colnames are Tumor sample barcode (the column exists in maf data to match with Donor_ids)


###############################################################
# candidate <- 'gc19_pc.cds::gencode::SGK1::ENSG00000118515.7' # 'gc19_pc.cds::gencode::ZNF608::ENSG00000168916.11'
# cohort <- 'Lymph-BNHL' #'Eso-AdenoCa' #'Liver-HCC'




##################################################################

# drivers = shared_drivers
# drivers <- lapply(drivers, function(s){
#   s[!grepl('enhancers', s)]
# })
# 
# 
# 
# df_expression_analysis <- c()
# for (i in 1:length(drivers)) {
#   cohort = names(drivers)[i]
#   print(cohort)
#   
#   cancer_specific_dat <- load_cancer_specific_data(path_donorInfo, cohort, all_mutData, all_sampleInfo)
#   
#   
#   cancerSp_mutData = cancer_specific_dat[[1]]
#   cancerSp_sampleInfo = cancer_specific_dat[[2]]
#   
#   cancer_specific_CN_RNA <- load_cohort_specific_CNV_RNA(complete_count_matrix, complete_cn, cancerSp_sampleInfo)
#   
#   
#   cohort_specific_CN <- cancer_specific_CN_RNA$cns_withRNAseqData
#   cohort_specific_expression <- cancer_specific_CN_RNA$restricted_count_matrix
#   n_RNAseq <- cancer_specific_CN_RNA$n_RNAseq
#   n_CN <- cancer_specific_CN_RNA$n_CN
#   
#   candidate_info <- c()
#   if(n_RNAseq != 0){
#     for (j in 1:length(drivers[[i]])) {
#       
#       candidate = drivers[[i]][j]
#       
#       if(length(candidate) != 0 && !is.na(candidate)){
#         print(candidate)
#         
#         candidateDat <- restrict_CN_RNA_forCandidate(cohort_specific_CN, cohort_specific_expression, candidate)
#         
#         candidate_rna_subset = candidateDat$candidate_rna_subset
#         candidate_cn_subset = candidateDat$candidate_cn_subset
#         N_with_both_RNA_CN <- candidateDat$N_with_both_RNA_CN
#         
#         if (N_with_both_RNA_CN != 0) {
#           model_data <- generate_data_glm(candidate, cancerSp_mutData, donorInfo, cohort,
#                                           candidate_rna_subset, candidate_cn_subset)
#           
#           final_df = model_data$final_df
#           N = model_data$totSamples
#           nMutated <- model_data$nMutated
#           nMutated_withRNAseq_CN <- model_data$nMutated_withRNAseq_CN
#           
#           if (nMutated_withRNAseq_CN >= 3) {
#             p_value_CN <- compute_pVal_quassiPois_CN(final_df, cohort, Test = 'LRT')
#             # p_value_CN <- p_value_CN$`Pr(>Chi)`[2]
#             p_value_SCNA <- compute_pVal_quassiPois_SCNA(final_df, cohort, Test = 'LRT')
#             # p_value_SCNA <- p_value_SCNA$`Pr(>Chi)`[2]
#           } else {
#             p_value_CN <- NA
#             p_value_SCNA <- NA
#           }
#           
#           candidate_info <- c(cohort, candidate, n_RNAseq, n_CN, N_with_both_RNA_CN, 
#                               N, nMutated, nMutated_withRNAseq_CN, p_value_CN,
#                               p_value_SCNA)
#           
#           df_expression_analysis <- rbind(df_expression_analysis, candidate_info)
#         }
#         
#         }
#         
#     }
#   }
#   
# }
# 
# colnames(df_expression_analysis) <- c('cohort', 'PCAWG_ID', 'total_RNAseq_data', 'total_CN_data', 
#                                       'nSamples_with_both_RNAandCN', 'N', 
#                                       '#Mutated_samples', '#Mutated_samples_withRNAseq_CN',
#                                       'p_value_CN',
#                                       'p_value_SCNA')
# df_expression_analysis <- as.data.frame(df_expression_analysis)
# 
# 
# df_expression_analysis <- left_join(df_expression_analysis, ann, by =c('PCAWG_ID' = 'PCAWG_IDs'))
# 
# dir.create('../external/BMR/output/Res_reviewerComments/RNAexpression/', recursive = T, showWarnings = F)
# fwrite(df_expression_analysis, file = '../external/BMR/output/Res_reviewerComments/RNAexpression/shared_drivers_expression.tsv', sep = '\t')
# 
# 
# drivers = eMET_unique
# drivers <- lapply(drivers, function(s){
#   s[!grepl('enhancers', s)]
# })
# 
# fwrite(df_expression_analysis, file = '../external/BMR/output/Res_reviewerComments/RNAexpression/eMET_drivers_expression.tsv', sep = '\t')

###################################################################################
drivers = fread('../external/BMR/output/Res_reviewerComments/driverDiscovery/suppTab_eMET_allCohorts.tsv')
enhancers <- drivers[which(grepl('enhancers', drivers$binID)),]
drivers <- as.data.frame(drivers[which(!grepl('enhancers', drivers$binID)),])

all_cohorts <- unique(drivers$cohort)

df_expression_analysis <- c()
for (i in 1:length(all_cohorts)) {
  cohort = all_cohorts[i]
  print(cohort)
  
  cancer_specific_dat <- load_cancer_specific_data(path_donorInfo, cohort, all_mutData, all_sampleInfo)
  
  
  cancerSp_mutData = cancer_specific_dat[[1]]
  cancerSp_sampleInfo = cancer_specific_dat[[2]]
  
  cancer_specific_CN_RNA <- load_cohort_specific_CNV_RNA(complete_count_matrix, complete_cn, cancerSp_sampleInfo)
  
  
  cohort_specific_CN <- cancer_specific_CN_RNA$cns_withRNAseqData
  cohort_specific_expression <- cancer_specific_CN_RNA$restricted_count_matrix
  n_RNAseq <- cancer_specific_CN_RNA$n_RNAseq
  n_CN <- cancer_specific_CN_RNA$n_CN
  
  candidate_info <- c()
  if(n_RNAseq != 0){
    cancerSp_drivers <- drivers[which(drivers$cohort == cohort),]
    cancerSp_drivers <- cancerSp_drivers$binID
    for (j in 1:length(cancerSp_drivers)) {
      
      candidate = cancerSp_drivers[j]
      
      if(length(candidate) != 0 && !is.na(candidate)){
        print(candidate)
        
        candidateDat <- restrict_CN_RNA_forCandidate(cohort_specific_CN, cohort_specific_expression, candidate)
        
        candidate_rna_subset = candidateDat$candidate_rna_subset
        candidate_cn_subset = candidateDat$candidate_cn_subset
        N_with_both_RNA_CN <- candidateDat$N_with_both_RNA_CN
        
        if (N_with_both_RNA_CN != 0) {
          model_data <- generate_data_glm(candidate, cancerSp_mutData, donorInfo, cohort,
                                          candidate_rna_subset, candidate_cn_subset)
          
          final_df = model_data$final_df
          N = model_data$totSamples
          nMutated <- model_data$nMutated
          nMutated_withRNAseq_CN <- model_data$nMutated_withRNAseq_CN
          
          if (nMutated_withRNAseq_CN >= 3) {
            p_value_CN <- compute_pVal_quassiPois_CN(final_df, cohort, Test = 'LRT')
            # p_value_CN <- p_value_CN$`Pr(>Chi)`[2]
            p_value_SCNA <- compute_pVal_quassiPois_SCNA(final_df, cohort, Test = 'LRT')
            # p_value_SCNA <- p_value_SCNA$`Pr(>Chi)`[2]
          } else {
            p_value_CN <- NA
            p_value_SCNA <- NA
          }
          
          candidate_info <- c(cohort, candidate, n_RNAseq, n_CN, N_with_both_RNA_CN, 
                              N, nMutated, nMutated_withRNAseq_CN, p_value_CN,
                              p_value_SCNA)
          
          df_expression_analysis <- rbind(df_expression_analysis, candidate_info)
        }
        
      }
      
    }
  }
  
}

colnames(df_expression_analysis) <- c('cohort', 'PCAWG_ID', 'total_RNAseq_data', 'total_CN_data', 
                                      'nSamples_with_both_RNAandCN', 'N', 
                                      '#Mutated_samples', '#Mutated_samples_withRNAseq_CN',
                                      'p_value_CN',
                                      'p_value_SCNA')
df_expression_analysis <- as.data.frame(df_expression_analysis)


df_expression_analysis <- left_join(df_expression_analysis, ann, by =c('PCAWG_ID' = 'PCAWG_IDs'))

dir.create('../external/BMR/output/Res_reviewerComments/RNAexpression/', recursive = T, showWarnings = F)
df_expression_analysis <- left_join(df_expression_analysis, drivers, by = c('cohort', 'PCAWG_ID' = 'binID'))

fwrite(df_expression_analysis, file = '../external/BMR/output/Res_reviewerComments/RNAexpression/supp_drivers_expression.tsv', sep = '\t')


