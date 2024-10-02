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
  print(paste0('#samples with RNAseq data: ', sum(colnames(count_matrix) %in% (sampleInfo$aliquot_id)))) 
  
  # Load gene-level CNVs
  print(paste0('#samples with CN data: ', sum(colnames(cn) %in% (sampleInfo$Tumor_Sample_Barcode))))  
  
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
  
  list(restricted_count_matrix, cns_withRNAseqData)
}


restrict_CN_RNA_forCandidate <- function(cohort_specific_CN, cohort_specific_expression, candidate){
  pattern <- strsplit(candidate, '::')[[1]][3]
  
  
  candidate_cn <- cohort_specific_CN[which(rownames(cohort_specific_CN) == pattern),]
  candidate_cn <- candidate_cn[,which(!is.na(candidate_cn[1,]))]
  candidate_rna <- cohort_specific_expression[grep(paste0(':', pattern, ':'), rownames(cohort_specific_expression)),]
  candidate_rna <- candidate_rna[,which(!is.na(candidate_rna[1,]))]
  
  # Get the common column names
  common_cols <- intersect(colnames(candidate_cn), colnames(candidate_rna))
  print(paste0('#samples with both CN and RNA data for "', candidate, '" is ', length(common_cols)))
  # Subset and sort columns for both dataframes
  candidate_cn_subset <- candidate_cn[, common_cols[order(common_cols)]]
  dim(candidate_cn_subset)
  candidate_rna_subset <- candidate_rna[, common_cols[order(common_cols)]]
  dim(candidate_rna_subset)
  
  list(candidate_rna_subset, candidate_cn_subset)
}

generate_data_glm <- function(candidate, mutData, donorInfo, cohort, candidate_rna_subset, candidate_cn_subset){
  mutated_IDs <- mutData[which(mutData$PCAWG_ID == candidate),]
  Donors <- colnames(candidate_cn_subset)
  MUT <- ifelse(Donors %in% mutated_IDs$donor_id, 'mutated', 'unmutated')
  
  totSamples <- nrow(donorInfo[which(donorInfo$cohort1 == cohort ),])
  
  print(paste0('total number of samples in the cohort: ', totSamples))
  print(paste0('number of mutated samples in the cohort: ', length(unique(mutated_IDs$donor_id))))
  print(paste0(sum(MUT == 'mutated'), ' form ', length(MUT), ' samples have RNAseq and CN data'))
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
  final_df
}


compute_pVal_quassiPois_CN <- function(final_df, Test = 'LRT'){ #Test = 'LRT' or 'F'
  # Fit the quasi-Poisson GLM
  glm_model <- glm(FPKM_UQ ~ MUT + CN , 
                   data = final_df, 
                   family = quasipoisson())
  
  # Perform likelihood ratio test to assess the significance of MUT
  # Reduced model without MUT
  glm_model_reduced <- glm(FPKM_UQ ~ CN , 
                           data = final_df, 
                           family = quasipoisson())
  
  # Likelihood ratio test (using anova)
  lrt <- anova(glm_model_reduced, glm_model, test = Test)
  
}


compute_pVal_quassiPois_SCNA <- function(final_df, Test = 'LRT'){ #Test = 'LRT' or 'F'
  # Fit the quasi-Poisson GLM
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

candidate <- 'gc19_pc.cds::gencode::SGK1::ENSG00000118515.7'
cohort <- 'Lymph-BNHL' #'Eso-AdenoCa' #'Liver-HCC'
cancer_specific_dat <- load_cancer_specific_data(path_donorInfo, cohort, all_mutData, all_sampleInfo)


cancerSp_mutData = cancer_specific_dat[[1]]
cancerSp_sampleInfo = cancer_specific_dat[[2]]

cancer_specific_CN_RNA <- load_cohort_specific_CNV_RNA(complete_count_matrix, complete_cn, cancerSp_sampleInfo)


cohort_specific_CN <- cancer_specific_CN_RNA[[2]]
cohort_specific_expression <- cancer_specific_CN_RNA[[1]]

candidateDat <- restrict_CN_RNA_forCandidate(cohort_specific_CN, cohort_specific_expression, candidate)

candidate_rna_subset = candidateDat[[1]]
candidate_cn_subset = candidateDat[[2]]

model_data <- generate_data_glm(candidate, cancerSp_mutData, donorInfo, cohort,
                                candidate_rna_subset, candidate_cn_subset)

final_df = model_data
Test = 'LRT'
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
x <- summary(glm_model)
x[["coefficients"]]