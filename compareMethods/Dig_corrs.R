rm(list = ls())

library(data.table)
ann <- fread('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', sep = ',')
ann <- ann[(ann$in_CGC == TRUE) | (ann$in_CGC_literature == TRUE) | 
             (ann$in_CGC_new == TRUE) | (ann$in_oncoKB == TRUE) | 
             (ann$in_pcawg == TRUE), ]
drivers <- ann$PCAWG_IDs



# Define base directories
base_dir <- '../../Dig/output/elemDriver/'
base_dir_eMET <- '../external/BMR/output/reviewerComments/'

# List of element types to iterate over
elemTypes <- c("gc19_pc.3utr", "gc19_pc.promCore", "enhancers",  
               "gc19_pc.5utr", "gc19_pc.ss", "gc19_pc.cds")

# Initialize an empty list to store results for all element types
all_results <- list()

cohorts <- c("Biliary-AdenoCA"  ,            "Bladder-TCC",
"Bone-Leiomyo"      ,           "Bone-Osteosarc",
"Breast-AdenoCa"    ,           "CNS-GBM",
"CNS-Medullo"       ,           "ColoRect-AdenoCA",
"Eso-AdenoCa"      ,            "Head-SCC",
"Kidney-RCC"      ,             "Liver-HCC",
 "Lung-AdenoCA"    ,             "Lung-SCC",
"Lymph-BNHL"        ,           "Lymph-CLL",
"Ovary-AdenoCA"      ,          "Panc-AdenoCA",
"Panc-Endocrine"     ,   
"Pancan_SNV_MNV_INDEL"  ,       "Prost-AdenoCA",
"Skin-Melanoma"         ,       "Stomach-AdenoCA",
 "Uterus-AdenoCA")

# Loop over each element type
for (elem_type in elemTypes) {
  
  corrs <- c()
  eMET_corrs <- c()
  
  for (cohort in cohorts) {
    # Read mutated elements and filter for non-zero mutations
    if (cohort == "Pancan_SNV_MNV_INDEL") {
      mutated_elems <- fread('../../iDriver/extdata/procInput/BMRs/observed/Pancan-no-skin-melanoma-lymph/test_y.tsv')
    } else {
      mutated_elems <- fread(paste0('../../iDriver/extdata/procInput/BMRs/observed/', cohort, '/test_y.tsv'))
    }
    
    mutated_elems <- as.data.frame(mutated_elems[which(mutated_elems$nMut != 0), ])
    mutated_elems <- mutated_elems[!(mutated_elems$binID %in% drivers), ]
    
    # Read full data frame and filter by mutated elements
    full_df <- fread(paste0(base_dir, cohort, '.results.txt'))
    full_df <- full_df[which(full_df$ELT %in% mutated_elems$binID), ]
    
    # Filter by the current element type
    df <- full_df[grepl(elem_type, full_df$ELT), ]
    
    # Calculate observed and predicted rates and their correlation
    obsRate <- df$R_OBS / df$R_SIZE
    predRate <- df$MU / df$R_SIZE
    # predRate <- (df$ALPHA * df$THETA) / df$R_SIZE
    corr <- cor(obsRate, predRate, method = 'spearman')
    corrs <- c(corrs, corr)
    
    # Read eMET results and extract correlation for the element type
    if (cohort == "Pancan_SNV_MNV_INDEL") {
      path_eMET <- paste0(base_dir_eMET,'/Pancan-no-skin-melanoma-lymph/eMET/eMET_assessment.tsv')
    } else {
      path_eMET <- paste0(base_dir_eMET, cohort, '/eMET/eMET_assessment.tsv')
    }
    
    eMET_results <- data.frame(fread(path_eMET))
    row_idx <- which(eMET_results$V1 == 'corr_eMET')
    col_idx <- which(colnames(eMET_results) == elem_type)
    eMET_corrs <- c(eMET_corrs, eMET_results[row_idx, col_idx])
  }
  
  # Store results in a data frame for the current element type
  df <- data.frame('cohort' = cohorts, 'corr_Dig' = corrs, 'corr_eMET' = eMET_corrs, element_type = elem_type)
  all_results[[elem_type]] <- df
}

# Combine all results into a single data frame
final_results <- do.call(rbind, all_results)
final_results$cohort <- gsub('Pancan_SNV_MNV_INDEL', 'Pancan-no-skin-melanoma-lymph', final_results$cohort)
# View the final results
print(final_results)
rownames(final_results) <- 1:nrow(final_results)

dir.create('../../make_features/external/BMR/output/Res_reviewerComments/compare2Dig/', recursive = T, showWarnings = F)
fwrite(final_results, file = '../../make_features/external/BMR/output/Res_reviewerComments/compare2Dig/BMR_comparison_Mu__woDrivers.tsv',
       sep = '\t')
