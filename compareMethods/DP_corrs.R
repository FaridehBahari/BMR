rm(list = ls())

library(data.table)

# Define base directories
base_dir <- '../output/elemDriver/'
base_dir_eMET <- '../../make_features/external/BMR/output/reviewerComments/'

# List of element types to iterate over
elemTypes <- c("gc19_pc.3utr", "gc19_pc.promCore", "enhancers",  
               "gc19_pc.5utr", "gc19_pc.ss", "gc19_pc.cds")

ann <- fread('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', sep = ',')
ann <- ann[(ann$in_CGC == TRUE) | (ann$in_CGC_literature == TRUE) | 
             (ann$in_CGC_new == TRUE) | (ann$in_oncoKB == TRUE) | 
             (ann$in_pcawg == TRUE), ]
drivers <- ann$PCAWG_IDs

# Initialize an empty list to store results for all element types
all_results <- list()

# Get all cohorts 
cohorts <- c( 'Pancan-no-skin-melanoma-lymph', "Liver-HCC", "Bladder-TCC" ,
              "ColoRect-AdenoCA" , "Lymph-BNHL",
              "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
              "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
              "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",
              "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
              "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc")

# Loop over each element type
for (elem_type in elemTypes) {
  
  corrs <- c()
  eMET_corrs <- c()
  
  for (cohort in cohorts) {
    print(cohort)
    
    # Read mutated elements and filter for non-zero mutations
    mutated_elems <- fread(paste0('../../iDriver/extdata/procInput/BMRs/observed/', cohort, '/test_y.tsv'))
    
    mutated_elems <- as.data.frame(mutated_elems[which(mutated_elems$nMut != 0), ])
    mutated_elems <- mutated_elems[!(mutated_elems$binID %in% drivers), ]
    
    # Read full data frame and filter by mutated elements
    full_df <- fread(paste0('../../DriverPower/output/', cohort, '/non0_', cohort, '.result.tsv')) #fread(paste0('../../DriverPower/output/', cohort, '/', cohort, '.result.tsv'))
    full_df <- full_df[which(full_df$binID %in% mutated_elems$binID), ]
    
    # Filter by the current element type
    df <- full_df[grepl(elem_type, full_df$binID), ]
    
    # Calculate observed and predicted rates and their correlation
    obsRate <- df$nMut / df$length
    predRate <- df$nPred / df$length
    # predRate <- (df$ALPHA * df$THETA) / df$R_SIZE
    corr <- cor(obsRate, predRate, method = 'spearman')
    corrs <- c(corrs, corr)
    
    # Read eMET results and extract correlation for the element type
    path_eMET <- paste0(base_dir_eMET, cohort, '/eMET/eMET_assessment.tsv')
    
    eMET_results <- data.frame(fread(path_eMET))
    row_idx <- which(eMET_results$V1 == 'corr_eMET')
    col_idx <- which(colnames(eMET_results) == elem_type)
    eMET_corrs <- c(eMET_corrs, eMET_results[row_idx, col_idx])
  }
  
  # Store results in a data frame for the current element type
  df <- data.frame('cohort' = cohorts, 'corr_DP.non0' = corrs, 'corr_eMET' = eMET_corrs, element_type = elem_type)
  all_results[[elem_type]] <- df
  print('-------------------')
}

# Combine all results into a single data frame
final_results <- do.call(rbind, all_results)

# View the final results
print(final_results)
rownames(final_results) <- 1:nrow(final_results)

dir.create('../../make_features/external/BMR/output/Res_reviewerComments/compare2DP/', recursive = T, showWarnings = F)
fwrite(final_results, file = '../../make_features/external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_non0_woDrivers.tsv', sep = '\t')
