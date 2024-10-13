rm(list = ls())

library(ComplexHeatmap)
library(circlize)
library(data.table)

import_allCohorts <- function(paths, ass_type, cohorts){
  cohorts <- ifelse(cohorts == 'Pancan-no-skin-melanoma-lymph', 'Pan-cancer', cohorts)
  SNVs <- lapply(paths, fread)
  SNVs <- do.call(rbind, lapply(SNVs, function(s){
    s[grepl(ass_type, s$V1),]
  }))
  
  
  method <- unique(ifelse(grepl('eMET', SNVs$V1), 'eMET', 'Intergenic'))
  SNVs <- SNVs[,-c('V1')]
  if (method == 'Intergenic') {
    SNVs <- SNVs[,-c('train')]
  }
  
  SNVs <- t(SNVs)
  colnames(SNVs) <- cohorts
  
  rownames(SNVs) <- paste0(rownames(SNVs), '__', method)
  SNVs
}

save_heatmap_SNV_nonSNV <- function(data, variant, path_save = '../external/BMR/plots/'){
  
  mat_data <- as.matrix(data)
  elems <- unlist(lapply(strsplit(rownames(mat_data), '__'), function(s){
    s[1]
  }))
  methods <- unlist(lapply(strsplit(rownames(mat_data), '__'), function(s){
    s[2]
  }))
  
  # Define colors for the heatmap
  col_fun <- colorRamp2(c(0.5, 0.7, 0.9), c('blue','white','red')) 
  
  
  method_annotation <- HeatmapAnnotation(
    Method = methods,
    which = "row",
    col = list(Method = c( Intergenic = '#386cb0',
                           eMET = '#e41a1c'))
  )
  
  element_names <- c('enhancers' = 'Enhancers',
                     'gc19_pc.3utr' = '3\' UTR',
                     'gc19_pc.5utr' = '5\' UTR',
                     'gc19_pc.cds' = 'CDS',
                     'gc19_pc.promCore' = 'Core Promoter',
                     'gc19_pc.ss' = 'Splice site'
                     # , 'lncrna.ncrna' = 'lncRNA',
                     # 'lncrna.promCore' = 'lncRNA Promoter'
  )
  
  
  # Map the element types to their corresponding names using element_names
  element_type_labels <- element_names[elems]
  
  element_type_annotation <- HeatmapAnnotation(
    `element type` = element_type_labels,        # Use descriptive element names
    which = "row",
    col = list(`element type` = c("3' UTR" = "#14b556", 
                                  "CDS" = "#CC99BB",
                                  "5' UTR" = "#AAAA44", 
                                  "Core Promoter" = "#77AADD",
                                  "Splice site" = "#11266d", 
                                  "Enhancers" = '#b03914'))
  )
  
  cell_function <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat_data[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
  }
  
  # Create the heatmap
  ht <- Heatmap(mat_data, 
                name = "Correlation", 
                col = col_fun, 
                row_split = elems,
                row_labels = methods,
                show_row_names = F, 
                show_column_names = TRUE, 
                left_annotation = element_type_annotation,  
                right_annotation = method_annotation,
                column_names_rot = 80,
                cluster_rows = FALSE,                   
                cluster_columns = FALSE,
                cell_fun = cell_function,
                row_title_rot = 0,
                heatmap_legend_param = list(title = "Correlation", 
                                            legend_direction = "horizontal")
  )
  
  draw(ht, 
       annotation_legend_side = "bottom",    
       heatmap_legend_side = "bottom") 
  
  png(paste0(path_save, variant, "_heatmap.png"),
      width = 10, height = 6.5, units = "in", res = 300)
  draw(ht, 
       annotation_legend_side = "bottom",    # Place annotation legends at the bottom
       heatmap_legend_side = "bottom")       # Place heatmap legend at the bottom
  dev.off()
  
}


######################################
ass_type <- 'corr'
cohorts <- c( "Pancan-no-skin-melanoma-lymph","Liver-HCC", "Bladder-TCC" ,"ColoRect-AdenoCA" , "Lymph-BNHL",
                        "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
                        "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
                        "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",
                        "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
                        "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc") #

files_snv_eMET <- paste0('../external/BMR/output/Res_reviewerComments/SNV_nonSNV/', cohorts, '/SNV/eMET/eMET_assessment.tsv')
files_snv_intergenic <- paste0('../external/BMR/output/Res_reviewerComments/SNV_nonSNV/', cohorts, '/SNV/GBM/GBM_assessments.tsv')





SNVs_eMET <- import_allCohorts(files_snv_eMET, ass_type, cohorts)
SNVs_intergenic <- import_allCohorts(files_snv_intergenic, ass_type, cohorts)

data <- rbind(SNVs_eMET, SNVs_intergenic)

files_nonsnv_eMET <- paste0('../external/BMR/output/Res_reviewerComments/SNV_nonSNV/', cohorts, '/nonSNV/eMET/eMET_assessment.tsv')
files_nonsnv_intergenic <- paste0('../external/BMR/output/Res_reviewerComments/SNV_nonSNV/', cohorts, '/nonSNV/GBM/GBM_assessments.tsv')

nonSNVs_eMET <- import_allCohorts(files_nonsnv_eMET, ass_type, cohorts)
nonSNVs_intergenic <- import_allCohorts(files_nonsnv_intergenic, ass_type, cohorts)

data2 <- rbind(nonSNVs_eMET, nonSNVs_intergenic)

#######################################



save_heatmap_SNV_nonSNV(data, variant = 'SNV')
save_heatmap_SNV_nonSNV(data2, variant = 'nonSNV')
