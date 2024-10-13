rm(list = ls())

library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(tidyr)
library(circlize)


ad <- fread('../external/BMR/output/Res_reviewerComments/compare2AD/BMR_comparison_nSamples_woDrivers.tsv')
dp <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_non0_woDrivers.tsv')
all_methods <- left_join(dp, ad, by = c("cohort" ,"corr_eMET" ,"element_type"))
dig <- fread('../external/BMR/output/Res_reviewerComments/compare2Dig/BMR_comparison_Mu__woDrivers.tsv')
dig <- dig[, c("cohort","corr_Dig","element_type"),]
all_methods <- left_join(all_methods, dig, by = c("cohort" ,"element_type"))
baseline <- fread('../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/assessments/all_cohorts_assessments.tsv')


baseline <- baseline[which(!baseline$element.type %in% c('lncrna.ncrna', 'lncrna.promCore')),]

colnames(baseline) <- c( "baseline_corr", "element_type" , "cohort")
all_methods <- left_join(all_methods, baseline, by = c("cohort" ,"element_type"))

colnames(all_methods) <- c('cohort', 'DriverPower', 'eMET', 'element_type', 'ActiveDriverWGS', 'Dig', 'Baseline')

cohorts_toExclude <- c('Thy-AdenoCA', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC', 'CNS-PiloAstro', 'Kidney-ChRCC')
all_methods <- all_methods[which(!all_methods$cohort %in% cohorts_toExclude),]

all_methods$Baseline <- as.numeric(all_methods$Baseline)
all_methods$eMET <- as.numeric(all_methods$eMET)
all_methods$ActiveDriverWGS <- as.numeric(all_methods$ActiveDriverWGS)
all_methods$Dig <- as.numeric(all_methods$Dig)
all_methods$DriverPower <- as.numeric(all_methods$DriverPower)



# Melt the data to long format for ggplot
long_data <- melt(all_methods, 
                  id.vars = c("cohort", "element_type"), 
                  measure.vars = c("Baseline", "eMET", 'ActiveDriverWGS', 'Dig', 'DriverPower'),
                  variable.name = "Method", 
                  value.name = "Value")

# Sample Data
data <- long_data
data$METHOD_elem <- paste0(data$Method, '::', data$element_type)
data <- data[,c('cohort', 'METHOD_elem', 'Value')]
data$cohort <- ifelse(data$cohort == 'Pancan-no-skin-melanoma-lymph', 'Pan-cancer', data$cohort)
data
# Data reshaping
heatmap_data <- reshape2::dcast(data, METHOD_elem ~ cohort, value.var = "Value", fun.aggregate = mean)

elems <- unlist(lapply(strsplit(heatmap_data$METHOD_elem, '::'), function(s){
  s[2]
}))

heatmap_data$element_type <- unlist(lapply(strsplit(heatmap_data$METHOD_elem, '::'), function(s){
  s[2]
}))

heatmap_data$method <- unlist(lapply(strsplit(heatmap_data$METHOD_elem, '::'), function(s){
  s[1]
}))

idx <- apply(heatmap_data, 1, function(s){
  sum(is.na(s))==0}
)
heatmap_data <- heatmap_data[,-1]
heatmap_data <- heatmap_data[idx,]


element_names <- c('enhancers' = 'Enhancers',
                   'gc19_pc.3utr' = '3\' UTR',
                   'gc19_pc.5utr' = '5\' UTR',
                   'gc19_pc.cds' = 'CDS',
                   'gc19_pc.promCore' = 'Core Promoter',
                   'gc19_pc.ss' = 'Splice site'
                   # , 'lncrna.ncrna' = 'lncRNA',
                   # 'lncrna.promCore' = 'lncRNA Promoter'
)


# Annotations
# Map the element types to their corresponding names using element_names
element_type_labels <- element_names[heatmap_data$element_type]

# Update the element type annotation using element names as labels
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



element_types <- heatmap_data$element_type
element_factor <- element_types



# Define a function to print values inside the cells
cell_function <- function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
}


method_annotation <- HeatmapAnnotation(
  Method = heatmap_data$method,
  which = "row",
  col = list(Method = c( DriverPower = '#386cb0',
                         DP_non0 = 'black',
                         DP_with0 = 'orange',
                         eMET = '#e41a1c', 
                         ActiveDriverWGS = '#a6dba0',
                         Baseline = 'skyblue',
                         Dig = "yellow"))
)

# Map the element types to their corresponding names using element_names
element_factor <- factor(heatmap_data$element_type, 
                         levels = names(element_names), 
                         labels = element_names)


define_cohort_order <- function(path_donorInfo, path_pcawg_supp, based_on){
  donorInfo <- fread(path_donorInfo)
  donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),]
  
  sample_sheet <- fread(path_pcawg_supp)
  cohorts_nMut <- donorInfo[,c('cohort1', 'freq')]
  cohorts_nMut <- cohorts_nMut[!duplicated(cohorts_nMut)]
  
  included_cohorts <- c("Liver-HCC", "Bladder-TCC" ,"ColoRect-AdenoCA" , "Lymph-BNHL",
                        "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
                        "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
                        "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",
                        "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
                        "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc")
  
  
  if (based_on == 'N') {
    cohorts_nSample <- data.frame(table(donorInfo$cohort1))
    cohorts_nSample <- cohorts_nSample[which(cohorts_nSample$Var1 %in% included_cohorts),]
    cohorts_nSample <- cohorts_nSample[order(cohorts_nSample$Freq, decreasing = T),]
    #""
    vec <- c(2253, cohorts_nSample$Freq)
    cohorts_order <- c("Pan-cancer", as.character(cohorts_nSample$Var1))
    
  } else if (based_on == 'nMut') {
    # Summing up frequencies per cohort
    cohorts_nMut <- data.frame(cohorts_nMut %>%
                                 group_by(cohort1) %>%
                                 summarise(total_freq = sum(freq)))
    cohorts_nMut <- cohorts_nMut[which(cohorts_nMut$cohort1 %in% included_cohorts),]
    cohorts_nMut <- cohorts_nMut[order(cohorts_nMut$total_freq, decreasing = T),]
    cohorts_order <- c("Pan-cancer", as.character(cohorts_nMut$cohort1))
    vec <- c(21183078, cohorts_nMut$total_freq)
    
  }
  
  list(cohorts_order, vec)
}



# Define the custom column order
path_donorInfo <- '../../iDriver/extdata/procInput/iDriverInputs/donorInfo.tsv'
path_pcawg_supp <- '../../pcawg_sample_sheet.tsv'

donorInfo <- fread(path_donorInfo)
donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),]


# Step 1: Calculate the number of mutations and donors for each cohort
cohort_summary <- donorInfo %>%
  group_by(cohort1) %>%
  summarise(nMut =  sum(freq),  # Count of mutations (assuming each row represents a mutation)
            nDonors = length(unique(D_id)))  # Count of unique donors

cohort_summary <- rbind(c('Pan-cancer', 21183078, 2253), cohort_summary)
cohort_summary$nMut <- as.numeric(cohort_summary$nMut)
cohort_summary$nDonors <- as.numeric(cohort_summary$nDonors)
cohort_summary$perDonorMuts <- cohort_summary$nMut/cohort_summary$nDonors


included_cohorts <- c('Pan-cancer', "Liver-HCC", "Bladder-TCC" ,"ColoRect-AdenoCA" , "Lymph-BNHL",
                      "Uterus-AdenoCA" , "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
                      "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
                      "Breast-AdenoCa" , "Biliary-AdenoCA", "Eso-AdenoCa",
                      "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
                      "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc")


cohort_summary <- cohort_summary[which(cohort_summary$cohort1 %in% included_cohorts),]


# Ensure that cohorts are in the correct order as per your heatmap
cohort_summary <- cohort_summary[order(cohort_summary$perDonorMuts, decreasing = T),]

nMut_scaled <- scales::rescale(cohort_summary$nMut)
nDonors_scaled <- scales::rescale(cohort_summary$nDonors)
perDonorMuts_scaled <- scales::rescale(cohort_summary$perDonorMuts)

METHODs = heatmap_data$method
# Reorder columns of the matrix according to custom_order
heatmap_data <- heatmap_data[, cohort_summary$cohort1]



# Create a top annotation with gradient-filled boxes for nMut and nDonor
top_annotation <- HeatmapAnnotation(
  nMut = anno_simple(nMut_scaled, 
                     col = colorRamp2(c(1.000000000, 0.116781036 ,
                                        0.071816182, 0.042586189, 0.031538380, 0.027397082,
                                        0.005543080, 0.004752904, 0.002738043, 0),
                                      rev(c('#fcfbfd','#efedf5','#dadaeb','#bcbddc',
                                            '#9e9ac8','#807dba','#6a51a3','#4a1486', '#3f007d', 'black'))), 
                     border = TRUE),
  nDonors = anno_simple(nDonors_scaled, 
                        col = colorRamp2(c(1.000000000, 0.032286996, 0.093721973,
                                           0.039013453, 0.014349776, 0.008968610,
                                           0.006726457, 0.026008969, 0.052914798, 0), 
                                         rev(c('#f7fcfd','#e5f5f9','#ccece6',
                                               '#99d8c9','#66c2a4','#41ae76',
                                               '#238b45','#005824', '#00441b', 'black'))),
                        border = TRUE),
  perDonor_nMut = anno_simple(perDonorMuts_scaled, 
                        col = colorRamp2(c(1, 0.8, 0.7, 0.6, 0.4,0.3, 0.2, 0.1, 0), 
                                         rev(c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb',
                                               '#41b6c4','#1d91c0','#225ea8','#253494',
                                               '#081d58'))),
                        border = TRUE),
  annotation_width = unit(c(1, 1), "cm"),  # Adjust width as needed
  which = "col"  # Specify that this annotation is for columns
)



log_perDonorMuts <- round(log10(cohort_summary$perDonorMuts), 2)

top_annotation <- HeatmapAnnotation(
  nMut = anno_simple(nMut_scaled, 
                     col = colorRamp2(c(1.000000000, 0.116781036 ,
                                        0.071816182, 0.042586189, 0.031538380, 0.027397082,
                                        0.005543080, 0.004752904, 0.002738043, 0),
                                      rev(c('#fcfbfd','#efedf5','#dadaeb','#bcbddc',
                                            '#9e9ac8','#807dba','#6a51a3','#4a1486', '#3f007d', 'black'))), 
                     border = TRUE),
  nDonors = anno_simple(nDonors_scaled, 
                        col = colorRamp2(c(1.000000000, 0.032286996, 0.093721973,
                                           0.039013453, 0.014349776, 0.008968610,
                                           0.006726457, 0.026008969, 0.052914798, 0), 
                                         rev(c('#f7fcfd','#e5f5f9','#ccece6',
                                               '#99d8c9','#66c2a4','#41ae76',
                                               '#238b45','#005824', '#00441b', 'black'))),
                        border = TRUE),
  perDonorMuts = anno_simple(perDonorMuts_scaled, 
                             col = colorRamp2(c(0, .05, .1, .15, .2, .25, .3, .4, .5, .6,  1), 
                                              c('white', '#ffffd9','#edf8b1','#f7fcb9', '#c7e9b4','#addd8e', '#7fcdbb',
                                                      '#41b6c4','#1d91c0','#225ea8','#253494')),
                             border = TRUE),
  log_perDonorMuts = anno_text(log_perDonorMuts, 
                               location = 2,  # Adjust text position
                               just = "center", # Center the text
                               gp = gpar(fontsize = 10),  # Customize font size
                               rot = 0),  # Set rotation to 0 for horizontal text
  annotation_width = unit(c(1, 1), "cm"),  # Adjust width as needed
  which = "col"  # Specify that this annotation is for columns
)


# Step 3: Update your heatmap drawing code to include the top annotation
ht <- Heatmap(as.matrix(heatmap_data), 
              name = "Value",                         
              row_split = element_factor,             
              row_labels = METHODs,       
              show_row_names = F,                    
              show_column_names = TRUE,              
              column_names_rot = 80,                 
              row_title_rot = 0,                      
              cluster_rows = FALSE,                   
              cluster_columns = FALSE,                
              cell_fun = cell_function,              
              left_annotation = element_type_annotation,  
              right_annotation = method_annotation,
              top_annotation = top_annotation,  # Add the top annotation here
              heatmap_legend_param = list(title = "Correlation", 
                                          legend_direction = "horizontal"))


# Draw the heatmap with the top annotation
# draw(ht, 
#      annotation_legend_side = "bottom",    
#      heatmap_legend_side = "bottom")       

# Save the heatmap as PNG
png("../external/BMR/plots/heatmap_cancerTypes_tools.png",
    width = 10, height = 9, units = "in", res = 300)
draw(ht, 
     annotation_legend_side = "bottom",    # Place annotation legends at the bottom
     heatmap_legend_side = "bottom")       # Place heatmap legend at the bottom
dev.off()
