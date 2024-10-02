################################## stack_barPlot (baseline: filter just block)###########################################
rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)

ad <- fread('../external/BMR/output/Res_reviewerComments/compare2AD/BMR_comparison_nSamples.tsv')
dp <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_non0.tsv')
all_methods <- left_join(dp, ad, by = c("cohort" ,"corr_eMET" ,"element_type"))
dig <- fread('../external/BMR/output/Res_reviewerComments/compare2Dig/BMR_comparison_Mu.tsv')
dig <- dig[, c("cohort","corr_Dig","element_type"),]
all_methods <- left_join(all_methods, dig, by = c("cohort" ,"element_type"))


cohorts <- c('Biliary-AdenoCA', 'Bladder-TCC', 'Bone-Cart', 'Bone-Epith',
             'Bone-Leiomyo', 'Bone-Osteosarc', 'Breast-AdenoCa', 'Breast-DCIS',
             'Breast-LobularCa', 'Cervix-AdenoCA', 'Cervix-SCC', 'CNS-GBM',
             "CNS-Medullo", "CNS-Oligo", "CNS-PiloAstro", "ColoRect-AdenoCA",             
             "Eso-AdenoCa", "Head-SCC", "Kidney-ChRCC", "Kidney-RCC",                    
             "Liver-HCC", "Lung-AdenoCA" , "Lung-SCC", "Lymph-BNHL",                    
             "Lymph-CLL", "Lymph-NOS", "Myeloid-AML", "Myeloid-MPN",                  
             "Ovary-AdenoCA", "Panc-AdenoCA", "Panc-Endocrine",               
             "Pancan-no-skin-melanoma-lymph", "Prost-AdenoCA", "Skin-Melanoma",                
             "Stomach-AdenoCA", "Thy-AdenoCA", "Uterus-AdenoCA")

baseline <- data.frame('baseline_corr' = NA, 'element.type' = NA, 'cohort' = NA)
for (cohort in cohorts) {
  
  x <- fread(paste0('../external/BMR/output/Res_reviewerComments/baseLine_localRates/', cohort, '_baselineAssessment.tsv'))
  x$cohort <- cohort
  baseline <- rbind(baseline, x)
}


baseline <- baseline[which(!is.na(baseline$cohort)),]
baseline <- baseline[which(!baseline$element.type %in% c('lncrna.ncrna', 'lncrna.promCore')),]

colnames(baseline) <- c( "baseline_corr", "element_type" , "cohort" )


all_methods <- left_join(all_methods, baseline, by = c("cohort" ,"element_type"))
all_methods$DP <- all_methods$corr_eMET - all_methods$corr_DP.non0
all_methods$AD <- all_methods$corr_eMET - all_methods$corr_AD
all_methods$Dig <- all_methods$corr_eMET - all_methods$corr_Dig
all_methods$baseline <- all_methods$corr_eMET - all_methods$baseline_corr

all_methods <- all_methods[,c('cohort', 'element_type', "DP", "AD" , "Dig", "baseline")]


elem <- 'gc19_pc.3utr'
data <- all_methods[grep(elem, all_methods$element_type),]


data_long <- gather(data, group, value, DP:baseline) %>%
  arrange(factor(cohort, levels = (data$cohort))) %>% 
  mutate(cohort=factor(cohort, levels=unique(cohort)))

library(ggplot2)

# plot
ggplot(data_long, aes(fill=group, y=value, x=cohort)) + 
  geom_bar(position="stack", stat="identity")

create_bench_colours <- function() {
  
  c( DP = '#386cb0',
     eMET = '#e41a1c', 
     AD = '#a6dba0',
     baseline = 'skyblue',
     Dig = "yellow" 
  )
  
}

for (elem in unique(all_methods$element_type)) {
  print(elem)
  data <- all_methods[grep(elem, all_methods$element_type),]
  
  
  data_long <- gather(data, group, value, DP:baseline) %>%
    arrange(factor(cohort, levels = (data$cohort))) %>% 
    mutate(cohort=factor(cohort, levels=unique(cohort)))
  
  print(ggplot(data_long, aes(fill=group, y=value, x=cohort)) + 
          geom_bar(position="stack", stat="identity")+
          scale_fill_manual(values = create_bench_colours()) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                # legend.title = element_blank(),
                text = element_text(size = 15),
                axis.text.x = element_text(angle = 65, hjust = 1),
                axis.line = element_line(colour = "black")))
}
##########################stack_barPlot (baseline: filter all blocks of an elem) ######################################
rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)

ad <- fread('../external/BMR/output/Res_reviewerComments/compare2AD/BMR_comparison_nSamples.tsv')
dp <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_non0.tsv')
all_methods <- left_join(dp, ad, by = c("cohort" ,"corr_eMET" ,"element_type"))
dig <- fread('../external/BMR/output/Res_reviewerComments/compare2Dig/BMR_comparison_Mu.tsv')
dig <- dig[, c("cohort","corr_Dig","element_type"),]
all_methods <- left_join(all_methods, dig, by = c("cohort" ,"element_type"))


cohorts <- c('Biliary-AdenoCA', 'Bladder-TCC', 'Bone-Cart', 'Bone-Epith',
             'Bone-Leiomyo', 'Bone-Osteosarc', 'Breast-AdenoCa', 'Breast-DCIS',
             'Breast-LobularCa', 'Cervix-AdenoCA', 'Cervix-SCC', 'CNS-GBM',
             "CNS-Medullo", "CNS-Oligo", "CNS-PiloAstro", "ColoRect-AdenoCA",             
             "Eso-AdenoCa", "Head-SCC", "Kidney-ChRCC", "Kidney-RCC",                    
             "Liver-HCC", "Lung-AdenoCA" , "Lung-SCC", "Lymph-BNHL",                    
             "Lymph-CLL", "Lymph-NOS", "Myeloid-AML", "Myeloid-MPN",                  
             "Ovary-AdenoCA", "Panc-AdenoCA", "Panc-Endocrine",               
             "Pancan-no-skin-melanoma-lymph", "Prost-AdenoCA", "Skin-Melanoma",                
             "Stomach-AdenoCA", "Thy-AdenoCA", "Uterus-AdenoCA")

baseline <- fread('../external/BMR/output/Res_reviewerComments/baseLine_localRates/input2_filter_allBlocks_of_the_elem/all_cohorts_assessments.tsv')


baseline <- baseline[which(!is.na(baseline$cohort)),]
baseline <- baseline[which(!baseline$element.type %in% c('lncrna.ncrna', 'lncrna.promCore')),]

colnames(baseline) <- c( "baseline_corr", "element_type" , "cohort")


all_methods <- left_join(all_methods, baseline, by = c("cohort" ,"element_type"))
all_methods$DP <- all_methods$corr_eMET - all_methods$corr_DP.non0
all_methods$AD <- all_methods$corr_eMET - all_methods$corr_AD
all_methods$Dig <- all_methods$corr_eMET - all_methods$corr_Dig
all_methods$baseline <- all_methods$corr_eMET - all_methods$baseline_corr

all_methods <- all_methods[,c('cohort', 'element_type', "DP", "AD" , "Dig", "baseline")]


# elem <- 'gc19_pc.3utr'
# data <- all_methods[grep(elem, all_methods$element_type),]
# 
# 
# data_long <- gather(data, group, value, DP:baseline) %>%
#   arrange(factor(cohort, levels = (data$cohort))) %>% 
#   mutate(cohort=factor(cohort, levels=unique(cohort)))
# 
library(ggplot2)

# # plot
# ggplot(data_long, aes(fill=group, y=value, x=cohort)) + 
#   geom_bar(position="stack", stat="identity")

create_bench_colours <- function() {
  
  c( DP = '#386cb0',
     eMET = '#e41a1c', 
     AD = '#a6dba0',
     baseline = 'skyblue',
     Dig = "yellow" 
  )
  
}

for (elem in unique(all_methods$element_type)) {
  print(elem)
  data <- all_methods[grep(elem, all_methods$element_type),]
  
  
  data_long <- gather(data, group, value, DP:baseline) %>%
    arrange(factor(cohort, levels = (data$cohort))) %>% 
    mutate(cohort=factor(cohort, levels=unique(cohort)))
  
  print(ggplot(data_long, aes(fill=group, y=value, x=cohort)) + 
          geom_bar(position="stack", stat="identity")+
          scale_fill_manual(values = create_bench_colours()) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                # legend.title = element_blank(),
                text = element_text(size = 15),
                axis.text.x = element_text(angle = 65, hjust = 1),
                axis.line = element_line(colour = "black")))
}
######################### Heatmap (baseline: filter all blocks of an elem) ###########################
rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)

ad <- fread('../external/BMR/output/Res_reviewerComments/compare2AD/BMR_comparison_nSamples.tsv')
dp <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_non0.tsv')
all_methods <- left_join(dp, ad, by = c("cohort" ,"corr_eMET" ,"element_type"))
dig <- fread('../external/BMR/output/Res_reviewerComments/compare2Dig/BMR_comparison_Mu.tsv')
dig <- dig[, c("cohort","corr_Dig","element_type"),]
all_methods <- left_join(all_methods, dig, by = c("cohort" ,"element_type"))
baseline <- fread('../external/BMR/output/Res_reviewerComments/baseLine_localRates/input2_filter_allBlocks_of_the_elem/all_cohorts_assessments.tsv')


baseline <- baseline[which(!baseline$element.type %in% c('lncrna.ncrna', 'lncrna.promCore')),]

colnames(baseline) <- c( "baseline_corr", "element_type" , "cohort")


all_methods <- left_join(all_methods, baseline, by = c("cohort" ,"element_type"))
dp0 <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_with0.tsv')
all_methods <- left_join(all_methods, dp0, by = c("cohort" ,"element_type"))
colnames(all_methods) <- c('cohort', 'DP_non0', 'eMET', 'element_type', 'AD', 'Dig', 'baseline', 'DP_with0', 'eMET2')
all_methods <- all_methods[,c('cohort', 'element_type', 'baseline', 'eMET',  'AD', 'Dig',  'DP_with0', 'DP_non0')]
# all_methods <- all_methods[which(all_methods$element_type != 'gc19_pc.ss'),]

cohorts_toExclude <- c('Thy-AdenoCA', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC', 'CNS-PiloAstro', 'Kidney-ChRCC')
all_methods <- all_methods[which(!all_methods$cohort %in% cohorts_toExclude),]

all_methods$baseline <- as.numeric(all_methods$baseline)
all_methods$eMET <- as.numeric(all_methods$eMET)
all_methods$AD <- as.numeric(all_methods$AD)
all_methods$Dig <- as.numeric(all_methods$Dig)
all_methods$DP_with0 <- as.numeric(all_methods$DP_with0)
all_methods$DP_non0 <- as.numeric(all_methods$DP_non0)


library(tidyr)
library(ggplot2)
library(pheatmap)

# Assuming your data is in 'all_methods' data.table

# Step 1: Split data by element_type
split_data <- split(all_methods, all_methods$element_type)


# Step 2: Loop through each element_type and create heatmaps
for (element in names(split_data)) {
  data_element <- split_data[[element]]
  
  data_element_long <- pivot_longer(data_element, 
                                    cols = c(baseline, eMET, AD, Dig, DP_with0, DP_non0), 
                                    names_to = "method", 
                                    values_to = "value")
  print(ggplot(data_element_long, aes(x = cohort, y = method, fill = value)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(value, 2)), 
              color = "black", size = 3.5)  +  # Add text labels

    scale_fill_gradientn(colors = c("#053061", "#2166AC", "#2166AC", "#4393C3", "#92C5DE" ,
                                    "#D1E5F0", "#E7E1EF", "#D4B9DA" , "#CE1256", "#980043" ),
                         breaks = c(-.1,0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                         labels = c(-.1, 0,.1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                         limits = c(-.1, 1),
                         na.value = "grey50",
                         guide = "colorbar") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15),) +
    #labs(x = "Feature Category", y = "Element type")
    labs(x = "", y = "")
  )
}


library(ggplot2)
library(dplyr)

# Replace NA values with 1000
data_element_long <- data_element_long %>%
  mutate(value = ifelse(is.na(value), 1000, value))

# Plot heatmap with ggplot2
ggplot(data_element_long, aes(x = cohort, y = method, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0.5, na.value = "grey", limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(fill = "Value", 
       x = "Cohort", 
       y = "Method", 
       title = "Heatmap of Methods by Cohort") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3)


##################################### box plots ############################################

rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)

ad <- fread('../external/BMR/output/Res_reviewerComments/compare2AD/BMR_comparison_nSamples.tsv')
dp <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_non0.tsv')
all_methods <- left_join(dp, ad, by = c("cohort" ,"corr_eMET" ,"element_type"))
dig <- fread('../external/BMR/output/Res_reviewerComments/compare2Dig/BMR_comparison_Mu.tsv')
dig <- dig[, c("cohort","corr_Dig","element_type"),]
all_methods <- left_join(all_methods, dig, by = c("cohort" ,"element_type"))
baseline <- fread('../external/BMR/output/Res_reviewerComments/baseLine_localRates/input2_filter_allBlocks_of_the_elem/all_cohorts_assessments.tsv')


baseline <- baseline[which(!baseline$element.type %in% c('lncrna.ncrna', 'lncrna.promCore')),]

colnames(baseline) <- c( "baseline_corr", "element_type" , "cohort")


all_methods <- left_join(all_methods, baseline, by = c("cohort" ,"element_type"))
dp0 <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_with0.tsv')
all_methods <- left_join(all_methods, dp0, by = c("cohort" ,"element_type"))
colnames(all_methods) <- c('cohort', 'DP_non0', 'eMET', 'element_type', 'AD', 'Dig', 'baseline', 'DP_with0', 'eMET2')
all_methods <- all_methods[,c('cohort', 'element_type', 'baseline', 'eMET',  'AD', 'Dig',  'DP_with0', 'DP_non0')]
# all_methods <- all_methods[which(all_methods$element_type != 'gc19_pc.ss'),]

cohorts_toExclude <- c('Thy-AdenoCA', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC', 'CNS-PiloAstro', 'Kidney-ChRCC')
all_methods <- all_methods[which(!all_methods$cohort %in% cohorts_toExclude),]

all_methods$baseline <- as.numeric(all_methods$baseline)
all_methods$eMET <- as.numeric(all_methods$eMET)
all_methods$AD <- as.numeric(all_methods$AD)
all_methods$Dig <- as.numeric(all_methods$Dig)
all_methods$DP_with0 <- as.numeric(all_methods$DP_with0)
all_methods$DP_non0 <- as.numeric(all_methods$DP_non0)


# Melt the data to long format for ggplot
long_data <- melt(all_methods, 
                  id.vars = c("cohort", "element_type"), 
                  measure.vars = c("baseline", "eMET", "AD", "Dig", "DP_with0", "DP_non0"),
                  variable.name = "Method", 
                  value.name = "Value")


# Reorder 'Method' based on the decreasing median value
long_data[, Method := factor(Method, levels = names(sort(tapply(Value, Method, median), decreasing = TRUE)))]

# Create the grouped box plot with unfilled boxes and datapoints
ggplot(long_data, aes(x = Method, y = Value, color = element_type)) +
  geom_boxplot(position = position_dodge(width = 0.8), fill = NA, outlier.shape = NA) +  # Unfilled boxes, remove outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.5, alpha = 0.6) +  # Show datapoints
  facet_wrap(~element_type) +  # Separate plots for each element_type group
  theme_minimal() +
  labs(title = "Grouped Box Plot of Methods by Element Type (Sorted)",
       x = "Method",
       y = "Values",
       color = "Element Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Create the grouped box plot with unfilled boxes and datapoints
ggplot(long_data, aes(x = Method, y = Value, color = Method)) +
  geom_boxplot(position = position_dodge(width = 0.75), fill = NA, outlier.shape = NA) +  # Unfilled boxes, remove outliers
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.2, alpha = 0.6) +  # Show datapoints
  facet_wrap(~element_type, nrow = 2, strip.position = "bottom") +  # Move element-type label to the bottom
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Correlation",
       color = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(color = "black", fill = "white"),  # Add a box around element-type labels
        strip.placement = "outside")  # Ensure the label box is outside the plot area


# Create the grouped box plot with unfilled boxes and datapoints
ggplot(long_data, aes(x = Method, y = Value)) +
  geom_boxplot(aes(color = Method), position = position_dodge(width = 0.75), fill = NA, outlier.shape = NA) +  # Color boxes by method
  geom_jitter(aes(color = cohort), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 1.2, alpha = 0.6) +  # Color points by cohort
  facet_wrap(~element_type, nrow = 2, strip.position = "bottom") +  # Move element-type label to the bottom
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(color = "black", fill = "white"),  # Add a box around element-type labels
        strip.placement = "outside")  # Ensure the label box is outside the plot area

###################################### final heatmap ############################################
# Load the library
rm(list = ls())

library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(tidyr)

ad <- fread('../external/BMR/output/Res_reviewerComments/compare2AD/BMR_comparison_nSamples.tsv')
dp <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_non0.tsv')
all_methods <- left_join(dp, ad, by = c("cohort" ,"corr_eMET" ,"element_type"))
dig <- fread('../external/BMR/output/Res_reviewerComments/compare2Dig/BMR_comparison_Mu.tsv')
dig <- dig[, c("cohort","corr_Dig","element_type"),]
all_methods <- left_join(all_methods, dig, by = c("cohort" ,"element_type"))
baseline <- fread('../external/BMR/output/Res_reviewerComments/baseLine_localRates/input2_filter_allBlocks_of_the_elem/all_cohorts_assessments.tsv')


baseline <- baseline[which(!baseline$element.type %in% c('lncrna.ncrna', 'lncrna.promCore')),]

colnames(baseline) <- c( "baseline_corr", "element_type" , "cohort")


all_methods <- left_join(all_methods, baseline, by = c("cohort" ,"element_type"))
dp0 <- fread('../external/BMR/output/Res_reviewerComments/compare2DP/BMR_comparison_with0.tsv')
all_methods <- left_join(all_methods, dp0, by = c("cohort" ,"element_type"))
colnames(all_methods) <- c('cohort', 'DP_non0', 'eMET', 'element_type', 'AD', 'Dig', 'baseline', 'DP_with0', 'eMET2')
# all_methods <- all_methods[,c('cohort', 'element_type', 'baseline', 'eMET',  'AD', 'Dig',  'DP_with0', 'DP_non0')]
all_methods <- all_methods[,c('cohort', 'element_type', 'baseline', 'eMET',  'AD', 'Dig', 'DP_non0')]
colnames(all_methods) <- c('cohort', 'element_type', 'baseline', 'eMET',  'ActiveDriverWGS', 'Dig', 'DriverPower')

cohorts_toExclude <- c('Thy-AdenoCA', 'Myeloid-MPN', 'CNS-Oligo', 'Cervix-SCC', 'CNS-PiloAstro', 'Kidney-ChRCC')
all_methods <- all_methods[which(!all_methods$cohort %in% cohorts_toExclude),]

all_methods$baseline <- as.numeric(all_methods$baseline)
all_methods$eMET <- as.numeric(all_methods$eMET)
all_methods$ActiveDriverWGS <- as.numeric(all_methods$ActiveDriverWGS)
all_methods$Dig <- as.numeric(all_methods$Dig)
all_methods$DriverPower <- as.numeric(all_methods$DriverPower)



# Melt the data to long format for ggplot
long_data <- melt(all_methods, 
                  id.vars = c("cohort", "element_type"), 
                  measure.vars = c("baseline", "eMET", 'ActiveDriverWGS', 'Dig', 'DriverPower'),
                  variable.name = "Method", 
                  value.name = "Value")

# Sample Data
data <- long_data
data$METHOD_elem <- paste0(data$Method, '::', data$element_type)
data <- data[,c('cohort', 'METHOD_elem', 'Value')]
# Data reshaping
heatmap_data <- reshape2::dcast(data, METHOD_elem ~ cohort, value.var = "Value")
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
# Annotations
element_type_annotation <- HeatmapAnnotation(
  `element type` = heatmap_data$element_type,
  which = "row",
  col = list(`element type` = c("gc19_pc.3utr" = "#14b556", "gc19_pc.cds" = "#CC99BB",
                              "gc19_pc.5utr" = "#AAAA44", "gc19_pc.promCore" = "#77AADD",
                              "gc19_pc.ss" = "#11266d", "enhancers" = '#b03914'))
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
                         baseline = 'skyblue',
                         Dig = "yellow"))
)



# Generate the heatmap with custom cell function and row split by element-type
Heatmap(as.matrix(heatmap_data[, 1:24]), 
        name = "Value",                         # Name of the color legend
        row_split = element_factor,             # Slice rows by element-type
        row_labels = heatmap_data$method,       # Show method names (rows)
        show_row_names = F,                 # Hide row names (optional if adding method names in the cells)
        show_column_names = TRUE,               # Show cohort names (columns)
        row_title_rot = 0,                      # Rotate the row titles
        cluster_rows = FALSE,                   # Disable row clustering
        cluster_columns = FALSE,                # Disable column clustering
        cell_fun = cell_function,               # Function to print values inside the cells
        left_annotation = element_type_annotation,  # Element type annotation
        right_annotation = method_annotation,
        heatmap_legend_param = list(title = "Score", legend_direction = "horizontal"))


# Define the custom column order
path_donorInfo <- '../../iDriver/extdata/procInput/iDriverInputs/donorInfo.tsv'
path_pcawg_supp <- '../../pcawg_sample_sheet.tsv'

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
  
  if (based_on == 'nMuts') {
    cohorts_nSample <- data.frame(table(donorInfo$cohort1))
    cohorts_nSample <- cohorts_nSample[which(cohorts_nSample$Var1 %in% included_cohorts),]
    cohorts_nSample <- cohorts_nSample[order(cohorts_nSample$Freq, decreasing = T),]
    #"Pancan-no-skin-melanoma-lymph"
    
    cohorts_order <- c("Pancan-no-skin-melanoma-lymph", as.character(cohorts_nSample$Var1))
  } else if (based_on == 'signatures') {
    
  }
  
  cohorts_order
}

custom_order <- define_cohort_order(path_donorInfo, path_pcawg_supp, based_on = 'nMuts')

METHODs = heatmap_data$method
# Reorder columns of the matrix according to custom_order
heatmap_data <- heatmap_data[, custom_order]

# Rotate column names by 65 degrees, set element-type annotations and apply the custom column order
Heatmap(as.matrix(heatmap_data), 
        name = "Value",                         
        row_split = element_factor,             
        row_labels = METHODs,       
        show_row_names = F,                    
        show_column_names = TRUE,              
        column_names_rot = 80,                 # Rotate column names by 65 degrees
        row_title_rot = 0,                      
        cluster_rows = FALSE,                   
        cluster_columns = FALSE,                
        cell_fun = cell_function,              
        left_annotation = element_type_annotation,  # Ensure element-type annotation is shown
        right_annotation = method_annotation,
        heatmap_legend_param = list(title = "Correlation", legend_direction = "horizontal"))

############################## tissueSp Ftrs vs all Ftrs #################################
rm(list = ls())

library(data.table)
library(ggplot2)

create_method_colours <- function() {
  
  c( `NN (MSE loss)` = '#e41a1c' , `NN (Poisson loss)` = '#984ea3', GLM = 'grey',
     Intergenic = '#386cb0', XGBoost = '#386cb0', `all features` = '#386cb0',
     `tissue-specific features` = 'lightblue',
     `Variable-size intergenic` = '#386cb0', `Variable size` = '#386cb0',
     eMET = '#e41a1c', `element-specific` = '#cecdc9', 
     RF =  '#4daf4a', 
     `All features` = '#386cb0',
     PCA = '#e7d4e8',
     AE = '#a6dba0',
     # `Mutated` = '#386cb0', 
     # `Mutated and unmutated` = 'lightblue',
     # `Long mutated` = '#e41a1c',
     # `Long mutated and unmutated` = 'pink',
     `#mutations >= 1` = '#386cb0', 
     All = 'lightblue',
     `#mutations >= 1 & Length > 100` = '#e41a1c',
     `Length > 100` = 'pink',
     `10k` = '#d01c8b', `50k` =  '#f1b6da',
     `100k` = "#b8e186", `1M` = "#4dac26" 
  )
  
}


included_cohorts <- c("Liver-HCC", "ColoRect-AdenoCA" ,
                      "Uterus-AdenoCA" , "Kidney-RCC", "Lung-SCC",
                      "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
                      "Breast-AdenoCa",
                      "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
                      "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc")

files <- paste0('../external/BMR/output/Res_reviewerComments/tissueSpFtrs/', included_cohorts, '/GBM/GBM_assessments.tsv')

ass <- lapply(files, fread)

# Add a cohort column to each assessment
ass <- Map(function(dt, cohort) {
  dt[, cohort := cohort]
  return(dt)
}, ass, included_cohorts)

# Check one of the modified assessments
all_assessments <- do.call(rbind, ass)
all_assessments$usedFeatures <- 'tissue-specific features'
ass <- all_assessments[which(all_assessments$V1 == 'corr_GBM'),]
ass


files <- paste0('../external/BMR/output/reviewerComments/', included_cohorts, '/GBM/GBM_assessments.tsv')

ass <- lapply(files, fread)

# Add a cohort column to each assessment
ass <- Map(function(dt, cohort) {
  dt[, cohort := cohort]
  return(dt)
}, ass, included_cohorts)

# Check one of the modified assessments
all_assessments_base <- do.call(rbind, ass)
all_assessments_base$usedFeatures <- 'all features'

ass <- rbind(all_assessments_base, all_assessments)

ass_type = 'corr'
ass_elem <- ass[grepl(ass_type, ass$V1), ]

library(tidyr)

# Reshape data from wide to long format
ass_elem_long <- melt(ass_elem, 
                      id.vars = c("V1", "train", "cohort", "usedFeatures"), 
                      measure.vars = c("enhancers", "gc19_pc.3utr", "gc19_pc.5utr", "gc19_pc.cds", 
                                       "gc19_pc.promCore", "gc19_pc.ss"), 
                      variable.name = "element_type", 
                      value.name = 'ass')
if (ass_type == 'corr') {
  Y_label = 'correlation'
}

element_names <- c('enhancers' = 'Enhancers',
                   'gc19_pc.3utr' = '3\' UTR',
                   'gc19_pc.5utr' = '5\' UTR',
                   'gc19_pc.cds' = 'CDS',
                   'gc19_pc.promCore' = 'Core Promoter',
                   'gc19_pc.ss' = 'Splice site'
                   # , 'lncrna.ncrna' = 'lncRNA',
                   # 'lncrna.promCore' = 'lncRNA Promoter'
)

ass_elem_long$element_type <- factor(ass_elem_long$element_type, 
                                     levels = names(element_names), 
                                     labels = element_names)

# Create a grouped bar plot comparing usedFeatures for each element-type
ggplot(ass_elem_long, aes(x = cohort, y = ass, fill = usedFeatures)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~element_type, nrow = 2, strip.position = "bottom", labeller = label_value) +
  labs(x = "Cohort", y = ass_type, fill = NULL) +
  
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), 
        text = element_text(size = 15),
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(color = "black", fill = "lightgrey"),  # Element names in boxes
        strip.text = element_text(face = "bold", size = 10)) +
  scale_fill_manual(values = create_method_colours()) +  # Set the interior colors
  scale_color_manual(values = "black") +  # Set the border color
  labs(y = Y_label, x = '')


