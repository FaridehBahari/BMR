rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)

if (.Platform$OS.type == "windows") {
  # setwd('A:/myThesis/make_features/BMR/')
  setwd("C:/Active/projects/make_features/BMR/")
}
# source('plots/functions.R')

create_method_colours <- function() {
  c(GBM = '#daa520', RF = "#117744", GLM = "#aaaeba",
    NN_MSEloss = "#DDDD77", nn_mseLoss = "#DDDD77", 
    NN_PoisLoss = "#11266d", nn_poisLoss = "#11266d",
    # `1k` = '#edf8fb', `10k` = '#b3cde3', `50k` = '#8c96c6',
    # `100k` = '#8856a7', `1M` = '#810f7c'
    `1k` = '#fafad2', `10k` = '#f0e68c', `50k` = '#daa520',
    `100k` = '#b8860b', `1M` = '#85754e'
    )
}


create_element_colours <- function() {
  c(enhancers = '#b03914', gc19_pc.3utr = "#14b556", gc19_pc.5utr = "#AAAA44",
    gc19_pc.cds = "#CC99BB", gc19_pc.promCore = "#77AADD", 
    gc19_pc.ss = "#11266d", lncrna.ncrna = '#daa520', lncrna.promCore = "#DDAA77"
  )
}


c("#114477","#77AADD",   "#AAAA44", "#77CCCC", "#771122",
"#DDAA77","#88CCAA",
"#117744","#CC99BB", "#DDDD77", 
"#4477AA","#44AA77",
"#AA7744", "#774411",
"#117777","#AA4455",

"#AA4488", "#771155",
"#777711",


"#44AAAA", 

"#DD7788", "#52000D","#202020"
)
# (fill="#69b3a2", alpha=0.5)


# https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=5


extract_binSize_from_path <- function(path){
  spl <- unlist(strsplit(path, "/"))
  model <- spl[length(spl)-2]
  model
}

extract_ass_methods <- function(path_ass_method, ass_type, element){
  ass <- lapply(path_ass_method, fread)
  
  ass_ele <- lapply(ass, function(x){
    x = data.frame(x)
    row_idx = grep( ass_type, x$V1)
    col_idx = grep( element, colnames(x))
    x[row_idx, col_idx]
  })
  ass_ele_method <- unlist(ass_ele)
  ass_ele_method
}

extract_model_from_path <- function(path){
  spl <- unlist(strsplit(path, "/"))
  model <- spl[length(spl)-1]
  model
}


grouped_barPlot_model_elements <- function(path_ass_full, ass_type, compare){
  elements <- c('enhancers', 'gc19_pc.3utr', 'gc19_pc.5utr',
                'gc19_pc.cds', 'gc19_pc.promCore',
                'gc19_pc.ss', 'lncrna.ncrna', 'lncrna.promCore')
  
  elem_dfs <- list()
  
  for (element in elements) {
    
    assessment <- extract_ass_methods(path_ass_full, ass_type, element)
    
    
    model <- c()
    for (path in path_ass_full) {
      
      
      if(compare == 'per_binSize'){
        m = extract_binSize_from_path(path)
      } else if (compare == 'per_models_varSize') {
        m = extract_model_from_path(path)
      }
      
      model = c(model, m)
      
    }
    
    df <- data.frame(cbind(model, assessment, element))
    elem_dfs <- append(elem_dfs, list(df))
  }
  
  elem_dfs <- bind_rows(elem_dfs)
  
  # Order the data first by element and then by decreasing assessment
  elem_dfs <- elem_dfs %>%
    arrange(element, desc(assessment))
  
  # Create groups for laying out the barplot correctly
  elem_dfs <- elem_dfs %>%
    group_by(element) %>%
    mutate(id = row_number())
  
  # Create the circular or radial barplot
  
  # Convert 'assessment' to numeric
  elem_dfs$assessment <- as.numeric(elem_dfs$assessment)
  
  ggplot(elem_dfs, aes(x = element, y = assessment, fill = model, group = interaction(model, id))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    # coord_polar(start = 0) +  # This turns the barplot into a circular barplot
    ylim(0, 1) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),
          legend.position = "right",
          panel.grid = element_blank()) +
    scale_fill_manual(values = create_method_colours()) +
    labs(title = "Barplot of Model Assessments by Element",
         y = ass_type,
         fill = "Model")
}

plot_validation_boxplot <- function(directory_paths, metric, compare) {
  all_data <- list()
  
  for (directory_path in directory_paths) {
    files <- list.files(directory_path)
    
    for (file in files) {
      if (grepl("_assessment.tsv$", file)) {
        metric_data <- read.table(file.path(directory_path, file), sep='\t',
                                  header=FALSE, row.names=1) %>%
          t() %>%
          as.data.frame()
        
        if(compare == 'per_binSize'){
          model = extract_binSize_from_path(directory_path)
        } else if (compare == 'per_models_varSize') {
          model = extract_model_from_path(directory_path)
        }
        
        metric_name <- grep(metric, colnames(metric_data), value = TRUE)
        metric_data <- metric_data %>%
          mutate(Model = model, Value = as.numeric(metric_data[, metric_name]))
        
        all_data[[length(all_data) + 1]] <- select(metric_data, Model, Value)
      }
    }
  }
  
  if (length(all_data) > 0) {
    combined_data <- bind_rows(all_data)
    
    ordered <- unique(combined_data %>%
                        arrange(desc(Value)) %>%
                        pull(Model))
    
    
    # Use the ordered metrics as a factor level for model
    combined_data$Model <- factor(combined_data$Model, levels = ordered)
    
    
    
    # Use ggplot2 to plot the grouped box plot
    ggplot(combined_data, aes(x = Model, y = Value, fill = Model)) +
      geom_boxplot(position = position_dodge(width = 0.6), width = 0.5) +  # Adjust the width as needed
      labs(y = metric, title = " performance of different models on intergenic validation sets using variable-size bins ") +
      
      theme(legend.position="right")+
      scale_fill_manual(values = create_method_colours()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    print(paste("No data found for metric '", metric, "' in the specified directories."))
  }
}




DownSampling_linePlot_perModel_allElems <- function(path_ass_DS, ass_type){
  model <- extract_model_from_path(path_ass_DS[1])
  elements <- c('enhancers', 'gc19_pc.3utr', 'gc19_pc.5utr',
                'gc19_pc.cds', 'gc19_pc.promCore',
                'gc19_pc.ss', 'lncrna.ncrna', 'lncrna.promCore')
  
  elem_dfs <- list()
  
  for (element in elements) {
    
    assessment <- extract_ass_methods(path_ass_DS, ass_type, element)
    
    
    sampleSize <- c()
    for (path in path_ass_DS) {
      m = extract_binSize_from_path(path)
      sampleSize = c(sampleSize, m)
      
    }
    
    df <- data.frame(cbind(sampleSize, assessment, element))
    elem_dfs <- append(elem_dfs, list(df))
  }
  
  elem_dfs <- bind_rows(elem_dfs)
  
  
  library(ggplot2)
  
  # Convert the sampleSize column to a factor to ensure correct ordering in the plot
  elem_dfs$sampleSize <- factor(elem_dfs$sampleSize, levels = c("FullSet", "DS1M", "DS800k", "DS600k", "DS300k", "DS100k", "DS50k"))
  
  # Create the line plot
  ggplot(data = elem_dfs, aes(x = sampleSize, y = assessment, 
                              group = element, color = element)) +
    geom_line() +
    geom_point() +
    labs(title = "Reduction of Sample Size Monitoring",
         x = "Sample Size",
         y = ass_type) +
    theme_minimal()
  
  
  
  # Create a new dataframe containing only the "FullSet" and "DS50k" assessments for element
  full_ds50 <- subset(elem_dfs, sampleSize %in% c("FullSet", "DS600k", "DS50k"))
  full_ds50$assessment <- as.character(round(as.numeric(full_ds50$assessment), 3))
  
  
  # Plot the line plot with "FullSet" and "DS50k" assessments for each element
  ggplot(elem_dfs, aes(x = sampleSize, y = assessment, group = element, color = element)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = create_element_colours()) +
    geom_text(data = full_ds50, aes(label = assessment), hjust = -0.2, vjust = 0) +
    labs(title = paste0("Reduction of Sample Size Monitoring in Different Elements using ", model, " model"),
         x = "Sample Size",
         y = ass_type) +
    theme_minimal() +
    theme(axis.title.y = element_text(),
          axis.text.y = element_blank(),  # Remove y-axis values
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))  # Remove y-axis ticks
  
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#         Compare different models on variable-size intergenic bins
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
############## elements:
path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/NN_PoisLoss/NN_PoisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/NN_MSEloss/NN_MSEloss_assessments.tsv")


ass_type <- 'corr'
grouped_barPlot_model_elements(path_ass_full, ass_type, 'per_models_varSize')



############## validation Sets:

directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/rep_train_test/')


plot_validation_boxplot(directory_paths, 'corr', 'per_models_varSize') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#     bin size effect figures (fixed-concatenated intervals) ...GBM applied
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_binSizes <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/GBM_assessments.tsv'
                       # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1k/GBM/GBM_assessments.tsv'
)

ass_type <- 'corr'
grouped_barPlot_model_elements(path_ass_binSizes, ass_type, 'per_binSize')


binEffect_directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/rep_train_test/'
                               # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1k/GBM/rep_train_test/'
)

plot_validation_boxplot(binEffect_directory_paths, 'corr', 'per_binSize') 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#     bin size effect figures (fixed-window intervals)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_binSizes <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/1M/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/100k/GBM/GBM_assessments.tsv',
                       # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/50k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/10k/GBM/GBM_assessments.tsv'
                       # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/1k/GBM/GBM_assessments.tsv'
)

ass_type <- 'corr'
grouped_barPlot_model_elements(path_ass_binSizes, ass_type, 'per_binSize')


binEffect_directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/',
                               # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/rep_train_test/'
                               # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1k/GBM/rep_train_test/'
)

plot_validation_boxplot(binEffect_directory_paths, 'corr', 'per_binSize') 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                       dimension reduction figures
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/NN_PoisLoss/NN_PoisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/NN_MSEloss/NN_MSEloss_assessments.tsv")

path_ass_dimReds_pca <- c("../external/BMR/output/dimReduction_effect/PCA/GBM/GBM_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/RF/RF_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/NN_PoisLoss/NN_PoisLoss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/NN_MSEloss/NN_MSEloss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/GLM/GLM_PCA_rm_nonMutated_assessments.tsv")

path_ass_dimReds_ae <- c("../external/BMR/output/dimReduction_effect/AE/GBM/GBM_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/AE/RF/RF_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/AE/NN_PoisLoss/NN_PoisLoss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/AE/NN_MSEloss/NN_MSEloss_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/GLM/GLM_PCA_rm_nonMutated_assessments.tsv")


ass_type <- 'corr'
element <- 'enhancer'
original_data <- extract_ass_methods(path_ass_full, ass_type, element)
PCA <- extract_ass_methods(path_ass_dimReds_pca, ass_type, element)
AE <- extract_ass_methods(path_ass_dimReds_ae, ass_type, element)

df <- data.frame(cbind(PCA, AE, original_data))


df$model <- c("GBM", "RF", "NN_PoisLoss", "NN_MSEloss", "GLM")
df[which(df$model == 'GLM'), which(colnames(df) == 'original_data')] = NA

library(dplyr)
library(tidyr)
library(ggplot2)

# pivot_longer to create the long format dataframe as you've done before
df_long <- tidyr::pivot_longer(df, cols = c("PCA", "AE", "original_data"), names_to = "method")

# Create a new column for the combination of model and method, if necessary
df_long$model_method <- paste(df_long$model, df_long$method, sep = "_")

# Arrange the models by decreasing PCA value
pca_ordered <- df_long %>%
  filter(method == "PCA") %>%
  arrange(desc(value)) %>%
  pull(model)

# Use the ordered PCA as a factor level for model
df_long$model <- factor(df_long$model, levels = pca_ordered)

# Set the levels for the method so that the bars will be ordered as original, PCA, AE within each group
df_long$method <- factor(df_long$method, levels = c("original_data", "PCA", "AE"))

# Create the ggplot 
ggplot(df_long, aes(x = model, y = value, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black") +
  labs(title = "The Effect of Dimension Reduction Using AE and PCA on each model",
       x = "Model",
       y = paste0(element, ' ', ass_type )) +
  scale_fill_manual(values = c(
    'original_data' = '#446655', # You need to specify colors for 'original_data', PCA and AE
    'PCA' = '#9a9a53',
    'AE' = '#e7e79f'
    # Add other colors if needed
  )) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#########################################################
library(tidyverse)

# Function to extract method and model from directory path
extract_info_from_path <- function(directory_path) {
  parts <- strsplit(directory_path, '/')[[1]]
  method <- parts[which(parts == 'dimReduction_effect') + 1]
  model <- parts[which(parts == 'dimReduction_effect') + 2]
  return(c(method = method, model = model))
}

# Function to plot grouped box plot for a metric
plot_metric_boxplot <- function(directory_paths, metric) {
  all_data <- list()
  
  for (directory_path in directory_paths) {
    files <- list.files(directory_path)
    
    for (file in files) {
      if (grepl("_assessment.tsv$", file)) {
        metric_data <- read.table(file.path(directory_path, file), sep='\t',
                                  header=FALSE, row.names=1) %>%
          t() %>%
          as.data.frame()
        
        info <- extract_info_from_path(directory_path)
        
        metric_name <- grep(metric, colnames(metric_data), value = TRUE)
        metric_data <- metric_data %>%
          mutate(Method = info['method'], Model = info['model'], Value = as.numeric(metric_data[, metric_name]))
        
        all_data[[length(all_data) + 1]] <- select(metric_data, Model, Method, Value)
      }
    }
  }
  
  if (length(all_data) > 0) {
    combined_data <- bind_rows(all_data)
    
    pca_ordered <- unique(combined_data %>%
      filter(Method == "PCA") %>%
      arrange(desc(Value)) %>%
      pull(Model))
    
    
    # Use the ordered PCA as a factor level for model
    combined_data$Model <- factor(combined_data$Model, levels = pca_ordered)
    
    # Set the levels for the method so that the bars will be ordered as original, PCA, AE within each group
    combined_data$Method <- factor(combined_data$Method, levels = c("original", "PCA", "AE"))
    
    
    
    # Use ggplot2 to plot the grouped box plot
    ggplot(combined_data, aes(x = Model, y = Value, fill = Method)) +
      geom_boxplot(position = position_dodge(width = 0.6), width = 0.5) +  # Adjust the width as needed
      labs(y = metric, title = paste0(metric, " Distribution Comparison between Models using different dimension reduction Methods")) +
      
      theme(legend.position="top")+
      scale_fill_manual(values = c(
        'original' = '#446655', # You need to specify colors for 'original_data', PCA and AE
        'PCA' = '#9a9a53',
        'AE' = '#e7e79f'
        # Add other colors if needed
      )) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    print(paste("No data found for metric '", metric, "' in the specified directories."))
  }
}

directory_paths_PCA <- c(
  "../external/BMR/output/dimReduction_effect/PCA/GBM/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/RF/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/NN_PoisLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/NN_MSEloss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/GLM/rep_train_test/"
)

metric <- 'corr'

directory_paths_AE <- c(
  "../external/BMR/output/dimReduction_effect/AE/GBM/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/RF/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/NN_PoisLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/NN_MSEloss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/GLM/rep_train_test/"
)

directory_paths_orig <- c(
  "../external/BMR/output/dimReduction_effect/original/GBM/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/original/RF/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/original/NN_PoisLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/original/NN_MSEloss/rep_train_test/"
)

plot_metric_boxplot(c(directory_paths_PCA, directory_paths_AE, directory_paths_orig), metric)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                       Down sampling figures
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_DS <- c('../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/GBM/GBM_assessments.tsv')

ass_type = 'corr'

DownSampling_linePlot_perModel_allElems(path_ass_DS, ass_type)

# geom_errorbar( aes(ymin = len-sd, ymax = len+sd),width = 0.2)
# https://r-graph-gallery.com/104-plot-lines-with-error-envelopes-ggplot2.html


###################################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)

if (.Platform$OS.type == "windows") {
  # setwd('A:/myThesis/make_features/BMR/')
  setwd("C:/Active/projects/make_features/BMR/")
}
# source('plots/functions.R')


prepare_df_stat <- function(paths_eval_DS, statistic){ # statistic can be 'Mean' or 'Variance'
  
  df <- do.call(rbind, lapply(paths_eval_DS, fread))
  
  print(1)
  
  sample_size <- c()
  model = c()
  for (path_eval_DS in paths_eval_DS) {
    M = extract_model_from_path(path_eval_DS)
    model = c(model, M)
    sample_size <- c(sample_size, extract_binSize_from_path(path_eval_DS))
    
  }
  print(2)
  
  df_mean <- df[which(df$V1 == statistic),]
  df_mean$model <- model
  df_mean$sampleSize <- sample_size
  df_mean <- data.frame(df_mean)
  df_ass_mean <- df_mean[,c('V1', ass_type, "model", "sampleSize")]
  df_ass_mean <- melt(df_ass_mean, id.vars = c(ass_type,  "model", "sampleSize"))
  
  df_ass_mean
}

prepare_data_DS_val <- function(paths_eval_DS, ass_type){
  
  
  df_mean <- data.frame(prepare_df_stat(paths_eval_DS, 'Mean'))
  df_var <- data.frame(prepare_df_stat(paths_eval_DS, 'Variance'))
  
  print(3)
  
  df_mean$M_minus_sd <- df_mean[,ass_type] - sd(df_var[,ass_type])
  df_mean$M_plus_sd <- df_mean[,ass_type] + sd(df_var[,ass_type])
  df_mean$assessment <- ass_type
  df_mean$sampleSize <- factor(df_mean$sampleSize, levels = c("FullSet", "DS1M", 
                                                              "DS800k", "DS600k", 
                                                              "DS300k", "DS100k", 
                                                              "DS50k"))
  
  print(4)
  df_mean
}


DownSampling_eval <- function(paths_eval_DS, ass_type){
  
  df_mean <- prepare_data_DS_val(paths_eval_DS, ass_type)
  
  ggplot() +
    geom_line(data = df_mean, aes(x = sampleSize, y = get(ass_type), color = model,
                                  group = model)) +
    # geom_point(data = df_mean, aes(x = sampleSize, y = corr, color = model), size = 3) +
    geom_errorbar(data = df_mean, aes(x = sampleSize, ymin = M_minus_sd, ymax = M_plus_sd, 
                                      color = model), width = 0.2) +
    facet_wrap(~assessment, scales = "free_y") +
    scale_color_manual(values = create_method_colours()) +
    labs(x = "Model", y = ass_type, title = "Line Plot with Error Bars") +
    theme_minimal() +
    theme(axis.title.y = element_text(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))  
  
}
#######################################################


paths_eval_DS = c('../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/nn_mseLoss/model_metrics_summary.tsv',
                  #'../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/nn_mseLoss/model_metrics_summary.tsv',
                  #'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/RF/model_metrics_summary.tsv', 
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/nn_mseLoss/model_metrics_summary.tsv',
                  #'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/nn_mseLoss/model_metrics_summary.tsv',
                  #'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/nn_mseLoss/model_metrics_summary.tsv',
                  #'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/nn_mseLoss/model_metrics_summary.tsv',
                  #'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/nn_mseLoss/model_metrics_summary.tsv'
                  #'../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/RF/model_metrics_summary.tsv',
)



ass_type = 'corr'
DownSampling_eval(paths_eval_DS, ass_type)
