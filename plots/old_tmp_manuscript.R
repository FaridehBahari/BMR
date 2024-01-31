rm(list = ls())
if (.Platform$OS.type == "windows") {
  setwd('C:/Active/projects/make_features/BMR/')
}
source('plots/functions.R')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#         dimension reduction figures
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#         different models on callable variable size elements
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@