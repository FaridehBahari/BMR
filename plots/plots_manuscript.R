# (fill="#69b3a2", alpha=0.5)
# https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=5

rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(tidyverse)

setwd('C:/Active/projects/make_features/BMR/')

create_method_colours <- function() {
  c(GBM = '#daa520', GBM0 = '#f2dca5', 
    GBM_longer100 = '#838077', GBM0_longer100 = '#cecdc9',
    RF = "#117744", GLM = "#aaaeba", TL = "#77AADD", `element specific` = '#446c8b',
    NN_MSEloss = "#771155", nn_mseLoss = "#771155", 
    NN_PoisLoss = "#11266d", nn_poisLoss = "#11266d",
    # `1k` = '#edf8fb', `10k` = '#b3cde3', `50k` = '#8c96c6',
    # `100k` = '#8856a7', `1M` = '#810f7c'
    `1k` = '#F0FFFF', `10k` = '#ADD8E6', `50k` =  '#87CEEB',
    `100k` = '#5F9EA0', `1M` =  '#7393B3'
  )
}


create_element_colours <- function() {
  c(enhancers = '#b03914', Enhancers = '#b03914',
    `3\' UTR` = "#14b556", gc19_pc.3utr = "#14b556", 
    `5\' UTR` = "#AAAA44", gc19_pc.5utr = "#AAAA44",
    CDS = "#CC99BB", gc19_pc.cds = "#CC99BB",
    `Core Promoter` = "#77AADD", gc19_pc.promCore = "#77AADD", 
    `Splice site` = "#11266d", gc19_pc.ss = "#11266d", 
    lncRNA = '#daa520', lncrna.ncrna = '#daa520', 
    `lncRNA Promoter` = "#DDAA77", lncrna.promCore = "#DDAA77"
  )
}


extract_binSize_from_path <- function(path){
  spl <- unlist(strsplit(path, "/"))
  model <- spl[length(spl)-2]
  model
}

define_element_type <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[1]}))
  GenomicElement
  
}


extract_ass_methods <- function(path_ass_method, ass_type, element){
  ass <- lapply(path_ass_method, fread)
  
  ass_ele <- lapply(ass, function(x){
    x = data.frame(x)
    row_idx = grep( paste0(ass_type, '_'), x$V1)
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


grouped_barPlot_model_elements <- function(path_ass_full, ass_type, compare,
                                           path_save = '../external/BMR/plots/',
                                           save_name = '', model_names = NULL){
  dir.create(path_save, showWarnings = F, recursive = T)
  element_names <- c('enhancers' = 'Enhancers',
                     'gc19_pc.3utr' = '3\' UTR',
                     'gc19_pc.5utr' = '5\' UTR',
                     'gc19_pc.cds' = 'CDS',
                     'gc19_pc.promCore' = 'Core Promoter',
                     'gc19_pc.ss' = 'Splice site'
                     # , 'lncrna.ncrna' = 'lncRNA',
                     # 'lncrna.promCore' = 'lncRNA Promoter'
  )
  
  elements <- names(element_names)
  
  
  
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = 'MSE'
  } 
  
  
  elem_dfs <- list()
  
  for (element in elements) {
    
    assessment <- extract_ass_methods(path_ass_full, ass_type, element)
    
    
    model <- c()
    for (path in path_ass_full) {
      
      
      if(compare == 'per_binSize'){
        m = extract_binSize_from_path(path)
        if (m == 'var_size') {
          m = 'GBM'
        }
      } else if (compare == 'per_models_varSize') {
        m = extract_model_from_path(path)
      } else {
        m = model_names
      }
      
      model = c(model, m)
      
    }
    
    df <- data.frame(cbind(model, assessment, element = element_names[element]))
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
  
  ggplot(elem_dfs, aes(x = element, y = assessment,
                       fill = model, group = interaction(model, id))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    # coord_polar(start = 0) +  # This turns the barplot into a circular barplot
    ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    labs(title = paste0("Model Assessment by Element", "\n", compare),
         y = Y_lable,
         fill = "Model")
  
  ggsave(paste0(ass_type, "grouped_barPlot_model_elements_", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


grouped_barPlot_model_elements_MSE <- function(path_ass_full, ass_type = 'mse', compare,
                                               path_save = '../external/BMR/plots/',
                                               save_name = '', model_names = NULL){
  dir.create(path_save, showWarnings = F, recursive = T)
  element_names <- c('enhancers' = 'Enhancers',
                     'gc19_pc.3utr' = '3\' UTR',
                     'gc19_pc.5utr' = '5\' UTR',
                     'gc19_pc.cds' = 'CDS',
                     'gc19_pc.promCore' = 'Core Promoter',
                     'gc19_pc.ss' = 'Splice site'
                     # , 'lncrna.ncrna' = 'lncRNA',
                     # 'lncrna.promCore' = 'lncRNA Promoter'
  )
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = 'MSE'
  } 
  
  elem_dfs <- list()
  elements <- names(element_names)
  for (element in elements) {
    
    assessment <- extract_ass_methods(path_ass_full, ass_type, element)
    
    
    model <- c()
    for (path in path_ass_full) {
      
      if(compare == 'per_binSize'){
        m = extract_binSize_from_path(path)
        if (m == 'var_size') {
          m = 'GBM'
        }
      } else if (compare == 'per_models_varSize') {
        m = extract_model_from_path(path)
      } else {
        m = model_names
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
  
  ggplot(elem_dfs, aes(x = element, y = assessment,
                       fill = model, group = interaction(model, id))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    labs(title = paste0("Model Assessment by Element", "\n", compare),
         y = Y_lable,
         fill = "Model")
  
  ggsave(paste0(ass_type, "grouped_barPlot_model_elements_", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


plot_validation_boxplot <- function(directory_paths, metric, compare,
                                    path_save = '../external/BMR/plots/',
                                    save_name = '') {
  
  dir.create(path_save, showWarnings = F, recursive = T)
  
  if (metric == 'corr') {
    Y_lable = 'Correlation'
  } else if (metric == 'acc') {
    Y_lable = 'Accuracy'
  } else if (metric == 'mse') {
    Y_lable = 'MSE'
  } 
  
  
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
          if (model == 'var_size') {
            model = 'GBM'
          }
        } else if (compare == 'per_models_varSize') {
          model = extract_model_from_path(directory_path)
        }
        
        metric_name <- grep(paste0(metric, '_'), colnames(metric_data), value = TRUE)
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
      labs(y = Y_lable, title = " performance of different models on intergenic validation sets using variable-size bins ") +
      
      theme(legend.position="right")+
      scale_fill_manual(values = create_method_colours()) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 14),
            axis.line = element_line(colour = "black"))
    
    ggsave(paste0(metric, "_validation_boxplot", save_name,".png"),
           device = "png", width = 10, 
           bg = 'white',
           path = path_save)
  } else {
    print(paste("No data found for metric '", metric, "' in the specified directories."))
  }
}


DownSampling_linePlot_perModel_allElems <- function(path_ass_DS, ass_type,
                                                    path_save = '../external/BMR/plots/',
                                                    save_name = ''){
  dir.create(path_save, showWarnings = FALSE, recursive = TRUE)
  
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = 'MSE'
  }
  
  model <- extract_model_from_path(path_ass_DS[1])
  
  # Mapping dictionary for element names
  element_names <- c('enhancers' = 'Enhancers',
                     'gc19_pc.3utr' = '3\' UTR',
                     'gc19_pc.5utr' = '5\' UTR',
                     'gc19_pc.cds' = 'CDS',
                     'gc19_pc.promCore' = 'Core Promoter',
                     'gc19_pc.ss' = 'Splice site'
                     #,
                     # 'lncrna.ncrna' = 'lncRNA',
                     # 'lncrna.promCore' = 'lncRNA Promoter'
  )
  
  elements <- names(element_names)
  
  elem_dfs <- list()
  
  for (element in elements) {
    
    original_element_name <- element
    
    assessment <- extract_ass_methods(path_ass_DS, ass_type, original_element_name)
    
    
    sampleSize <- c()
    for (path in path_ass_DS) {
      m = extract_binSize_from_path(path)
      sampleSize = c(sampleSize, m)
      
    }
    
    df <- data.frame(cbind(sampleSize, assessment, element = element_names[original_element_name]))
    elem_dfs <- append(elem_dfs, list(df))
  }
  
  elem_dfs <- bind_rows(elem_dfs)
  
  
  # Convert the sampleSize column to a factor to ensure correct ordering in the plot
  elem_dfs$sampleSize <- factor(elem_dfs$sampleSize, levels = c("FullSet", "DS1M", "DS800k", "DS600k", "DS300k", "DS100k", "DS50k"))
  
  # Create the line plot
  plot <- ggplot(data = elem_dfs, aes(x = sampleSize, y = assessment, 
                                      group = element, color = element)) +
    geom_line() +
    geom_point() +
    labs(title = "Reduction of Sample Size Monitoring",
         x = "Sample Size",
         y = Y_lable) +
    theme_minimal()
  
  
  
  # Create a new dataframe containing only the "FullSet" and "DS50k" assessments for element
  # full_ds50 <- subset(elem_dfs, sampleSize %in% c("FullSet", "DS600k", "DS50k"))
  full_ds50 <- elem_dfs
  full_ds50$assessment <- as.character(round(as.numeric(full_ds50$assessment), 3))
  
  
  # Plot the line plot with "FullSet" and "DS50k" assessments for each element
  plot <- plot +
    geom_text(data = full_ds50, aes(label = assessment), hjust = -0.2, vjust = 0) +
    labs(title = paste0("Reduction of Sample Size Monitoring in Different Elements using ", model, " model"),
         x = "Sample Size",
         y = Y_lable) +
    scale_color_manual(values = create_element_colours()) +
    theme_minimal() +
    theme(axis.title.y = element_text(),
          axis.text.y = element_blank(),  # Remove y-axis values
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 14),
          axis.line = element_line(colour = "black"))  # Remove y-axis ticks
  
  ggsave(paste0(ass_type, "_DownSampling_linePlot_perModel_allElems", 
                save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
  
  print(plot)
  
}


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


DownSampling_eval <- function(paths_eval_DS, ass_type,
                              path_save = '../external/BMR/plots/',
                              save_name = ''){
  
  dir.create(path_save, showWarnings = F, recursive = T)
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = 'MSE'
  }
  df_mean <- prepare_data_DS_val(paths_eval_DS, ass_type)
  
  ggplot() +
    geom_line(data = df_mean, aes(x = sampleSize, y = get(ass_type), color = model,
                                  group = model)) +
    # geom_point(data = df_mean, aes(x = sampleSize, y = corr, color = model), size = 3) +
    geom_errorbar(data = df_mean, aes(x = sampleSize, ymin = M_minus_sd, ymax = M_plus_sd, 
                                      color = model), width = 0.2) +
    facet_wrap(~assessment, scales = "free_y") +
    scale_color_manual(values = create_method_colours()) +
    labs(x = "Model", y = Y_lable, title = "Line Plot with Error Bars") +
    theme_minimal() +
    theme(axis.title.y = element_text(),
          axis.text.y = element_text(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 14),
          axis.line = element_line(colour = "black"))  
  
  ggsave(paste0(ass_type, "_DownSampling_eval", 
                save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
  
}


extract_dimReducMethod_from_path <- function(path){
  spl <- unlist(strsplit(path, "/"))
  model <- spl[length(spl)-2]
  if (!model %in% c('PCA', 'AE') ) {
    model <- 'original'
  }
  model
}


grouped_barPlot_dimReduction <- function(path_ass_full, ass_type, 
                                         path_save = '../external/BMR/plots/',
                                         save_name = ''){
  
  dir.create(path_save, showWarnings = F, recursive = T)
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = 'MSE'
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
  
  elements <- names(element_names)
  
  
  elem_dfs <- list()
  
  for (element in elements) {
    
    assessment <- extract_ass_methods(path_ass_full, ass_type, element)
    
    
    model <- c()
    dim_red <- c()
    for (path in path_ass_full) {
      
      m = extract_model_from_path(path)
      model = c(model, m)
      
      dimReduc <- extract_dimReducMethod_from_path(path)
      dim_red <- c(dim_red, dimReduc)
      
    }
    
    df <- data.frame(cbind(model, dim_red, assessment, element = element_names[element]))
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
  
  # Reorder elements based on the model variable
  elem_dfs <- elem_dfs %>%
    arrange(model, element)
  
  # Plot the grouped bar plot
  ggplot(elem_dfs, aes(x = element, y = assessment, fill = model, alpha = dim_red)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    scale_alpha_manual(values = c(original = 1, PCA = 0.6, AE = 0.2))+
    labs(#title = "Model Assessment by Element",
      y = Y_lable,
      fill = "Model")
  
  ggsave(paste0(ass_type, "grouped_barPlot_dimReduction_", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}

extract_info_from_path <- function(directory_path) {
  parts <- strsplit(directory_path, '/')[[1]]
  method <- parts[length(parts)-2]
  if (!method %in% c('PCA', 'AE') ) {
    method <- 'original'
  }
  model <- parts[length(parts)-1]
  return(c(method = method, model = model))
}




grouped_boxplot_deimReduction <- function(directory_paths, metric, 
                                          path_save = '../external/BMR/plots/',
                                          save_name = '') {
  
  dir.create(path_save, showWarnings = F, recursive = T)
  if (metric == 'corr') {
    Y_lable = 'Correlation'
  } else if (metric == 'acc') {
    Y_lable = 'Accuracy'
  } else if (metric == 'mse') {
    Y_lable = 'MSE'
  }
  
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
          mutate(Method = info['method'], Model = info['model'],
                 Value = as.numeric(metric_data[, metric_name]))
        
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
    ggplot(combined_data, aes(x = Model, y = Value, fill = Model, alpha = Method)) +
      geom_boxplot(position = position_dodge(width = 0.6), width = 0.5) +  # Adjust the width as needed
      labs(y = Y_lable, title = paste0(metric, " Distribution Comparison between Models using different dimension reduction Methods")) +
      
      theme(legend.position="top",
            text = element_text(size = 14))+
      scale_fill_manual(values = create_method_colours()) +
      scale_alpha_manual(values = c(original = 1, PCA = 0.4, AE = 0.1))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    print(paste("No data found for metric '", metric, "' in the specified directories."))
  }
  
  ggsave(paste0(ass_type, "grouped_boxPlot_dimReduction_", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


grouped_barPlot_model_benchmark <- function(driver_based_on, coding_nonCoding,
                                            path_save = '../external/BMR/plots/',
                                            save_name = ''){
  
  path_bench_res <- paste0('../external/BMR/benchmark_binomTest/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_', driver_based_on, 'fdr.csv')
  dir.create(path_save, showWarnings = F, recursive = T)
  df <- fread(path_bench_res)
  df <- df[which(df$method %in% c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL" )),]
  
  # Get column indices where "CDS" is present in column names
  cds_cols <- grep("CDS", names(df))
  
  # Split dataframe based on column indices
  cds_df <- df[, c(2, cds_cols), with = FALSE]
  nc_df <- df[, -cds_cols, with = FALSE]
  
  colnames(cds_df) <- gsub('CDS_', '', toupper(colnames(cds_df)))
  colnames(nc_df) <- gsub('NC_', '', toupper(colnames(nc_df)))
  
  if (coding_nonCoding == 'coding') {
    df = cds_df
  } else if(coding_nonCoding == 'non-coding')  {
    df = nc_df
  }
  
  df = df[, -c('N_ELEMNTS', 'NTPS', 'NHITS')]
  df <- melt(setDT(df), id.vars = "METHOD", c("PRECISIONS", "RECALLS", "F1",
                                              "AUC", "AUPR"))
  colnames(df) <- c('model', 'metrics', 'value')
  
  # Order the data first by element and then by decreasing assessment
  df <- df %>%
    arrange(metrics, desc(model ))
  
  # Create groups for laying out the barplot correctly
  dfs <- df %>%
    group_by(metrics) %>%
    mutate(id = metrics)
  
  ggplot(dfs, aes(x = metrics, y = value, fill = model, group = interaction(model, id))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    # coord_polar(start = 0) +  # This turns the barplot into a circular barplot
    ylim(min(0, dfs$value), max(1, dfs$value)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    labs(title = paste0(coding_nonCoding, ':\ngold standard based on ', based_on),
         y = 'value',
         x = 'metrics',
         fill = "Model")
  
  ggsave(paste0(coding_nonCoding,"_", driver_based_on, "_grouped_barPlot_model_benchmark", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


precisinRecall_dotPlot <- function(driver_based_on, coding_nonCoding,
                                   path_save = '../external/BMR/plots/',
                                   save_name = ''){
  
  path_bench_res <- paste0('../external/BMR/benchmark_binomTest/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_', driver_based_on, 'fdr.csv')
  dir.create(path_save, showWarnings = F, recursive = T)
  df <- fread(path_bench_res)
  df <- df[which(df$method %in% c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL" )),]
  
  # Get column indices where "CDS" is present in column names
  cds_cols <- grep("CDS", names(df))
  
  # Split dataframe based on column indices
  cds_df <- df[, c(2, cds_cols), with = FALSE]
  nc_df <- df[, -cds_cols, with = FALSE]
  
  
  colnames(cds_df) <- gsub('CDS_', '', toupper(colnames(cds_df)))
  colnames(nc_df) <- gsub('NC_', '', toupper(colnames(nc_df)))
  
  if (coding_nonCoding == 'coding') {
    df = cds_df
  } else if(coding_nonCoding == 'non-coding')  {
    df = nc_df
  }
  
  dfs = df
  
  dfs$model <- factor(df$METHOD)
  coding_nonCoding = 'coding'
  ggplot(dfs, aes(x = RECALLS, y = PRECISIONS)) +
    geom_point(aes(colour = model)) +
    # coord_polar(start = 0) +  # This turns the barplot into a circular barplot
    ylim(min(0, dfs$PRECISIONS), max(1, dfs$PRECISIONS)) + 
    xlim(min(0, dfs$RECALLS), max(1, dfs$RECALLS)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_colour_manual(values = create_method_colours()) +
    labs(title = paste0(coding_nonCoding, ':\ngold standard based on ', based_on),
         x = 'Recall', y = 'Precision'
    )
  
  ggsave(paste0(coding_nonCoding,"_", driver_based_on, "_PrecisionRecall", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


top_100hit_barplot <- function(based_on, coding_nonCoding,
                               path_save = '../external/BMR/plots/',
                               save_name = ''){
  
  path_bench_res <- paste0('../external/BMR/benchmark_binomTest/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_',
                           based_on, 'fixedNumberOfElems.csv')
  
  dir.create(path_save, showWarnings = F, recursive = T)
  df <- fread(path_bench_res)
  df <- df[which(df$method %in% c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL" )),]
  
  # Get column indices where "CDS" is present in column names
  cds_cols <- grep("CDS", names(df))
  
  # Split dataframe based on column indices
  cds_df <- df[, c(2, cds_cols), with = FALSE]
  nc_df <- df[, -cds_cols, with = FALSE]
  
  colnames(cds_df) <- gsub('CDS_', '', toupper(colnames(cds_df)))
  colnames(nc_df) <- gsub('NC_', '', toupper(colnames(nc_df)))
  
  if (coding_nonCoding == 'coding') {
    df = cds_df
  } else if(coding_nonCoding == 'non-coding')  {
    df = nc_df
  }
  
  
  # Order the data first by element and then by decreasing assessment
  df <- df %>%
    arrange(METHOD , desc( NTPS))
  
  
  ggplot(df, aes(x = METHOD, y = NTPS, fill = METHOD)) +
    geom_bar(stat = "identity", width = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    labs(title = paste0(coding_nonCoding, ':\ngold standard based on ', based_on),
         y = 'number of true positives',
         x = 'model')
  
  ggsave(paste0(coding_nonCoding,"_", based_on, "_nTPs_barPlot", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


histogram_perElement_length <- function(path_responseTab_PCAWG,
                                        path_save = '../external/BMR/plots/',
                                        save_name = ''){
  
  y <- fread(path_responseTab_PCAWG)
  
  y$`element type` <- define_element_type(y$binID)
  
  y <- y[which(!y$`element type` %in% c('lncrna.ncrna', 'lncrna.promCore')), ]
  
  # Mapping dictionary for element types
  element_types <- c('enhancers' = 'Enhancers',
                     'gc19_pc.3utr' = '3\' UTR',
                     'gc19_pc.5utr' = '5\' UTR',
                     'gc19_pc.cds' = 'CDS',
                     'gc19_pc.promCore' = 'Core Promoter',
                     'gc19_pc.ss' = 'Splice site')
  
  # Replace the element type names in the data frame with modified names for plotting
  y$`element type` <- element_types[y$`element type`]
  
  # Plot histogram
  ggplot(y, aes(x = length, fill = `element type`)) +
    geom_histogram(binwidth = 20) +
    facet_wrap(~`element type`, scales = "free") +
    theme_minimal() +
    theme(text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_element_colours()) +  # Use the function to generate colors
    labs(title = "Length Distribution for Each Element Type",
         x = "Length",
         y = "Frequency")
  
  ggsave(paste0(save_name,"perElement_length_histogram.png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


# Function to read .tsv files and extract 'corr' column
read_corr_from_tsv <- function(directory) {
  tsv_file <- file.path(directory, "model_metrics_summary.tsv")
  if (file.exists(tsv_file)) {
    data <- data.frame(fread(tsv_file))
    corr <- data$corr
    # print(length(corr)) 
    return(corr)
  } else {
    warning(paste("No file found in directory:", directory))
    return(NULL)
  }
}



importance_plot <- function(parent_directory, path_GBM, LOO = FALSE,
                            path_save = '../external/BMR/plots/',
                            save_name = ''){
  
  # List all subdirectories
  subdirectories <- list.dirs(parent_directory, recursive = FALSE)
  
  # Read corr values from each subdirectory
  
  df <- data.frame()
  for (directory in subdirectories) {
    metric = read_corr_from_tsv(directory)
    df <- rbind(df, (metric))
  }
  colnames(df) <- c('Mean', 'variance')
  
  corr_data <- df$Mean
  
  # Create a data frame with subdirectory names and corr values
  plot_data <- data.frame(`Feature Group` = sub("GBM_", "", basename(subdirectories)),
                          Correlation = corr_data)
  plot_data$Feature.Group = sapply(plot_data$Feature.Group, function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
  plot_data <- plot_data[which(!plot_data$Feature.Group %in% c('DNA_methylation', 'APOBEC')),]
  plot_data$Feature.Group <- gsub('_', " ", plot_data$Feature.Group)
  
  full_corr_gbm <- fread(path_GBM)$corr[1]
  
  if (LOO) {
    
    plot_data$Correlation <- full_corr_gbm - plot_data$Correlation
    y_axis <- 'Difference of correlation'
  } else {
    y_axis <- 'Correlation'
  }
  
  plot_data$Feature.Group <- factor(plot_data$Feature.Group, 
                                    levels = plot_data$Feature.Group[order(-plot_data$Correlation)])
  
  # Plot using ggplot2
  p <- ggplot(plot_data, aes(x = Feature.Group, y = Correlation, fill = Feature.Group)) +
    geom_bar(stat = "identity") +
    labs(y = "Correlation", fill = NULL) +
    scale_fill_manual(values = c(`Epigenetic mark` = '#1b9e77',
                                  `RNA expression` = '#d95f02',
                                  `Replication timing`= '#7570b3',
                                  HiC = '#e7298a',
                                  `DNA accessibility` = '#66a61e',
                                  `NucleotideContext` = '#e6ab02',
                                 Conservation = '#a6761d')) +
    
    #c('#081d58','#253494','#225ea8','#1d91c0','#41b6c4','#c7e9b4','#edf8b1') 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.title = element_blank()) +
    xlab("")
  
  if (!LOO) {
    p <- p + geom_hline(yintercept = full_corr_gbm, linetype = "dashed", linewidth = 1.3)
  }
  
  ggsave(paste0("FeatureImportance", save_name,".png"),
         plot = p, device = "png", width = 10, 
         bg = 'white',
         path = path_save)
  
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#####     fig1: Compare different models on variable-size intergenic bins####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
############## elements:

path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/nn_poisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/nn_mseLoss_assessments.tsv")


ass_types <- c('corr', 'acc', 'mse')
for (ass_type in ass_types) {
  if (ass_type != 'mse') {
    grouped_barPlot_model_elements(path_ass_full, ass_type, 'per_models_varSize', save_name = 'new')
  } else {
    grouped_barPlot_model_elements_MSE(path_ass_full, compare = 'per_models_varSize', save_name = 'new')
  }
}

############## validation Sets:

directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/rep_train_test/')

ass_types <- c('corr', 'acc', 'mse')
for (ass_type in ass_types) {
  plot_validation_boxplot(directory_paths, ass_type, 'per_models_varSize', save_name = 'new') 
  
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
######     fig1_supp: mutated-unmutated bins ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100/rep_train_test/')

ass_types <- c('corr', 'acc', 'mse')


path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0/GBM0_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100/GBM_longer100_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100/GBM0_longer100_assessments.tsv")


for (ass_type in ass_types) {
  plot_validation_boxplot(directory_paths, ass_type, 'per_models_varSize', save_name = 'unMutatedANDlength')
  grouped_barPlot_model_elements(path_ass_full, ass_type, 'per_models_varSize', save_name = 'unMutatedANDlength')
  
}
grouped_barPlot_model_elements_MSE(path_ass_full, compare = 'per_models_varSize', save_name = 'unMutatedANDlength')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig2.Down sampling figures ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_DS <- c('../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/GBM/GBM_assessments.tsv')

for (ass_type in ass_types) {
  DownSampling_linePlot_perModel_allElems(path_ass_DS, ass_type)
  
}



paths_eval_DS = c('../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/nn_mseLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/nn_mseLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/nn_mseLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/nn_mseLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/nn_mseLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/nn_mseLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/RF/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/GBM/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/nn_poisLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/nn_mseLoss/model_metrics_summary.tsv',
                  '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/RF/model_metrics_summary.tsv'
)



for (ass_type in ass_types) {
  DownSampling_eval(paths_eval_DS, ass_type)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig3.Transfer learning figures ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_full <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv',
                   '../external/BMR/output/element-specific/GBM/GBM_ensemble_bootstraps100_assessment.tsv',
                   '../external/BMR/output/TL/GBM/GBM_ensemble_bootstraps100_assessment.tsv')

ass_types = c('corr', 'acc', 'mse')
compare <- ''
model_names <- c('GBM', 'element specific', 'TL')

for (ass_type in ass_types) {
  if (ass_type != 'mse') {
    grouped_barPlot_model_elements(path_ass_full, ass_type, compare,
                                   save_name = 'TLvsGBM', model_names = model_names)
  } else {
    grouped_barPlot_model_elements_MSE(path_ass_full, ass_type, compare,
                                       save_name = 'TLvsGBM', model_names = model_names)
  }
}



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig.4   dimension reduction figures ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

path_ass_orig <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/nn_poisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/nn_mseLoss_assessments.tsv")

path_ass_dimReds_pca <- c("../external/BMR/output/dimReduction_effect/PCA/GBM/GBM_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/RF/RF_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/nn_poisLoss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/nn_mseLoss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/GLM/GLM_assessments.tsv")

path_ass_dimReds_ae <- c("../external/BMR/output/dimReduction_effect/AE/GBM/GBM_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/RF/RF_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/nn_poisLoss/nn_poisLoss_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/nn_mseLoss/nn_mseLoss_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/GLM/GLM_assessments.tsv")


path_ass_full <- c(path_ass_orig, path_ass_dimReds_pca, path_ass_dimReds_ae)

ass_types <- c('corr', 'acc', 'mse')
for (ass_type in ass_types) {
  grouped_barPlot_dimReduction(path_ass_full, ass_type)
}




########################
directory_paths_PCA <- c(
  "../external/BMR/output/dimReduction_effect/PCA/GBM/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/RF/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/GLM/rep_train_test/"
)


directory_paths_AE <- c(
  "../external/BMR/output/dimReduction_effect/AE/GBM/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/RF/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/nn_poisLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/nn_mseLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/AE/GLM/rep_train_test/"
)

directory_paths_orig <- c(
  "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/",
  "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/rep_train_test/",
  "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/rep_train_test/",
  "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/rep_train_test/"
)

directory_paths <- c(directory_paths_PCA, directory_paths_AE, directory_paths_orig)

for (ass_type in ass_types) {
  grouped_boxplot_deimReduction(directory_paths, ass_type)
}



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig.bin size effect figures (fixed-concatenated intervals) ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_binSizes <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/GBM_assessments.tsv'
                       # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1k/GBM/GBM_assessments.tsv'
)

ass_types <- c('corr', 'acc', 'mse')
save_name = 'concatenated'
for (ass_type in ass_types) {
  
  if (ass_type != 'mse') {
    grouped_barPlot_model_elements(path_ass_binSizes, ass_type, compare = 'per_binSize', save_name = save_name)
  } else {
    grouped_barPlot_model_elements_MSE(path_ass_binSizes, ass_type, 
                                       compare = 'per_binSize', save_name = save_name)
  }
}


binEffect_directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/rep_train_test/'
                               # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1k/GBM/rep_train_test/'
)


for (ass_type in ass_types) {
  plot_validation_boxplot(binEffect_directory_paths, ass_type, 'per_binSize', save_name = save_name) 
}

#####     fig. bin size effect figures (fixed-window intervals) #####     

path_ass_binSizes <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/1M/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/100k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/50k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/10k/GBM/GBM_assessments.tsv'
)

save_name = 'splitted'
for (ass_type in ass_types) {
  
  if (ass_type != 'mse') {
    grouped_barPlot_model_elements(path_ass_binSizes, ass_type, 'per_binSize', save_name = save_name)
  } else {
    grouped_barPlot_model_elements_MSE(path_ass_binSizes, 'per_binSize', save_name = save_name)
  }
  
}


binEffect_directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/1M/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/100k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/50k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/10k/GBM/rep_train_test/'
)


for (ass_type in ass_types) {
  plot_validation_boxplot(binEffect_directory_paths, ass_type, 
                          'per_binSize', save_name = save_name) 
}




######################### driver identification plots ##########################
save_name = 'binomTest'
coding_nonCodings <- c('coding', 'non-coding')
based_ons <- c('all', 'in_CGC_new', 'in_pcawg', 'in_oncoKB', 'any')
for (coding_nonCoding in coding_nonCodings) {
  for (based_on in based_ons) {
    grouped_barPlot_model_benchmark(based_on, coding_nonCoding, save_name = save_name)
    precisinRecall_dotPlot(based_on, coding_nonCoding, save_name = save_name)
    top_100hit_barplot(based_on, coding_nonCoding, save_name = save_name)
  }
}

####################### PCAWG length distribution plots ########################
path_responseTab_PCAWG <- '../external/BMR/rawInput/responseTabs/Pan_Cancer/reg_elems.tsv'
histogram_perElement_length(path_responseTab_PCAWG)

######################## feature importance plots ############################

parent_directory <- "../external/BMR/output/featureImportance/"
path_GBM <- '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/model_metrics_summary.tsv'
importance_plot(parent_directory, path_GBM)



parent_directory <- "../external/BMR/output/featureImportance_LOO/"
path_GBM <- '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/model_metrics_summary.tsv'
importance_plot(parent_directory, path_GBM = path_GBM, save_name = 'LOO', LOO = T)
############################################################################

# elementSp vs TL























#################################################################
# #length distribution
# 
# library(ggplot2)
# 
# 
# # Read the data
# 
# y <- fread('../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv')
# elem_length1 <- y$length
# 
# y <- fread('../external/BMR/rawInput/responseTabs_bedtools/Pan_Cancer/var_bins.tsv')
# elem_length2 <- y$length
# 
# y <- fread('../external/BMR/rawInput/responseTabs/Pan_Cancer/reg_elems.tsv')
# elem_length3 <- y$length
# 
# y <- fread('../../BMR_proj/external/rawInput/Pan_Cancer_train_y.tsv')
# elem_length4 <- y$length
# 
# # Combine all element lengths into a single vector
# all_elem_lengths <- c(elem_length1, elem_length2, elem_length3, elem_length4)
# 
# # Determine the maximum value among all element lengths
# max_length <- max(all_elem_lengths)
# 
# 
# # Create violin plots for each vector separately
# plot1 <- ggplot(data = NULL, aes(x = "", y = elem_length1)) +
#   geom_violin(fill = "skyblue", color = "black") +
#   labs(title = "", x = "non-PCAWG callable variable bins", y = "") +
#   scale_y_log10(limits = c(1, max_length))+
#   theme_minimal()
# 
# plot2 <- ggplot(data = NULL, aes(x = "", y = elem_length2)) +
#   geom_violin(fill = "skyblue", color = "black") +
#   labs(title = "", x = "Full set variable bins", y = "") +
#   scale_y_log10(limits = c(1, max_length))+
#   theme_minimal()
# 
# plot3 <- ggplot(data = NULL, aes(x = "", y = elem_length3)) +
#   geom_violin(fill = "skyblue", color = "black") +
#   labs(title = "", x = "PCAWG elements", y = "Length") +
#   scale_y_log10(limits = c(1, max_length))+
#   theme_minimal()
# 
# plot4 <- ggplot(data = NULL, aes(x = "", y = elem_length4)) +
#   geom_violin(fill = "skyblue", color = "black") +
#   labs(title = "", x = "DP intergenic bins", y = "") +
#   scale_y_log10(limits = c(1, max_length))+
#   theme_minimal()
# 
# # Arrange the plots in a grid
# grid.arrange(plot3, plot1, plot2, plot4, ncol = 4)
# 





# 
# # Read the data
# 
library(ggplot2)

# Read the data and extract element lengths with non-zero nMut values
y1 <- fread('../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv')
elem_length1 <- y1$length[which(y1$nMut != 0)]

y2 <- fread('../external/BMR/rawInput/responseTabs_bedtools/Pan_Cancer/var_bins.tsv')
elem_length2 <- y2$length[which(y2$nMut != 0)]

y3 <- fread('../external/BMR/rawInput/responseTabs/Pan_Cancer/reg_elems.tsv')
elem_length3 <- y3$length[which(y3$nMut != 0)]

y4 <- fread('../../BMR_proj/external/rawInput/Pan_Cancer_train_y.tsv')
elem_length4 <- y4$length[which(y4$nMut != 0)]

# Combine all element lengths into a single vector
all_elem_lengths <- c(elem_length1, elem_length2, elem_length3, elem_length4)

# Determine the maximum value among all element lengths
max_length <- max(all_elem_lengths)

# Create violin plots for each vector separately, setting the same y-axis limits
plot1 <- ggplot(data = NULL, aes(x = "", y = elem_length1)) +
  geom_violin(fill = '#7393B3', color = "black") +
  labs(title = "", x = "non-PCAWG callable variable bins", y = "") +
  scale_y_log10(limits = c(1, max_length)) +  # Set y-axis limits
  theme_minimal()

plot2 <- ggplot(data = NULL, aes(x = "", y = elem_length2)) +
  geom_violin(fill = "#FFDAB9", color = "black") +
  labs(title = "", x = "Full set variable bins", y = "") +
  scale_y_log10(limits = c(1, max_length)) +  # Set y-axis limits
  theme_minimal()

plot3 <- ggplot(data = NULL, aes(x = "", y = elem_length3)) +
  geom_violin(fill = "#e77349", color = "black") +
  labs(title = "", x = "PCAWG elements", y = "Length") +
  scale_y_log10(limits = c(1, max_length)) +  # Set y-axis limits
  theme_minimal()

plot4 <- ggplot(data = NULL, aes(x = "", y = elem_length4)) +
  geom_violin(fill = "#FFDAB9", color = "black") +
  labs(title = "", x = "DP intergenic bins", y = "") +
  scale_y_log10(limits = c(1, max_length)) +  # Set y-axis limits
  theme_minimal()

# Arrange the plots in a grid
grid.arrange(plot3, plot1, plot2, plot4, ncol = 4)
# 


grid.arrange(plot3, plot1, ncol = 2)

# path_tetsY <- '../external/BMR/rawInput/responseTabs/Pan_Cancer/reg_elems.tsv'
# length_distribution_perElement <- function(path_tetsY){
#   
#   y <- fread(path_tetsY)
#   y$element_type <- define_element_type(y$binID)
#   elements <- unique(y$element_type)
#   max_length <- max(y$length)
#   
#   elem = 'enhancers'
#   len <- y[which(y$element_type == elem),]$length
#   
#   p1 <- ggplot(data = NULL, aes(x = "", y = len)) +
#     geom_violin(fill ="#b03914" , color = "black") +
#     labs(title = "", x = elem, y = "") +
#     scale_y_log10(limits = c(1, max_length)) +  
#     theme_minimal()
#   
#   elem = 'gc19_pc.cds'
#   len <- y[which(y$element_type == elem),]$length
#   
#   p2 <- ggplot(data = NULL, aes(x = "", y = len)) +
#     geom_violin(fill = "#CC99BB" , color = "black") +
#     labs(title = "", x = 'CDS', y = "") +
#     scale_y_log10(limits = c(1, max_length)) +  
#     theme_minimal()
#   
#   elem = 'gc19_pc.ss'
#   len <- y[which(y$element_type == elem),]$length
#   
#   p3 <- ggplot(data = NULL, aes(x = "", y = len)) +
#     geom_violin(fill = "#11266d"   , color = "black") +
#     labs(title = "", x = 'splice site', y = "") +
#     scale_y_log10(limits = c(1, max_length)) +  
#     theme_minimal()
#   
#   elem = 'gc19_pc.5utr'
#   len <- y[which(y$element_type == elem),]$length
#   
#   p4 <- ggplot(data = NULL, aes(x = "", y = len)) +
#     geom_violin(fill = "#AAAA44", color = "black") +
#     labs(title = "", x = '5UTR', y = "") +
#     scale_y_log10(limits = c(1, max_length)) +  
#     theme_minimal()
#   
#   elem = 'gc19_pc.3utr'
#   len <- y[which(y$element_type == elem),]$length
#   
#   p5 <- ggplot(data = NULL, aes(x = "", y = len)) +
#     geom_violin(fill = "#14b556" , color = "black") +
#     labs(title = "", x = '3UTR', y = "") +
#     scale_y_log10(limits = c(1, max_length)) +  
#     theme_minimal()
#   
#   elem = 'gc19_pc.promCore'
#   len <- y[which(y$element_type == elem),]$length
#   
#   p6 <- ggplot(data = NULL, aes(x = "", y = len)) +
#     geom_violin(fill = "#FFDAB9", color = "black") +
#     labs(title = "", x = 'core promoter', y = "") +
#     scale_y_log10(limits = c(1, max_length)) +  
#     theme_minimal()
#   
#   elem = 'lncrna.ncrna'
#   len <- y[which(y$element_type == elem),]$length
#   
#   # p7 <- ggplot(data = NULL, aes(x = "", y = len)) +
#   #   geom_violin(fill = "#FFDAB9", color = "black") +
#   #   labs(title = "", x = 'lncRNA', y = "") +
#   #   scale_y_log10(limits = c(1, max_length)) +  
#   #   theme_minimal()
#   # 
#   # elem = 'lncrna.promCore'
#   # len <- y[which(y$element_type == elem),]$length
#   # 
#   # p8 <- ggplot(data = NULL, aes(x = "", y = len)) +
#   #   geom_violin(fill = "#FFDAB9", color = "black") +
#   #   labs(title = "", x = 'lncRNA promoter', y = "") +
#   #   scale_y_log10(limits = c(1, max_length)) +
#   #   theme_minimal()
#   
#   
#   grid.arrange(p1, p2, p3, p4, p5, p6, 
#                # p7, p8, 
#                ncol = 2)
# }
# 


# Load necessary libraries
library(ggplot2)
library(reshape2)

# Your data
data <- data.frame(
  Element = c("lncrna.ncrna", "lncrna.promCore", "gc19_pc.ss", "enhancers", "gc19_pc.cds", "gc19_pc.promCore", "gc19_pc.5utr", "gc19_pc.3utr"),
  Ratio = c(0.8030734636896665, 0.6287743442778786, 0.2988726618976751, 0.5485366783685965, 0.6551608689501593, 0.6339806785697778, 0.3374977825852253, 0.6047107207493484,
            0.4685042593494515, 0.2802511309420116, 0.7783057631595223, 0.5532080267569323, 0.3803545917843505, 0.4845719063059615, 0.7607092444284935, 0.5178973122878467,
            0.8395220670591235, 0.648191657968322, 0.3570648648928321, 0.5063690916847651, 0.7241949462797541, 0.623272801289832, 0.3418936567919746, 0.6262748972663467,
            0.7036465906981347, 0.49817102987443995, 0.8542001855971655, 0.6285330456499335, 0.682057440661196, 0.7497892902441917, 0.8417247397451391, 0.6535073110879858,
            0.9472286800869942, 0.8809027605832732, 0.9191374533271238, 0.9137284342594494, 0.8916824121964794, 0.8881900581257646, 0.9261607895366636, 0.9059664874794969,
            0.7727922566719677, 0.5898240983728147, 0.5142062218301741, 0.5860233680938737, 0.7173929786744792, 0.566321736071273, 0.5678976000557882, 0.6678578377194958,
            0.38139557954689574, 0.16584506287594716, 0.33894087398252604, 0.1731014177445291, 0.1544577900050936, 0.15563884122730207, 0.3235945016598348, 0.17792215622915933),
  Group = rep(c("HiC", "nucleotideContext", "Replication_timing", "DNA_accessibility", "Epigenetic_mark", "RNA_expression", "conservation"), each = 8)
)

data$Element <- factor(heatmap_data$Element, levels = c("lncrna.ncrna", "lncrna.promCore", 
                                                        "gc19_pc.ss", "enhancers", "gc19_pc.cds",
                                                        "gc19_pc.promCore", "gc19_pc.5utr", 
                                                        "gc19_pc.3utr"))





library(RColorBrewer)

# Plot the heatmap
ggplot(data, aes(x = Group, y = Element, fill = Ratio)) +
  geom_tile() +
  geom_text(aes(label = round(Ratio, 2)), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient(low =  "darkblue", high = "#D53E4F", 
                      breaks = c(.1, .2,  .4, .6,  .8, .9),
                      labels = brewer.pal(6, "Spectral"),
                      limits = c(0, 1),
                      # na.value = "grey50",
                      guide = "colorbar")  +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Ratios", x = "Group", y = "Element")




colnames(data) <- c("Element type", "Importance per feature category \n (mean correlation ratio)",   "Feature Category" )

# Plot the heatmap
ggplot(data, aes(x = `Feature Category`, y = `Element type`, fill = `Importance per feature category \n (mean correlation ratio)`)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(`Importance per feature category \n (mean correlation ratio)`, 2)), 
            color = "black", size = 3.5)  +  # Add text labels
  scale_fill_gradientn(colors = c( 'black', "#43589F", '#7393B3', "#3288BD", "#66C2A5", "#ABDDA4", 
                                   "#E6F598", "#FDAE61", "#D53E4F",
                                   "#9E0142"),# rev(brewer.pal(6, "Spectral")),
                       breaks = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                       labels = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                       limits = c(0, 1),
                       na.value = "grey50",
                       guide = "colorbar") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),) +
  labs(x = "Feature Category", y = "Element type")

# "#5E4FA2" "#3288BD" "#66C2A5" "#ABDDA4" "#E6F598" "#FFFFBF" "#FEE08B" "#FDAE61" "#F46D43" "#D53E4F" "#9E0142"
