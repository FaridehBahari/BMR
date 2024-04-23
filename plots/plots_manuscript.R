# (fill="#69b3a2", alpha=0.5)
# https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=5

rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)


create_method_colours <- function() {
  c(GBM = '#daa520', GBM0 = '#f2dca5', 
    GBM_longer100 = '#838077', GBM0_longer100 = '#cecdc9',
    RF = "#117744", GLM = "#aaaeba",
    NN_MSEloss = "#771155", nn_mseLoss = "#771155", 
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


grouped_barPlot_model_elements <- function(path_ass_full, ass_type, compare,
                                           path_save = '../external/BMR/plots/',
                                           save_name = ''){
  dir.create(path_save, showWarnings = F, recursive = T)
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
    ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    labs(title = "Model Assessment by Element",
         y = ass_type,
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
  dir.create(path_save, showWarnings = F, recursive = T)
  
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
          text = element_text(size = 14),
          axis.line = element_line(colour = "black"))  # Remove y-axis ticks
  
  ggsave(paste0(ass_type, "_DownSampling_linePlot_perModel_allElems", 
                save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
  
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
  
  elements <- c('enhancers', 'gc19_pc.3utr', 'gc19_pc.5utr',
                'gc19_pc.cds', 'gc19_pc.promCore',
                'gc19_pc.ss', 'lncrna.ncrna', 'lncrna.promCore')
  
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
    
    df <- data.frame(cbind(model, dim_red, assessment, element))
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
      y = ass_type,
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
      labs(y = metric, title = paste0(metric, " Distribution Comparison between Models using different dimension reduction Methods")) +
      
      theme(legend.position="top")+
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


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#####     fig1: Compare different models on variable-size intergenic bins####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
############## elements:

path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/nn_poisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/nn_mseLoss_assessments.tsv")


ass_type <- 'corr'
grouped_barPlot_model_elements(path_ass_full, ass_type, 'per_models_varSize', save_name = 'new')



############## validation Sets:

directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/rep_train_test/')


plot_validation_boxplot(directory_paths, 'corr', 'per_models_varSize', save_name = 'new') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
######     fig2: mutated-unmutated bins ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100/rep_train_test/')

plot_validation_boxplot(directory_paths, 'corr', 'per_models_varSize', save_name = 'unMutatedANDlength')

path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0/GBM0_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100/GBM_longer100_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100/GBM0_longer100_assessments.tsv")


ass_type <- 'corr'
grouped_barPlot_model_elements(path_ass_full, ass_type, 'per_models_varSize', save_name = 'unMutatedANDlength')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig3.Down sampling figures ##### 
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



ass_type = 'corr'
DownSampling_eval(paths_eval_DS, ass_type)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#####fig.4   dimension reduction figures ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

path_ass_orig <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/nn_poisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/nn_mseLoss_assessments.tsv")

path_ass_dimReds_pca <- c("../external/BMR/output/dimReduction_effect/PCA/GBM/GBM_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/RF/RF_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/nn_poisLoss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/nn_mseLoss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/GLM/GLM_PCA_rm_nonMutated_assessments.tsv")

path_ass_dimReds_ae <- c("../external/BMR/output/dimReduction_effect/AE/GBM/GBM_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/RF/RF_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/nn_poisLoss/nn_poisLoss_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/nn_mseLoss/nn_mseLoss_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/GLM/GLM_PCA_rm_nonMutated_assessments.tsv")


path_ass_full <- c(path_ass_orig, path_ass_dimReds_pca, path_ass_dimReds_ae)

ass_type <- 'corr'

grouped_barPlot_dimReduction(path_ass_full, ass_type)

########################
directory_paths_PCA <- c(
  "../external/BMR/output/dimReduction_effect/PCA/GBM/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/RF/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/rep_train_test/",
  "../external/BMR/output/dimReduction_effect/PCA/GLM/rep_train_test/"
)

metric <- 'corr'

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
grouped_boxplot_deimReduction(directory_paths, metric)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig5.bin size effect figures (fixed-concatenated intervals) ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_binSizes <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/GBM_assessments.tsv'
                       # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1k/GBM/GBM_assessments.tsv'
)

ass_type <- 'corr'
save_name = 'concatenated'
grouped_barPlot_model_elements(path_ass_binSizes, ass_type, 'per_binSize', save_name = save_name)


binEffect_directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/rep_train_test/'
                               # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1k/GBM/rep_train_test/'
)

plot_validation_boxplot(binEffect_directory_paths, 'corr', 'per_binSize', save_name = save_name) 


#####     fig5. bin size effect figures (fixed-window intervals) #####     

path_ass_binSizes <- c(#'../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/1M/GBM/GBM_assessments.tsv',
                       # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/100k/GBM/GBM_assessments.tsv',
                       # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/50k/GBM/GBM_assessments.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/10k/GBM/GBM_assessments.tsv'
)

save_name = 'splitted'
ass_type <- 'corr'
grouped_barPlot_model_elements(path_ass_binSizes, ass_type, 'per_binSize', save_name = save_name)


binEffect_directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/rep_train_test/',
                               # '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/rep_train_test/',
                               '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/rep_train_test/'
)

plot_validation_boxplot(binEffect_directory_paths, 'corr', 'per_binSize', save_name = save_name) 





