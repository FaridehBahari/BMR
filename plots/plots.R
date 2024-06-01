# (fill="#69b3a2", alpha=0.5)
# https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=5

rm(list = ls())

library(dplyr)
library(reshape2)
library(gridExtra)
library(ggrepel)
library(scales)
library(ggplot2)
library(cowplot)
library(data.table)
library(grid)



# setwd('C:/Active/projects/make_features/BMR/')

create_method_colours <- function() {
  
  c( `NN (MSE loss)` = '#e41a1c' , `NN (Poisson loss)` = '#984ea3', GLM = 'grey',
     Intergenic = '#386cb0', XGBoost = '#386cb0', 
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




create_element_colours <- function() {
  
  c(enhancers = '#b03914', Enhancers = '#b03914',
    `3\' UTR` = "#14b556", gc19_pc.3utr = "#14b556", 
    `5\' UTR` = "#AAAA44", gc19_pc.5utr = "#AAAA44",
    CDS = "#CC99BB", gc19_pc.cds = "#CC99BB",
    `Core Promoter` = "#77AADD", gc19_pc.promCore = "#77AADD", 
    `Splice site` = "#11266d", gc19_pc.ss = "#11266d", 
    lncRNA = '#daa520', lncrna.ncrna = '#daa520', 
    `lncRNA Promoter` = "#DDAA77", lncrna.promCore = "#DDAA77",
    Intergenic = "#AFBCCB",
    `Intergenic with overlap` = "#E8CEC1"
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

extract_dimReducMethod_from_path <- function(path){
  spl <- unlist(strsplit(path, "/"))
  model <- spl[length(spl)-2]
  if (!model %in% c('PCA', 'AE') ) {
    model <- 'All features'
  }
  model
}


plot_validation_boxplot <- function(directory_paths, metric, compare, GBM_name,
                                    path_save = '../external/BMR/plots/',
                                    save_name = '') {
  
  dir.create(path_save, showWarnings = F, recursive = T)
  
  if (metric == 'corr') {
    Y_lable = 'Correlation'
  } else if (metric == 'acc') {
    Y_lable = 'Accuracy'
  } else if (metric == 'mse') {
    Y_lable = '-log10(MSE)'
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
          if (model == 'var_size' | grepl('GBM', model)) {
            model = GBM_name
          } 
          
        } else if (compare == 'per_models_varSize') {
          model = extract_model_from_path(directory_path)
          
          if ( grepl('GBM', model)) {
            model = GBM_name
          } 
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
    
    if (metric == 'mse') {
      combined_data$Value <- -log10(combined_data$Value)
    } 
    
    # Replace values in the Model column
    combined_data <- combined_data %>%
      mutate(Model = case_when(
        Model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
        Model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
        TRUE ~ Model  # Keep the original value if it doesn't match any condition
      ))
    ordered <- unique(combined_data %>%
                        arrange(desc(Value)) %>%
                        pull(Model))
    
    
    
    # Use the ordered metrics as a factor level for model
    combined_data$Model <- factor(combined_data$Model, levels = ordered)
    
    
    
    # Use ggplot2 to plot the grouped box plot
    ggplot(combined_data, aes(x = Model, y = Value, fill = Model)) +
      geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, show.legend = F) +  # Adjust the width as needed
      labs(y = Y_lable, x = '') +
      #, title = " performance of different models on intergenic validation sets using variable-size bins ") +
      
      theme(legend.position="right")+
      scale_fill_manual(values = create_method_colours()) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            # legend.title = element_blank(),
            text = element_text(size = 15),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.line = element_line(colour = "black"))
    
    ggsave(paste0(metric, "_validation_boxplot", save_name,".png"),
           device = "png", width = 4, height = 5,
           bg = 'white',
           path = path_save)
  } else {
    print(paste("No data found for metric '", metric, "' in the specified directories."))
  }
}


# Helper function to extract the legend from a ggplot object
extract_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[which(sapply(g$grobs, function(x) x$name) == "guide-box")][[1]]
  return(legend)
}

grouped_barPlot_model_elements <- function(path_ass_full, ass_type, compare, GBM_name = '',
                                           path_save = '../external/BMR/plots/',
                                           save_name = '', model_names = NULL, manual_limit = NULL) {
  dir.create(path_save, showWarnings = FALSE, recursive = TRUE)
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
    Y_label = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_label = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_label = '-log10(MSE)'
  } 
  
  elem_dfs <- list()
  
  for (element in elements) {
    assessment <- extract_ass_methods(path_ass_full, ass_type, element)
    model <- c()
    for (path in path_ass_full) {
      if (compare == 'per_binSize') {
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
  elem_dfs <- elem_dfs %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ GBM_name,
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  elem_dfs <- elem_dfs %>%
    arrange(element, desc(assessment))
  elem_dfs <- elem_dfs[!duplicated(elem_dfs),]
  elem_dfs <- elem_dfs %>%
    group_by(element) %>%
    mutate(id = row_number())
  
  elem_dfs$assessment <- as.numeric(elem_dfs$assessment)
  
  if (ass_type == 'mse') {
    elem_dfs$assessment <- -log10(elem_dfs$assessment)
  }
  
  # Create the bar plot without legend
  plot <- ggplot(elem_dfs, aes(x = element, y = assessment,
                               fill = model, color = model, group = interaction(model, id))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black", show.legend = FALSE) +  # Specify the border color
     
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          text = element_text(size = 15),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +  # Set the interior colors
    scale_color_manual(values = "black") +  # Set the border color
    labs(y = Y_label, x = '')  # Labels
  
  
  if (manual_limit) {
    plot <- plot + ylim(9, max(1, elem_dfs$assessment))
  } else {
    plot <- plot + ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment))
  }
  
  ggsave(paste0(ass_type, "grouped_barPlot_model_elements_", save_name, ".png"),
         plot = plot,
         device = "png", width = 5, height = 5,
         bg = 'white',
         path = path_save)
  
  # Extract and save the legend separately
  legend_plot <- ggplot(elem_dfs, aes(x = element, y = assessment,
                                      fill = model, color = model, group = interaction(model, id))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black", show.legend = TRUE) +  # Specify the border color
    ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment)) + 
    theme_minimal() +
    theme(legend.position = "right") +
    scale_fill_manual(values = create_method_colours()) +  # Set the interior colors
    scale_color_manual(values = "black")  # Set the border color
  
  legend <- extract_legend(legend_plot)
  
  if (!is.null(legend)) {
    legend_plot <- ggdraw() + draw_grob(legend)
    ggsave(paste0("models_legend", save_name, ".png"),
           plot = legend_plot,
           device = "png", width = 2, height = 3,
           bg = 'white',
           path = path_save)
  } else {
    message("No legend found in the plot.")
  }
}


grouped_barPlot_model_elements_MSE <- function(path_ass_full, ass_type = 'mse',
                                               compare, GBM_name = '',
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
    
    df <- data.frame(cbind(model, assessment, element = element_names[element]))
    elem_dfs <- append(elem_dfs, list(df))
  }
  
  elem_dfs <- bind_rows(elem_dfs)
  elem_dfs <- elem_dfs[!duplicated(elem_dfs),]
  elem_dfs <- elem_dfs %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ GBM_name,
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  
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
    geom_bar(stat = "identity", position = "dodge", width = 0.5, 
             color = "black", show.legend = F) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(),
          legend.title = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    scale_color_manual(values = "black") +
    scale_y_continuous(trans = trans_new("pseudo_log", 
                                         transform = function(x) {-log10(x)}, 
                                         inverse = function(x) 10^(-x)),
                       breaks = c(10^{-12}, 10^{-9}, 10^{-6}, 10^{-3}),
                       labels = c(expression(10^-12), expression(10^-9), expression(10^-6), expression(10^-3))) +
    labs(y = Y_lable, fill = "Model")
  
  
  ggsave(paste0(ass_type, "grouped_barPlot_model_elements_", save_name,".png"),
         device = "png", width = 6, height = 6,
         bg = 'white',
         path = path_save)
}

grouped_barPlot_TL_elements_MSE <- function(path_ass_full, ass_type = 'mse',
                                               compare, GBM_name = '',
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
    
    df <- data.frame(cbind(model, assessment, element = element_names[element]))
    elem_dfs <- append(elem_dfs, list(df))
  }
  
  elem_dfs <- bind_rows(elem_dfs)
  elem_dfs <- elem_dfs[!duplicated(elem_dfs),]
  elem_dfs <- elem_dfs %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ GBM_name,
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  
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
    geom_bar(stat = "identity", position = "dodge", width = 0.5, 
             color = "black", show.legend = F) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),
          legend.position = "right",
          text = element_text(size = 14),
          panel.grid = element_blank(),
          legend.title = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +
    scale_color_manual(values = "black") +
    labs(y = Y_lable, fill = "Model")
  
  
  ggsave(paste0(ass_type, "grouped_barPlot_model_elements_", save_name,".png"),
         device = "png", width = 6, height = 5,
         bg = 'white',
         path = path_save)
}



DownSampling_linePlot_perModel_allElems <- function(path_ass_DS, ass_type, 
                                                    path_save = '../external/BMR/plots/',
                                                    save_name = ''){
  
  # Directory creation
  dir.create(path_save, showWarnings = FALSE, recursive = TRUE)
  
  # Define Y label based on assessment type
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = 'MSE'
  }
  
  # Extract model from path
  model <- extract_model_from_path(path_ass_DS[1])
  
  # Mapping dictionary for element names
  element_names <- c('enhancers' = 'Enhancers',
                     'gc19_pc.3utr' = '3\' UTR',
                     'gc19_pc.5utr' = '5\' UTR',
                     'gc19_pc.cds' = 'CDS',
                     'gc19_pc.promCore' = 'Core Promoter',
                     'gc19_pc.ss' = 'Splice site')
  
  elements <- names(element_names)
  
  elem_dfs <- list()
  
  # Loop through elements and create dataframes
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
  colnames(elem_dfs) <- c('sampleSize', 'assessment', 'element type')
  # Convert the sampleSize column to a factor to ensure correct ordering in the plot
  elem_dfs$sampleSize <- factor(elem_dfs$sampleSize, levels = c("FullSet", "DS1M", "DS800k", "DS600k", "DS300k", "DS100k", "DS50k"))
  elem_dfs$assessment <- as.numeric(elem_dfs$assessment)
  
  
  # Create the line plot
  plot <- ggplot(data = elem_dfs, aes(x = sampleSize, y = assessment, 
                                      group = `element type`, color = `element type`)) +
    geom_line() +
    geom_point( show.legend = FALSE) +
    labs(#title = "Reduction of Sample Size Monitoring",
      x = "",
      y = '') +
    theme_minimal()
  
  # Create a new dataframe containing only the "FullSet" and "DS50k" assessments for element
  full_ds50 <- elem_dfs
  
  # Conditionally format the assessment values based on ass_type
  if (ass_type == 'mse') {
    full_ds50$assessment <- scientific(as.numeric(full_ds50$assessment), digits = 3)
  } else {
    full_ds50$assessment <- as.character(round(as.numeric(full_ds50$assessment), 3))
  }
  
  full_ds50$assessment <- as.numeric(full_ds50$assessment)
  # Plot the line plot with "FullSet" and "DS50k" assessments for each element
  plot <- plot +
    geom_text_repel(data = full_ds50, aes(label = assessment), hjust = -0.2, vjust = 0)+
    labs(#title = paste0("Reduction of Sample Size Monitoring in Different Elements using ", model, " model"),
      x = "",
      y = Y_lable) +
    scale_color_manual(values = create_element_colours()) +
    theme_minimal() +
    theme(plot.margin = unit(c(2, 1, 1, 1), "cm"),  # Adjust the plot margins
          axis.title.y = element_text(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),  # Remove y-axis values
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          text = element_text(size = 14),
          axis.line = element_line(colour = "black"))+
    guides(color = guide_legend(override.aes = list(label = "")))  
  
  
  ggsave(paste0(ass_type, "_DownSampling_linePlot_perModel_allElems", 
                save_name,".png"),
         device = "png", width = 8,  height = 5.5,
         bg = 'white',
         path = path_save)
  
  print(plot)
  
}


prepare_df_stat <- function(paths_eval_DS, statistic){ # statistic can be 'Mean' or 'Variance'
  
  df <- do.call(rbind, lapply(paths_eval_DS, fread))
  
  sample_size <- c()
  model = c()
  for (path_eval_DS in paths_eval_DS) {
    M = extract_model_from_path(path_eval_DS)
    model = c(model, M)
    sample_size <- c(sample_size, extract_binSize_from_path(path_eval_DS))
    
  }
  
  df_mean <- df[which(df$V1 == statistic),]
  df_mean$model <- model
  df_mean$sampleSize <- sample_size
  df_mean <- data.frame(df_mean)
  df_ass_mean <- df_mean[,c('V1', ass_type, "model", "sampleSize")]
  # df_ass_mean <- melt(df_ass_mean, id.vars = c(ass_type,  "model", "sampleSize"))
  
  df_ass_mean
}


prepare_data_DS_val <- function(paths_eval_DS, ass_type){
  
  
  df_mean <- data.frame(prepare_df_stat(paths_eval_DS, 'Mean'))
  df_var <- data.frame(prepare_df_stat(paths_eval_DS, 'Variance'))
  
  
  df_mean$M_minus_sd <- df_mean[,ass_type] - sqrt(df_var[,ass_type])
  df_mean$M_plus_sd <- df_mean[,ass_type] + sqrt(df_var[,ass_type])
  df_mean$assessment <- ass_type
  df_mean$sampleSize <- factor(df_mean$sampleSize, levels = c("FullSet", "DS1M", 
                                                              "DS800k", "DS600k", 
                                                              "DS300k", "DS100k", 
                                                              "DS50k"))
  
  df_mean
}


DownSampling_eval <- function(paths_eval_DS, ass_type,GBM_name,
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
  
  df_mean <- df_mean %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ GBM_name,
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  ggplot() +
    geom_line(data = df_mean, aes(x = sampleSize, y = get(ass_type), color = model,
                                  group = model)) +
    # geom_point(data = df_mean, aes(x = sampleSize, y = corr, color = model), size = 3) +
    geom_errorbar(data = df_mean, aes(x = sampleSize, ymin = M_minus_sd, ymax = M_plus_sd,
                                      color = model), width = 0.2) +
    #facet_wrap(~assessment, scales = "free_y") +
    scale_color_manual(values = create_method_colours()) +
    labs(x = '', y = Y_lable, title = "") +
    theme_minimal() +
    
    theme(axis.title.y = element_text(),
          legend.title = element_blank(),
          axis.title.x = element_text(),
          axis.text.y = element_text(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          text = element_text(size = 15),
          axis.line = element_line(colour = "black"))
  
  ggsave(paste0(ass_type, "_DownSampling_eval",
                save_name,".png"),
         device = "png", width = 8, height = 5.5,
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
    y_axis <- 'Decrease in Average Correlation'
    font_size = 12
  } else {
    y_axis <- 'Correlation'
    font_size = 14
  }
  
  plot_data$Feature.Group <- gsub("NucleotideContext", "Nucleotide Content", plot_data$Feature.Group)
  plot_data$Feature.Group <- factor(plot_data$Feature.Group, 
                                    levels = plot_data$Feature.Group[order(-plot_data$Correlation)])
  
  # Plot using ggplot2
  p <- ggplot(plot_data, aes(x = Feature.Group, y = Correlation, fill = Feature.Group)) +
    geom_bar(stat = "identity", show.legend=FALSE) +
    labs(y = "Correlation", fill = NULL) +
    
    scale_fill_manual(values = c(`Epigenetic mark` = "#66C2A5",
                                 `RNA expression` = "#3288BD",
                                 `Replication timing`= "#43589F",
                                 HiC = "#D53E4F",
                                 `DNA accessibility` = "#ABDDA4",
                                 `Nucleotide Content` = "#D4B9DA",
                                 Conservation = "#D1E5F0"))  +
    #   'black', "#43589F", '#7393B3', "#3288BD", "#66C2A5", "#ABDDA4", 
    # "#E6F598", "#FDAE61", "#D53E4F",
    # "#9E0142"
    # c('#081d58','#253494','#225ea8','#1d91c0','#41b6c4','#c7e9b4','#edf8b1')
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = font_size),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black"),
          # legend.title = element_blank()
          ) +
    xlab("")+
    ylab(y_axis)
  
  if (!LOO) {
    p <- p + geom_hline(yintercept = full_corr_gbm, linetype = "dashed", linewidth = 1.3)+
      annotate("text", x = Inf, y = full_corr_gbm, label = "Intergenic",
               hjust = 1, vjust = -0.5, size = 3.5)
  }
  
  ggsave(paste0("FeatureImportance", save_name,".png"),
         plot = p, device = "png", width = 4.5, height = 4.25,
         bg = 'white',
         path = path_save)
  
}

perElem_importance_heatmap <- function(path_ratio_perElemGroups,
                                       path_save = '../external/BMR/plots/',
                                       save_name = ''){
  dir.create(path_save, showWarnings = F, recursive = T)
  
  data <- fread(path_ratio_perElemGroups)
  data <- data[, c(2, 3,4)]
  data <- data[which(!data$Element %in% c('lncrna.ncrna', 'lncrna.promCore')),]
  
  data$feature_category = sapply(data$feature_category, function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
  data <- data[which(!data$feature_category %in% c('DNA_methylation', 'APOBEC')),]
  data$feature_category <- gsub('_', " ", data$feature_category)
  
  
  colnames(data) <- c("Element type",
                      "Importance per feature category \n (correlation ratio)",
                      "Feature Category")
  
  element_names <- c('Enhancers', '3\' UTR', '5\' UTR', 'CDS', 'Core Promoter', 'Splice site')
  names(element_names) <- c('enhancers', 'gc19_pc.3utr', 'gc19_pc.5utr', 'gc19_pc.cds', 'gc19_pc.promCore', 'gc19_pc.ss')
  
  data$`Element type` <- element_names[data$`Element type`]
  
  
  data$`Element type` <- factor(data$`Element type`)
  
  
  
  # Plot the heatmap
  ggplot(data, aes(x = `Feature Category`, y = `Element type`, fill = `Importance per feature category \n (correlation ratio)`)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(`Importance per feature category \n (correlation ratio)`, 2)), 
              color = "black", size = 3.5)  +  # Add text labels
    # scale_fill_gradientn(colors = c( 'black', "#43589F", '#7393B3', "#3288BD", "#66C2A5", "#ABDDA4", 
    #                                  "#E6F598", "#FDAE61", "#D53E4F",
    #                                  "#9E0142"),# rev(brewer.pal(6, "Spectral")),
    #                      breaks = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    #                      labels = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
    #                      limits = c(0, 1),
    #                      na.value = "grey50",
    #                      guide = "colorbar") +
    scale_fill_gradientn(colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE" ,
                                    "#D1E5F0", "#E7E1EF", "#D4B9DA" , "#CE1256", "#980043" ),
                         breaks = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                         labels = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                         limits = c(0, 1),
                         na.value = "grey50",
                         guide = "colorbar") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15),) +
    #labs(x = "Feature Category", y = "Element type")
    labs(x = "", y = "")
  
  ggsave(paste0("perElement_variableImportance_heatmap", 
                save_name,".png"),
         device = "png", width = 8, height = 5.5,
         bg = 'white',
         path = path_save)
  
}




perElem_importance_heatmap <- function(path_ratio_perElemGroups,
                                       path_save = '../external/BMR/plots/',
                                       save_name = '') {
  dir.create(path_save, showWarnings = F, recursive = T)
  
  data <- fread(path_ratio_perElemGroups)
  data <- data[, c(2, 3,4)]
  data <- data[which(!data$Element %in% c('lncrna.ncrna', 'lncrna.promCore')),]
  
  data$feature_category = sapply(data$feature_category, function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
  data <- data[which(!data$feature_category %in% c('DNA_methylation', 'APOBEC')),]
  data$feature_category <- gsub('_', " ", data$feature_category)
  
  colnames(data) <- c("Element type",
                      "Importance per feature category \n (correlation ratio)",
                      "Feature Category")
  
  element_names <- c('Enhancers', '3\' UTR', '5\' UTR', 'CDS', 'Core Promoter', 'Splice site')
  names(element_names) <- c('enhancers', 'gc19_pc.3utr', 'gc19_pc.5utr', 'gc19_pc.cds', 'gc19_pc.promCore', 'gc19_pc.ss')
  
  data$`Element type` <- element_names[data$`Element type`]
  
  data$`Element type` <- factor(data$`Element type`)
  
  # Create the heatmap
  heatmap <- ggplot(data, aes(x = `Feature Category`, y = `Element type`, fill = `Importance per feature category \n (correlation ratio)`)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(`Importance per feature category \n (correlation ratio)`, 2)), 
              color = "black", size = 3.5) +
    scale_fill_gradientn(colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE" ,
                                    "#D1E5F0", "#E7E1EF", "#D4B9DA" , "#CE1256", "#980043"),
                         breaks = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                         labels = c(0, .1, .2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                         limits = c(0, 1),
                         na.value = "grey50",
                         guide = "colorbar") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 15)) +
    labs(x = "", y = "")
  
  # Extract the legend
  legend <- get_legend(heatmap)
  
  # Remove the legend from the heatmap
  heatmap <- heatmap + theme(legend.position = "none")
  
  # Save the heatmap
  ggsave(paste0("perElement_variableImportance_heatmap", save_name, ".png"),
         plot = heatmap,
         device = "png", width = 5, height = 4,
         bg = 'white',
         path = path_save)
  
  # Save the legend
  legend_plot <- ggdraw() + draw_grob(legend)
  ggsave(paste0("perElement_variableImportance_legend", save_name, ".png"),
         plot = legend_plot,
         device = "png", width = 4, height = 4,
         bg = 'white',
         path = path_save)
}


save_nTP_nHit_plot <- function(path_driver_performance, GBM_name,
                               path_save = '../external/BMR/plots/',
                               save_name = ''){
  dir.create(path_save, showWarnings = F, recursive = T)
  
  df <- fread(path_driver_performance)
  data <- df[which(df$method %in% c('GBM', 'eMET')), c('CDS_nTPs', 'CDS_nHits',
                                                       'NC_nTPs', 'NC_nHits')]
  
  data$CDS_nHits <- data$CDS_nHits -data$CDS_nTPs
  data$NC_nHits <- data$NC_nHits - data$NC_nTPs
  data$method = c(GBM_name, 'eMET') 
  
  # Reshape data for plotting
  data_long <- tidyr::pivot_longer(data, cols = -method, names_to = "Variable", values_to = "Value")
  
  data_long$facet = unlist(lapply(strsplit( data_long$Variable, "_"), function(s){
    
    if (s[1] == 'NC') {
      'Non-coding'
    } else {
      s[1]
    }
    
  }))
  data_long$Variable <- gsub('CDS_', '', data_long$Variable)
  data_long$Variable <- gsub('NC_', '', data_long$Variable)
  
  colnames(data_long) <- c('method', 'Variable', 'Count', 'facet')
  
  data_long <- data_long %>%
    mutate(Variable = case_when(
      Variable == 'nTPs' ~ '#Reported True Positives',
      Variable == 'nHits' ~ '#Unreported Candidates'
    ))
  data_long$Variable <- factor(data_long$Variable, levels = c('#Unreported Candidates',
                                                              '#Reported True Positives'))
  ggplot(data_long, aes(x = method, y = Count, fill = Variable)) + 
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~ facet) +
    theme_minimal() +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5), 
              color = "white", size = 4) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x = element_blank(),  # Remove x-axis label
          legend.position = "top",
          text = element_text(size = 15),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = c(`#Reported True Positives` = "#b30000", 
                                 `#Unreported Candidates` = "#999999")) +
    labs(fill = NULL)  # Remove legend title
  
  ggsave(paste0("nTPs_nHits", 
                save_name,".png"),
         device = "png", width = 5.2, height = 4.8,
         bg = 'white',
         path = path_save)
}



length_distribution_perElement <- function(path_tetsY, path_intergenic,
                                           path_save = '../external/BMR/plots/',
                                           save_name = ''){
  dir.create(path_save, showWarnings = F, recursive = T)
  
  
  # y2 <- fread(path_intergenic2)
  # y2 <- y2[which(y2$nMut != 0),]
  # y2$element_type <- 'Intergenic with overlap'
  
  
  y1 <- fread(path_intergenic)
  y1 <- y1[which(y1$nMut != 0),]
  
  y1$element_type <- 'Intergenic'
  
  y <- fread(path_tetsY)
  y$element_type <- define_element_type(y$binID)
  y <- y[which(!y$element_type %in% c('lncrna.promCore', 'lncrna.ncrna')),]
  
  element_names <- c('Enhancers', '3\' UTR', '5\' UTR', 'CDS', 'Core Promoter', 'Splice site')
  names(element_names) <- c('enhancers', 'gc19_pc.3utr', 'gc19_pc.5utr', 'gc19_pc.cds', 'gc19_pc.promCore', 'gc19_pc.ss')
  
  y$element_type <- element_names[y$element_type]
  
  y <- rbind(y, y1)
  # y <- rbind(y, y2)
  
  
  ggplot(y, aes(x = element_type, y = length, fill = element_type)) +
    geom_violin(trim = FALSE) +
    labs(x = "", # Removed the title, set x-label to empty
         y = "Length") +
    theme_minimal() +
    scale_fill_manual(values = create_element_colours()) + # Changed to scale_fill_manual
    scale_y_log10(limits = c(1, max(y$length, na.rm = TRUE))) + # Safe handling of NA values
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),
          text = element_text(size = 14),
          axis.line = element_line(colour = "black"))
  
  ggsave(paste0("lengthDistribution", 
                save_name,".png"),
         device = "png", width = 8, height = 6,
         bg = 'white',
         path = path_save)
}


# grouped_barPlot_dimReduction <- function(path_ass_full, ass_type, GBM_name,
#                                          path_save = '../external/BMR/plots/',
#                                          save_name = ''){
#   
#   dir.create(path_save, showWarnings = F, recursive = T)
#   if (ass_type == 'corr') {
#     Y_lable = 'Correlation'
#   } else if (ass_type == 'acc') {
#     Y_lable = 'Accuracy'
#   } else if (ass_type == 'mse') {
#     Y_lable = 'MSE'
#   }
#   
#   element_names <- c('enhancers' = 'Enhancers',
#                      'gc19_pc.3utr' = '3\' UTR',
#                      'gc19_pc.5utr' = '5\' UTR',
#                      'gc19_pc.cds' = 'CDS',
#                      'gc19_pc.promCore' = 'Core Promoter',
#                      'gc19_pc.ss' = 'Splice site'
#                      # , 'lncrna.ncrna' = 'lncRNA',
#                      # 'lncrna.promCore' = 'lncRNA Promoter'
#   )
#   
#   elements <- names(element_names)
#   
#   
#   elem_dfs <- list()
#   
#   for (element in elements) {
#     
#     assessment <- extract_ass_methods(path_ass_full, ass_type, element)
#     
#     
#     model <- c()
#     dim_red <- c()
#     for (path in path_ass_full) {
#       
#       m = extract_model_from_path(path)
#       model = c(model, m)
#       print(model)
#       dimReduc <- extract_dimReducMethod_from_path(path)
#       dim_red <- c(dim_red, dimReduc)
#       
#     }
#     
#     df <- data.frame(cbind(model, dim_red, assessment, element = element_names[element]))
#     elem_dfs <- append(elem_dfs, list(df))
#   }
#   
#   elem_dfs <- bind_rows(elem_dfs)
#   elem_dfs$dim_red <- factor(elem_dfs$dim_red, levels = c('Original', 'PCA', 'AE'))
#   
#   # Order the data first by element and then by decreasing assessment
#   elem_dfs <- elem_dfs %>%
#     arrange(element, desc(as.numeric(assessment)))
#   
#   # Create groups for laying out the barplot correctly
#   elem_dfs <- elem_dfs %>%
#     group_by(element) %>%
#     mutate(id = row_number())
#   
#   # Create the circular or radial barplot
#   
#   # Convert 'assessment' to numeric
#   elem_dfs$assessment <- as.numeric(elem_dfs$assessment)
#   
#   # Reorder elements based on the model variable
#   elem_dfs <- elem_dfs %>%
#     arrange(model, element)
#   
#   elem_dfs <- elem_dfs %>%
#     mutate(model = case_when(
#       grepl('GBM', model) ~ GBM_name,
#       model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
#       model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
#       TRUE ~ model  # Keep the original value if it doesn't match any condition
#     ))
#   
#   # Plot the grouped bar plot
#   ggplot(elem_dfs, aes(x = model, y = assessment, fill = dim_red)) +
#     geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black") +
#     ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment)) +
#     theme_classic() +  # Using theme_classic() as base theme
#     facet_wrap(~element, strip.position = "bottom") +
#     scale_fill_manual(values = create_method_colours()) +
#     scale_alpha_manual(values = c(Original = 1, PCA = 0.6, AE = 0.2)) +
#     labs(y = Y_lable, x = '') +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#       axis.title.x = element_blank(),
#       legend.position = "right",
#       text = element_text(size = 14),
#       panel.grid = element_blank(),
#       axis.line = element_line(colour = "black"),
#       strip.placement = "outside"
#     )
#   
#   
#   ggsave(paste0(ass_type, "grouped_barPlot_dimReduction_", save_name,".png"),
#          device = "png", width = 10, 
#          bg = 'white',
#          path = path_save)
# }

grouped_barPlot_dimReduction <- function(path_ass_full, paths_eval_DimRed, ass_type, GBM_name,
                                         path_save = '../external/BMR/plots/',
                                         save_name = ''){
  
  dir.create(path_save, showWarnings = F, recursive = T)
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = '-log10(MSE)'
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
    arrange(element, desc(as.numeric(assessment)))
  
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
  
  elem_dfs <- elem_dfs %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ GBM_name,
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  intergenic_mean <- data.frame(prepare_df_stat(paths_eval_DimRed, 'Mean'))
  
  colnames(intergenic_mean) <- c("assessment",       ass_type , "model", "dim_red")
  selected_cols <- c("model", "dim_red", ass_type )
  intergenic_mean <- intergenic_mean[,selected_cols]
  colnames(intergenic_mean) <- c("model", "dim_red", "assessment" )
  
  intergenic_mean$dim_red <- ifelse(intergenic_mean$dim_red == 'var_size', 
                                    'All features', intergenic_mean$dim_red )
  
  intergenic_mean <- intergenic_mean %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ 'XGBoost',
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  intergenic_mean$element <- 'Intergenic'
  intergenic_mean$id <- 20
  
  elem_dfs <- rbind(elem_dfs, intergenic_mean)
  elem_dfs$dim_red <- factor(elem_dfs$dim_red, levels = c('All features', 'PCA', 'AE'))
  
  if (ass_type == 'mse') {
    elem_dfs$assessment <- -log10(elem_dfs$assessment)
  }
  
  
  # Plot the grouped bar plot
  ggplot(elem_dfs, aes(x = model, y = assessment, fill = dim_red)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black") +
    ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment)) +
    theme_classic() +  # Using theme_classic() as base theme
    facet_wrap(~element, strip.position = "bottom", ncol = 4) +
    scale_fill_manual(values = create_method_colours()) +
    labs(y = Y_lable, x = '') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(), 
      text = element_text(size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.placement = "outside"
    )
  
  
  ggsave(paste0(ass_type, "NEWgrouped_barPlot_dimReduction_", save_name,".png"),
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}

plot_validation_boxplot_lengthMut <- function(directory_paths, metric, compare,
                                              path_save = '../external/BMR/plots/',
                                              save_name = '') {
  
  dir.create(path_save, showWarnings = F, recursive = T)
  
  if (metric == 'corr') {
    Y_lable = 'Correlation'
  } else if (metric == 'acc') {
    Y_lable = 'Accuracy'
  } else if (metric == 'mse') {
    Y_lable = '-log10(MSE)'
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
    model_names <- c('All',
                     '#mutations >= 1 & Length > 100',
                     'Length > 100', 
                     '#mutations >= 1')
    names(model_names) <- c('GBM0_rep', 'GBM_longer100_rep', 'GBM0_longer100_rep', 'GBM_rep')
    # names(model_names) <- c('GBM0', 'GBM_longer100', 'GBM0_longer100', 'GBM')
    
    combined_data$Model <- model_names[combined_data$Model]
    
    if (ass_type == 'mse') {
      combined_data$Value <- -log10(combined_data$Value)
    }
    
    ordered <- unique(combined_data %>%
                        arrange(desc(Value)) %>%
                        pull(Model))
    # Use the ordered metrics as a factor level for model
    combined_data$Model <- factor(combined_data$Model, levels = ordered)
    
    # Use ggplot2 to plot the grouped box plot
    ggplot(combined_data, aes(x = Model, y = Value, fill = Model)) +
      geom_boxplot(position = position_dodge(width = 0.6), width = 0.5, show.legend = F) +  # Adjust the width as needed
      labs(y = Y_lable, x = '') +
      #, title = " performance of different models on intergenic validation sets using variable-size bins ") +
      
      theme(legend.position="right")+
      scale_fill_manual(values = create_method_colours()) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.title = element_blank(),
            text = element_text(size = 16),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.line = element_line(colour = "black"))
    
    ggsave(paste0(metric, "_validation_boxplot", save_name,".png"),
           device = "png", width = 6,
           bg = 'white',
           path = path_save)
  } else {
    print(paste("No data found for metric '", metric, "' in the specified directories."))
  }
}


grouped_barPlot_model_elements_lengthMut <- function(path_ass_full, ass_type, compare,
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
    Y_lable = '-log10(MSE)'
  } 
  
  
  elem_dfs <- list()
  
  for (element in elements) {
    
    assessment <- extract_ass_methods(path_ass_full, ass_type, element)
    
    
    model <- c()
    for (path in path_ass_full) {
      
      
      if(compare == 'per_binSize'){
        m = extract_binSize_from_path(path)
        
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
  
  model_names <- c('All',
                   '#mutations >= 1 & Length > 100',
                   'Length > 100', 
                   '#mutations >= 1')
  names(model_names) <- c('GBM0', 'GBM_longer100', 'GBM0_longer100', 'GBM')
  # names(model_names) <- c('GBM0_rep', 'GBM_longer100_rep', 'GBM0_longer100_rep', 'GBM_rep')
  
  elem_dfs$model <- model_names[elem_dfs$model]
  
  
  
  # Order the data first by element and then by decreasing assessment
  elem_dfs <- elem_dfs %>%
    arrange(element, desc(assessment))
  elem_dfs <- elem_dfs[!duplicated(elem_dfs),]
  
  # Create groups for laying out the barplot correctly
  elem_dfs <- elem_dfs %>%
    group_by(element) %>%
    mutate(id = row_number())
  
  # Create the circular or radial barplot
  
  # Convert 'assessment' to numeric
  elem_dfs$assessment <- as.numeric(elem_dfs$assessment)
  
  ggplot(elem_dfs, aes(x = element, y = assessment,
                       fill = model, color = model, group = interaction(model, id))) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black") +  # Specify the border color
    ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "right",
          text = element_text(size = 14),
          legend.title = element_blank(),
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = create_method_colours()) +  # Set the interior colors
    scale_color_manual(values = "black") +  # Set the border color
    labs(y = Y_lable, x = '')  # Labels
  
  
  ggsave(paste0(ass_type, "grouped_barPlot_model_elements_", save_name,".png"),
         device = "png", width = 7, 
         bg = 'white',
         path = path_save)
}

Bin_model_comparson_plot <- function(path_ass_full, paths_interegnics, ass_type, GBM_name,
                                     path_save = '../external/BMR/plots/',
                                     save_name = ''){
  
  dir.create(path_save, showWarnings = F, recursive = T)
  if (ass_type == 'corr') {
    Y_lable = 'Correlation'
  } else if (ass_type == 'acc') {
    Y_lable = 'Accuracy'
  } else if (ass_type == 'mse') {
    Y_lable = '-log10(MSE)'
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
    Bins <- c()
    for (path in path_ass_full) {
      
      m = extract_model_from_path(path)
      model = c(model, m)
      
      binSize <- extract_binSize_from_path(path)
      Bins <- c(Bins, binSize)
      
    }
    
    df <- data.frame(cbind(model, Bins, assessment, element = element_names[element]))
    elem_dfs <- append(elem_dfs, list(df))
  }
  
  elem_dfs <- bind_rows(elem_dfs)
  
  elem_dfs <- elem_dfs %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ GBM_name,
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  elem_dfs <- elem_dfs %>%
    mutate(Bins = case_when(
      (Bins == 'var_size') ~ 'Variable size',
      TRUE ~ Bins  # Keep the original value if it doesn't match any condition
    ))
  
  # Order the data first by element and then by decreasing assessment
  elem_dfs <- elem_dfs %>%
    arrange(element, desc(as.numeric(assessment)))
  
  # Create groups for laying out the barplot correctly
  elem_dfs <- elem_dfs %>%
    group_by(element) %>%
    mutate(id = row_number())
  
  intergenic_mean <- data.frame(prepare_df_stat(paths_interegnics, 'Mean'))
  
  colnames(intergenic_mean) <- c("assessment",       ass_type , "model", "Bins")
  selected_cols <- c("model", "Bins", ass_type )
  intergenic_mean <- intergenic_mean[,selected_cols]
  colnames(intergenic_mean) <- c("model", "Bins", "assessment" )
  
  intergenic_mean$Bins <- ifelse(intergenic_mean$Bins == 'var_size', 
                                 'Variable size', intergenic_mean$Bins )
  
  intergenic_mean <- intergenic_mean %>%
    mutate(model = case_when(
      grepl('GBM', model) ~ 'XGBoost',
      model %in% c('nn_poisLoss', 'NN_PoisLoss') ~ 'NN (Poisson loss)',
      model %in% c('nn_mseLoss', 'NN_MSELoss') ~ 'NN (MSE loss)',
      TRUE ~ model  # Keep the original value if it doesn't match any condition
    ))
  
  intergenic_mean$element <- 'Intergenic'
  intergenic_mean$id <- 20
  
  # Convert 'assessment' to numeric
  elem_dfs$assessment <- as.numeric(elem_dfs$assessment)
  
  elem_dfs <- rbind(elem_dfs, intergenic_mean)
  elem_dfs$Bins <- factor(elem_dfs$Bins, levels = c("Variable size", "10k",
                                                    "50k", "100k", "1M"))
  
  if (ass_type == 'mse') {
    elem_dfs$assessment <- -log10(elem_dfs$assessment)
  }
  
  
  
  
  # Reorder elements based on the model variable
  elem_dfs <- elem_dfs %>%
    arrange(model, element)
  
  # Plot the grouped bar plot
  p = ggplot(elem_dfs, aes(x = model, y = assessment, fill = Bins)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black") +
    
    theme_classic() +  # Using theme_classic() as base theme
    facet_wrap(~element, strip.position = "bottom", ncol = 4) +
    scale_fill_manual(values = create_method_colours()) +
    labs(y = Y_lable, x = '') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(), 
      text = element_text(size = 14),
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.placement = "outside"
    )
  
  p = p + ylim(min(0, elem_dfs$assessment), max(1, elem_dfs$assessment))
  ggsave(paste0(ass_type, "BINsize", save_name,".png"),
         plot = p,
         device = "png", width = 10, 
         bg = 'white',
         path = path_save)
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### Fig S1. length distribution  ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_tetsY <- '../external/BMR/rawInput/responseTabs/Pan_Cancer/reg_elems.tsv'
path_intergenic <- '../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv'

length_distribution_perElement(path_tetsY, path_intergenic, save_name = 'FigS1')


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### Fig. 2a   models and bin sizes ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_full <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/nn_poisLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/nn_mseLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/GBM_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/RF/RF_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_poisLoss/nn_poisLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_mseLoss/nn_mseLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/GBM_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/RF/RF_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_poisLoss/nn_poisLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_mseLoss/nn_mseLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/GBM_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/RF/RF_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_poisLoss/nn_poisLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_mseLoss/nn_mseLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/GBM_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/RF/RF_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_poisLoss/nn_poisLoss_assessments.tsv',
                   '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_mseLoss/nn_mseLoss_assessments.tsv'
)

ass_types = c('corr', 'mse')


paths_interegnics <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/GBM/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/RF/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_poisLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/1M/nn_mseLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/GBM/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/RF/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_poisLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/100k/nn_mseLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/GBM/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/RF/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_poisLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/50k/nn_mseLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/RF/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_poisLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/nn_mseLoss/model_metrics_summary.tsv',
                       '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/10k/GBM/model_metrics_summary.tsv')




for (ass_type in ass_types) {
  if (ass_type == 'mse') {
    save_name = 'FigS2_a'
  } else if (ass_type == 'corr') {
    save_name = 'Fig2_a'
  }
  Bin_model_comparson_plot(path_ass_full, paths_interegnics, ass_type,
                           save_name = save_name, GBM_name = 'XGBoost')
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### Fig. 2bc Down sampling figures ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_DS <- c('../external/BMR/output/with_RepliSeq_HiC/DownSampling/FullSet/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS1M/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS800k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS600k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS300k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS100k/GBM/GBM_assessments.tsv',
                 '../external/BMR/output/with_RepliSeq_HiC/DownSampling/DS50k/GBM/GBM_assessments.tsv')

for (ass_type in ass_types) {
  if (ass_type == 'mse') {
    save_name = 'FigS2_c'
  } else if (ass_type == 'corr') {
    save_name = 'Fig2_c'
  }
  DownSampling_linePlot_perModel_allElems(path_ass_DS, ass_type, save_name =save_name)
  
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


ass_types <- c('corr', 'mse')
for (ass_type in ass_types) {
  if (ass_type == 'mse') {
    save_name = 'FigS2_b'
  } else if (ass_type == 'corr') {
    save_name = 'Fig2_b'
  }
  DownSampling_eval(paths_eval_DS, ass_type, GBM_name = 'XGBoost',
                    save_name =save_name)
}


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#####     FigS3: Compare different models on variable-size intergenic bins####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
############## validation Sets:

directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/rep_train_test/')

ass_types <- c('corr', 'mse')
for (ass_type in ass_types) {
  
  plot_validation_boxplot(directory_paths, ass_type, 'per_models_varSize', 
                          GBM_name = 'XGBoost', save_name = 'FigS3') 
  
}
############## elements:

path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/nn_poisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/nn_mseLoss_assessments.tsv")


ass_types <- c('corr', 'mse')
for (ass_type in ass_types) {
  
  if (ass_type != 'mse') {
    manual_limit = F
  } else {
    manual_limit == T
  }
  
  grouped_barPlot_model_elements(path_ass_full, ass_type, 'per_models_varSize',
                                 GBM_name = 'XGBoost', manual_limit = manual_limit, save_name = 'FigS3')
}










#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig S4. mutated-unmutated bins ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#HERE#
directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_rep/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_rep/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100_rep/rep_train_test/',
                     '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100_rep/rep_train_test/')

ass_types <- c('corr', 'mse')


path_ass_full <- c("../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0/GBM0_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100/GBM_longer100_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100/GBM0_longer100_assessments.tsv")


for (ass_type in ass_types) {
  
  
  
  plot_validation_boxplot_lengthMut(directory_paths, ass_type, 'per_models_varSize', save_name = 'unMutatedANDlength_FigS4')
  grouped_barPlot_model_elements_lengthMut(path_ass_full, ass_type, 'per_models_varSize', save_name = 'unMutatedANDlength_FigS4')
  
}

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the data from CSV files
df1 <- fread('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_longer100_rep/rep_train_test/_length_intervals_assessment.tsv')
df2 <- fread('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_rep/rep_train_test/_length_intervals_assessment.tsv')
df3 <- fread('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM_longer100_rep/rep_train_test/_length_intervals_assessment.tsv')
df4 <- fread('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM0_rep/rep_train_test/_length_intervals_assessment.tsv')
# Add a column to each dataframe to identify the source
df1$source <- 'Length > 100' #'GBM0_longer100_rep'
df2$source <-  '#mutations >= 1' #'GBM_rep'
df3$source <- '#mutations >= 1 & Length > 100' #'GBM_longer100_rep'
df4$source <- 'All' # 'GBM0'

# Combine the dataframes
df_combined <- bind_rows(df1, df2, df3, df4)

# Handle missing values
df_combined <- df_combined %>%
  mutate(corr = ifelse(is.na(corr), 0, corr))

df_combined$V1 <- factor(df_combined$V1, levels = c( "0-20", "20-50", "50-100",
                                                     "100-200", "200-400","400-1500",
                                                     "1500<"))

# Create the bar plot
plot <- ggplot(df_combined, aes(x = V1, y = corr, fill = source)) +
  geom_bar(stat = 'identity', position = 'dodge', width = 0.5, color = "black", show.legend = F) +
  labs(#title = 'Comparison of Correlation Values Across Length Intervals',
    x = '',
    y = 'Mean correlation') +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = create_method_colours()) 

# Save the plot as a PNG file
ggsave('../external/BMR/plots/length_intervals_comparison_FigS4.png', bg = 'white',plot, width = 4)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig3a.Transfer learning figures ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_ass_full <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv',
                   '../external/BMR/output/element-specific/GBM/GBM_ensemble_bootstraps100_assessment.tsv',
                   '../external/BMR/output/TL/eMET/eMET_ensemble_bootstraps100_assessment.tsv')

ass_types = c('corr', 'mse')
compare <- ''
model_names <- c('Intergenic', 'element-specific', 'eMET')

for (ass_type in ass_types) {
  manual_limit = F
  if (ass_type == 'mse') {
    
    grouped_barPlot_TL_elements_MSE(path_ass_full, ass_type, compare, GBM_name = 'Intergenic',
                                    save_name = 'eMETvsGBM_MSE_realValues', model_names = model_names)
    
  } else {
    grouped_barPlot_model_elements(path_ass_full, ass_type, compare, GBM_name = 'Intergenic',
                                   save_name = 'eMETvsGBM_Fig3a', model_names = model_names, 
                                   manual_limit = manual_limit)
  }
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig3b.Transfer learning drivers ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
path_driver_performance <- '../external/BMR/benchmark_.05/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_in_oncoKBfdr.csv'
save_nTP_nHit_plot(path_driver_performance, GBM_name = 'Intergenic', save_name = 'Fig3b')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig.4   dimension reduction figures ####
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


path_ass_orig <- c("../external/BMR/output/TL/eMET/eMET_ensemble_bootstraps100_assessment.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/RF_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/nn_poisLoss_assessments.tsv",
                   "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/nn_mseLoss_assessments.tsv")

path_ass_dimReds_pca <- c("../external/BMR/output/dimReduction_effect/PCA/eMET/eMET_ensemble_bootstraps100_assessment.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/GBM/GBM_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/RF/RF_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/nn_poisLoss_assessments.tsv",
                          "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/nn_mseLoss_assessments.tsv")

path_ass_dimReds_ae <- c("../external/BMR/output/dimReduction_effect/AE/eMET/eMET_ensemble_bootstraps100_assessment.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/GBM/GBM_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/RF/RF_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/nn_poisLoss/nn_poisLoss_assessments.tsv",
                         "../external/BMR/output/dimReduction_effect/AE/nn_mseLoss/nn_mseLoss_assessments.tsv")


path_ass_full <- c(path_ass_orig, path_ass_dimReds_pca, path_ass_dimReds_ae)
paths_eval_DimRed <- c("../external/BMR/output/dimReduction_effect/AE/GBM/model_metrics_summary.tsv",
                       "../external/BMR/output/dimReduction_effect/AE/RF/model_metrics_summary.tsv",
                       "../external/BMR/output/dimReduction_effect/AE/nn_poisLoss/model_metrics_summary.tsv",
                       "../external/BMR/output/dimReduction_effect/AE/nn_mseLoss/model_metrics_summary.tsv",
                       "../external/BMR/output/dimReduction_effect/PCA/GBM/model_metrics_summary.tsv",
                       "../external/BMR/output/dimReduction_effect/PCA/RF/model_metrics_summary.tsv",
                       "../external/BMR/output/dimReduction_effect/PCA/nn_poisLoss/model_metrics_summary.tsv",
                       "../external/BMR/output/dimReduction_effect/PCA/nn_mseLoss/model_metrics_summary.tsv",
                       "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/model_metrics_summary.tsv",
                       "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/RF/model_metrics_summary.tsv",
                       "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_mseLoss/model_metrics_summary.tsv",
                       "../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/nn_poisLoss/model_metrics_summary.tsv")

ass_types <- c('corr', 'mse')
for (ass_type in ass_types) {
  
  if (ass_type == 'mse') {
    save_name = ''
  } else if (ass_type == 'corr') {
    save_name = 'Fig4'
  }
  grouped_barPlot_dimReduction(path_ass_full, paths_eval_DimRed, ass_type,
                               GBM_name = 'XGBoost', save_name = save_name)
  
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##### fig5.feature importance plots ##### 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


parent_directory <- "../external/BMR/output/featureImportance/"
path_GBM <- '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/model_metrics_summary.tsv'
importance_plot(parent_directory, path_GBM, save_name = 'Fig5a')



parent_directory <- "../external/BMR/output/featureImportance_LOO/"
path_GBM <- '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/model_metrics_summary.tsv'
importance_plot(parent_directory, path_GBM = path_GBM, save_name = 'LOO_FigS5', LOO = T)


path_ratio_perElemGroups <- '../external/BMR/output/eMET_GroupImportance/importanceRatios.csv'
perElem_importance_heatmap(path_ratio_perElemGroups, save_name = 'Fig5b')



###################################################################################


# Define parameters for the binomial distribution
n <- 50  # number of trials
p <- .5  # probability of success

# Generate the binomial distribution data
x <- 0:n
y <- dbinom(x, n, p)

# Create a data frame for plotting
data <- data.frame(x, y)

# Define the critical value for the right tail (e.g., k = 35 for demonstration)
critical_value <- 35

# Plot the binomial distribution curve
p <- ggplot(data, aes(x, y)) +
  geom_col(fill = "lightblue") +
  theme_minimal() +
  xlab("#Mutations") +
  ylab("Probability") +
  theme(
        legend.position = "right",
        text = element_text(size = 14),
        panel.grid = element_blank(),
        legend.title = element_blank(), 
        axis.line = element_line(colour = "black"))

# Add the vertical red line for the critical value
p <- p + geom_vline(xintercept = critical_value, color = "red", linetype = "dashed", size = 1) +
  geom_col(data = subset(data, x >= critical_value), aes(x, y), fill = "red", alpha = 0.5)

# Display the plot
print(p)


















# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ##### fig. bin size effect figures (fixed-window intervals) ##### 
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 
# path_ass_binSizes <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_assessments.tsv',
#                        '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/1M/GBM/GBM_assessments.tsv',
#                        '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bebedtools_splitted/100k/GBM/GBM_assessments.tsv',
#                        '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/50k/GBM/GBM_assessments.tsv',
#                        '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/10k/GBM/GBM_assessments.tsv'
# )
# 
# ass_types <- c('corr', 'mse')
# save_name = 'splitted'
# for (ass_type in ass_types) {
#   
#   if (ass_type != 'mse') {
#     grouped_barPlot_model_elements(path_ass_binSizes, ass_type, 
#                                    compare = 'per_binSize', GBM_name = 'Variable-size intergenic',
#                                    save_name = save_name)
#   } else {
#     grouped_barPlot_model_elements_MSE(path_ass_binSizes, compare = 'per_binSize', 
#                                        GBM_name = 'Variable-size intergenic',
#                                        save_name = save_name)
#   }
#   
# }
# 
# 
# binEffect_directory_paths <- c('../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/rep_train_test/',
#                                '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/1M/GBM/rep_train_test/',
#                                '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splitted/100k/GBM/GBM_assessments.tsv',
#                                '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/50k/GBM/rep_train_test/',
#                                '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/bedtools_splited/10k/GBM/rep_train_test/'
# )
# 
# 
# for (ass_type in ass_types) {
#   plot_validation_boxplot(binEffect_directory_paths, ass_type,
#                           GBM_name = 'Variable-size intergenic',
#                           compare = 'per_binSize', save_name = save_name) 
# }
# 
# 
# grouped_barPlot_model_benchmark <- function(driver_based_on, coding_nonCoding,
#                                             path_save = '../external/BMR/plots/',
#                                             save_name = ''){
#   
#   path_bench_res <- paste0('../external/BMR/benchmark_binomTest/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_', driver_based_on, 'fdr.csv')
#   dir.create(path_save, showWarnings = F, recursive = T)
#   df <- fread(path_bench_res)
#   df <- df[which(df$method %in% c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL" )),]
#   
#   # Get column indices where "CDS" is present in column names
#   cds_cols <- grep("CDS", names(df))
#   
#   # Split dataframe based on column indices
#   cds_df <- df[, c(2, cds_cols), with = FALSE]
#   nc_df <- df[, -cds_cols, with = FALSE]
#   
#   colnames(cds_df) <- gsub('CDS_', '', toupper(colnames(cds_df)))
#   colnames(nc_df) <- gsub('NC_', '', toupper(colnames(nc_df)))
#   
#   if (coding_nonCoding == 'coding') {
#     df = cds_df
#   } else if(coding_nonCoding == 'non-coding')  {
#     df = nc_df
#   }
#   
#   df = df[, -c('N_ELEMNTS', 'NTPS', 'NHITS')]
#   df <- melt(setDT(df), id.vars = "METHOD", c("PRECISIONS", "RECALLS", "F1",
#                                               "AUC", "AUPR"))
#   colnames(df) <- c('model', 'metrics', 'value')
#   
#   # Order the data first by element and then by decreasing assessment
#   df <- df %>%
#     arrange(metrics, desc(model ))
#   
#   # Create groups for laying out the barplot correctly
#   dfs <- df %>%
#     group_by(metrics) %>%
#     mutate(id = metrics)
#   
#   ggplot(dfs, aes(x = metrics, y = value, fill = model, group = interaction(model, id))) +
#     geom_bar(stat = "identity", position = "dodge", width = 0.5) +
#     # coord_polar(start = 0) +  # This turns the barplot into a circular barplot
#     ylim(min(0, dfs$value), max(1, dfs$value)) + 
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#           axis.title.x = element_blank(),
#           legend.position = "right",
#           text = element_text(size = 14),
#           panel.grid = element_blank(), 
#           axis.line = element_line(colour = "black")) +
#     scale_fill_manual(values = create_method_colours()) +
#     labs(title = paste0(coding_nonCoding, ':\ngold standard based on ', based_on),
#          y = 'value',
#          x = 'metrics',
#          fill = "Model")
#   
#   ggsave(paste0(coding_nonCoding,"_", driver_based_on, "_grouped_barPlot_model_benchmark", save_name,".png"),
#          device = "png", width = 10, 
#          bg = 'white',
#          path = path_save)
# }
# 
# 
# precisinRecall_dotPlot <- function(driver_based_on, coding_nonCoding,
#                                    path_save = '../external/BMR/plots/',
#                                    save_name = ''){
#   
#   path_bench_res <- paste0('../external/BMR/benchmark_binomTest/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_', driver_based_on, 'fdr.csv')
#   dir.create(path_save, showWarnings = F, recursive = T)
#   df <- fread(path_bench_res)
#   df <- df[which(df$method %in% c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL" )),]
#   
#   # Get column indices where "CDS" is present in column names
#   cds_cols <- grep("CDS", names(df))
#   
#   # Split dataframe based on column indices
#   cds_df <- df[, c(2, cds_cols), with = FALSE]
#   nc_df <- df[, -cds_cols, with = FALSE]
#   
#   
#   colnames(cds_df) <- gsub('CDS_', '', toupper(colnames(cds_df)))
#   colnames(nc_df) <- gsub('NC_', '', toupper(colnames(nc_df)))
#   
#   if (coding_nonCoding == 'coding') {
#     df = cds_df
#   } else if(coding_nonCoding == 'non-coding')  {
#     df = nc_df
#   }
#   
#   dfs = df
#   
#   dfs$model <- factor(df$METHOD)
#   coding_nonCoding = 'coding'
#   ggplot(dfs, aes(x = RECALLS, y = PRECISIONS)) +
#     geom_point(aes(colour = model)) +
#     # coord_polar(start = 0) +  # This turns the barplot into a circular barplot
#     ylim(min(0, dfs$PRECISIONS), max(1, dfs$PRECISIONS)) + 
#     xlim(min(0, dfs$RECALLS), max(1, dfs$RECALLS)) + 
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#           
#           legend.position = "right",
#           text = element_text(size = 14),
#           panel.grid = element_blank(), 
#           axis.line = element_line(colour = "black")) +
#     scale_colour_manual(values = create_method_colours()) +
#     labs(title = paste0(coding_nonCoding, ':\ngold standard based on ', based_on),
#          x = 'Recall', y = 'Precision'
#     )
#   
#   ggsave(paste0(coding_nonCoding,"_", driver_based_on, "_PrecisionRecall", save_name,".png"),
#          device = "png", width = 10, 
#          bg = 'white',
#          path = path_save)
# }
# 
# 
# top_100hit_barplot <- function(based_on, coding_nonCoding,
#                                path_save = '../external/BMR/plots/',
#                                save_name = ''){
#   
#   path_bench_res <- paste0('../external/BMR/benchmark_binomTest/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_',
#                            based_on, 'fixedNumberOfElems.csv')
#   
#   dir.create(path_save, showWarnings = F, recursive = T)
#   df <- fread(path_bench_res)
#   df <- df[which(df$method %in% c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL" )),]
#   
#   # Get column indices where "CDS" is present in column names
#   cds_cols <- grep("CDS", names(df))
#   
#   # Split dataframe based on column indices
#   cds_df <- df[, c(2, cds_cols), with = FALSE]
#   nc_df <- df[, -cds_cols, with = FALSE]
#   
#   colnames(cds_df) <- gsub('CDS_', '', toupper(colnames(cds_df)))
#   colnames(nc_df) <- gsub('NC_', '', toupper(colnames(nc_df)))
#   
#   if (coding_nonCoding == 'coding') {
#     df = cds_df
#   } else if(coding_nonCoding == 'non-coding')  {
#     df = nc_df
#   }
#   
#   
#   # Order the data first by element and then by decreasing assessment
#   df <- df %>%
#     arrange(METHOD , desc( NTPS))
#   
#   
#   ggplot(df, aes(x = METHOD, y = NTPS, fill = METHOD)) +
#     geom_bar(stat = "identity", width = 0.5) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#           legend.position = "right",
#           text = element_text(size = 14),
#           panel.grid = element_blank(), 
#           axis.line = element_line(colour = "black")) +
#     scale_fill_manual(values = create_method_colours()) +
#     labs(title = paste0(coding_nonCoding, ':\ngold standard based on ', based_on),
#          y = 'number of true positives',
#          x = 'model')
#   
#   ggsave(paste0(coding_nonCoding,"_", based_on, "_nTPs_barPlot", save_name,".png"),
#          device = "png", width = 10, 
#          bg = 'white',
#          path = path_save)
# }
# 
# 
# histogram_perElement_length <- function(path_responseTab_PCAWG,
#                                         path_save = '../external/BMR/plots/',
#                                         save_name = ''){
#   
#   y <- fread(path_responseTab_PCAWG)
#   
#   y$`element type` <- define_element_type(y$binID)
#   
#   y <- y[which(!y$`element type` %in% c('lncrna.ncrna', 'lncrna.promCore')), ]
#   
#   # Mapping dictionary for element types
#   element_types <- c('enhancers' = 'Enhancers',
#                      'gc19_pc.3utr' = '3\' UTR',
#                      'gc19_pc.5utr' = '5\' UTR',
#                      'gc19_pc.cds' = 'CDS',
#                      'gc19_pc.promCore' = 'Core Promoter',
#                      'gc19_pc.ss' = 'Splice site')
#   
#   # Replace the element type names in the data frame with modified names for plotting
#   y$`element type` <- element_types[y$`element type`]
#   
#   # Plot histogram
#   ggplot(y, aes(x = length, fill = `element type`)) +
#     geom_histogram(binwidth = 20) +
#     facet_wrap(~`element type`, scales = "free") +
#     theme_minimal() +
#     theme(text = element_text(size = 14),
#           panel.grid = element_blank(), 
#           axis.line = element_line(colour = "black")) +
#     scale_fill_manual(values = create_element_colours()) +  # Use the function to generate colors
#     labs(title = "Length Distribution for Each Element Type",
#          x = "Length",
#          y = "Frequency")
#   
#   ggsave(paste0(save_name,"perElement_length_histogram.png"),
#          device = "png", width = 10, 
#          bg = 'white',
#          path = path_save)
# }
# 
# 
# ######################### driver identification plots ##########################
# save_name = 'binomTest'
# coding_nonCodings <- c('coding', 'non-coding')
# based_ons <- c('all', 'in_CGC_new', 'in_pcawg', 'in_oncoKB', 'any')
# for (coding_nonCoding in coding_nonCodings) {
#   for (based_on in based_ons) {
#     grouped_barPlot_model_benchmark(based_on, coding_nonCoding, save_name = save_name)
#     precisinRecall_dotPlot(based_on, coding_nonCoding, save_name = save_name)
#     top_100hit_barplot(based_on, coding_nonCoding, save_name = save_name)
#   }
# }
# 
# ####################### PCAWG length distribution plots ########################
# path_responseTab_PCAWG <- '../external/BMR/rawInput/responseTabs/Pan_Cancer/reg_elems.tsv'
# histogram_perElement_length(path_responseTab_PCAWG)