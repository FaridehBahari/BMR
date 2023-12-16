library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library("ie2misc")

save_dotPlots <- function(path_pred, path_val, path_elem, color,
                          path_save, element){
  if (element == "intergenic") {
    elem = "^[b]+[0-9]+$"
  } else elem = element
  
  dir.create(paste0(path_save, "scatterPlot/", element, "/"),
             showWarnings = F, recursive = T)
  pred <- fread(path_pred)
  vals = fread(path_val)
  elems = fread(path_elem)
  elems$obsRates <- elems$nMut / (elems$length * elems$N)
  elems <- elems[which(elems$nMut != 0), c('binID', 'obsRates')]
  vals <- vals[which(vals$nMut != 0),c('binID', 'obsRates')]
  obs <- rbind(elems, vals)
  df <- left_join(obs, pred, by = 'binID')
  df <- df[which(df$obsRates != 0),]
  df <- df[grep(elem, df$binID),]
  
  model_name = unlist(strsplit(path_pred, "/"))[length(unlist(strsplit(path_pred, "/"))) - 1]
  
  # Extract the directory
  directory <- dirname(path_pred)
  #spear_corr_df = cor(df$obsRates, df$predRate, method = "spearman")
  MAE = mae(df$obsRates, df$predRate)
  path_assess_table = paste0(directory, "/", model_name, '_assessments.tsv')
  assessTab = data.frame(fread(path_assess_table))
  acc = assessTab[1,element]
  corr = assessTab[2,element]
  MSE = assessTab[3,element]
  
  
  max_limit = max(c(df$obsRates, df$predRate))
  min_limit = min(c(df$obsRates, df$predRate))

  ggplot(data=df, aes(x=predRate, y=obsRates))+
    # stat_poly_line() +
    # stat_poly_eq() +
    geom_point(colour = color, alpha = 0.2)+
    # scale_x_log10()+
    # scale_y_log10()+
    theme_bw()+
    labs(title = paste0(model_name, "_", element))+
    # stat_cor(method="spearman")+ 
    annotate("text", x = (min_limit+1e-7)*100, y = max_limit, label = paste0("\n", "acc:", acc, "\n",
    "corr:", corr, "\n",
    "MSE:", MSE, "\n",
    "MAE:", MAE)) +
    coord_cartesian(xlim=c(min_limit, max_limit), ylim=c(min_limit, max_limit))
    #create_method_colours(model_name)

  ggsave(paste0(element, "_", model_name, "_","obs_pred_dotPlot.png"),
         device = "png", width = 8, height = 8,
         path = paste0(path_save, "scatterPlot/", element, "/"))
}

################# box plots for comparing observed vs different method predicted rates #######
save_boxPlot <- function(element, path_preds, path_elem, path_val, path_save, save_name, include_obsRates = TRUE){
  if (element == "intergenic") {
    elem = "^[b]+[0-9]+$"
  } else elem = element
  
  all_colors <- determine_color_elem(element)
  
  
  dir.create(paste0(path_save, "boxPlot/"), showWarnings = F, recursive = T)
  
  vals = fread(path_val)
  elems = fread(path_elem)
  elems$obsRates <- elems$nMut / (elems$length * elems$N)
  elems <- elems[which(elems$nMut != 0), c('binID', 'obsRates')]
  colnames(elems) = c('binID', 'Rate')
  vals <- vals[, c('binID', 'obsRates')]
  colnames(vals) = c('binID', 'Rate')
  
  
  dfs <- fread(path_preds[1])
  colnames(dfs) = c('binID', 'Rate')
  dfs$model_name <- unlist(strsplit(path_preds[1], "/"))[length(unlist(strsplit(path_preds[1], "/"))) - 1]
  # dfs = left_join(obs, dfs, by = 'binID')
  # dfs$mse = (dfs$obsRates - dfs$predRate)^2
  # dfs$meanSE = mean(dfs$mse)
  
  obs <- rbind(elems, vals)
  obs <- obs[which(obs$Rate != 0),]
  obs$model_name <- 'observedRate'
  
  for(path_pred in path_preds[2:length(path_preds)]){
    model_name = unlist(strsplit(path_pred, "/"))[length(unlist(strsplit(path_pred, "/"))) - 1]
    tmp_pred = fread(path_pred)
    colnames(tmp_pred) = c('binID', 'Rate')
    tmp_pred$model_name = model_name
    tmp_pred = tmp_pred[which(tmp_pred$binID %in% obs$binID), ]
    dfs = rbind(dfs, tmp_pred)
    
  }
  
  if (include_obsRates){
    
    dfs = rbind(dfs, obs)
  }
  
  dfs <- dfs[grep(elem, dfs$binID),]
  
  model_names = unique(dfs$model_name)
  ggplot(data=dfs)+
    geom_boxplot(aes(x=model_name, y=Rate, color = model_name))+
    scale_y_log10()+
    theme_bw()+
    labs(title = element)+
    theme(axis.text.x = element_text(angle = 90,hjust = 1 ),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    create_method_colours(model_names, all_colors)
  
  ggsave(paste0(element,save_name, "boxPlot.png"),
         device = "png",
         path = paste0(path_save, "/boxPlot/"))
  
}


prepare_obs_pred_df <- function(element, path_preds, path_elem, path_val){
  
  if (element == "intergenic") {
    elem = "^[b]+[0-9]+$"
  } else elem = element
  
  vals = fread(path_val)
  elems = fread(path_elem)
  elems$obsRates <- elems$nMut / (elems$length * elems$N)
  elems <- elems[which(elems$nMut != 0), c('binID', 'obsRates')]
  
  vals <- vals[, c('binID', 'obsRates')]
  obs <- rbind(elems, vals)
  obs <- obs[which(obs$obsRates != 0),]
  
  preds <- do.call(rbind, lapply(path_preds, function(x){
    df <- fread(x)
    colnames(df) = c('binID', 'predRates')
    df$model_name <- unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/"))) - 1]
    df
  }))
  
  dfs = left_join(obs, preds, by = 'binID')
  dfs$sqError = (dfs$obsRates - dfs$predRates)^2
  dfs$absError = abs((dfs$obsRates - dfs$predRates))
  dfs <- dfs[grep(elem, dfs$binID),]
  
  
  dfs
}




# save_barPlot <- function(element, path_preds, path_elem, path_val, path_save, save_name, measurement){
#   
#   dir.create(paste0(path_save, "barPlot/"), showWarnings = F, recursive = T)
#   x = prepare_obs_pred_df(element, path_preds, path_elem, path_val)
#   models <- unique(x$model_name)
#   
#   measures <- c()
#   for (i in 1:length(models)) {
#     measure_df <- data.frame(x[which(x$model_name == models[i]),])
#     msr <- mean(measure_df[,c(measurement)])
#     measures <- c(measures, msr)
#   }
#   
#   dat = data.frame(cbind(models, measures))
#   dat$measures = as.numeric(dat$measures)
#   dat <- dat[order(dat$measures, decreasing = TRUE),]
#   
#   dat$model_name <- factor(dat$models, levels = unique(dat$models))
#   
#   p = ggplot(data=dat, aes(x=measures, y=model_name, fill=model_name))+
#     geom_bar(stat = "identity", width=0.7)+
#     # scale_y_log10()+
#     theme_bw()+
#     labs(title = element, x = "squared error", y = "model")+
#     theme(axis.text.x = element_text(angle = 90,hjust = 1 , size = 4),
#           axis.text.y = element_text(size = 4),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           legend.text = element_text(size = 3))
#   p + create_method_colours(dat$model_name)
#   
#   
#   ggsave(paste0(element,save_name, "barPlot.png"),
#          device = "png",
#          path = paste0(path_save, "/barPlot/"))
#   
# }


save_barPlot <- function(element, path_assessments, path_save, save_name, measurement){
  all_colors <- determine_color_elem(element)
  dir.create(paste0(path_save, "barPlot/"), showWarnings = F, recursive = T)
  x = prepare_assessments(path_assessments, element, measurement)
  models <- unique(x$models)
  paste0(path_save, '/', measurement, '/')
  
  x$measures = as.numeric(x$measures)
  orders = order(x$measures, decreasing = TRUE)
  x <- x[orders,]
  all_colors <- all_colors[orders]
  x$models <- factor(x$models, levels = unique(x$models))
  
  p = ggplot(data=x, aes(x=measures, y=models, fill=models))+
    geom_bar(stat = "identity", width=0.7)+
    # scale_y_log10()+
    theme_bw()+
    labs(title = element, x = measurement, y = "model")+
    theme(axis.text.x = element_text(angle = 90,hjust = 1 , size = 5),
          axis.text.y = element_text(size = 5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 4))
  p = p + create_method_colours(x$models, all_colors)
  
  
  ggsave(paste0(element,save_name,"_", measurement, "_barPlot.png"),
         device = "png",
         width = 10,
         height = 6,
         path = paste0(path_save, "/barPlot/"))
  
}




# save_barPlot2 <- function(element, path_preds, path_elem, path_val, path_save,
#                          save_name, measurement, patterns = NULL) {
#   dir.create(paste0(path_save, "barPlot/"), showWarnings = FALSE, recursive = TRUE)
#   x = prepare_obs_pred_df(element, path_preds, path_elem, path_val)
#   models <- unique(x$model_name)
#   
#   measures <- c()
#   for (i in 1:length(models)) {
#     measure_df <- data.frame(x[which(x$model_name == models[i]),])
#     msr <- mean(measure_df[,c(measurement)])
#     measures <- c(measures, msr)
#   }
#   
#   dat = data.frame(cbind(models, measures))
#   
#   if (length(patterns) != 0) {
#     # Initialize the "group" column with "models"
#     dat$group <- dat$models
#     
#     # Loop through the patterns and group_names
#     for (i in 1:length(patterns)) {
#       pattern <- patterns[i]
#       
#       # Assign the group based on the pattern
#       dat$group <- ifelse(grepl(pattern, dat$models, ignore.case = TRUE), 
#                           paste0(pattern, "_models"), dat$group) 
#     } 
#   }
#   
#   dat$measures = as.numeric(dat$measures)
#   dat <- dat[order(dat$measures, decreasing = TRUE),]
#   
#   dat$model_name <- factor(dat$models, levels = unique(dat$models))
#   
#   
#   p <- ggplot(data = dat, aes(x = group, y = measures, fill = model_name)) +
#     geom_bar(stat = "identity", position = "dodge", width = 0.7) +
#     theme_bw() +
#     labs(title = element, y = "squared error", x = "model") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           legend.text = element_text(size = 6))
#   
#   ggsave(paste0(element, save_name, "groupedBarPlot.png"),
#          plot = p,
#          device = "png",
#          path = paste0(path_save, "/groupedBarPlot/"))
# }


create_method_colours <- function(METHODs, all_colors){
  # all_colors <- c('black', '#01962348', '#800000','#008000', '#000080', '#808000', '#800080', '#008080', '#808080',
  #                 '#2E3E51', '#E04C41', '#236CD6', '#43AC63', '#F2B94A', '#9659B2', '#D15357', '#3DB182',
  #                 "#6E016B", "#E6AB02", "#7FC97F",
  #                 "#11266d", "#8DD3C7",
  #                 "#666666", "#FB8072",
  #                 "#386CB0", "#B2182B", "#2e6f12")
  
  METHODs <- factor(METHODs, levels = unique(METHODs))
  
  method_color_df <- data.frame(METHODs, all_colors[1:length(METHODs)])
  method_colors <- method_color_df[,2]
  
  names(method_colors) <- levels(METHODs)
  colScale <- scale_fill_manual(name = "method", values = method_colors)
  return(colScale)
}


extract_modelsInfo <- function(path_assessments){
  main_method = unlist(lapply(path_assessments, function(x){
    main_model = unlist(strsplit(x, "/"))[4]
    main_model
  }))
  
  LOSS= unlist(lapply(path_assessments, function(x){
    loss = tolower(unlist(strsplit(x, "/"))[5])
    loss
  }))
  
  activationFunction = unlist(lapply(path_assessments, function(x){
    actFunc = ifelse(unlist(strsplit(x, "/"))[6] == 'PReLU', 'PReLU', 'leaky_ReLU')
    actFunc
  }))
  
  model_name <- unlist(lapply(path_assessments, function(x){
    model.name = unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))-1]
    model.name
  }))
  
  model_info <- data.frame(cbind(model_name, main_method, activationFunction, LOSS))
  model_info
}


prepare_assessments <- function(path_assessments, element, measurement){
  model_info = extract_modelsInfo(path_assessments)
  
  x = do.call(rbind, lapply(path_assessments, function(x){
    df <- fread(x)
    df$model_name <- unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/"))) - 1]
    data.frame(df)
  }))
  models <- unique(x$model_name)
  
  measures <- c()
  for (i in 1:length(models)) {
    measure_df <- data.frame(x[which(x$model_name == models[i]),element])
    
    if (measurement == 'acc') {
      msr = measure_df[1,]
    } else if  (measurement == 'corr') {
      msr = measure_df[2,]
    } else if  (measurement == 'mse') {
      msr = measure_df[3,]
    } 
    
    measures <- c(measures, msr)
  }
  
  dat = data.frame(cbind(models, measures))
  # dat$S_C <- ifelse(grepl('siam', dat$models), 'Siam', 'Classic')
  dat = left_join(dat, model_info, by = c('models' = 'model_name'))
  dat
}



save_groupedBarPlots_SC <- function(element, path_assessments, path_save,
                          save_name, measurement, loss, patterns = NULL) {
  dir.create(paste0(path_save, "barPlot/"), showWarnings = FALSE,
             recursive = TRUE)
  dat = prepare_assessments(path_assessments, element, measurement)
  dat = dat[which(dat$LOSS == loss),]
  
  if (length(patterns) != 0) {
    # Initialize the "group" column with "models"
    dat$group <- dat$models
    
    # Loop through the patterns and group_names
    for (i in 1:length(patterns)) {
      pattern <- patterns[i]
      
      # Assign the group based on the pattern
      dat$group <- ifelse(grepl(pattern, dat$models, ignore.case = TRUE), 
                          paste0(pattern, "_models"), dat$group) 
    } 
  }
  
  dat$measures = as.numeric(dat$measures)
  dat$model_name <- factor(dat$models, levels = unique(dat$models))
  
  
  p <- ggplot(data = dat, aes(x = group, y = measures, fill = main_method)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    theme_bw() +
    labs(title = paste0(element, '_',loss), y = measurement, x = "model") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 6))
  p = p + scale_fill_manual(values = c("#11266d", "#E6AB02"))
  
  ggsave(paste0(element, save_name, "_",loss, "_", measurement,
                "_","groupedBarPlot.png"),
         plot = p,
         device = "png",
         path = paste0(path_save, "/groupedBarPlot/"))
}


save_groupedBarPlots_SC_actFunc <- function(element, path_assessments, path_save,
                                    save_name, measurement, loss, patterns = NULL) {
  dir.create(paste0(path_save, "barPlot/"), showWarnings = FALSE,
             recursive = TRUE)
  dat = prepare_assessments(path_assessments, element, measurement)
  dat = dat[which(dat$LOSS == loss),]
  
  if (length(patterns) != 0) {
    # Initialize the "group" column with "models"
    dat$group <- dat$models
    
    # Loop through the patterns and group_names
    for (i in 1:length(patterns)) {
      pattern <- patterns[i]
      
      # Assign the group based on the pattern
      dat$group <- ifelse(grepl(pattern, dat$models, ignore.case = TRUE), 
                          paste0(pattern, "_models"), dat$group) 
    } 
  }
  
  dat$measures = as.numeric(dat$measures)
  dat$model_name <- factor(dat$models, levels = unique(dat$models))
  
  dat$cols = paste0(dat$activationFunction, "_", dat$main_method)
  dat <- dat[order(dat$measures, decreasing = TRUE),]
  
  
  p <- ggplot(data = dat, aes(x = group, y = measures, fill = cols)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    theme_bw() +
    labs(title = paste0(element, '_',loss), y = measurement, x = "model") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 6))
  p = p + scale_fill_manual(values = c('#586919', '#d67406',
                                       #'#728920', '#aaaeba',
                                       '#070e40', '#b89f06'))
  
  ggsave(paste0(element, save_name, measurement, "_", loss,
                "_","groupedBarPlot.png"),
         plot = p,
         device = "png",
         path = paste0(path_save, "/groupedBarPlot/"))
}



create_method_colours <- function(METHODs, all_colors){
  # all_colors <- c('black', '#01234567', '#800000','#008000', '#000080', '#808000', '#800080', '#008080', '#808080',
  #                 '#2E3E51', '#E04C41', '#236CD6', '#43AC63', '#F2B94A', '#9659B2', '#D15357', '#3DB182',
  #                 "#6E016B", "#E6AB02", "#7FC97F",
  #                 "#11266d", "#8DD3C7",
  #                 "#666666", "#FB8072",
  #                 "#386CB0", "#B2182B", "#2e6f12")
  
  METHODs <- factor(METHODs, levels = unique(METHODs))
  
  method_color_df <- data.frame(METHODs, all_colors[1:length(METHODs)])
  method_colors <- method_color_df[,2]
  
  names(method_colors) <- levels(METHODs)
  colScale <- scale_fill_manual(name = "method", values = method_colors)
  return(colScale)
}

determine_color_elem <- function(element){
  elements <- c("intergenic", "gc19_pc.cds", "enhancers", "gc19_pc.3utr", "gc19_pc.5utr",
                "gc19_pc.promCore", "gc19_pc.ss", "lncrna.ncrna","lncrna.promCore")
  colors <- list(c('#513d4a', '#8e6b82', '#cc99bb', '#dbb7cf', '#ead6e3'),
                 c('#2f4458', '#476684', '#77aadd', '#9fc3e7', '#c8ddf1'),
                 c('#2f5151', '#538e8e', '#77cccc', '#ade0e0', '#d6efef'),
                 c('#446655', '#5f8e76', '#88ccaa', '#abdbc3', '#cfeadd'), # green
                 c('#58582f', '#9a9a53', '#dddd77', '#e7e79f', '#f1f1c8'),
                 c('#58442f', '#9a7653', '#ddaa77', '#e7c39f', '#f1ddc8'),
                 c('#684419', '#925f23', '#d18932', '#deac6f', '#eccfad'),
                 c('#582f36', '#844751', '#c66b7a', '#e3929f', '#eebbc3'),
                 c('#3b0811', '#530b17', '#771122', '#ad707a', '#d6b7bc')
                 )
  
  for (i in 1:length(elements)) {
    if (elements[i] == element) {
      color = colors[[i]]
    }
  }
  color
}

