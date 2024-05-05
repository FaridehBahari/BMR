

path_bench_res <- '../../make_features/external/BMR/benchmark/Pancan-no-skin-melanoma-lymph/tables/table_GoldStd_basedon_allfdr.csv'

grouped_barPlot_model_elements <- function(path_ass_full, ass_type, compare,
                                           path_save = '../external/BMR/plots/',
                                           save_name = ''){
  dir.create(path_save, showWarnings = F, recursive = T)
  df <- fread(path_bench_res)
  df <- df[which(df$method %in% c("GBM", "RF", "nn_mseLoss", "nn_poisLoss", "TL" )),]
  
  
  # Get column indices where "CDS" is present in column names
  cds_cols <- grep("CDS", names(df))
  
  # Split dataframe based on column indices
  cds_df <- df[, c(1, cds_cols)]
  nc_df <- df[, -cds_cols]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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

