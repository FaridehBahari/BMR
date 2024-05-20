######### ActivePathways ###############
# https://cran.r-project.org/web/packages/ActivePathways/vignettes/ActivePathways-vignette.html
rm(list = ls())
# install.packages("ActivePathways")
library(ActivePathways)

scores <- read.table(
  system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package = 'ActivePathways'), 
  header = TRUE, sep = '\t', row.names = 'Gene')
scores <- as.matrix(scores)
head(scores)

scores[is.na(scores)] <- 1


gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
ActivePathways(scores, gmt.file)

nrow(ActivePathways(scores, gmt.file, significant = 0.05))

nrow(ActivePathways(scores, gmt.file, significant = 0.1))

#### GMT objects >>> The GMT is structured as a list of terms (e.g., molecular pathways, biological processes, etc.). In the GMT object, each term is a list containing an id, a name, and the list of genes associated with this term.
gmt <- read.GMT(gmt.file)
names(gmt[[1]])

# Look at the genes annotated to the first term
gmt[[1]]$genes

# Get the full name of Reactome pathway 2424491
gmt$`REAC:2424491`$name

# processing step for GMT files >>>>  the removal of gene sets that are too large or small
gmt <- Filter(function(term) length(term$genes) >= 10, gmt)
gmt <- Filter(function(term) length(term$genes) <= 500, gmt)

ActivePathways(scores, gmt)

ActivePathways(scores, gmt.file, geneset.filter = c(10, 500))

# To be continued...
# Backgrou



##############################################################################
library(data.table)


path_inference <- '../external/BMR/output/TL/GBM/inference/GBM_inference_binomTest.tsv'
df <- fread(path_inference)
df <- df[!grep('lncrna', df$binID),]
df <- df[!grep('enhancers', df$binID),]


get_info <- function(binID){
  gene_name <- paste0(unlist(strsplit(binID, '::'))[3], '::', unlist(strsplit(binID, '::'))[4])
  element <- unlist(strsplit(binID, '::'))[1]
  
  return(list(gene_name, element))
}
df <- data.frame(df)
x = apply(df, 1, function(s){
  binID = s[1]
  
  info = get_info(binID)
  cbind(info[[1]], info[[2]])
})

x= t(x)
p_value <- df$raw_p_binom

x <- cbind(x, p_value)
head(x)

# First, ensure it's a data.table
x = data.frame(x)
colnames(x) <- c('gene', 'element_type', 'p_value')
library(tidyr)

# Assuming your data frame is named 'x'
wide_df <- spread(x, element_type, p_value)
head(wide_df)

rownames(wide_df) <- wide_df$gene

wide_df <- wide_df[,c(2:ncol(wide_df))]


wide_df[is.na(wide_df)] <- 1

# Assuming your data frame is named 'wide_df'
wide_df_numeric <- wide_df

# Convert all values to numeric
wide_df_numeric[] <- lapply(wide_df_numeric, as.numeric)

# Convert to matrix
scores <- as.matrix(wide_df_numeric)

# Print the resulting matrix
head(scores)

