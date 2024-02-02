rm(list = ls())

source("scripts/functions.R")
library(rtracklayer)
library(dplyr)
library(genomation)
library(data.table)
library(stats)
library(ggplot2)
# arguments:
path_to_trainGE <- "../extdata/input/PCAWG/train_elements.tsv.gz"

# make df from train genomic elements that contains each bin starts 
trainGE_gr <- tsv.gz2GRanges(path_to_trainGE)

binstarts <- as.data.frame(trainGE_gr@ranges@start)
bins <- trainGE_gr$name
chr <- as.character(trainGE_gr@seqnames)
df <- as.data.frame(cbind(chr, bins, binstarts))

load(file ="../extdata/input/PCAWG/complete_dataTable_train.RData")
data <- as.data.frame(t(complete_dataTable))

load("../extdata/output/RegLmntDriver/poisson_model_KidneyRCC_26Features.RData")

pdct <- predict(poisson_model, as.data.frame(data), type = "response") # !!!!!!!!! the order of pdct and freq_table_TrainGEs bins are not the same:
# pdct order is based on the output of the generate_complete_dataTable function
fit <- fitted(poisson_model) #should be the same as pdct
bins <- rownames(data)
data <- cbind(data, pdct, fit, bins)

target <- df$bins
data_ordered <- data[match(target, data$bins),]

all_data <- cbind(df, data_ordered[,27:30])
colnames(all_data) <- c("chr", "bins", "binStart", "length", "denominator", "obsMut", "pdctMut")

BMR_table <- all_data %>% group_by(chr)
chr_groups <- group_split(BMR_table)

chrs <- c()
for (i in 1:22) {
  chrs <- c(chrs, unique(chr_groups[[i]]$chr))
}

chr <- which(chrs == "chr22")


chr_1 <- chr_groups[[chr]]
chr_1 <- chr_1[order(chr_1$binStart),]

adj_obsMut <- chr_1$obsMut / chr_1$length
chr_1 <- cbind(chr_1, adj_obsMut)
table_Mut_chr1 <- as.data.frame(sort(table(chr_1$adj_obsMut)))

MyData_thrshld <- function(chr_1, threshold){
  highly_Mut_bins <- as.data.frame(chr_1[chr_1$adj_obsMut>threshold,]) # the mean of adj_obsMut is 0.000187
  data <- highly_Mut_bins[,c(3,8)]
  head(data)
  
  kmeans_clstrng <- kmeans(data$binStart, 5)
  data_clstrs <- as.factor(kmeans_clstrng$cluster)
  myData <- cbind(highly_Mut_bins, data_clstrs)
  
  myData
}

threshold <- 0
head(MyData_thrshld(chr_1, threshold))

# within-class scatter matrix


clusteringScore <-function(x, gr) {
  Sw = sum(aggregate(x, list(gr), function(z) { 
    vZ = ifelse(length(z) == 1, 0, var(z))
    vZ*length(z) /length(x) 
  })[, 2])
  Sm = var(x)
  
  Sm / Sw
}

obsStat_function <- function(chr_1, threshold){
  
  myData <- MyData_thrshld(chr_1, threshold)
  
  x <- myData[,3]
  gr = myData[,9]
  
  obsStat = clusteringScore(x, gr)
  
  obsStat
}

obsStat_function(chr_1, threshold)

### calculating the expected values by random sampling from all binStarts (with large n)

expStats <- function(chr_1, threshold){
  
  myData <- MyData_thrshld(chr_1, threshold)
  
  stats=c()
  for(i in 1:1000) {
    
    idx <- sample(1:nrow(chr_1), nrow(myData))
    
    chr_idx <- chr_1[idx, ]
    
    kmeans_clstrng_all <- kmeans(chr_idx$binStart, 5)
    
    data_clstrs_all <- as.factor(kmeans_clstrng_all$cluster)
    
    myDataExp <- cbind(chr_idx, data_clstrs_all)
    
    gr <- myDataExp$data_clstrs_all
    x <- myDataExp$binStart
    
    stats = c(stats, clusteringScore(x, gr) )
  }
  
  stats
}

# head(expStats(chr_1, threshold))



# plotting:
hist(chr_1$adj_obsMut)
hist(chr_1$adj_obsMut, breaks = 80)
abline(v=mean(chr_1$adj_obsMut), col = 3)

# sensitivity Analysis:
thresholds <- c(0, .00017, .0003, .0005, .001)


pValue <- c()

for(i in 1:length(thresholds)){
  threshold <- thresholds[i]
  obsStat <- obsStat_function(chr_1, threshold)
  stats <- expStats(chr_1, threshold)
  p <- mean(stats>= obsStat)
  
  pValue <- cbind(pValue, p)
}

as.data.frame(rbind(thresholds, pValue))

obsStat <- obsStat_function(chr_1, threshold)
stats <- expStats(chr_1, threshold)
mean(stats>= obsStat)


############### correlation analysis for features and #obsMuts ############
###########################################################################
cor_Features_obsMut <- sapply(data_ordered[,c(1:26)],
                              function(x) cor(as.numeric(x), data_ordered$Response))

ordered_corrs <- sort(abs(cor_Features_obsMut))

cor.test(data_ordered$Response, data_ordered$wgEncodeUwRepliSeqHepg2WaveSignalRep1)

## plotting the correlations:
pairs(data_ordered[,c(1, 29)], pch = 19)

idx_bins_chr1 <- chr_1$bins

pairs(data_ordered[idx_bins_chr1,c(8, 29)], pch = 19)
