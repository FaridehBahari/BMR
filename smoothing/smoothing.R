rm(list = ls())
# setwd('A:/myThesis/make_features/BMR/')
setwd('C:/Active/projects/make_features/BMR/')
# library(rtracklayer)
library(dplyr)
# library(genomation)
library(data.table)
library(stats)
# library(ggplot2)

path_to_trainGE <- "../external/database/bins/proccessed/callable_intergenic_intervals_wo_pcawg.bed6"
path_intergenic_response <- '../external/BMR/rawInput/responseTabs/Pan_Cancer/var_bins.tsv'
path_pred_train <- '../external/BMR/output/bin_size_effect/var_size_nonMutInclude_longer50/GBM/GBM_predTrain.tsv'

path_testGE <- '../external/database/bins/raw/PCAWG_test_genomic_elements.bed6'
path_pcawg_response <- '../external/BMR/rawInput/responseTabs/Pan_Cancer/reg_elems.tsv'
path_predTest <- '../external/BMR/output/bin_size_effect/var_size_nonMutInclude_longer50/GBM/GBM_predTest.tsv'


# load intergenics 
intergenic <- fread(path_to_trainGE)
intergenic <- intergenic[,1:4]
colnames(intergenic) <- c('chr', 'start', 'end', 'binID')
Y_intergenic <- fread(path_intergenic_response)
Intergenics <- left_join(Y_intergenic, intergenic, by = 'binID')
pred_intergenic <- fread(path_pred_train)
Intergenics <- left_join(Intergenics, pred_intergenic, by = 'binID')


regLms <- fread(path_testGE)
regLms <- regLms[,1:4]
colnames(regLms) <- c('chr', 'start', 'end', 'binID')
Y_regLms <- fread(path_pcawg_response)
regLms <- left_join(Y_regLms, regLms, by = 'binID')
pred_regLms <- fread(path_predTest)
regLms <- left_join(regLms, pred_regLms, by = 'binID')

all_data <- rbind(Intergenics, regLms)
all_data$obsRate <- all_data$nMut/(all_data$length * all_data$N)

BMR_table <- all_data %>% group_by(chr)
chr_groups <- group_split(BMR_table)

chrs <- c()
for (i in 1:length(chr_groups)) {
  chrs <- c(chrs, unique(chr_groups[[i]]$chr))
}

chr <- which(chrs == "chr22")


given <- chr_groups[[chr]]
given <- given[order(given$start),]
plot(given$start, given$predRate)
plot(given$start, given$obsRate)

############ check for dependency of bins #############
#######################################################
##1) acf
acf(given$nMut, lag.max = max(given$end))
acf(given$obsRate, lag.max = max(given$end))

acf(na.omit(given$predRate), lag.max = max(given$end))

d <- given$predRate - given$obsRate
acf(na.omit(d), lag.max = max(given$end))

##2)By using clustering scores >>>>> checkBinsDep.R script in supplementary


################ Smoothing >>>>>>>>>>> Nadaraya and Watson’s #############
##########################################################
given <- given[which(!is.na(given$predRate)),]
x <- given$start
y <- given$predRate


# Kernel Regression (Nadaraya and Watson’s Kernel Smoother):
smoothed <- ksmooth(x, y, kernel = c("box", "normal"), bandwidth = 1000000,
          range.x = range(x))


given$kernelReg <- smoothed$y


idx <- which(!is.na(given$kernelReg))
new_given <- given[idx,]
cor(new_given$obsRate, new_given$kernelReg, method = 'spearman') # bandW = 50k: 0.115137
cor(new_given$obsRate, new_given$predRate, method = 'spearman') #  bandW = 50k: 0.386814
############################################################


# define span (the function is from http://users.stat.umn.edu/~helwig/notes/smooth-notes.html)

loess.gcv <- function(x, y){
  nobs <- length(y)
  xs <- sort(x, index.return = TRUE)
  x <- xs$x
  y <- y[xs$ix]
  tune.loess <- function(s){
    lo <- loess(y ~ x, span = s)
    mean((lo$fitted - y)^2) / (1 - lo$trace.hat/nobs)^2
  }
  os <- optimize(tune.loess, interval = c(.01, 99))$minimum
  lo <- loess(y ~ x, span = os)
  list(x = x, y = lo$fitted, df = lo$trace.hat, span = os)
}


gcv <- loess.gcv(x,y) 

fit_LOESS <- loess(given$predRate ~ given$start, span = 10/nrow(given)) # for chr5 gcv$span = 0.01042671
plot(given$start, given$predRate)
lines(given$start, fit_LOESS$fitted, col = 5)
quantile(fit_LOESS$fitted)

idx = 1:1000
plot(given$binStart[idx], given$pdctMut[idx])
lines(given$binStart[idx], fit_LOESS$fitted[idx], col = 2, lwd=2)


smoothed_model <- predict(fit_LOESS) # the same as fit_LOESS$fitted

############ for test elements:###############
##############################################
path_to_testGEs <- "../extdata/input/PCAWG/PCAWG_test_genomic_elements.bed12.gz"
path_to_MAF <- "../extdata/input/PCAWG/Kidney_RCC_MAF.RData"

load("../extdata/input/PCAWG/test_features_restricted.RData") # a subset of the second HDF5 file (~1.5 Gb)
test_features <- restricted_features

tests <- MutAnnotation_by_RegElmnts(path_to_MAF, path_to_testGEs)
tests$all_GEs_grl
gr_tests_blocks <- unlist(tests$all_GEs_grl)

GE_name <- names(gr_tests_blocks)
GE_start <- start(gr_tests_blocks)
chr <- seqnames(gr_tests_blocks)
df_tests <- as.data.frame(cbind(chr, GE_name, GE_start))
