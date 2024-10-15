# rm(list = ls())
# 
# library(data.table)
# 
# ################## cohorts #############################
# included_cohorts <- c("Pancan-no-skin-melanoma-lymph", "Liver-HCC", "ColoRect-AdenoCA" ,
#                       "Uterus-AdenoCA" , "Kidney-RCC", "Lung-SCC", "Biliary-AdenoCA",
#                       "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
#                       "Breast-AdenoCa", "Bladder-TCC", "Eso-AdenoCa", 
#                       "Lymph-BNHL", "Lymph-CLL",
#                       "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
#                       "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc")
# 
# ############# AD #####################
# files <- paste0('../../ActiveDriverWGSR/output/', included_cohorts, '/AD_result.tsv')
# 
# ass <- lapply(files, fread)
# 
# # Add a cohort column to each assessment
# responses <- lapply(paste0('../../iDriver/extdata/procInput/BMRs/observed/', included_cohorts, '/test_y.tsv'), fread)
# responses <- lapply(responses, function(s){
#   s <- s[which(s$nMut !=0),]
# })
# 
# # Add a cohort column to each assessment
# ass <- Map(function(dt, cohort) {
#   dt[, cohort := cohort]
#   return(dt)
# }, ass, included_cohorts)
# 
# 
# 
# ass <- Map(function(ass_dt, resp_dt) {
#   # Filter rows in ass_dt where ELT is in binID of resp_dt
#   ass_dt <- ass_dt[ELT %in% resp_dt$binID]
#   return(ass_dt)
# }, ass, responses)
# 
# 
# drivers <- lapply(ass, function(s){
#   print(nrow(s))
#   s$fdr = p.adjust(s$pp_element, method = 'fdr')
#   s = s[which(s$fdr < 0.05), c('id', 'pp_element', 'fdr','cohort')]
# })
# 
# drivers <- do.call(rbind, drivers)
# colnames(drivers) <- c('PCAWG_ID', 'p_value', 'fdr','cohort')
# drivers$tools <- 'ActiveDriverWGS'
# dir.create('../external/BMR/output/Res_reviewerComments/driverDiscovery/', recursive = T, showWarnings = F)
# fwrite(drivers, file = '../external/BMR/output/Res_reviewerComments/driverDiscovery/AD_allCohorts.tsv', sep = '\t')
# 
# ######################### DP ######################################
# 
# files <- paste0('../../DriverPower/output/', included_cohorts, '/non0_', included_cohorts, '.result.tsv')
# 
# ass <- lapply(files, fread)
# 
# # Add a cohort column to each assessment
# ass <- Map(function(dt, cohort) {
#   dt[, cohort := cohort]
#   return(dt)
# }, ass, included_cohorts)
# 
# drivers <- lapply(ass, function(s){
#   s = s[which(s$nMut !=0), ]
#   print(nrow(s))
#   s$fdr = p.adjust(s$raw_p, method = 'fdr')
#   s = s[which(s$fdr < 0.05), c('binID', 'raw_p', 'fdr','cohort')]
# })
# 
# drivers <- do.call(rbind, drivers)
# colnames(drivers) <- c('PCAWG_ID', 'p_value', 'fdr','cohort')
# drivers$tools <- 'DriverPower'
# dir.create('../external/BMR/output/Res_reviewerComments/driverDiscovery/', recursive = T, showWarnings = F)
# fwrite(drivers, file = '../external/BMR/output/Res_reviewerComments/driverDiscovery/DP_allCohorts.tsv', sep = '\t')
# 
# ######################### Dig ######################################
# 
# cohort_name <- ifelse(included_cohorts == 'Pancan-no-skin-melanoma-lymph', 'Pancan_SNV_MNV_INDEL', included_cohorts)
# files <- paste0('../../Dig/output/elemDriver/', cohort_name, '.results.txt')
# ass <- lapply(files, fread)
# responses <- lapply(paste0('../../iDriver/extdata/procInput/BMRs/observed/', included_cohorts, '/test_y.tsv'), fread)
# responses <- lapply(responses, function(s){
#   s <- s[which(s$nMut !=0),]
# })
# 
# # Add a cohort column to each assessment
# ass <- Map(function(dt, cohort) {
#   dt[, cohort := cohort]
#   return(dt)
# }, ass, included_cohorts)
# 
# 
# 
# ass <- Map(function(ass_dt, resp_dt) {
#   # Filter rows in ass_dt where ELT is in binID of resp_dt
#   ass_dt <- ass_dt[ELT %in% resp_dt$binID]
#   return(ass_dt)
# }, ass, responses)
# 
# 
# 
# drivers <- lapply(ass, function(s){
#   print(nrow(s))
#   s$fdr = p.adjust(s$PVAL_MUT_BURDEN, method = 'fdr')
#   s = s[which(s$fdr < 0.05), c('ELT', "PVAL_MUT_BURDEN", 'fdr', "cohort")]
# })
# 
# drivers <- do.call(rbind, drivers)
# colnames(drivers) <- c('PCAWG_ID', 'p_value', 'fdr','cohort')
# drivers$tools <- 'Dig'
# fwrite(drivers, file = '../external/BMR/output/Res_reviewerComments/driverDiscovery/Dig_allCohorts.tsv', sep = '\t')
# 
# ######################### eMET ######################################
# 
# files <- paste0('../external/BMR/output/reviewerComments/', included_cohorts, '/eMET/inference/eMET_inference.tsv')
# 
# ass <- lapply(files, fread)
# 
# # Add a cohort column to each assessment
# ass <- Map(function(dt, cohort) {
#   dt[, cohort := cohort]
#   return(dt)
# }, ass, included_cohorts)
# 
# drivers <- lapply(ass, function(s){
#   s = s[which(s$nMut !=0), ]
#   print(nrow(s))
#   s$fdr = p.adjust(s$p_value, method = 'fdr')
#   s = s[which(s$fdr < 0.05), c('binID', 'p_value', 'fdr','cohort')]
# })
# 
# drivers <- do.call(rbind, drivers)
# colnames(drivers) <- c('PCAWG_ID', 'p_value', 'fdr','cohort')
# drivers$tools <- 'eMET'
# dir.create('../external/BMR/output/Res_reviewerComments/driverDiscovery/', recursive = T, showWarnings = F)
# fwrite(drivers, file = '../external/BMR/output/Res_reviewerComments/driverDiscovery/eMET_allCohorts.tsv', sep = '\t')


rm(list = ls())

library(data.table)

################## cohorts #############################
included_cohorts <- c("Pancan-no-skin-melanoma-lymph", "Liver-HCC", "ColoRect-AdenoCA" ,
                      "Uterus-AdenoCA" , "Kidney-RCC", "Lung-SCC", "Biliary-AdenoCA",
                      "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
                      "Breast-AdenoCa", "Bladder-TCC", "Eso-AdenoCa", 
                      "Lymph-BNHL", "Lymph-CLL",
                      "CNS-GBM", "Panc-AdenoCA" , "Lung-AdenoCA" ,"Prost-AdenoCA",
                      "Ovary-AdenoCA" , "Bone-Leiomyo", "CNS-Medullo","Bone-Osteosarc")

# ############# AD #####################
files <- paste0('../../ActiveDriverWGSR/output/', included_cohorts, '/AD_result.tsv')

ass <- lapply(files, fread)

# Add a cohort column to each assessment
responses <- lapply(paste0('../../iDriver/extdata/procInput/BMRs/observed/', included_cohorts, '/test_y.tsv'), fread)
responses <- lapply(responses, function(s){
  s <- s[which(s$nMut !=0),]
})

# Add a cohort column to each assessment
ass <- Map(function(dt, cohort) {
  dt[, cohort := cohort]
  return(dt)
}, ass, included_cohorts)



ass <- Map(function(ass_dt, resp_dt) {
  # Filter rows in ass_dt where ELT is in binID of resp_dt
  ass_dt <- ass_dt[id %in% resp_dt$binID]
  return(ass_dt)
}, ass, responses)


drivers <- lapply(ass, function(s){
  print(nrow(s))
  s$fdr = p.adjust(s$pp_element, method = 'fdr')
  s = s[, c('id', 'pp_element', 'fdr','cohort')]
  colnames(s) <- c('binID', "p_value", 'fdr', "cohort")
  fwrite(s, paste0('../../ActiveDriverWGSR/output/', unique(s$cohort), '/final_AD_result.tsv'), sep = '\t')
})



# ######################### Dig ######################################
define_element_type <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[1]}))
  GenomicElement
  
}

cohort_name <- ifelse(included_cohorts == 'Pancan-no-skin-melanoma-lymph', 'Pancan_SNV_MNV_INDEL', included_cohorts)
files <- paste0('../../Dig/output/elemDriver/', cohort_name, '.results.txt')
ass <- lapply(files, fread)
responses <- lapply(paste0('../../iDriver/extdata/procInput/BMRs/observed/', included_cohorts, '/test_y.tsv'), fread)
responses <- lapply(responses, function(s){
  s <- s[which(s$nMut !=0),]
})

# Add a cohort column to each assessment
ass <- Map(function(dt, cohort) {
  dt[, cohort := cohort]
  return(dt)
}, ass, included_cohorts)



ass <- Map(function(ass_dt, resp_dt) {
  # Filter rows in ass_dt where ELT is in binID of resp_dt
  ass_dt <- ass_dt[ELT %in% resp_dt$binID]
  return(ass_dt)
}, ass, responses)



drivers <- lapply(ass, function(s){
  print(nrow(s))
  s$fdr = p.adjust(s$PVAL_MUT_BURDEN, method = 'fdr')
  s = s[, c('ELT', "PVAL_MUT_BURDEN", 'fdr', "cohort")]
  colnames(s) <- c('binID', "p_value", 'fdr', "cohort")
  fwrite(s, paste0('../../Dig/output/elemDriver/', unique(s$cohort), '_final_Dig_result.tsv'), sep = '\t')
})

######################### DP ######################################

files <- paste0('../../DriverPower/output/', included_cohorts, '/non0_', included_cohorts, '.result.tsv')

ass <- lapply(files, fread)

# Add a cohort column to each assessment
ass <- Map(function(dt, cohort) {
  dt[, cohort := cohort]
  return(dt)
}, ass, included_cohorts)

drivers <- lapply(ass, function(s){
  s = s[which(s$nMut !=0), ]
  s = s[!(grepl('lncrna', s$binID)),]
  print(nrow(s))
  s$fdr = p.adjust(s$raw_p, method = 'fdr')
  s = s[, c('binID', 'raw_p', 'fdr','cohort')]
  colnames(s) <- c('binID', 'p_value', 'fdr','cohort')
  
  fwrite(s, paste0('../../DriverPower/output/', unique(s$cohort), '/final_DP_result.tsv'), sep = '\t')
  s = s[which(s$fdr < 0.05), c('binID', 'p_value', 'fdr','cohort')]
})

drivers <- do.call(rbind, drivers)
colnames(drivers) <- c('PCAWG_ID', 'p_value', 'fdr','cohort')
drivers$tools <- 'DriverPower'
dir.create('../external/BMR/output/Res_reviewerComments/driverDiscovery/', recursive = T, showWarnings = F)
fwrite(drivers, file = '../external/BMR/output/Res_reviewerComments/driverDiscovery/DP_allCohorts.tsv', sep = '\t')

######################### eMET ######################################

files <- paste0('../external/BMR/output/reviewerComments/', included_cohorts, '/eMET/inference/eMET_inference.tsv')

ass <- lapply(files, fread)

# Add a cohort column to each assessment
ass <- Map(function(dt, cohort) {
  dt[, cohort := cohort]
  return(dt)
}, ass, included_cohorts)

drivers <- lapply(ass, function(s){
  s = s[which(s$nMut !=0), ]
  print(nrow(s))
  s$fdr = p.adjust(s$p_value, method = 'fdr')
  s = s[which(s$fdr < 0.05), c('binID', 'p_value', 'fdr','cohort')]
})

drivers <- do.call(rbind, drivers)
colnames(drivers) <- c('PCAWG_ID', 'p_value', 'fdr','cohort')
drivers$tools <- 'eMET'
dir.create('../external/BMR/output/Res_reviewerComments/driverDiscovery/', recursive = T, showWarnings = F)
fwrite(drivers, file = '../external/BMR/output/Res_reviewerComments/driverDiscovery/eMET_allCohorts.tsv', sep = '\t')
