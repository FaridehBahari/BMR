rm(list = ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)

path_bed <- '../external/database/bins/proccessed/PCAWG_callable.bed6'
win_size <- 100000


gr_elements <- import.bed(path_bed)
gr_elements$blockID <- paste0('b', 1:length(gr_elements))

# gr_extended <- resize(gr_elements, width = width(gr_elements) + 100000, fix = "center")
# gr_split <- split(gr_extended, gr_extended$name)
gr_split <- resize(gr_elements, width = width(gr_elements) + win_size, fix = "center")
gr_split <- split(gr_split, gr_split$name)


library(doParallel)
library(foreach)


num_cores <- 40  
registerDoParallel(num_cores)

t0 = Sys.time()
# Step 3: Parallelize the loop using foreach
filtered_bg <- foreach(i = 1:length(gr_split), .combine = c, .packages = "GenomicRanges") %dopar% {
  if (i %% 5000 == 0) {
    print(i) # i = 121790
  }
  
  gr_elem <- gr_elements[which(gr_elements$name == names(gr_split)[i])]
  bg_containElem <- gr_split[[i]]
  gr_elem$length <- width(gr_elem)
  
  filtered_bg_grs <- GRanges()
  # Remove overlaps within each element
  for (j in 1:length(bg_containElem)) {
    filtered_bg_gr <- GenomicRanges::setdiff(bg_containElem[j], gr_elem)
    filtered_bg_gr$name <- unique(gr_elem$name)
    ID = bg_containElem[j]$blockID
    filtered_bg_gr$blockID <- ID
    filtered_bg_gr$length = gr_elem[which(gr_elem$blockID == ID)]$length
    filtered_bg_grs <- c(filtered_bg_grs, filtered_bg_gr)
  }
  
  return(filtered_bg_grs)
}


print(filtered_bg)

print(Sys.time() - t0)

# filtered_bg <- GRanges()
# for (i in 1:length(gr_split)) {
#   if (i%%5000 == 0) {
#     print(i)
#   }
#   gr_elem <- gr_elements[which(gr_elements$name == names(gr_split)[i])]
#   
#   # Remove overlaps within each element
#   filtered_bg_gr = GenomicRanges::setdiff(gr_split[[i]], gr_elem)
#   filtered_bg_gr$name <- unique(gr_elem$name)
#   filtered_bg <- c(filtered_bg, filtered_bg_gr)
# }



select_cohort <- function(path_donorInfo, 
                          cohort, exclude_lymph_melanoma = TRUE,
                          exclude_hyper_mutated = TRUE){
  
  donorInfo <- fread(path_donorInfo)
  
  if (exclude_lymph_melanoma) {
    exceptions <- c("Skin-Melanoma", "SKCM-US",
                    "Lymph-NOS", 
                    "Lymph-CLL", "CLLE-ES",
                    "Lymph-BNHL", "MALY-DE", "DLBC-US")
    donorInfo <- donorInfo[-which(donorInfo$cohort1 %in% exceptions),] # 2279 donors
  }
  
  if (exclude_hyper_mutated) {
    donorInfo <- donorInfo[which(!donorInfo$HyperMut_donor),] 
  }
  
  if (!cohort %in% c('Pan_Cancer', 'Pancan-no-skin-melanoma-lymph')) {
    donorInfo <- donorInfo[which(donorInfo$cohort1 == cohort),]
  } 
  
  donorInfo <- donorInfo[,c("D_id","freq" )]
  colnames(donorInfo) <- c('donor_id', 'totalMutNrs')
  donorInfo
  
  
  
}


map_muts <- function(gr, testGE_gr){
  
  
  names(testGE_gr) = testGE_gr$name
  
  GenomicElement <- unlist(lapply(strsplit(names(testGE_gr), "[::]"), function(x){x[1]}))
  
  lenElement <- width(testGE_gr)
  
  mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
  
  callable_GE <- split(testGE_gr, names(testGE_gr))
  lenElement_new <- sum(width(callable_GE))
  
  ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
  
  idx_gr <- queryHits(ov)
  idx_callable_GE <- subjectHits(ov)
  
  gr_mappedGE <- gr[idx_gr]
  
  GE_col=unique(mcols(callable_GE, level = "within")[,"GenomicElement"])[idx_callable_GE]
  EL_col = sum(mcols(callable_GE, level = "within")[,"lenElement"])[idx_callable_GE]
  name_col = names(callable_GE)[idx_callable_GE]
  
  mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE),
                                  GenomicElement=GE_col,
                                  elemenntLength = EL_col,
                                  name = name_col)
  
  idx_completely_black <- which(lenElement_new == 0)
  
  
  
  if(length(idx_completely_black) != 0){
    lenElement_new <- lenElement_new[-idx_completely_black]
    callable_GE <- callable_GE[-idx_completely_black]
  }
  
  
  list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE,
       lenAllTests = lenElement_new)
  
}



create_response <- function(all_elements, n_donors){
  
  obs_Muts <- all_elements$MAF_GE_mapped
  df <- as.data.frame(DataFrame(obs_Muts))
  df <- df[,c('name', 'D_id')]
  
  # Convert df to a data.table (if it's not already one)
  setDT(df)
  
  # Create the result table using data.table syntax
  resTab <- df[, .N, by = .(name, D_id)][
    N != 0, .(nMut = sum(N), nSample = .N), by = .(binID = name)
  ]
  
  idx_lenElement <- resTab$binID
  all_lmnt_length <- all_elements$lenAllTests
  ElementLength <- all_lmnt_length[idx_lenElement]
  resTab$length <- ElementLength
  
  resTab$N <- n_donors
  
  all_binIds <- names(all_elements[["lenAllTests"]])
  nonMut_ids <- all_binIds[which(!all_binIds %in% resTab$binID)]
  res_nonMutated <- data.frame(cbind('binID' = nonMut_ids, 
                                     'nMut' = rep(0, length(nonMut_ids)),
                                     'nSample' = rep(0, length(nonMut_ids)),
                                     'length' = all_lmnt_length[!names(all_lmnt_length) %in% resTab$binID], 
                                     'N' = rep(n_donors, length(nonMut_ids))))
  
  resTab <- rbind(resTab, res_nonMutated)
  
  return(resTab)
}


############################### create block-wise background response tables ####################
dir.create('../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/', recursive=T)
load('../../iDriver/extdata/procInput/mut_PCAWG/all_samples.RData')
gr <- gr[!duplicated(gr)]

path_donorInfo <- '../../iDriver/extdata/procInput/iDriverInputs/donorInfo.tsv'


cohorts <- c('Pancan-no-skin-melanoma-lymph', "Liver-HCC" ,
  "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
  "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
  "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
  "Breast-AdenoCa","Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",
  "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
  "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",
  "CNS-Medullo", "CNS-Oligo", "Cervix-SCC",
  "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Osteosarc",
  "Breast-LobularCa",  "Lymph-NOS", "Myeloid-AML", "Bone-Epith",
  "Cervix-AdenoCA","Breast-DCIS","Bone-Cart" #,"Myeloid-MDS"
) 


for (cohort in cohorts) {
  
  print(cohort)
  
  exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
  
  
  
  
  donorInfo <- select_cohort(path_donorInfo,
                             cohort,
                             exclude_lymph_melanoma,
                             exclude_hyper_mutated = T)
  
  gr_cohort <- gr[which(gr$D_id %in% donorInfo$donor_id)]
  
  
  
  bg_bed <- filtered_bg
  bg_bed$PCAWG_ID <- bg_bed$name
  bg_bed$name <- bg_bed$blockID
  
  t0 <- Sys.time()
  all_elements <- map_muts(gr_cohort, bg_bed)
  print(paste0('map_muts total time: ', Sys.time() - t0))
  
  n_donors <- length(unique(gr_cohort$D_id))
  
  t0 <- Sys.time()
  resTab <- create_response(all_elements, n_donors)
  print(paste0('create_response total time: ', Sys.time() - t0))
  
  resTab
  
  df <- data.frame('PCAWG_ID' = bg_bed$PCAWG_ID,
                   'elem_ID' = bg_bed$name, 
                   'PCAWG_block_length' = bg_bed$length)
  
  df <- df[!duplicated(df),]
  
  pred_tab <- left_join(resTab, df, by = c('binID'= 'elem_ID'))
  pred_tab <- pred_tab[!duplicated(pred_tab),]
  # pred_tab <- as.data.frame(left_join(resTab, df, by = c('binID'= 'elem_ID')))
  # pred_tab <- pred_tab[,which(!colnames(pred_tab) %in% c('binID', 'PCAWG_block_length'))]
  # pred_tab <- pred_tab[!duplicated(pred_tab),]
  
  
  fwrite(pred_tab, file = paste0('../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/',cohort, '.tsv'), sep = '\t')
}

################################# estimate mutation rate based on weighted means of block lengths ###############################################
rm(list = ls())
library(data.table)
library(dplyr)

# setwd('C:/Active/projects/make_features/BMR/')

dir.create('../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/assessments/', showWarnings = F, recursive = T)
dir.create('../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/predRates/', showWarnings = F, recursive = T)
ann <- fread('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', sep = ',')


ann <- ann[(ann$in_CGC == TRUE) | (ann$in_CGC_literature == TRUE) | 
             (ann$in_CGC_new == TRUE) | (ann$in_oncoKB == TRUE) | 
             (ann$in_pcawg == TRUE), ]
drivers <- ann$PCAWG_IDs
cohorts <- c('Pancan-no-skin-melanoma-lymph', "Liver-HCC" ,
             "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
             "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
             "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
             "Breast-AdenoCa","Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",
             "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
             "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",
             "CNS-Medullo", "CNS-Oligo", "Cervix-SCC",
             "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Osteosarc",
             "Breast-LobularCa",  "Lymph-NOS", "Myeloid-AML", "Bone-Epith",
             "Cervix-AdenoCA","Breast-DCIS","Bone-Cart" #,"Myeloid-MDS"
) 

assessments <- c()
for (cohort in cohorts) {
  
  print(cohort)
  
  test_y <- fread(paste0('../../iDriver/extdata/procInput/BMRs/observed/', cohort,
                         '/test_y.tsv'))
  test_y$obsRate <- test_y$nMut/test_y$length
  
  final_resTab <- fread(paste0('../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/',cohort, '.tsv'))
  
  obs_tab <- test_y[,c('binID', 'obsRate')]
  df <- left_join(final_resTab, obs_tab, by = c('PCAWG_ID' = 'binID'))
  df$nMut = as.numeric( df$nMut)
  df$length = as.numeric( df$length)
  df$PCAWG_block_length = as.numeric(df$PCAWG_block_length)
  
  
  # Efficient calculation using dplyr
  results <- df %>%
    group_by(PCAWG_ID) %>%
    summarise(
      total_length = sum(PCAWG_block_length),
      rates = nMut / length,
      pred_elem = sum((PCAWG_block_length / sum(PCAWG_block_length)) * rates) / n(),
      obs_elem = unique(obsRate),
      .groups = "drop") %>%
    ungroup()
  
  results2 <- results[,c('PCAWG_ID', 'total_length', 'pred_elem', 'obs_elem')]
  results2 <- results2[!duplicated(results2),]
  results2 <- results2[which(results2$obs_elem != 0),]
  
  define_element_type <- function(binID_vector){
    
    s <- strsplit(binID_vector, "[::]")
    GenomicElement <- unlist(lapply(s, function(x){x[1]}))
    GenomicElement
    
  }
  
  results2$elemType <- define_element_type(results2$PCAWG_ID)
  fwrite(results2, paste0('../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/predRates/',cohort, '.tsv'), sep = '\t')
  
  elems <- unique(results2$elemType)
  
  # remove known drivers from assessments
  results2 <- results2[!(results2$PCAWG_ID %in% drivers), ]
  results2 <- as.data.frame(results2)
  corrs <- c()
  for (elemType in elems) {
    elem <- results2[which(results2$elemType == elemType),]
    corrs <- c(corrs, cor(elem$obs_elem, elem$pred_elem, method = 'spearman'))
  }
  
  ass <- data.frame('baseline_corr' = corrs, 'element-type' = elems)
  ass$cohort <- cohort
  assessments <- rbind(assessments, ass)
  
}


fwrite(assessments, file = '../external/BMR/output/Res_reviewerComments/baseLine_localRates/final/assessments/all_cohorts_assessments.tsv',
       , sep = '\t')


