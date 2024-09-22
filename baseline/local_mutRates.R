######################## strategy number 1 ###################################
# just extend each block to 100kb window
rm(list = ls())
# setwd('C:/Active/projects/make_features/BMR/')
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)


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

create_downstream_background = function(gr_elements, win_size, this_genome) {
  
  
  # background is plus/minus window around every component of element
  # make sure flanking coordinates stay within chromosomal boundaries
  # add one extra nucleotide of slack on both ends - trinuc tabulation needs that
  
  bg_starts = GenomicRanges::start(gr_elements) - win_size
  bg_ends = GenomicRanges::start(gr_elements) - 1
  
  # max pos is end of chromosome minus one
  max_chr_pos = GenomeInfoDb::seqlengths(this_genome)[as.character(GenomicRanges::seqnames(gr_elements))]
  max_chr_pos = max_chr_pos - 1
  
  # min pos is 2nd pos of chromosome
  bg_starts[bg_starts < 2] = 2
  bg_ends[bg_ends > max_chr_pos] = max_chr_pos[bg_ends > max_chr_pos]
  gr_background = GenomicRanges::GRanges(GenomicRanges::seqnames(gr_elements), 
                                         IRanges::IRanges(bg_starts, bg_ends))
  
  mcols(gr_background) <- gr_elements$name
  names(mcols(gr_background)) <- gsub("^X$", "name", names(mcols(gr_background)))
  gr_background
}


create_upstream_background = function(gr_elements, win_size, this_genome) {
  
  
  # background is plus/minus window around every component of element
  # make sure flanking coordinates stay within chromosomal boundaries
  # add one extra nucleotide of slack on both ends - trinuc tabulation needs that
  
  bg_starts = GenomicRanges::end(gr_elements) + 1 
  bg_ends = GenomicRanges::end(gr_elements) + win_size 
  
  # max pos is end of chromosome minus one
  max_chr_pos = GenomeInfoDb::seqlengths(this_genome)[as.character(GenomicRanges::seqnames(gr_elements))]
  max_chr_pos = max_chr_pos - 1
  
  # min pos is 2nd pos of chromosome
  bg_starts[bg_starts < 2] = 2
  bg_ends[bg_ends > max_chr_pos] = max_chr_pos[bg_ends > max_chr_pos]
  gr_background = GenomicRanges::GRanges(GenomicRanges::seqnames(gr_elements), 
                                         IRanges::IRanges(bg_starts, bg_ends))
  
  mcols(gr_background) <- gr_elements$name
  names(mcols(gr_background)) <- gsub("^X$", "name", names(mcols(gr_background)))
  gr_background
  
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

#################### prepare input ###########
path_bed <- '../external/database/bins/proccessed/PCAWG_callable.bed6'
win_size <- 50000
this_genome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens

gr_elements <- import.bed(path_bed)
gr_elements$PCAWG_ID <- gr_elements$name
gr_elements$name <- paste0('ele', 1:length(gr_elements) )



#### create background bed ####
gr_background_ds = create_downstream_background(gr_elements, win_size, this_genome)
gr_background_us = create_upstream_background(gr_elements, win_size, this_genome)

bg_bed <- c(gr_background_ds, gr_background_us)
######### annotate muts ##########
dir.create('../external/BMR/output/Res_reviewerComments/baseLine_localRates/input/', recursive=T)
load('../../iDriver/extdata/procInput/mut_PCAWG/all_samples.RData')

path_donorInfo <- '../../iDriver/extdata/procInput/iDriverInputs/donorInfo.tsv'

# cohort <- 'Pancan-no-skin-melanoma-lymph'

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
  all_elements <- map_muts(gr_cohort, bg_bed)
  
  n_donors <- length(unique(gr_cohort$D_id))
  resTab <- create_response(all_elements, n_donors)
  
  ele_PCAWGs = data.frame('background' = gr_elements$name,
                          'PCAWG_ID' = gr_elements$PCAWG_ID, 
                          'length_PCAWG_ID' = width(gr_elements))
  
  final_resTab = left_join(resTab, ele_PCAWGs, by = c('binID' = 'background'))
  colnames(final_resTab)[1] = 'background_id'
  fwrite(final_resTab, file = paste0('../external/BMR/output/Res_reviewerComments/baseLine_localRates/input/',
                                     cohort, '_responseTable_all100k.tsv'), sep = '\t')
  
  ##### assesments bg1 ####
  test_y <- fread(paste0('../../iDriver/extdata/procInput/BMRs/observed/', cohort,
                         '/test_y.tsv'))
  test_y$obsRate <- test_y$nMut/test_y$length
  
  # obs <- test_y[which(test_y$binID == 'gc19_pc.3utr::gencode::AL627309.1::ENSG00000237683.5'),]
  # x = df[which(df$PCAWG_ID == 'gc19_pc.3utr::gencode::AL627309.1::ENSG00000237683.5'),]
  # 
  # total_length = sum(x$length_PCAWG_ID)
  # weights <- x$length_PCAWG_ID/total_length
  # 
  # rates <- x$nMut/x$length
  # pred_elem <- sum(weights*rates)/length(rates)
  # obs_elem = obs$nMut/obs$length
  
  
  obs_tab <- test_y[,c('binID', 'obsRate')]
  df <- left_join(final_resTab, obs_tab, by = c('PCAWG_ID' = 'binID'))
  df$nMut = as.numeric( df$nMut)
  df$length = as.numeric( df$length)
  df$length_PCAWG_ID = as.numeric(df$length_PCAWG_ID)
  
  # Efficient calculation using dplyr
  results <- df %>%
    group_by(PCAWG_ID) %>%
    summarise(
      total_length = sum(length_PCAWG_ID),
      rates = nMut / length,
      pred_elem = sum((length_PCAWG_ID / sum(length_PCAWG_ID)) * rates) / n(),
      obs_elem = unique(obsRate)
    ) %>%
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
  
  elems <- unique(results2$elemType)
  
  corrs <- c()
  for (elemType in elems) {
    elem <- results2[which(results2$elemType == elemType),]
    corrs <- c(corrs, cor(elem$obs_elem, elem$pred_elem, method = 'spearman'))
  }
  
  ass <- data.frame('baseline_corr' = corrs, 'element-type' = elems)
  fwrite(ass, file = paste0('../external/BMR/output/Res_reviewerComments/baseLine_localRates/',
                            cohort, '_baselineAssessment.tsv'), sep = '\t')
  
}

