######################## strategy number 1 ###################################
# just extend each block to 100kb window
rm(list = ls())
setwd('C:/Active/projects/make_features/BMR/')
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
###############################

path_bed <- '../external/database/bins/proccessed/PCAWG_callable.bed6'
win_size <- 50000
this_genome <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens

gr_elements <- import.bed(path_bed)
gr_elements$PCAWG_ID <- gr_elements$name
gr_elements$name <- paste0('ele', 1:length(gr_elements) )


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



#### bg number1 ####
gr_background_ds = create_downstream_background(gr_elements, win_size, this_genome)
gr_background_us = create_upstream_background(gr_elements, win_size, this_genome)

bg_bed <- c(gr_background_ds, gr_background_us)


######### annotate muts ##########

load('../../iDriver/extdata/procInput/mut_PCAWG/pan_cancer_FILTEREDsamples.RData')

all_elements <- map_muts(gr, bg_bed)

n_donors <- length(unique(gr$D_id))
resTab <- create_response(all_elements, n_donors)

ele_PCAWGs = data.frame('background' = gr_elements$name,
                        'PCAWG_ID' = gr_elements$PCAWG_ID, 
                        'length_PCAWG_ID' = width(gr_elements))

final_resTab = left_join(resTab, ele_PCAWGs, by = c('binID' = 'background'))
colnames(final_resTab)[1] = 'background_id'
fwrite(final_resTab, file = '../external/BMR/output/Res_reviewerComments/baseLine_localRates/input/responseTable_all100k.tsv', sep = '\t')

##### assesments bg1 ####
final_resTab <- fread('../external/BMR/output/Res_reviewerComments/baseLine_localRates/input/responseTable_all100k.tsv')

test_y <- fread('../../iDriver/extdata/procInput/BMRs/observed/Pancan-no-skin-melanoma-lymph/test_y.tsv')
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

ass <- data.frame('bg1' = corrs, 'element-type' = elems)



get_identified_drivers <- function(path_ann_pcawg_IDs, based_on) {
  
  # Load the CSV file into a DataFrame
  df <- fread(path_ann_pcawg_IDs)
  
  if (based_on == "all") {
    # Filter rows where at least one of the specified columns is TRUE
    filtered_df <- df[which(df$in_CGC | df$in_CGC_literature | df$in_CGC_new | df$in_oncoKB | df$in_pcawg), ]
  } else if (based_on == "in_pcawg") {
    # Filter rows where in_pcawg is TRUE
    filtered_df <- df[df$in_pcawg, ]
  }
  
  # Select the 'PCAWG_IDs' column from the filtered DataFrame
  drivers <- filtered_df$PCAWG_IDs
  
  return(drivers)
}


drivers <- get_identified_drivers('../external/BMR/procInput/ann_PCAWG_ID_complement.csv', 'all')

results2 <- results2[which(!results2$PCAWG_ID %in% drivers),]

####################### edited Strategy number 1 (bg1.1) ######################
# need to be corrected
rm(list = ls())
setwd('C:/Active/projects/make_features/BMR/')
library(GenomicRanges)
library(rtracklayer)


path_bed <- '../external/database/bins/proccessed/PCAWG_callable.bed6'
win_size <- 100000


gr_elements <- import.bed(path_bed)

# gr_extended <- resize(gr_elements, width = width(gr_elements) + 100000, fix = "center")
# gr_split <- split(gr_extended, gr_extended$name)
gr_split <- resize(gr_elements, width = width(gr_elements) + 100000, fix = "center")
gr_split <- split(gr_split, gr_split$name)

library(doParallel)
library(foreach)


num_cores <- 20  # Use one less than the total number of cores
registerDoParallel(num_cores)

t0 = Sys.time()
# Step 3: Parallelize the loop using foreach
filtered_bg <- foreach(i = 1:length(gr_split), .combine = c, .packages = "GenomicRanges") %dopar% {
  if (i %% 5000 == 0) {
    print(i)
  }
  
  gr_elem <- gr_elements[which(gr_elements$name == names(gr_split)[i])]
  
  # Remove overlaps within each element
  filtered_bg_gr <- GenomicRanges::setdiff(gr_split[[i]], gr_elem)
  filtered_bg_gr$name <- unique(gr_elem$name)
  
  return(filtered_bg_gr)
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

load('../../iDriver/extdata/procInput/mut_PCAWG/pan_cancer_FILTEREDsamples.RData')
bg_bed <- filtered_bg
bg_bed$PCAWG_ID <- bg_bed$name
bg_bed$name <- paste0('ele', 1:length(bg_bed))

t0 <- Sys.time()
all_elements <- map_muts(gr, bg_bed)
print(paste0('map_muts total time: ', Sys.time() - t0))

n_donors <- length(unique(gr$D_id))

t0 <- Sys.time()
resTab <- create_response(all_elements, n_donors)
print(paste0('create_response total time: ', Sys.time() - t0))

resTab

df <- data.frame('PCAWG_ID' = bg_bed$PCAWG_ID,
                 'elem_ID' = bg_bed$name)
df2 = data.frame('PCAWG_ID' = gr_elements$name, 'PCAWG_block_length' = width(gr_elements))
df <- left_join(df, df2, by = 'PCAWG_ID')

pred_tab <- left_join(resTab, df, by = c('binID'= 'elem_ID'))

fwrite(pred_tab, file = '../external/BMR/output/Res_reviewerComments/baseLine_localRates/input/background_responseTab_mainElement_filtered_s1.1.tsv', sep = '\t')

########################## Strategy number 2 #############
# # extend each block and remove all pcawg elements from the extenxions:
# 
# cut -f1,2 ../../../../Projects/bahari_work/database/fasta/hg19.fa.fai > hg19.genome
# 
# bedtools slop -i ../external/database/bins/proccessed/PCAWG_callable.bed6 -g ../../../../Projects/bahari_work/database/fasta/hg19.fa -b 50000 > ../external/BMR/output/Res_reviewerComments/baseLine_localRates/callable_pcawg.bed
# 
# 
# bedtools subtract -a ../external/BMR/output/Res_reviewerComments/baseLine_localRates/extended_callable_pcawg.bed -b ../
#   external/database/bins/raw/PCAWG_test_genomic_elements.bed6 > ../external/BMR/output/Res_reviewerComments/baseLine_localRates/extended_pcawgFiltered.bed

rm(list = ls())
library(GenomicRanges)
library(rtracklayer)

path_bed <- '../external/BMR/output/Res_reviewerComments/baseLine_localRates/input/extended_pcawgFiltered.bed'

bg_bed2 <-  import.bed(path_bed)
bg_bed2$PCAWG_ID <- bg_bed2$name
bg_bed2$name <- paste0('ele', 1:length(bg_bed2))

load('../../iDriver/extdata/procInput/mut_PCAWG/pan_cancer_FILTEREDsamples.RData')
all_elements2 <- map_muts(gr, bg_bed2)

resTab2 <- create_response(all_elements2, n_donors)
