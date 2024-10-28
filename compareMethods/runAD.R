rm(list = ls())

library(ActiveDriverWGS)
library(data.table)
library(dplyr)
library(GenomicRanges)

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





path_donorInfo <- '../../iDriver/extdata/procInput/iDriverInputs/donorInfo.tsv'
path_to_gr='../../iDriver/extdata/procInput/mut_PCAWG/all_samples.RData'
path_to_elements <- "../external/database/bins/proccessed/PCAWG_callable.bed12"


cohorts <- c('Pancan-no-skin-melanoma-lymph', "Liver-HCC", 
             "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
             "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
             "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
             "Breast-AdenoCa","Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",         
             "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
             "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",    
             "Lymph-NOS", "CNS-Medullo", "Myeloid-AML", "CNS-Oligo", "Cervix-SCC",
             "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
             "Cervix-AdenoCA", "Breast-DCIS", "Bone-Cart", "Myeloid-MDS" )

cohort = "Liver-HCC"
exclude_lymph_melanoma = ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)



donorInfo <- select_cohort(path_donorInfo, 
                           cohort,
                           exclude_lymph_melanoma,
                           exclude_hyper_mutated = T)
load(path_to_gr)

gr_cohort <- gr[which(gr$D_id %in% donorInfo$donor_id)]


# create the mutData long format for the cohort of interest and scores of interest
df <- data.frame('chr' = as.character(seqnames(gr_cohort)),
                 'pos1'= start(gr_cohort),
                 'pos2' = end(gr_cohort),
                 'ref' = gr_cohort$ref_al,
                 'alt' = gr_cohort$alt_al, 
                 'patient' = gr_cohort$D_id)
head(df)

df_unique <- df[!duplicated(df),]

PCAWG_elems = prepare_elements_from_BED12(path_to_elements)


df_unique <- df_unique[which(grepl("[ATGC\\-]", c(df_unique$ref))),]
df_unique <- df_unique[which(grepl("[ATGC\\-]", c(df_unique$alt))),]

dir.create('../../ActiveDriverWGSR/output/recovery/', showWarnings = F, recursive = T)
results = ActiveDriverWGS(mutations = df_unique,
                          elements = PCAWG_elems,
                          recovery.dir = '../../ActiveDriverWGSR/output/recovery/',
                          mc.cores = 50)

dir.create(paste0('../../ActiveDriverWGSR/output/', cohort, '/'), showWarnings = F, recursive = T)
fwrite(results, file = paste0("../../ActiveDriverWGS/output/", cohort, "/results.tsv"), sep = '\t')