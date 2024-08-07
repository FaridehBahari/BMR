rm(list = ls())

library(stringr)
library(RColorBrewer)
library(rtracklayer)

setwd('A:/myThesis/activeProj/iDriver/')
setwd('C:/Active/projects/iDriver/iDriver/')

map_muts <- function(gr, path_to_GEs, element,callable_GE_train, 
                     count_blocks=TRUE, test=NULL, train=NULL){
  
  if (element == "test") {
    
    load(path_to_GEs)
    
    
    testGE_gr <- unlist(testGE)
    
    g <- names(testGE_gr)
    s <- strsplit(g, "[::]")
    GenomicElement <- unlist(lapply(s, function(x){x[1]}))
    
    lenElement <- width(testGE_gr)
    
    mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
    
    callable_GE <- relist(testGE_gr, testGE)
    lenElement_new <- sum(width(callable_GE))
    
  } else if (element == "train") {
    
    #callable_GE <- make_callable_trainGE(path_to_GEs, path_to_callable)
    callable_GE <- callable_GE_train
    lenElement_new <- sum(width(callable_GE))
  }
  
  
  
  ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
  
  # if we run this function for each element, length(unique(queryHits(ov))) is nMut and
  # length(unique(subjectHits(ov))) is number of mutated blocks. >>> ov <- findOverlaps(gr, testGE$`gc19_pc.cds::gencode::ARID1A::ENSG00000117713.13`, ignore.strand = TRUE)
  
  idx_gr <- queryHits(ov)
  idx_callable_GE <- subjectHits(ov)
  
  gr_mappedGE <- gr[idx_gr]
  
  GE_Mutated <- callable_GE[idx_callable_GE]  
  
  if (element == "test") {
    mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE), 
                                    GenomicElement=unique(mcols(GE_Mutated, level = "within")[,"GenomicElement"]), 
                                    elemenntLength = sum(mcols(GE_Mutated, level = "within")[,"lenElement"]),
                                    name = names(GE_Mutated))
  } else if (element == "train") {
    mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE), 
                                    name = names(GE_Mutated))
  }
  
  idx_completely_black <- which(lenElement_new == 0)
  
  
  
  if(length(idx_completely_black) != 0){
    lenElement_new <- lenElement_new[-idx_completely_black]
    callable_GE <- callable_GE[-idx_completely_black]
  }
  
  if (count_blocks) {
    
    GE_blocks <- unlist(callable_GE)
    
    ov_block <- findOverlaps(gr, GE_blocks, ignore.strand = TRUE)
    idx_gr_block <- queryHits(ov_block)
    idx_GE_block <- subjectHits(ov_block)
    
    gr_blocks_mappedGE <- gr[idx_gr_block]
    GE_blocks_mutated <- GE_blocks[idx_GE_block]
    
    mcols(gr_blocks_mappedGE) <- DataFrame(mcols(gr_blocks_mappedGE), 
                                           binID=names(GE_blocks_mutated),
                                           block_ID=paste0(names(GE_blocks_mutated), "**",1:length(gr_blocks_mappedGE)),
                                           block_length=width(GE_blocks_mutated))
    
    list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE, 
         lenAllTests = lenElement_new, gr_blocks = gr_blocks_mappedGE)
  } else {
    list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE, 
         lenAllTests = lenElement_new)
  }
  
}

generate_mut_df <- function(gr, testGEs, sigHits){
  
  all_elements <- map_muts(gr, path_to_GEs, path_to_callable, element = 'test',
                           callable_GE_train = '', count_blocks = T)
  gr_annotated <- all_elements$MAF_GE_mapped
  
  gr_annotated <- gr_annotated[which(gr_annotated$name %in% sigHits)]
  
  mutation <- as.data.frame(cbind(gr_annotated$name,
                                  as.character(gr_annotated$var_type),
                                  as.numeric(start(gr_annotated))))
  colnames(mutation) <- c("binID", "var_type", "X_position")
  
  mutation
}

########### Step_2: Visualize the mutations in the element of interest as a lolipop plot #######
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

lollipop_full_length <- function(gene, mutation, testGE){
  
  
  GeneMutation <- mutation[mutation$binID == gene, ]
  
  coordinats <- subset(as.data.frame(table(as.character(GeneMutation$var_type),
                                           GeneMutation$X_position)), Freq != 0)
  
  colnames(coordinats) <- c("Variant_type", "X", "freq")
  coordinats$Variant_type <- gsub("SNP", "SNV", coordinats$Variant_type)
  
  y_ylim <- if(max(coordinats$freq) > 5) {
    y_ylim = max(coordinats$freq)
  }else{
    y_ylim = 5}
  
  element_coordinates <- testGE[[gene]]
  L <- sum(width(element_coordinates))
  
  S <- start(element_coordinates) # - (start(element_coordinates)[1])    #exon start
  E <- end(element_coordinates) #- start(element_coordinates)[1]        #exon End
  c <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"),
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), 
         brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), 
         brewer.pal(12, "Set3"), brewer.pal(9, "Set1")) #pallet
  
  gene_info = define_element_type(gene)
  
  # Set the font size for axis labels, axis numbers, and main title
  par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2)
  
  # Define the plot
  plot(x = 1,                 
       xlab = "Position of Mutations", 
       ylab = "#Mutations",
       xlim = c(min(as.numeric(as.character(coordinats$X))), 
                max(as.numeric(as.character(coordinats$X)))), 
       ylim = c(-1, y_ylim),
       frame = FALSE,
       type = "n",
       main = paste0( gene_info[[2]][1], " ", gene_info[[1]][1])
       # main = paste0(gene_info[[1]][1], " of ", gene_info[[2]][1])
  )
  
  # Add a gray rectangle at the bottom of the plot
  rect(min(as.numeric(as.character(coordinats$X))) - 100, 
       -.3, 
       max(as.numeric(as.character(coordinats$X))),
       0, density = NA, angle = 45, col = "grey", border = NA)
  
  # Add blue rectangles for element coordinates
  for(i in 1:length(element_coordinates)){
    rect(S[i], -.4, E[i], 0.2, density = NA, angle = 45, col = '#386cb0', border = NA)
  }
  
  # Add segments for mutation positions
  for(i in 1:nrow(coordinats)){
    segments(as.numeric(as.character(coordinats$X[i])),
             0, as.numeric(as.character(coordinats$X[i])),
             coordinats$freq[i], lwd = .2)
  }
  
  # Add points for mutations
  for (i in 1:nrow(coordinats)) {
    if (coordinats$Variant_type[i] == 'SNV') {
      point_type = 16
      col = "#4dac26" 
    } else if (coordinats$Variant_type[i] == 'INS' | coordinats$Variant_type[i] == 'DEL') {
      point_type = 17
      col = 'red' 
    } else if (coordinats$Variant_type[i] == 'DNP' | coordinats$Variant_type[i] == 'TNP' | coordinats$Variant_type[i] == 'ONP') {
      point_type = 8
      col = '#984ea3'
    }
    
    points(x = as.numeric(as.character(coordinats$X[i])), 
           y = coordinats$freq[i], 
           pch = point_type, 
           col = col)
  }
  
  nMut <- nrow(GeneMutation)
  
  # Define legend labels, point types, and colors
  legend_labels <- c("SNV", "INS/DEL") #, "DNV/TNV/ONV")
  legend_pch <- c(16, 17, 8)
  legend_colors <- c("#4dac26", "red") #, "#984ea3")
  
  # Add the legend
  legend(x = "topright",
         legend = legend_labels,
         pch = legend_pch,
         col = legend_colors,
         title = paste0("Element length = ", L),
         ncol = 1, horiz = FALSE, cex = 1.2, pt.cex = 1.3)
  
}

################### Step3: save lolipop plots for a group of genes ##############
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
define_element_type <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[1]}))
  GenomicElements <- unlist(lapply(strsplit(GenomicElement, "[.]"), function(x){x[length(x)]}))
  gene_name <- unlist(lapply(s, function(x){x[5]}))
  
  GEs <- c()
  for (GenomicElement in GenomicElements) {
    if (GenomicElement == 'enhancers') {
      GenomicElement = 'enhancer'
    } else if (GenomicElement == '3utr') {
      GenomicElement = "3'UTR"
    } else if (GenomicElement == '5utr') {
      GenomicElement = "3'UTR"
    } else if (GenomicElement == 'cds') {
      GenomicElement = 'CDS'
    } else if (GenomicElement == 'promCore') {
      GenomicElement = 'Core Promoter'
    } else if (GenomicElement == 'ss') {
      GenomicElement = 'Splice site'
    } 
    
    GEs <- c(GEs, GenomicElement)
  }
  
  
  list(GEs, gene_name)
  
}


save_lolipops <- function(mutation, testGE,
                          path_out = '../extdata/lolipops/',
                          save_name = ''){
  selected_genes <- unique(mutation$binID)
  dir.create(path_out, showWarnings = F, recursive = T)
  
  info <- define_element_type(selected_genes)
  
  elemType <- info[[1]]
  gene_name <- info[[2]]
  
  for (i in 1:length(selected_genes)) {
    print(i)
    file_name <- paste0(gene_name[i], '_',elemType[i])
    png( paste0(path_out, file_name, save_name, ".png"))
    
    lollipop_full_length(selected_genes[i], mutation, testGE)
    dev.off()
  }
  
}
################################################################################

path_gr <- '../extdata/procInput/mut_PCAWG/pan_cancer_FILTEREDsamples.RData'
path_to_GEs = '../extdata/procInput/callable_testGE.RData' 

load(path_gr)
load(path_to_GEs)
sigHits <- c('gc19_pc.3utr::gencode::ADH1B::ENSG00000196616.8', 
             'gc19_pc.3utr::gencode::BHMT::ENSG00000145692.10',
             'gc19_pc.3utr::gencode::BRINP2::ENSG00000198797.6',
             'gc19_pc.3utr::gencode::GLYR1::ENSG00000140632.12',
             'gc19_pc.3utr::gencode::IFI16::ENSG00000163565.14',
             'gc19_pc.3utr::gencode::PDPR::ENSG00000090857.9',
             'gc19_pc.3utr::gencode::TCP10::ENSG00000203690.7',
             'gc19_pc.promCore::gencode::COL4A2::ENSG00000134871.13',
             'gc19_pc.promCore::gencode::OR5T3::ENSG00000172489.5',
             'gc19_pc.promCore::gencode::SRPRB::ENSG00000144867.7',
             'gc19_pc.promCore::gencode::TPTE::ENSG00000166157.12',
             'gc19_pc.promCore::gencode::ZSCAN5B::ENSG00000197213.5')

mutation <- generate_mut_df(gr, testGE, sigHits)

save_lolipops(mutation, testGE)


################################################################################
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

path_to_GEs = '../extdata/procInput/callable_testGE.RData' 
load(path_to_GEs)


ex <-  testGE$`gc19_pc.3utr::gencode::IFI16::ENSG00000163565.14`
BSgenomeViews(Hsapiens, ex)

#GRanges object with 1 range and 1 metadata column:
# seqnames              ranges strand |       hit
# <Rle>           <IRanges>  <Rle> | <logical>
#   [1]     chr1 159024692-159024945      + |      TRUE
# -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

as.character(DNAStringSet(BSgenomeViews(Hsapiens, ex)))
# "AATCTGGATGTCATTGACGATAATGTTTATGGAGATAAGGTCTAAGTGCCTAAAAAAATGTACATATACCTGGTTGAAATACAACACTATACATACACACCACCATATATACTAGCTGTTAATCCTATGGAATGGGGTATTGGGAGTGCTTTTTTAATTTTTCATAGTTTTTTTTTAATAAAATGGCATATTTTGCATCTACAACTTCTATAATTTGAAAAAATAAATAAACATTATCTTTTTTGTGAAAGGAA"
# 159024859
159024692  - 159024859 +1



#################################################################################
tpte_wt <-  testGE$`gc19_pc.promCore::gencode::TPTE::ENSG00000166157.12`
BSgenomeViews(Hsapiens, tpte_wt)

mutation[grepl('TPTE', mutation$binID), ]

as.character(reverseComplement(DNAStringSet(BSgenomeViews(Hsapiens, tpte_wt))))
# [1] "CGCCTCCAAGACCCCCCCCCAACAAAAAGGAGCGTCCCCCACCCCTACCCCCGCCCGGAGGACTTAGGGCCTGG"                                                                                                                                                                                                                                                                                            
# [2] "CTCGGGCGCGGAGCTAAGTGTAGGCGCCGGGGGTCCCTAGAGCCGCCGGGGCGCAGCGAGTCCGGCGCTGGGTAACTGTTGGGTCAGAAACTGTTCAGGTAGCAGCTGTTGTGCCCTCCCTTGGCCCCGCCGCTCGGAGACGCCCCGCCCCCTGCCTTGAACGGCCGCCCGGCCCCGCCCCAGCGCCCACGTGACTAGCATAGGCGCGCCCCCGTTCCGCCCGCCGCCGCAGACTCCGCCTCCGGGACGCGAGCGAGCGGCGAGCGCGCGCACTACCAGTTCTTGCTCGGCGACTCCCGCGCACGCGCGCGCCGTGCCACCCTCCCCGCACCCCTCCTCCCGCCATCCGGCTTAACGT"

10990771 - 10990763 +1

#################################################################################

gr_del <- GenomicRanges::GRanges(seqnames = 'chr1',
                                 IRanges::IRanges(start = 159024692,
                                                  width = 1),strand = '+')
BSgenomeViews(Hsapiens, gr_del)


gr_del <- GenomicRanges::GRanges(seqnames = 'chr21',
                                 IRanges::IRanges(start = 10990771,
                                                  width = 1),strand = '+')
BSgenomeViews(Hsapiens, gr_del)
