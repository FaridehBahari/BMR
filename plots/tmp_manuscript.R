rm(list = ls())

if (.Platform$OS.type == "windows") {
  setwd('A:/myThesis/BMR_proj/RankBased_BMR/')
}
source('plots/functions.R')
path_val <- "../external/rawInput/validation_sets_10folds/Pan_Cancer_validate_y_fold_1.tsv"
path_save <- "../external/plots_manuscript/"
path_elem <- "../external/rawInput/Pan_Cancer_test_y.tsv"
path_preds <- c("../external/output/GBM/GBM_predTest.tsv",
                "../external/output/RF/RF_predTest.tsv",
                "../external/output/GLM/GLM_predTest.tsv",
                "../external/output/classic_NN/poisLoss/PReLU/classicPois_sevenL_do3_bs512_PRelUactivation/classicPois_sevenL_do3_bs512_PRelUactivation_predTest.tsv")
############## boxPlots ###############
element <- 'intergenic'
save_name = 'withoutObsRates'
# save_boxPlot(element, path_preds, path_elem, path_val, path_save, save_name)
save_boxPlot(element, path_preds, path_elem, path_val, path_save, save_name, include_obsRates = FALSE)

######################## ordered barplot ###################################
path_assessments <- c("../external/output/GBM/GBM_assessments.tsv",
                "../external/output/RF/RF_assessments.tsv",
                "../external/output/classic_NN/poisLoss/PReLU/classicPois_sevenL_do3_bs512_PRelUactivation_bn/classicPois_sevenL_do3_bs512_PRelUactivation_bn_assessments.tsv",
                "../external/output/classic_NN/MSEloss/PReLU/classicMSE_sevenL_do3_bs512_PReLU/classicMSE_sevenL_do3_bs512_PReLU_assessments.tsv",
                "../external/output/GLM/GLM_assessments.tsv"
                )

elements <- c("intergenic", "gc19_pc.cds", "enhancers", "gc19_pc.3utr", "gc19_pc.5utr",
              "gc19_pc.promCore", "gc19_pc.ss", "lncrna.ncrna","lncrna.promCore")

save_name <- 'test_bar'
measurements = c('mse', 'corr', 'acc')
for (element in elements) {
  for (measurement in measurements){
    save_barPlot(element, path_assessments, path_save, save_name, measurement)
  }
}

# colors:
c('#2f5151', '#538e8e', '#77cccc', '#ade0e0', '#d6efef') #blue
c('#2f4458', '#476684', '#77aadd', '#9fc3e7', '#c8ddf1') #cian
c('#513d4a', '#8e6b82', '#cc99bb', '#dbb7cf', '#ead6e3') #purple
c('#446655', '#5f8e76', '#88ccaa', '#abdbc3', '#cfeadd') # green
c('#58582f', '#9a9a53', '#dddd77', '#e7e79f', '#f1f1c8') # olive
c('#58442f', '#9a7653', '#ddaa77', '#e7c39f', '#f1ddc8') # peach-brown
c('#582f36', '#844751', '#c66b7a', '#e3929f', '#eebbc3') # cherry


c('#3b0811', '#530b17', '#771122', '#ad707a', '#d6b7bc') # rouge
c('#684419', '#925f23', '#d18932', '#deac6f', '#eccfad') # orange