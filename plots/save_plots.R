rm(list = ls())

######################## define variables ###################################
if (.Platform$OS.type == "windows") {
  setwd('A:/myThesis/BMR_proj/RankBased_BMR/')
}
source('plots/functions.R')
path_preds <- c("../external/output/GBM/GBM_predTest.tsv",
                "../external/output/RF/RF_predTest.tsv",
                "../external/output/GLM/GLM_predTest.tsv",
                "../external/output/siam/mseLoss/siam_500_50fiveL_do3_bs512_n2_cd8_mse/siam_500_50fiveL_do3_bs512_n2_cd8_mse_predTest.tsv",
                "../external/output/siam/mseLoss/siam_threeL_do3_bs512_n2_cd8_mse/siam_threeL_do3_bs512_n2_cd8_mse_predTest.tsv",
                "../external/output/siam/mseLoss/siam_sevenL_do3_bs512_n2_cd8_mse/siam_sevenL_do3_bs512_n2_cd8_mse_predTest.tsv",
                "../external/output/siam/poisLoss/siam_500_50threeL_do3_bs512_n2_cd8_pois/siam_500_50threeL_do3_bs512_n2_cd8_pois_predTest.tsv",
                "../external/output/siam/poisLoss/siam_500_50fiveL_do3_bs512_n2_cd8_pois/siam_500_50fiveL_do3_bs512_n2_cd8_pois_predTest.tsv",
                "../external/output/siam/poisLoss/siam_500_50sevenL_do3_bs512_n2_cd8_pois/siam_500_50sevenL_do3_bs512_n2_cd8_pois_predTest.tsv",
                "../external/output/classic_NN/MSEloss/lossMSE_500_5threeL_do3_bs512/lossMSE_500_5threeL_do3_bs512_predTest.tsv",
                "../external/output/classic_NN/MSEloss/lossMSE_500_5fiveL_do3_bs512/lossMSE_500_5fiveL_do3_bs512_predTest.tsv",
                "../external/output/classic_NN/MSEloss/lossMSE_500_5sevenL_do3_bs512/lossMSE_500_5sevenL_do3_bs512_predTest.tsv",
                "../external/output/classic_NN/poisLoss/lossPois_500_5threeL_do3_bs512/lossPois_500_5threeL_do3_bs512_predTest.tsv",
                "../external/output/classic_NN/poisLoss/lossPois_500_5fiveL_do3_bs512/lossPois_500_5fiveL_do3_bs512_predTest.tsv",
                "../external/output/classic_NN/poisLoss/lossPois_500_5sevenL_do3_bs512/lossPois_500_5sevenL_do3_bs512_predTest.tsv",
                "../external/output/siam/mseLoss/PReLU/siam_threeL_do3_bs512_n2_cd8_mse_PReLU/siam_threeL_do3_bs512_n2_cd8_mse_PReLU_predTest.tsv",
                "../external/output/siam/mseLoss/PReLU/siam_fiveL_do3_bs512_n2_cd8_mse_PReLU/siam_fiveL_do3_bs512_n2_cd8_mse_PReLU_predTest.tsv",
                "../external/output/siam/mseLoss/PReLU/siam_sevenL_do3_bs512_n2_cd8_mse_PReLU/siam_sevenL_do3_bs512_n2_cd8_mse_PReLU_predTest.tsv",
                
                
                "../external/output/siam/poisLoss/PReLU/siam_threeL_do3_bs512_n2_cd8_poisLoss_PReLU/siam_threeL_do3_bs512_n2_cd8_poisLoss_PReLU_predTest.tsv",
                "../external/output/siam/poisLoss/PReLU/siam_fiveL_do3_bs512_n2_cd8_poisLoss_PReLU/siam_fiveL_do3_bs512_n2_cd8_poisLoss_PReLU_predTest.tsv",
                "../external/output/siam/poisLoss/PReLU/siam_sevenL_do3_bs512_n2_cd8_poisLoss_PReLU/siam_sevenL_do3_bs512_n2_cd8_poisLoss_PReLU_predTest.tsv",
                
                "../external/output/classic_NN/MSEloss/PReLU/classicMSE_threeL_do3_bs512_PReLU/classicMSE_threeL_do3_bs512_PReLU_predTest.tsv",
                "../external/output/classic_NN/MSEloss/PReLU/classicMSE_fiveL_do3_bs512_PReLU/classicMSE_fiveL_do3_bs512_PReLU_predTest.tsv",
                "../external/output/classic_NN/MSEloss/PReLU/classicMSE_sevenL_do3_bs512_PReLU/classicMSE_sevenL_do3_bs512_PReLU_predTest.tsv",
                
                "../external/output/classic_NN/poisLoss/PReLU/classicPois_yhreeL_do3_bs512_PRelUactivation/classicPois_yhreeL_do3_bs512_PRelUactivation_predTest.tsv",
                "../external/output/classic_NN/poisLoss/PReLU/classicPois_fiveL_do3_bs512_PRelUactivation/classicPois_fiveL_do3_bs512_PRelUactivation_predTest.tsv",
                "../external/output/classic_NN/poisLoss/PReLU/classicPois_sevenL_do3_bs512_PRelUactivation_bn/classicPois_sevenL_do3_bs512_PRelUactivation_bn_predTest.tsv"
                
                
)

path_assessments <- c(
  # "../external/output/GBM/GBM_assessments.tsv",
  #               "../external/output/RF/RF_assessments.tsv",
  #               "../external/output/GLM/GLM_assessments.tsv",
  "../external/output/siam/mseLoss/siam_500_50fiveL_do3_bs512_n2_cd8_mse/siam_500_50fiveL_do3_bs512_n2_cd8_mse_assessments.tsv",
  "../external/output/siam/mseLoss/siam_threeL_do3_bs512_n2_cd8_mse/siam_threeL_do3_bs512_n2_cd8_mse_assessments.tsv",
  "../external/output/siam/mseLoss/siam_sevenL_do3_bs512_n2_cd8_mse/siam_sevenL_do3_bs512_n2_cd8_mse_assessments.tsv",
  "../external/output/siam/poisLoss/siam_500_50threeL_do3_bs512_n2_cd8_pois/siam_500_50threeL_do3_bs512_n2_cd8_pois_assessments.tsv",
  "../external/output/siam/poisLoss/siam_500_50fiveL_do3_bs512_n2_cd8_pois/siam_500_50fiveL_do3_bs512_n2_cd8_pois_assessments.tsv",
  "../external/output/siam/poisLoss/siam_500_50sevenL_do3_bs512_n2_cd8_pois/siam_500_50sevenL_do3_bs512_n2_cd8_pois_assessments.tsv",
  "../external/output/classic_NN/MSEloss/lossMSE_500_5threeL_do3_bs512/lossMSE_500_5threeL_do3_bs512_assessments.tsv",
  "../external/output/classic_NN/MSEloss/lossMSE_500_5fiveL_do3_bs512/lossMSE_500_5fiveL_do3_bs512_assessments.tsv",
  "../external/output/classic_NN/MSEloss/lossMSE_500_5sevenL_do3_bs512/lossMSE_500_5sevenL_do3_bs512_assessments.tsv",
  "../external/output/classic_NN/poisLoss/lossPois_500_5threeL_do3_bs512/lossPois_500_5threeL_do3_bs512_assessments.tsv",
  "../external/output/classic_NN/poisLoss/lossPois_500_5fiveL_do3_bs512/lossPois_500_5fiveL_do3_bs512_assessments.tsv",
  "../external/output/classic_NN/poisLoss/lossPois_500_5sevenL_do3_bs512/lossPois_500_5sevenL_do3_bs512_assessments.tsv",
  "../external/output/siam/mseLoss/PReLU/siam_threeL_do3_bs512_n2_cd8_mse_PReLU/siam_threeL_do3_bs512_n2_cd8_mse_PReLU_assessments.tsv",
  "../external/output/siam/mseLoss/PReLU/siam_fiveL_do3_bs512_n2_cd8_mse_PReLU/siam_fiveL_do3_bs512_n2_cd8_mse_PReLU_assessments.tsv",
  "../external/output/siam/mseLoss/PReLU/siam_sevenL_do3_bs512_n2_cd8_mse_PReLU/siam_sevenL_do3_bs512_n2_cd8_mse_PReLU_assessments.tsv",
  
  
  "../external/output/siam/poisLoss/PReLU/siam_threeL_do3_bs512_n2_cd8_poisLoss_PReLU/siam_threeL_do3_bs512_n2_cd8_poisLoss_PReLU_assessments.tsv",
  "../external/output/siam/poisLoss/PReLU/siam_fiveL_do3_bs512_n2_cd8_poisLoss_PReLU/siam_fiveL_do3_bs512_n2_cd8_poisLoss_PReLU_assessments.tsv",
  "../external/output/siam/poisLoss/PReLU/siam_sevenL_do3_bs512_n2_cd8_poisLoss_PReLU/siam_sevenL_do3_bs512_n2_cd8_poisLoss_PReLU_assessments.tsv",
  
  "../external/output/classic_NN/MSEloss/PReLU/classicMSE_threeL_do3_bs512_PReLU/classicMSE_threeL_do3_bs512_PReLU_assessments.tsv",
  "../external/output/classic_NN/MSEloss/PReLU/classicMSE_fiveL_do3_bs512_PReLU/classicMSE_fiveL_do3_bs512_PReLU_assessments.tsv",
  "../external/output/classic_NN/MSEloss/PReLU/classicMSE_sevenL_do3_bs512_PReLU/classicMSE_sevenL_do3_bs512_PReLU_assessments.tsv",
  
  "../external/output/classic_NN/poisLoss/PReLU/classicPois_yhreeL_do3_bs512_PRelUactivation/classicPois_yhreeL_do3_bs512_PRelUactivation_assessments.tsv",
  "../external/output/classic_NN/poisLoss/PReLU/classicPois_fiveL_do3_bs512_PRelUactivation/classicPois_fiveL_do3_bs512_PRelUactivation_assessments.tsv",
  "../external/output/classic_NN/poisLoss/PReLU/classicPois_sevenL_do3_bs512_PRelUactivation_bn/classicPois_sevenL_do3_bs512_PRelUactivation_bn_assessments.tsv"
  
  
)



path_elem <- "../external/rawInput/Pan_Cancer_test_y.tsv"

path_val <- "../external/rawInput/validation_sets_10folds/Pan_Cancer_validate_y_fold_1.tsv"

path_save <- "../external/plots/"
elements <- c("intergenic" )#, "gc19_pc.cds", "enhancers", "gc19_pc.3utr", "gc19_pc.5utr",
#"gc19_pc.promCore", "gc19_pc.ss", "lncrna.ncrna",
#"lncrna.promCore")


colors <- c("#11266d", "#E6AB02",
            "#B2182B", "#2e6f12",
            "#6E016B", "#7FC97F",
            "#8DD3C7",
            "#666666", "#FB8072",
            "#386CB0",
            'black', '#01234567', '#800000','#008000', '#000080', '#808000', '#800080', '#008080', '#808080',
            '#2E3E51', '#E04C41', '#236CD6', '#43AC63', '#F2B94A', '#9659B2', '#D15357', '#3DB182',
            "#6E016B", "#E6AB02", "#7FC97F",
            "#11266d", "#8DD3C7",
            "#666666", "#FB8072",
            "#386CB0", "#B2182B", "#2e6f12")

save_name = 'all_compare'

######################## scatterplots ###################################
for(i in 1:length(path_preds)){
  for(element in elements){
    save_dotPlots(path_preds[i], path_val, path_elem, colors[i],
                  path_save, element)
  }
}
############## boxPlots ###############
save_boxPlot(element, path_preds, path_elem, path_val, path_save, save_name)

######################## ordered barplot ###################################

measurements = c('mse', 'corr', 'acc')
element = 'intergenic'
for (measurement in measurements){
  save_barPlot(element, path_assessments, path_save, save_name, measurement)
  
}


##########################################################
# patterns = c('mse', 'pois', 'siam')
# 
# save_name = 'siam_layers'
# save_name = 'classic_layers'
# save_name = 'Pois_Siam_Classic_layers'
# save_name = 'MSE_Siam_Classic_layers'
# patterns = c('threeL', 'fiveL', 'sevenL')
# save_barPlot2(element, path_preds, path_elem, path_val, path_save, save_name,
#               measurement, patterns)


########### two bars grouped bar plot ################

losses <- c('mseloss', 'poisloss')
measurements = c('acc', 'corr', 'mse')
patterns <- c('hreeL', 'fiveL', 'sevenL')
for (loss in losses) {
  for (measurement in measurements){
    save_groupedBarPlots_SC(element, path_assessments, path_save, save_name,
                            measurement, loss, patterns)
  }
}


##################### four bars grouped bar plot #################

losses <- c('mseloss', 'poisloss')
patterns <- c('hreeL', 'fiveL', 'sevenL')
for (loss in losses) {
  for (measurement in measurements){
    save_groupedBarPlots_SC_actFunc(element, path_assessments, path_save, save_name,
                                    measurement,loss, patterns)
  }
}


