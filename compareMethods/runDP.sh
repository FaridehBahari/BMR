#!/bin/bash

# List of cancer types
cohorts=("Uterus-AdenoCA" "Stomach-AdenoCA" "ColoRect-AdenoCA"
	 "Skin-Melanoma" "Prost-AdenoCA" "Panc-Endocrine"
	 "Panc-AdenoCA" "Ovary-AdenoCA" "Lymph-CLL" "Lymph-BNHL"
	 "Lung-SCC" "Lung-AdenoCA" "Liver-HCC" "Kidney-RCC"
	 "Head-SCC" "Eso-AdenoCa" "CNS-Medullo" 
	 "CNS-GBM" "Breast-AdenoCa" "Bone-Osteosarc" "Bone-Leiomyo"
	 "Bladder-TCC"  "Biliary-AdenoCA"

         "Kidney-ChRCC"   "CNS-PiloAstro"  "Myeloid-MPN"  "Thy-AdenoCA" "CNS-Oligo" "Cervix-SCC"  )

# Activate conda environment
conda activate DP2

# Loop through each cancer type
for cohort in "${cohorts[@]}"; do
# Create the output directory
mkdir -p output/$cohort

# train the model
nohup driverpower model --feature ../make_features/external/ftrMtrix/var_1372features.h5 \
--response ../iDriver/extdata/procInput/BMRs/observed/$cohort/train_y.tsv \
--method GBM --gbmParam ./driverpower/xgb_param.pkl \
--name $cohort --predict \
--modelDir ./output/$cohort/ > ../${cohort}_model_nohup

# use the model on the test data to save nPreds for each element
driverpower infer --feature ../make_features/external/ftrMtrix/pcawg_features_1372.h5 \
--response ../iDriver/extdata/procInput/BMRs/observed/$cohort/test_y.tsv \
--model ./output/$cohort/${cohort}.GBM.model.pkl \
--name $cohort --outDir ./output/$cohort/ > ../${cohort}_infer_nohup
done
  
