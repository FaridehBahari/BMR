import os
import subprocess

# List of cohorts
# 'Pancan-no-skin-melanoma-lymph', 
cohorts = [
    'Pan_Cancer', "Liver-HCC", 
    "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
    "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
    "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
    "Breast-AdenoCa", "Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",         
    "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
    "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",    
    "Lymph-NOS", "CNS-Medullo", "Myeloid-AML", "CNS-Oligo", "Cervix-SCC",
    "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
    "Cervix-AdenoCA", "Breast-DCIS", "Bone-Cart", "Myeloid-MDS"
]

# Path to the template ini file
template_ini_path = 'configs/rate_based/sim_setting_iDriver.ini'

# Read the template ini file
with open(template_ini_path, 'r') as file:
    ini_content = file.read()

# Loop over the cohorts and run the algorithm
for cohort in cohorts:
    print(cohort)
    # Replace the cohort name in the ini content
    modified_ini_content = ini_content.replace('Pancan-no-skin-melanoma-lymph', cohort)
    
    # Save the modified ini content to a new file
    modified_ini_path = f'configs/rate_based/sim_setting_iDriver_{cohort}.ini'
    with open(modified_ini_path, 'w') as file:
        file.write(modified_ini_content)
    
    # Run the algorithm with the modified ini file
    command = f'python -u RUN_BMR.py {modified_ini_path}'
    subprocess.run(command, shell=True)

# run the following command from BMR directory:
# nohup python -u run_BMR_allCohorts.py > ../run_allCohorts_GBM &
