import os

variable = """#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name={}_test
#SBATCH -c 1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16gb
set -e
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'
phenotype={}
"""

fixed = """# Create PRS score files for several p-value thresholds. Files numbered 0-4, corresponding values above.
for i in 6 28 29 # 6 - original betas, 28 - average of 20 sets, 29 - average of filtered sets
do
$plink2 \
--bfile data/ukb_merged/merged \
--remove data/ukb_populations/EUR_all.txt \
--score data/LDpred1/${phenotype}_final_LDpred.tsv 2 5 ${i} header no-mean-imputation \
--memory 10000 \
--out data/LDpred1/prs/${phenotype}_LDpred_scores_${i} 

done

rm ldpred_scripts/temp_scripts/{phenotype}.sh
"""


def main(phenotype):
    return variable.format(phenotype, phenotype) + fixed


if __name__ == "__main__":
    phenotypes = ['BMI','Height','RBC', 'Platelet', 'MCV', 'Monocyte', 'WBC', 'MCH', 'Eosinophil', 'Lymphocyte',
                  'DBP','SBP','Hb','Ht','MCHC','Neutrophil','Basophil']
    for phenotype in phenotypes:
        with open(f'ldpred_scripts/temp_scripts/{phenotype}.sh', 'w') as f:
            f.write(main(phenotype))
        os.system(f'sbatch ldpred_scripts/temp_scripts/{phenotype}.sh')
