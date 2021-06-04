import os

script_string = """#!/bin/bash
    #
    #SBATCH --account=mfplab
    #SBATCH --job-name={}
    #SBATCH -c 1
    #SBATCH --time=48:00:00
    #SBATCH --mem-per-cpu=8gb

    set -e

    plink='/rigel/mfplab/users/jm4454/plink/plink'
    plink2='/rigel/mfplab/users/jm4454/plink/plink2'

    phenotype={}

    for chromosome in $(seq 1 22);
    do
       $plink2 \ 
        --bfile data/ukb_merged/chr${chromosome} \
        --keep data/ukb_populations/${phenotype}.txt \
        --pheno data/phenotypes/full_phenotypes.pheno \
        --pheno-name $phenotype \
        --require-pheno $phenotype \
        --covar data/ukb_merged/covar_all_samples.covar \
        --covar-name age sex_covar age_sq age_sex age_sq_sex $(printf 'PC%i_AVG ' $(seq 1 20)) \
        --require-covar age sex_covar age_sq age_sex age_sq_sex $(printf 'PC%i_AVG ' $(seq 1 20)) \
        --vif 100000 \
        --memory 35000 \
        --glm hide-covar \
        --out data/gwas_results1/${phenotype}.chr${chromosome}
    done
    
    rm ${phenotype}.sh
    """
        
def main(phenotype):
    return script_string.format(phenotype, phenotype)


if __name__ == "__main__":
        
    phenotypes = ['Height', 'BMI' ,'RBC', 'Platelet', 'MCV', 'Monocyte', 'WBC', 'MCH', 'Eosinophil', 'Lymphocyte',
                  'DBP','SBP','Hb','Ht','MCHC','Neutrophil','Basophil']
    for phenotype in phenotypes:
        with open(f'{phenotype}.sh', 'w') as f:
            f.write(main(phenotype))
        os.system(f'sbatch {phenotype}.sh')
