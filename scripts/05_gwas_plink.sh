#!/bin/bash
#SBATCH -J GWAS
#SBATCH -o GWAS.o%j
#SBATCH -e GWAS.o%j
#SBATCH -p normal
#SBATCH -N 10
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

#Central directory version of plink not up to-date
plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC
do
  # Can't use Plink optimization for all phenotypes simultaneously because
  # each phenotype uses different samples.
  # The below hashtags refer to past GWAS, using genotype and UKBB
  for chromosome in $(seq 1 22);
  do
    $plink2\
      --bfile data/ukb_merged/chr${chromosome} \
      --keep data/ukb_populations/${phenotype}.txt \
      --pheno data/phenotypes/full_phenotypes.pheno \
      --pheno-name $phenotype \
      --require-pheno $phenotype \
      --covar data/ukb_merged/covar_all_samples.covar \
      --covar-name age sex_covar age_sq age_sex age_sq_sex $(printf "PC%i " $(seq 1 20)) \
      --require-covar age sex_covar age_sq age_sex age_sq_sex $(printf "PC%i " $(seq 1 20)) \
      --vif 100000 \
      --memory 35000 \
      --glm hide-covar \
      --out data/gwas_results/${phenotype}.chr${chromosome}
  done
done
