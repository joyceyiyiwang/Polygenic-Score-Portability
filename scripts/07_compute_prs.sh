#!/bin/bash
#SBATCH -J PRS
#SBATCH -o PRS.o%j
#SBATCH -e PRS.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 6:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e 

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC
do
    for threshold in 0 1 2 3 4
    do
    # Score each individual using the SNPs below a p-value threshold and GWAS betas
    $plink2 \
        --bfile data/ukb_merged/merged \
        --keep data/ukb_populations/${phenotype}.txt \
        --extract data/prs/${phenotype}_threshold_${threshold}.txt \
        --score data/gwas_results/${phenotype}_combined.glm.linear 3 6 9 header no-mean-imputation \
        --memory 10000 \
        --out data/prs/${phenotype}_${threshold}_scores
    done
done

