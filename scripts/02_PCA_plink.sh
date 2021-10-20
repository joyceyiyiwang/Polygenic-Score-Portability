#!/bin/bash
#SBATCH -J PCA
#SBATCH -o PCA.o%j
#SBATCH -e PCA.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 8
#SBATCH -t 5:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

# Fail if any command fails
set -e

#Central directory of plink software is not up to-date
plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Compute PCA on the 1000 Genomes project filtered and merged genotypes
$plink2 \
  --bfile data/kgp_merged/merged \
  --freq counts \
  --pca allele-wts 20 \
  --out data/kgp_merged/merged

#  Apply loadings from PCA on 1000 Genomes to 1000 Genomes (for consistency)
$plink2 \
  --bfile data/kgp_merged/merged \
  --read-freq data/kgp_merged/merged.acount \
  --score data/kgp_merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-25 \
  --out data/kgp_merged/projection

# Apply loadings from PCA computed on 1000 Genomes to the merged UK Biobank file
$plink2 \
  --bfile data/ukb_merged/merged \
  --read-freq data/kgp_merged/merged.acount \
  --score data/kgp_merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-25 \
  --out data/ukb_merged/projection
