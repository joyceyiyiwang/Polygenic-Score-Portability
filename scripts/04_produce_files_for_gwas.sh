#!/bin/bash
#SBATCH -J covariates_phenotypes
#SBATCH -o covariates_phenotypes.o%j
#SBATCH -e covariates_phenotypes.o%j
#SBATCH -p normal
#SBATCH -N 2
#SBATCH -n 6
#SBATCH -t 10:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

#module load anaconda/3-4.2.0
#module load R
#source activate prs
source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate prs1

# Creates covariate file with age, sex, age*sex, age^2, age^2 * sex and PC1, ..., PC20
Rscript 04a_create_covariates.R

# Creates data/phenotypes/full_phenotypes.pheno, combining all phenotypes for GWAS
# and inverse rank normal transforming phenotypes
Rscript 04b_create_phenotypes_file.R
