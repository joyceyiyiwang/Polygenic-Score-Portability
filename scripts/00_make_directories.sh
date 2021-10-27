#!/bin/bash
#SBATCH -J prepare_directories
#SBATCH -o prepare_directories.o%j
#SBATCH -e prepare_directories.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

mkdir data \
      data/ambiguous_indel_snps \
      data/ukb_filtered \
      data/ukb_merged \
      data/ukb_populations #will include separated pop groups
mkdir data/phenotypes \
      data/gwas_results \
      data/prs #only scores from C+T approach
mkdir data/LDpred \
      data/fst
mkdir img

