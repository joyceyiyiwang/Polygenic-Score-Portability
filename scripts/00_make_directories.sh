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
      data/intersecting_filtered \
      data/kgp_filtered \
      data/kgp_merged \
      data/kgp_meta \
      data/ukb_filtered \
      data/ukb_merged \
      data/ukb_meta \
      data/ukb_populations #will include separated pop groups
mkdir data/models \
      data/phenotypes \
      data/gwas_results \
      data/prs #only scores from C+T approach
mkdir data/kgp_populations \
      data/fst \
      data/LDpred \
      data/LDpred/prs #contains PRS for target individuals
mkdir data/LDpred/tmp-data #memory storage for LDpred code
mkdir data/LDpred/val_prs #contains PRS for validation set
mkdir data/prs_comparisons #contains agg data from C+T & LDpred
mkdir data/theory #LD & MAF patterns from target groups
mkdir data/theor_herit #LDSC results to calculate h2
mkdir data/theoretical #Destination of results from theor predictions
mkdir img

