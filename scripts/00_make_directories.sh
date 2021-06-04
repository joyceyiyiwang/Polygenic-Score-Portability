#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=prepare_directories
#SBATCH -c 1
#SBATCH --time=5:00
#SBATCH --mem-per-cpu=1gb

mkdir data/ambiguous_indel_snps \
      data/intersecting_filtered \
      data/kgp_filtered \
      data/kgp_merged \
      data/kgp_meta \
      data/ukb_filtered \
      data/ukb_merged \
      data/ukb_meta \
      data/ukb_populations \ #will include separated pop groups
      data/models \
      data/phenotypes \
      data/gwas_results \
      data/prs \ #only scores from C+T approach
      data/kgp_populations \
      data/fst \ 
      data/LDpred \ 
      data/LDpred/prs \ #contains PRS for target individuals
      data/LDpred/tmp-data \ #memory storage for LDpred code
      data/LDpred/val_prs \ #contains PRS for validation set
      data/prs_comparisons \ #contains agg data from C+T & LDpred
      data/theory \  #LD & MAF patterns from target groups
      data/theor_herit \  #LDSC results to calculate h2
      data/theoretical \ #Destination of results from theor predictions
      img \

