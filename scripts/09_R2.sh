#!/bin/bash
#SBATCH -J R2
#SBATCH -o R2.o%j
#SBATCH -e R2.o%j
#SBATCH -p normal
#SBATCH -N 2
#SBATCH -n 6
#SBATCH -t 10:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

Rscript 09a_run_R2_models.R
