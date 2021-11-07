#!/bin/bash
#SBATCH -J fst
#SBATCH -o fst.o%j
#SBATCH -e fst.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate prs1

Rscript 08a_edit_merged_fam.R
python 08b_fst_parallel.py

#Do afterwards
#Rscript 08c_final_fst_formatting.R
conda deactivate
