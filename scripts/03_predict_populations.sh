#!/bin/sh
#SBATCH -J classify_populations
#SBATCH -o classify_populations.o%j
#SBATCH -e classify_populations.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 10
#SBATCH -t 3:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

#module load anaconda/3-4.2.0
#source activate prs1
source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate prs1

# Train a classifier on 1000 Genomes to predict super population labels
# Save classifier and evaluate it on the UK Biobank PCs
python 03a_classify_ukb.py

# Create files with sample IDs for each group
# Also does train/test splitting for the "target" individuals.
python 03b_separate_populations.py

# Rename erroneously named file data/ukb_populations/EUR_all = all training EUR

# Create plots of PCAs of 1000 Genomes population classifier on 1KG and  UK Biobank
#Rscript 03c_plot_pca.R
