#!/bin/bash
#SBATCH -J Manhattan_plots
#SBATCH -o Manhattan_plots.o%j
#SBATCH -e Manhattan_plots.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 2:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate prs1

# Plot Manhattan plots
Rscript 05c_plot_ManhattanPlots.R
