#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=ldpred_adj
#SBATCH -c 1
#SBATCH --time=01-00:00:00
#SBATCH --mem-per-cpu=6gb

set -e

module load anaconda
module load R
source activate prs1

Rscript ldpred_scripts/04a_find_set.R
python ldpred_scripts/04b_test_ldpred_scores_batch.py
