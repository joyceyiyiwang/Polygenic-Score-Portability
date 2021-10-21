#!/bin/bash
#SBATCH -J CT
#SBATCH -o CT.o%j
#SBATCH -e CT.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 30:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e 

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate prs1

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
phenotype={}
thresholds=(5e-8 1e-5 1e-4 1e-3 1e-2)

# Clump GWAS results 
for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC
do
    for chromosome in $(seq 1 22);
    do
    # Convert the GWAS output file from the Plink 2 to Plink 1 
    python 06b_convert_plink2_glm_to_plink1.py \
    data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.linear \
    --output data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc

    # Clump GWAS results using LD from 1000 Genomes EUR super-population 
    $plink \
    --memory 35000 \
    --bfile data/LDpred/LD_WB_train_${chromosome} \
    --clump data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc \
    --clump-p1 0.01 \
    --clump-p2 1 \
    --clump-r2 0.5 \
    --clump-kb 250 \
    --out data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}
  	done

  	# Combine clumped SNPs across chromosomes
  	head -n 1 data/gwas_results/${phenotype}.chr1.${phenotype}.clumped > data/gwas_results/${phenotype}_combined.clumped 
  	tail -n +2 -q data/gwas_results/${phenotype}.chr*.${phenotype}.clumped >> data/gwas_results/${phenotype}_combined.clumped 

  	# Create files of SNPs meeting several p-value thresholds. Files numbered 0-4.
  	for threshold in 0 1 2 3 4
  	do
    	# Further filter clumped SNPs using p-value thresholds (removes multiallelic SNPs)
    	python 06c_filter_snps_for_prs.py \
      	data/gwas_results/${phenotype}_combined.clumped \
      	--threshold ${thresholds[$threshold]} \
      	--output data/prs/${phenotype}_threshold_${threshold}.txt
  	done

  	# Create combined GWAS result files for each phenotype
  	python 06d_combine_glm_threshold_4.py \
    	data/gwas_results/${phenotype}.chr*.${phenotype}.glm.linear \
    	--keep data/prs/${phenotype}_threshold_4.txt \
    	--output data/gwas_results/${phenotype}_combined.glm.linear
done

