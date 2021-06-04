import os

script_string = """#!/bin/bash
    #
    #SBATCH --account=mfplab
    #SBATCH --job-name=CT_{}
    #SBATCH -c 1
    #SBATCH --time=03:00:00
    #SBATCH --mem-per-cpu=8gb

    set -e 
    plink='/rigel/mfplab/users/jm4454/plink/plink'
    phenotype={}
    thresholds=(5e-8 1e-5 1e-4 1e-3 1e-2)
    
    # Clump GWAS results 
    for chromosome in $(seq 1 22);
    do
    # Convert the GWAS output file from the Plink 2 to Plink 1 
    python scripts2/06b_convert_plink2_glm_to_plink1.py \
    data/gwas_results1/${phenotype}.chr${chromosome}.${phenotype}.glm.linear \
    --output data/gwas_results1/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc

    # Clump GWAS results using LD from 1000 Genomes EUR super-population 
    $plink \
    --memory 35000 \
    --bfile data/LDpred1/LD_EUR_train_${chromosome} \
    --clump data/gwas_results1/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc \
    --clump-p1 0.01 \
    --clump-p2 1 \
    --clump-r2 0.5 \
    --clump-kb 250 \
    --out data/gwas_results1/${phenotype}.chr${chromosome}.${phenotype}
  	done

  	# Combine clumped SNPs across chromosomes
  	head -n 1 data/gwas_results1/${phenotype}.chr1.${phenotype}.clumped > 
    	data/gwas_results1/${phenotype}_combined.clumped 
  	tail -n +2 -q data/gwas_results1/${phenotype}.chr*.${phenotype}.clumped >>  
    	data/gwas_results1/${phenotype}_combined.clumped 

  	# Create files of SNPs meeting several p-value thresholds. Files numbered 0-4.
  	for threshold in 0 1 2 3 4
  	do
    	# Further filter clumped SNPs using p-value thresholds (removes multiallelic SNPs)
    	python scripts2/06c_filter_snps_for_prs.py \
      	data/gwas_results1/${phenotype}_combined.clumped \
      	--threshold ${thresholds[$threshold]} \
      	--output data/prs1/${phenotype}_threshold_${threshold}.txt
  	done

  	# Create combined GWAS result files for each phenotype
  	python scripts2/06d_combine_glm_threshold_4.py \
    	data/gwas_results1/${phenotype}.chr*.${phenotype}.glm.linear \
    	--keep data/prs1/${phenotype}_threshold_4.txt \
    	--output data/gwas_results1/${phenotype}_combined.glm.linear
        
	# Delete the script itself
    rm ${phenotype}_clumping.sh
	"""


def main(phenotype):
    return script_string.format(phenotype, phenotype)

if __name__ == "__main__":
        
    phenotypes = ['Height', 'BMI' ,'RBC', 'Platelet', 'MCV', 'Monocyte', 'WBC', 'MCH', 'Eosinophil', 'Lymphocyte',
                  'DBP','SBP','Hb','Ht','MCHC','Neutrophil','Basophil']
    for phenotype in phenotypes:
        with open(f'{phenotype}_clumping.sh', 'w') as f:
            f.write(main(phenotype))
        os.system(f'sbatch {phenotype}_clumping.sh')


