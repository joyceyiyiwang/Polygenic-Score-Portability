#!/bin/bash
#SBATCH -J make_filtered_genotypes_files
#SBATCH -o make_filtered_genotypes_files.o%j
#SBATCH -e make_filtered_genotypes_files.o%j
#SBATCH -p normal
#SBATCH -N 5
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink="/work2/06568/joyce_w/stampede2/software/plink/plink/plink"
plink2="/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2"

# The files must be copied to my directory in order for this script to work
gt='/work2/06568/joyce_w/stampede2/pgs_portability/data/genotype_calls'
meta='/work2/06568/joyce_w/stampede2/pgs_portability/data/ukb_meta'

#module load anaconda3
#source activate prs1
source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate prs1

for i in $(seq 1 22);
do
  # Identify indels and ambiguous variants and write them to a file
  python 01a_get_ambiguous_indel_snps.py \
    ${gt}/ukb_snp_chr${i}_v2.bim \
    -o data/ambiguous_indel_snps/chr${i}.snps

  # Convert from Plink 1 to Plink 2 and compute the SNPs to be dropped
  # Keep the unrelated individuals
  # Remove participants who opted out from the study
  $plink2 \
    --bfile ${gt}/ukb_snp_chr${i}_v2  \
    --bed ${gt}/ukb22418_c${i}_b0_v2.bed \
    --fam ${gt}/ukb22418_c${i}_b0_v2_s488244.fam \
    --keep data/extracted_phenotypes/ukb.unrelated.id.txt \
    --remove ${meta}/w61666_20210809.csv \
    --exclude data/ambiguous_indel_snps/chr${i}.snps \
    --maf 0.01 \
    --geno 0.01 \
    --indep-pairwise 1000 kb 1 0.2 \
    --make-pgen \
    --out data/ukb_filtered/chr${i}

  # Write a new file using only filtered SNPs for 1000 Genomes (Plink 1 format)
  $plink2 \
    --pfile data/ukb_filtered/chr${i} \
    --extract data/ukb_filtered/chr${i}.prune.in \
    --make-bed \
    --out data/ukb_filtered/chr${i}

  rm data/ukb_filtered/chr${i}.pgen
  rm data/ukb_filtered/chr${i}.pvar
  rm data/ukb_filtered/chr${i}.psam 

  # Append the output file path to a new file (for merging them all below)
  printf "data/ukb_filtered/chr%s\n" $i >> data/ukb_merged/ukb_merged_list${i}.txt

  $plink \
    --merge-list data/ukb_merged/ukb_merged_list${i}.txt \
    --make-bed \
    --out data/ukb_merged/chr${i}

  rm data/ukb_filtered/chr${i}.b*
  rm data/ukb_filtered/*.fam
  rm data/ukb_filtered/*.log

  # Append the output file path to a new file (for merging them all below)
  printf "data/ukb_merged/chr$%s\n" $i >> data/ukb_merged/ukb_merged_list.txt
done

# Correct the file names in data/ukb_merged/ukb_merged_list.txt before continuing

#$plink \
#  --merge-list data/ukb_merged/ukb_merged_list.txt \
#  --make-bed \
#  --out data/ukb_merged/merged

#$plink2 \
#  --bfile data/ukb_merged/merged \
#  --make-pgen \
#  --out data/ukb_merged/merged


