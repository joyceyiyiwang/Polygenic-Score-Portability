#!/bin/bash
#SBATCH -J make_filtered_genotypes_files
#SBATCH -o make_filtered_genotypes_files.o%j
#SBATCH -e make_filtered_genotypes_files.o%j
#SBATCH -p normal
#SBATCH -N 1 #5
#SBATCH -n 4 #16
#SBATCH -t 1:00:00 #24:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink="/work2/06568/joyce_w/stampede2/software/plink/plink"
plink2="/work2/06568/joyce_w/stampede2/software/plink/plink2"

wb='/corral-repl/utexas/Recombining-sex-chro/ukb/data/genotype_calls'
#nonwb='/rigel/mfplab/projects/ukb_hakhamanesh/imputed/plink_neale_nonWB'

#module load anaconda3
#source activate prs1
source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate prs1

#for i in $(seq 1 23);
for i in $(seq 22 22);
do
  # Identify indels and ambiguous variants and write them to a file
  python 01a_get_ambiguous_indel_snps.py \
    # Having an erro running the following line: 
    # FileNotFoundError: [Errno 2] No such file or directory: 
    # '/corral-repl/utexas/Recombining-sex-chro/ukb/data/genotype_calls/ukb_snp_chr22_v2.bim'
    ${wb}/ukb_snp_chr${i}_v2.bim \
    -o data/ambiguous_indel_snps/chr${i}.snps

  # Convert from Plink 1 to Plink 2 and compute the SNPs to be dropped
  $plink2 \
    --bfile ${wb}/ukb_snp_chr${i}_v2  \
    --keep ${wb}/ukb22418_chr${i}_b0_v2_s488244.fam \
    --remove /work2/06568/joyce_w/stampede2/pgs_portability/data/ukb_meta/excluded_samples.sam \
    --exclude data/ambiguous_indel_snps/chr${i}.snps \
    --maf 0.01 \
    --geno 0.01 \
    --indep-pairwise 1000 kb 1 0.2 \
    --make-pgen \
    --out data/ukb_filtered/wb_chr${i}

  # Write a new file using only filtered SNPs for 1000 Genomes (Plink 1 format)
  $plink2 \
    --pfile data/ukb_filtered/wb_chr${i} \
    --extract data/ukb_filtered/wb_chr${i}.prune.in \
    --make-bed \
    --out data/ukb_filtered/wb_chr${i}

  rm data/ukb_filtered/wb_chr${i}.pgen
  rm data/ukb_filtered/wb_chr${i}.pvar
  rm data/ukb_filtered/wb_chr${i}.psam 

#  $plink2 \
#    --bfile ${nonwb}/chr${i} \
#    --keep ${nonwb}/ukb_filter_indivs_nowhiteB_unrelated.txt \
#    --remove /rigel/mfplab/projects/prs-portability/data/ukb_meta/excluded_samples.sam \
#    --extract data/ukb_filtered/wb_chr${i}.prune.in \
#    --make-bed \
#    --out data/ukb_filtered/nonwb_chr${i}


  # Append the output file path to a new file (for merging them all below)
  printf "data/ukb_filtered/wb_chr%s\n" $i >> data/ukb_merged/ukb_merged_list${i}.txt
#  printf "data/ukb_filtered/nonwb_chr%s\n" $i >> data/ukb_merged/ukb_merged_list${i}.txt

  $plink \
    --merge-list data/ukb_merged/ukb_merged_list${i}.txt \
    --make-bed \
    --out data/ukb_merged/chr${i}

  rm data/ukb_filtered/wb_chr${i}.b*
#  rm data/ukb_filtered/nonwb_chr${i}.b*
  rm data/ukb_filtered/*.fam
  rm data/ukb_filtered/*.log

  printf "data/ukb_merged/chr$%s\n" $i >> data/ukb_merged/ukb_merged_list.txt
done


#$plink \
#  --merge-list data/ukb_merged/ukb_merged_list.txt \
#  --make-bed \
#  --out data/ukb_merged/merged

#$plink2 \
#  --bfile data/ukb_merged/merged \
#  --make-pgen \
#  --out data/ukb_merged/merged


