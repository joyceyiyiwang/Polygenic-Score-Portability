#!/bin/bash
#SBATCH -J make_filtered_1kg
#SBATCH -o make_filtered_1kg.o%j
#SBATCH -e make_filtered_1kg.o%j
#SBATCH -p normal
#SBATCH -N 5
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e
#wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz 
#wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples.20130502.ALL.ped 
#mv integrated_call_samples.20130502.ALL.ped data/kgp_meta
#mv ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz data/kgp_meta
#gzip -d data/kgp_meta/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'


for i in $(seq 1 22);
do
  # Using existing 1KG file which includes all MAF 1%, genotype 1%, and pruned intersecting SNPs with UKBB.
  # From existing variants, extract those that were filtered from UKBB genotype data.
  $plink2 \
  --vcf data/kgp_filtered/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  --max-alleles 2 \
  --extract data/ukb_merged/chr${i}.prune.in \
  --make-bed \
  --out data/kgp_filtered/chr${i}

  # Add the file name to a running list (for easier merging below)
  printf "data/kgp_filtered/chr%s\n" $i >> data/kgp_merged/kgp_merged_list.txt

done

for i in $'X';
do
  # Using existing 1KG file which includes all MAF 1%, genotype 1%, and pruned intersecting SNPs with UKBB.
  # From existing variants, extract those that were filtered from UKBB genotype data.
  $plink2 \
  --vcf data/kgp_filtered/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz \
  --max-alleles 2 \
  --extract data/ukb_merged/chr${i}.prune.in \
  --make-bed \
  --out data/kgp_filtered/chr${i}

  # Add the file name to a running list (for easier merging below)
  printf "data/kgp_filtered/chr%s\n" $i >> data/kgp_merged/kgp_merged_list.txt

done

for i in $'Y';
do
  # Using existing 1KG file which includes all MAF 1%, genotype 1%, and pruned intersecting SNPs with UKBB.
  # From existing variants, extract those that were filtered from UKBB genotype data.
  $plink2 \
  --vcf data/kgp_filtered/ALL.chr${i}.phase3_integrated_v1b.20130502.genotypes.vcf.gz \
  --max-alleles 2 \
  --extract data/ukb_merged/chr${i}.prune.in \
  --make-bed \
  --out data/kgp_filtered/chr${i}

  # Add the file name to a running list (for easier merging below)
  printf "data/kgp_filtered/chr%s\n" $i >> data/kgp_merged/kgp_merged_list.txt

done

$plink \
  --merge-list data/kgp_merged/kgp_merged_list.txt \
  --make-bed \
  --out data/kgp_merged/merged

rm data/kgp_filtered/chr*

#$plink2 \
#  --bfile data/kgp_merged/merged \
#  --make-pgen \
#  --out data/kgp_merged/merged
