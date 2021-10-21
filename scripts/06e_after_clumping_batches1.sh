#!/bin/bash
#SBATCH -J post-clump
#SBATCH -o post-clump.o%j
#SBATCH -e post-clump.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 5:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

module load anaconda/3-4.2.0
source activate prs1

plink="/work2/06568/joyce_w/stampede2/software/plink/plink/plink"

cat data/prs/*_threshold_4.txt | uniq > data/prs/all_prs_snps.txt

for chromosome in $(seq 1 22);
do
  $plink \
    --bfile data/ukb_merged/chr${chromosome} \
    --extract data/prs/all_prs_snps.txt \
    --remove data/ukb_populations/WB_train.txt \
    --memory 35000 \
    --make-bed \
    --out data/prs/chr${chromosome}_temp

  printf "data/prs/chr%s_temp\n" $chromosome >> data/prs/prs_merged_list.txt
done

# Combine extracted SNPs across chromosomes into a single Plink 1 file of all
#  SNPs that meet the least extreme level of significance for any of the traits.
$plink \
  --merge-list data/prs/prs_merged_list.txt \
  --make-bed \
  --memory 35000 \
  --out data/prs/CT_merged

rm data/prs/chr*_temp.*
