#!/bin/bash
#SBATCH -J EUR_train_geno
#SBATCH -o EUR_train_geno.o%j
#SBATCH -e EUR_train_geno.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 8
#SBATCH -t 5:00:00
#SBATCH -A Harpak-Lab-GWAS
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


set -e

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'


for chromosome in $(seq 1 22);
do
   $plink2 \
      --bfile data/ukb_merged/chr${chromosome} \
      --keep data/ukb_populations/LD_CEU_train.txt \
      --make-bed \
      --out data/LDpred/LD_CEU_train_${chromosome}

   printf "data/LDpred/LD_CEU_train_%s\n" $chromosome >> data/LDpred/LD_CEU_train_merged_list.txt
done

$plink \
  --merge-list data/LDpred/LD_CEU_train_merged_list.txt \
  --make-bed  \
  --out data/LDpred/LD_CEU_merged

#$plink \
#  --bfile data/ukb_merged/merged \
#  --keep data/ukb_populations/CEU_all.txt \
#  --make-bed \
#  --out data/ukb_merged/CEU_all
