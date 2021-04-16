#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=groups_fst
#SBATCH -c 2
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8gb

set -e

plink='/rigel/mfplab/users/jm4454/plink/plink'
#plink2='/moto/palab/users/jm4454/plink/plink2'


for j in $(seq 1 40);
do

$plink \
--bfile data/ukb_merged/merged \
--keep data/theory/fst_${j}.txt \
--freq \
--out data/theory/fst_${j}

rm data/theory/fst_${j}.nosex
rm data/theory/fst_${j}.log
 
for i in $(seq 1 22);
do
$plink \
--bfile data/ukb_merged/chr${i} \
--keep data/theory/fst_${j}.txt \
--r \
--ld-window-kb 100 \
--out data/theory/fst_${j}_${i}

rm data/theory/fst_${j}_${i}.nosex
rm data/theory/fst_${j}_${i}.log
done
done

