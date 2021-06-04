#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=freq_ld_groups
#SBATCH -c 1
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8gb

set -e

plink='/rigel/mfplab/users/jm4454/plink/plink'
#plink2='/moto/palab/users/jm4454/plink/plink2'

for k in fst wpc;
do

for j in $(seq 1 40);
do

$plink \
--bfile data/ukb_merged/merged \
--keep data/theory/${k}_${j}.txt \
--freq \
--out data/theory/${k}_${j}

rm data/theory/${k}_${j}.nosex
rm data/theory/${k}_${j}.log
 
for i in $(seq 1 22);
do
$plink \
--bfile data/ukb_merged/chr${i} \
--keep data/theory/${k}_${j}.txt \
--r \
--ld-window-kb 100 \
--out data/theory/${k}_${j}_${i}

rm data/theory/${k}_${j}_${i}.nosex
rm data/theory/${k}_${j}_${i}.log
done
done
done
