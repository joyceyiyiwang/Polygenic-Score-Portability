import os 

rscript = """
.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(tidyverse)
library(bigsnpr)
library(runonce)
library(Matrix)

group <- {}
print("LD Matrix")

# Creates .bk and .rd files
for(chr in 1:22){{
  snp_readBed(paste0("theor_herit/fst_",group,"_",chr,".bed"))
}}

NCORES <- 1
n_val <- 750


for (chr in 1:22) {{
  rds_file <- paste0("theor_herit/fst_",group,"_",chr,".rds")
  ukb <- snp_attach(rds_file)
  G <- ukb$genotypes
  CHR <- as.integer(ukb$map$chromosome)
  POS <- ukb$map$physical.pos
  
  
  POS2 <- bigsnpr::snp_asGeneticPos(CHR, POS, dir = "theor_herit/tmp-data",
                                    ncores = NCORES)
  
  ind.chr <- which(CHR == chr)
  save_run(snp_cor(G, ind.row=1:nrow(G),ind.col = ind.chr,
                   alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000,
                   ncores = NCORES),
           file = paste0("theor_herit/tmp-data/fst_",group,"_corr_chr", chr, ".rds"))
  
}}

set.seed(1)

##################
tmp <- tempfile(tmpdir = "theor_herit/tmp-data")
print("Corr")
for (chr in 1:22){{
  
#  rds_file <- paste0("theor_herit/fst_",group,"_",chr,".rds")
#  ukb <- snp_attach(rds_file)
#  G <- ukb$genotypes
#  CHR <- as.integer(ukb$map$chromosome)
#  POS <- ukb$map$physical.pos
  
  
  # Indices in 'sumstats'
 # ind.chr <- which(sumstats$chr == chr)
  # indices in 'G' / rsid replacing `NUM_ID`
 # ind.chr2 <- sumstats$rsid[ind.chr]
  # Indices in 'corr' / ind.chr replacing which(CHR == chr) 
 # ind.chr3 <- match(ind.chr2,sumstats$rsid[ind.chr])
  
  corr0 <- readRDS(paste0("theor_herit/tmp-data/fst_",group,"_corr_chr", chr, ".rds"))#[ind.chr3, ind.chr3]
  
  if (chr == 1){{
    ld <- colSums(corr0^2)
 #   corr <- as_SFBM(corr0, tmp)
  }} else {{
    ld <- c(ld, colSums(corr0^2))
  #  corr$add_columns(corr0, nrow(corr))
  }}
  
}}

#########################
print("LDSC")
phenotype <- c("Height","MCV","Platelet")

h2_est <- c()

for (pheno in phenotype){{
 
 # sumstat_files <- c(paste0("data/LDpred1/",pheno,"_merged_sumstats_ldpred.tsv"))
  #  df_beta <- read_tsv(sumstat_files)
	df_beta <- tibble()
	for(i in 1:22){{
	file <- paste0("data/gwas_results1/",pheno,".chr",i,".",pheno,".glm.linear")
	file_tsv <- read_tsv(file)
	df_beta <- df_beta %>% bind_rows(file_tsv)
}}
colnames(df_beta) <- c("chr","pos","rsid","a0","alt","a1","test","n_eff","beta","beta_se","t_stat","p","errcode")

    
    (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                    sample_size = n_eff, blocks = NULL,intercept=1)))
    (ld_h2_est <- ldsc[["h2"]])
    h2_est <- c(h2_est,ld_h2_est)
    
}}

final_table <- tibble(phenotype=phenotype,
                      h2 = h2_est, 
			group = rep(group,3))
final_name <- paste0("theor_herit/fst_",group,"_h2.tsv")
final_table %>% write_tsv(final_name)

"""

variable = """#!/bin/bash
#SBATCH --account=mfplab
#SBATCH --job-name=fst_{}_herit
#SBATCH -c 1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=12gb
set -e

module load anaconda
module load R
source activate prs1

plink="/rigel/mfplab/users/jm4454/plink/plink"
plink2="/rigel/mfplab/users/jm4454/plink/plink2"

for i in $(seq 1 22);
do
$plink \
--bfile data/ukb_merged/merged \
--chr $i \
--keep data/theory/fst_{}.txt \
--make-bed \
--out theor_herit/fst_{}_$i

rm theor_herit/fst_{}_$i.log
rm theor_herit/fst_{}_$i.nosex
done

Rscript transethnic_scripts/ldsc_herit_fst_{}.R
rm transethnic_scripts/ldsc_herit_fst_{}.R

for i in $(seq 1 22);
do
rm theor_herit/fst_{}_$i.b*
rm theor_herit/fst_{}_$i.f*
rm theor_herit/fst_{}_$i.r*
rm theor_herit/tmp-data/fst_{}_corr_chr$i.rds
done

rm transethnic_scripts/ldsc_herit_fst_{}.sh



"""

def main1(group):
    return rscript.format(group)

def main2(group):
    return variable.format(group,group,group,group,group,group,group,group,group,group,group,group)

if __name__ == "__main__":
    groups = [i for i in range(1,41)]
    for group in groups:
        with open(f'transethnic_scripts/ldsc_herit_fst_{group}.R','w') as f:
            f.write(main1(group))
        with open(f'transethnic_scripts/ldsc_herit_fst_{group}.sh', 'w') as f:
            f.write(main2(group))
        os.system(f'sbatch transethnic_scripts/ldsc_herit_fst_{group}.sh')

