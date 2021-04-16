import os 

rscript = """
library(tidyverse)

group <- {}
beta_maf <- c()
pheno <- "Height"


sumstats_file <- paste0("data/LDpred1/",pheno,"_final_LDpred.tsv")
sumstats <- read_tsv(sumstats_file, col_types=cols()) %>% select(rsid,LDpred)

file <- "data/theory/wpc_1.frq"
maf_0 <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>% select(SNP,MAF)

num <- 0
den <- 0

file <- paste0("data/theory/wpc_",group,".frq")
maf <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>%
    select(SNP,MAF)

ld_0 <- data.frame()
for (chr in 1:22){{
  file <- paste0("data/theory/wpc_1_",chr,".ld")
  delim <- read_delim(file, delim = ' ', trim_ws = T, col_types = cols())
  delim <- delim %>% select(-X8) %>% mutate(R2 = R^2) %>%
    select(SNP_A,SNP_B,R,R2)
  ld_0 <- bind_rows(ld_0,delim)
}}

causal_0 <- ld_0 %>% filter(R2 >= 0.45)
causal_snps <- unique(c(causal_0$SNP_A,causal_0$SNP_B))
snps <- setdiff(maf$SNP,causal_snps)
rm(causal_0)
  
for(snp in snps){{  
  maf1 <- maf_0 %>% filter(SNP==snp) %>% pull(MAF)
  maf2 <- maf %>% filter(SNP==snp) %>% pull(MAF)
  
  beta <- sumstats %>% filter(rsid == snp) %>% pull(LDpred) 
  num <- sum(num, ((maf1)*(1-maf1)*beta^2),na.rm=T)
  den <- sum(den,((maf2)*(1-maf2)*beta^2),na.rm=T)
}}

beta_maf <- c(beta_maf,num/den)
print(paste0("WPC sumstats,",group,",",beta_maf)
print(group)
print(num)
print(den)
print(beta_maf)



"""

variable = """#!/bin/bash
#SBATCH --account=mfplab
#SBATCH --job-name=ss_pc_{}
#SBATCH -c 1
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=12gb
set -e

module load anaconda
module load R
source activate prs1

Rscript transethnic_scripts/ss_pc_{}.R
rm transethnic_scripts/ss_pc_{}.R
rm transethnic_scripts/ss_pc_{}.sh
"""

def main1(group):
    return rscript.format(group)

def main2(group):
    return variable.format(group,group,group, group)

if __name__ == "__main__":
    groups = [i for i in range(2,41)]
    for group in groups:
        with open(f'transethnic_scripts/ss_pc_{group}.R','w') as f:
            f.write(main1(group))
        with open(f'transethnic_scripts/ss_pc_{group}.sh', 'w') as f:
            f.write(main2(group))
        os.system(f'sbatch transethnic_scripts/ss_pc_{group}.sh')

