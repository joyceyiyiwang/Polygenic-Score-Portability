import os 

rscript = """
library(tidyverse)

gen_group <- "{}"
pheno <- "{}"
group <- {}
beta_maf <- c()

print(paste0(gen_group," sumstats"))
print(group)
print(pheno)

sumstats_file <- paste0("data/LDpred1/",pheno,"_final_LDpred.tsv")
sumstats <- read_tsv(sumstats_file, col_types=cols()) %>% select(rsid,LDpred)

file <- paste0("data/theory/",gen_group,"_1.frq")
maf_0 <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>% select(SNP,MAF)

num <- 0
den <- 0

file <- paste0("data/theory/", gen_group,"_",group,".frq")
maf <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>%
    select(SNP,MAF)

ld_0 <- data.frame()
for (chr in 1:22){{
  file <- paste0("data/theory/", gen_group,"_1_",chr,".ld")
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
#for(snp in maf$SNP){{    
  maf1 <- maf_0 %>% filter(SNP==snp) %>% pull(MAF)
  maf2 <- maf %>% filter(SNP==snp) %>% pull(MAF)
  
  beta <- sumstats %>% filter(rsid == snp) %>% pull(LDpred) 
  num <- sum(num, ((maf1)*(1-maf1)*beta^2),na.rm=T)
  den <- sum(den,((maf2)*(1-maf2)*beta^2),na.rm=T)
}}

beta_maf <- c(beta_maf,num/den)
print(paste0(gen_group," sumstats,",group,",",beta_maf,",",pheno))
print(group)
print(num)
print(den)
print(beta_maf)

"""

variable = """#!/bin/bash
#SBATCH --account=mfplab
#SBATCH --job-name=ss_ldpred_{}_{}_{}
#SBATCH -c 1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=12gb
#SBATCH --output=R-%x.out
set -e

module load anaconda
module load R
source activate prs1

gen_group='{}'
phenotype='{}'
group='{}'

Rscript transethnic_scripts/ss_${{gen_group}}_${{phenotype}}_${{group}}.R
rm transethnic_scripts/ss_${{gen_group}}_${{phenotype}}_${{group}}.R
rm transethnic_scripts/ss_${{gen_group}}_${{phenotype}}_${{group}}.sh
"""

def main1(gen_group, phenotype, group):
    return rscript.format(gen_group, phenotype,group)

def main2(gen_group, phenotype, group):
    return variable.format(gen_group,phenotype,group,
                           gen_group, phenotype, group)

if __name__ == "__main__":
    groups = [i for i in range(2,41)]
    gen_groups = ["fst","wpc"]
    phenotypes = ["Platelet","Height","MCV"]
    for gen_group in gen_groups:
        for phenotype in phenotypes:
            for group in gen_groups:
                with open(f'transethnic_scripts/ss_{gen_group}_{phenotype}_{group}.R','w') as f:
                    f.write(main1(gen_group,phenotype, group))
                with open(f'transethnic_scripts/ss_{gen_group}_{phenotype}_{group}.sh', 'w') as f:
                    f.write(main2(gen_group, phenotype,group))
                    os.system(f'sbatch transethnic_scripts/ss_{gen_group}_{phenotype}_{group}.sh')

