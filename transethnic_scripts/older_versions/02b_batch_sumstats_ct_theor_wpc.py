import os 

rscript = """
library(tidyverse)
pheno <- "{}"
group <- {}

sumstats_file <- paste0("data/gwas_results1/",pheno,"_combined.glm.linear")
sumstats <- read_tsv(sumstats_file, col_types=cols()) %>% select(ID,BETA)

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

for(i in 1:4){{
  
snp_txt <- paste0("data/prs1/",pheno,"_threshold_",i,".txt")
select_snps <- read_tsv(snp_txt,col_names="rsid")
sumstats_sub <- sumstats %>% filter(ID %in% select_snps$rsid)

for(snp in snps){{
#for(snp in maf$SNP){{
  maf1 <- maf_0 %>% filter(SNP==snp) %>% pull(MAF)
  maf2 <- maf %>% filter(SNP==snp) %>% pull(MAF)

  beta <- sumstats_sub %>% filter(ID == snp) %>% pull(BETA)
  beta <- ifelse(length(beta)==0,0,beta)
  num <- sum(num, ((maf1)*(1-maf1)*beta^2),na.rm=T)
  den <- sum(den,((maf2)*(1-maf2)*beta^2),na.rm=T)
}}

beta_maf <- num/den

print(paste0("WPC sumstats,",group,",",pheno,",",i,",",beta_maf))
}}


"""

variable = """#!/bin/bash
#SBATCH --account=mfplab
#SBATCH --job-name=ss_wpc_ct_{}_{}
#SBATCH -c 1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=12gb
set -e

module load anaconda
module load R
source activate prs1

Rscript transethnic_scripts/ss_wpc_ct_{}_{}.R
rm transethnic_scripts/ss_wpc_ct_{}_{}.R
rm transethnic_scripts/ss_wpc_ct_{}_{}.sh
"""

def main1(phenotype,group):
    return rscript.format(phenotype,group)

def main2(phenotype,group):
    return variable.format(phenotype,group,phenotype, group, phenotype,group, phenotype,group)

if __name__ == "__main__":
    phenotypes = ["Height","Platelet","MCV"]
    groups = [i for i in range(2,41)]
    for phenotype in phenotypes:
        for group in groups:
            with open(f'transethnic_scripts/ss_wpc_ct_{phenotype}_{group}.R','w') as f:
                f.write(main1(phenotype,group))
            with open(f'transethnic_scripts/ss_wpc_ct_{phenotype}_{group}.sh', 'w') as f:
                f.write(main2(phenotype,group))
            os.system(f'sbatch transethnic_scripts/ss_wpc_ct_{phenotype}_{group}.sh')

