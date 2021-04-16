import os 

rscript = """
library(tidyverse)

group <- {}

ld_0 <- data.frame()
for (chr in 1:22){{
  file <- paste0("data/theory/fst_1_",chr,".ld")
  delim <- read_delim(file, delim = ' ', trim_ws = T, col_types = cols())
  delim <- delim %>% select(-X8) %>% mutate(R2 = R^2) %>%
    select(SNP_A,SNP_B,R,R2)
  ld_0 <- bind_rows(ld_0,delim)
}}

causal_0 <- ld_0 %>% filter(R2 >= 0.25)
file <- "data/theory/fst_1.frq"
maf_0 <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>% select(SNP,MAF)

rel_r2 <- c()

ld <- data.frame()
  for (chr in 1:22){{
    file <- paste0("data/theory/fst_",group,"_",chr,".ld")
    delim <- read_delim(file, delim = ' ', trim_ws = T, col_types = cols())
    delim <- delim %>% select(-X8) %>% mutate(R2 = R^2) %>%
      select(SNP_A,SNP_B,R,R2)
    ld <- bind_rows(ld,delim)
}}

file <- paste0("data/theory/fst_",group,".frq")
maf <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>%
    select(SNP,MAF)


causal_var <- unique(c(causal_0$SNP_A,causal_0$SNP_B))
causal_ct <- length(causal_var)
snps <- setdiff(maf$SNP,causal_var)

num <- 0
den <- 0

for (snp in snps){{
  
  snp1 <- ld_0 %>% filter(SNP_A == snp | SNP_B==snp) %>%
 	filter(SNP_A %in% causal_var | SNP_B %in% causal_var)
  snp2 <- ld %>% filter(SNP_A == snp | SNP_B == snp) %>%
	 filter(SNP_A %in% causal_var | SNP_B %in% causal_var)
  
  comb <- snp1 %>% full_join(snp2,by=c("SNP_A","SNP_B"))
  avgprod_r <- ifelse(is.na(sum(comb$R.x * comb$R.y,na.rm=T)),0,
                      sum(comb$R.x * comb$R.y))/causal_ct
  
  maf1 <- maf_0 %>% filter(SNP==snp) %>% pull(MAF)
  maf2 <- maf %>% filter(SNP==snp) %>% pull(MAF)

  num <- sum(num,(avgprod_r * sqrt((maf2*(1-maf2))/(maf1*(1-maf1)))),na.rm=T)  
  den <- sum(den,((snp1 %>% pull(R2) %>% sum(na.rm=T))/causal_ct),na.rm=T) 
}}
rel_r2 <- c(rel_r2,(num/den)^2)
print("Fst LD")
print(group)
print(num)
print(den)
print(rel_r2)

"""

variable = """#!/bin/bash
#SBATCH --account=mfplab
#SBATCH --job-name=fst_{}
#SBATCH -c 1
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=12gb
set -e

module load anaconda
module load R
source activate prs1

Rscript transethnic_scripts/fst_{}.R
rm transethnic_scripts/fst_{}.R
rm transethnic_scripts/fst_{}.sh
"""

def main1(group):
    return rscript.format(group)

def main2(group):
    return variable.format(group,group,group,group)

if __name__ == "__main__":
    groups = [i for i in range(2,41)]
    for group in groups:
        with open(f'transethnic_scripts/fst_{group}.R','w') as f:
            f.write(main1(group))
        with open(f'transethnic_scripts/fst_{group}.sh', 'w') as f:
            f.write(main2(group))
        os.system(f'sbatch transethnic_scripts/fst_{group}.sh')

