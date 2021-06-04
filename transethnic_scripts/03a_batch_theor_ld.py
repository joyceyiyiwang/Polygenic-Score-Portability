import os 

rscript = """
library(tidyverse)

gen_group <- "{}"
group <- {}
print(paste0(gen_group," ld"))
print(group)

ld_0 <- data.frame()
for (chr in 1:22){{
  file <- paste0("data/theory/",gen_group,"_1_",chr,".ld")
  delim <- read_delim(file, delim = ' ', trim_ws = T, col_types = cols())
  delim <- delim %>% select(-X8) %>% mutate(R2 = R^2) %>%
    select(SNP_A,SNP_B,R,R2)
  ld_0 <- bind_rows(ld_0,delim)
}}

causal_0 <- ld_0 %>% filter(R2 >= 0.45)
file <- paste0("data/theory/",gen_group,"_1.frq")
maf_0 <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>% select(SNP,MAF)

rel_r2 <- c()

ld <- data.frame()
  for (chr in 1:22){{
    file <- paste0("data/theory/",gen_group,"_",group,"_",chr,".ld")
    delim <- read_delim(file, delim = ' ', trim_ws = T, col_types = cols())
    delim <- delim %>% select(-X8) %>% mutate(R2 = R^2) %>%
      select(SNP_A,SNP_B,R,R2)
    ld <- bind_rows(ld,delim)
}}

file <- paste0("data/theory/",gen_group,"_",group,".frq")
maf <- read_delim(file, delim = ' ', trim_ws=T,col_types=cols()) %>%
    select(SNP,MAF)


causal_var <- unique(c(causal_0$SNP_A,causal_0$SNP_B))
causal_ct <- length(causal_var)
snps <- setdiff(maf$SNP,causal_var)

num <- 0
den <- 0

for(snp in snps){{
#for(snp in maf$SNP){{  
  
  snp1 <- ld_0 %>% filter(SNP_A == snp | SNP_B==snp) %>%
 	filter(SNP_A %in% causal_var | SNP_B %in% causal_var)
  snp2 <- ld %>% filter(SNP_A == snp | SNP_B == snp) %>%
	 filter(SNP_A %in% causal_var | SNP_B %in% causal_var)
  
  comb <- snp1 %>% full_join(snp2,by=c("SNP_A","SNP_B"))
  avgprod_r <- ifelse(is.na(sum(comb$R.x * comb$R.y,na.rm=T)),0,
                      sum(comb$R.x * comb$R.y))/nrow(comb)
  
  maf1 <- maf_0 %>% filter(SNP==snp) %>% pull(MAF)
  maf2 <- maf %>% filter(SNP==snp) %>% pull(MAF)

  num <- sum(num,(avgprod_r * sqrt((maf2*(1-maf2))/(maf1*(1-maf1)))),na.rm=T)  
  den <- sum(den,
           ((snp1 %>% pull(R2) %>% sum(na.rm=T))/nrow(snp1)),
            na.rm=T) 
}}
rel_r2 <- c(rel_r2,(num/den)^2)
print(paste0(gen_group," LD,",group,",",rel_r2))
print(group)
print(num)
print(den)
print(rel_r2)

"""

variable = """#!/bin/bash
#SBATCH --account=mfplab
#SBATCH --job-name=ld_{}_{}
#SBATCH -c 1
#SBATCH --time=01-00:00:00
#SBATCH --mem-per-cpu=12gb
#SBATCH --output=R-%x.out
set -e

module load anaconda
module load R
source activate prs1

gen_group='{}'
group='{}'

Rscript transethnic_scripts/${{gen_group}}_${{group}}.R
rm transethnic_scripts/${{gen_group}}_${{group}}.R
rm transethnic_scripts/${{gen_group}}_${{group}}.sh
"""

def main1(gen_group,group):
    return rscript.format(gen_group,group)

def main2(gen_group, group):
    return variable.format(gen_group,group,
                           gen_group,group)

if __name__ == "__main__":
    groups = [i for i in range(1,41)]
    gen_groups = ["fst","wpc"]
    for gen_group in gen_groups:
        for group in groups:
            with open(f'transethnic_scripts/{gen_group}_{group}.R','w') as f:
                f.write(main1(gen_group,group))
            with open(f'transethnic_scripts/{gen_group}_{group}.sh', 'w') as f:
                f.write(main2(gen_group,group))
            os.system(f'sbatch transethnic_scripts/{gen_group}_{group}.sh')

