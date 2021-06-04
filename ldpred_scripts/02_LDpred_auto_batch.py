import os 

rscript = """

.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(bigsnpr)
library(Matrix)
library(tidyverse)

phenotype <- "{}"
tsv_file <- paste0("data/LDpred1/",phenotype,"_merged_sumstats_ldpred.tsv")
sumstats <- read_tsv(tsv_file)

tmp <- tempfile(tmpdir = "data/LDpred1/tmp-data")

print("Calculating Correlation Matrices")

for (chr in 1:22){{
  
  rds_file <- paste0("data/LDpred1/LD_EUR_train_",chr,".rds")
  ukb <- snp_attach(rds_file)
  G <- ukb$genotypes
  CHR <- as.integer(ukb$map$chromosome)
  POS <- ukb$map$physical.pos
  
  
  # Indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  # indices in 'G' / rsid replacing `NUM_ID`
  ind.chr2 <- sumstats$rsid[ind.chr]
  # Indices in 'corr' / ind.chr replacing which(CHR == chr) 
  ind.chr3 <- match(ind.chr2,sumstats$rsid[ind.chr])
  
  corr0 <- readRDS(paste0("data/LDpred1/corr_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1){{
    df_beta <- sumstats[ind.chr, c("chr","beta", "beta_se", "n_eff", "rsid")]
    ld <- colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  }} else {{
    df_beta <- rbind(df_beta, sumstats[ind.chr, c("chr","beta", "beta_se", "n_eff", "rsid")])
    ld <- c(ld, colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }}
  
}}


#### Retrieve initial heritability estimate
print("Heritability estimates")
(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
(ld_h2_est <- ldsc[["h2"]])

h2_table <- tibble(
Trait=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte","WBC","Eosinophil","Lymphocyte",
        "DBP","SBP","Hb","Ht","MCHC","Neutrophil","Basophil"),
h2 = c(0.485,0.308,0.267,0.253,0.248,0.234,0.230,0.191,0.184,0.210,
        0.143,0.151,0.1054,0.0942,0.0532,0.164,0.0242))
(h2_est <- h2_table[h2_table$Trait==phenotype,2] %>% pull(h2))

####### LDpred Auto
print("LDpred - auto beginning")
#LDpred2 - auto
NCORES <- 1
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 20),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

beta_auto_df <- as.data.frame(beta_auto)
final_beta_auto <- rowMeans(beta_auto_df)

sumstats1 <- sumstats %>%
  cbind(beta_auto_df) %>%
  mutate(average_BETA=final_beta_auto) %>%
  select(-a0,-a1,-beta_se,-p,-n_eff)

export_file <- paste0("data/LDpred1/",phenotype,"_ldpred_auto_results.tsv")
sumstats1 %>% write_tsv(export_file)

p_auto <- sapply(multi_auto, function(auto) auto$p_est)
export_file <- paste0("data/LDpred1/",phenotype,"_ldpred_auto_p_est.tsv")
p_auto %>% as.data.frame() %>% write_tsv(export_file)

h2_auto <- sapply(multi_auto, function(auto) auto$h2_est)
export_file <- paste0("data/LDpred1/",phenotype,"_ldpred_auto_h2_est.tsv")
h2_auto %>% as.data.frame() %>% write_tsv(export_file)

"""

variable = """#!/bin/bash
#SBATCH --account=mfplab
#SBATCH --job-name={}
#SBATCH -c 1
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=12gb
set -e

module load anaconda
module load R
source activate prs1

pheno='{}'
Rscript ldpred_scripts/temp_scripts/${pheno}.R
rm ldpred_scripts/temp_scripts/${pheno}.R
rm ldpred_scripts/${pheno}.sh
"""

def main1(phenotype):
    return rscript.format(phenotype)

def main2(phenotype):
    return variable.format(phenotype,phenotype)

if __name__ == "__main__":
    phenotypes = ['BMI','Height' ,'RBC', 'Platelet', 'MCV', 'Monocyte', 'WBC', 'MCH', 'Eosinophil', 'Lymphocyte']
    for phenotype in phenotypes:
        with open(f'ldpred_scripts/temp_scripts/{phenotype}.R','w') as f:
            f.write(main1(phenotype))
        with open(f'ldpred_scripts/{phenotype}.sh', 'w') as f:
            f.write(main2(phenotype))
        os.system(f'sbatch ldpred_scripts/{phenotype}.sh')
