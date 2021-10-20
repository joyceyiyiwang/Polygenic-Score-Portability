
library(tidyverse)

#Import Covariates

# For the phenotypes, the first column is blank
sex_df <- read_table('/work2/06568/joyce_w/stampede2/pgs_portability/data/extracted_phenotypes/ukb.sex.txt',
col_names=T)
sex_df <- sex_df %>% separate(colnames(sex_df)[1], into = c("IID", "is_male"), sep = "\t")
sex_df$IID <- as.numeric(sex_df$IID)
sex_df$is_male <- as.numeric(sex_df$is_male)
age_df <- read_table('/work2/06568/joyce_w/stampede2/pgs_portability/data/extracted_phenotypes/ukb.age.txt',
col_names=T)
age_df <- age_df %>% separate(colnames(age_df)[1], into = c("IID", "age"), sep = "\t")
age_df$IID <- as.numeric(age_df$IID)
age_df$age <- as.numeric(age_df$age)
#pc_df <- read_tsv('data/ukb_merged/projection.sscore')
pc_df <- read_table('/work2/06568/joyce_w/stampede2/pgs_portability/data/extracted_phenotypes/ukb.pc.txt',
                     col_names=T)
cn = c()
for(i in 1:40){
  cn[i] = paste0("PC", i)
}
pc_df <- pc_df %>% separate(colnames(pc_df)[1], into = c("IID", cn), sep = "\t")
pc_df <- mutate_all(pc_df, function(x) as.numeric(x))
unrelated_df <- read_table('/work2/06568/joyce_w/stampede2/pgs_portability/data/extracted_phenotypes/ukb.unrelated.txt',
                     col_names=T)
unrelated_df <- unrelated_df %>% separate(colnames(unrelated_df)[1], into = c("IID", "unrelated"), sep = "\t")
unrelated_df$IID <- as.numeric(unrelated_df$IID)
unrelated_df$unrelated <- as.numeric(unrelated_df$unrelated)
unrelated_df <- subset(unrelated_df, unrelated == 1)
id <- unrelated_df$IID
pc_df <- subset(pc_df, IID %in% id)

df <- sex_df %>% full_join(age_df,by="IID") %>% full_join(pc_df,by="IID")
df <- df %>%
	filter(!is.na(age)) %>%
# Code sex as 0 = missing, 1 = male, 2 = female, as in plink .sample files
	mutate(sex_covar=(2-is_male)) %>%
	mutate(sex_covar=replace(sex_covar,is.na(sex_covar),0)) %>%
	mutate(age_sq = age^2) %>%
	mutate(age_sex = age*sex_covar) %>%
	mutate(age_sq_sex = age_sq * sex_covar) %>%
	select(IID,sex_covar,age,age_sq,age_sex,age_sq_sex,contains("PC")) %>%
#  na.omit()
  subset(IID %in% id)


# Import population files

#AFR_df <- read_delim('data/ukb_populations/AFR_all.txt', delim = ' ', trim_ws = T)
#AMR_df <- read_delim('data/ukb_populations/AMR_all.txt', delim = ' ', trim_ws = T)
#EAS_df <- read_delim('data/ukb_populations/EAS_all.txt', delim = ' ', trim_ws = T)
#SAS_df <- read_delim('data/ukb_populations/SAS_all.txt', delim = ' ', trim_ws = T)
CHB_df <- read_delim('data/ukb_populations/CHB_all.txt', delim = ' ', trim_ws = T)
JPT_df <- read_delim('data/ukb_populations/JPT_all.txt', delim = ' ', trim_ws = T)
CHS_df <- read_delim('data/ukb_populations/CHS_all.txt', delim = ' ', trim_ws = T)
CDX_df <- read_delim('data/ukb_populations/CDX_all.txt', delim = ' ', trim_ws = T)
KHV_df <- read_delim('data/ukb_populations/KHV_all.txt', delim = ' ', trim_ws = T)
#CHD_df <- read_delim('data/ukb_populations/CHD_all.txt', delim = ' ', trim_ws = T)
TSI_df <- read_delim('data/ukb_populations/TSI_all.txt', delim = ' ', trim_ws = T)
GBR_df <- read_delim('data/ukb_populations/GBR_all.txt', delim = ' ', trim_ws = T)
FIN_df <- read_delim('data/ukb_populations/FIN_all.txt', delim = ' ', trim_ws = T)
IBS_df <- read_delim('data/ukb_populations/IBS_all.txt', delim = ' ', trim_ws = T)
YRI_df <- read_delim('data/ukb_populations/YRI_all.txt', delim = ' ', trim_ws = T)
LWK_df <- read_delim('data/ukb_populations/LWK_all.txt', delim = ' ', trim_ws = T)
GWD_df <- read_delim('data/ukb_populations/GWD_all.txt', delim = ' ', trim_ws = T)
MSL_df <- read_delim('data/ukb_populations/MSL_all.txt', delim = ' ', trim_ws = T)
ESN_df <- read_delim('data/ukb_populations/ESN_all.txt', delim = ' ', trim_ws = T)
ASW_df <- read_delim('data/ukb_populations/ASW_all.txt', delim = ' ', trim_ws = T)
ACB_df <- read_delim('data/ukb_populations/ACB_all.txt', delim = ' ', trim_ws = T)
MXL_df <- read_delim('data/ukb_populations/MXL_all.txt', delim = ' ', trim_ws = T)
PUR_df <- read_delim('data/ukb_populations/PUR_all.txt', delim = ' ', trim_ws = T)
CLM_df <- read_delim('data/ukb_populations/CLM_all.txt', delim = ' ', trim_ws = T)
PEL_df <- read_delim('data/ukb_populations/PEL_all.txt', delim = ' ', trim_ws = T)
GIH_df <- read_delim('data/ukb_populations/GIH_all.txt', delim = ' ', trim_ws = T)
PJL_df <- read_delim('data/ukb_populations/PJL_all.txt', delim = ' ', trim_ws = T)
BEB_df <- read_delim('data/ukb_populations/BEB_all.txt', delim = ' ', trim_ws = T)
STU_df <- read_delim('data/ukb_populations/STU_all.txt', delim = ' ', trim_ws = T)
ITU_df <- read_delim('data/ukb_populations/ITU_all.txt', delim = ' ', trim_ws = T)

eur_train <- read_delim('data/ukb_populations/CEU_all.txt', delim = ' ', trim_ws = T)
eur_test <- read_delim('data/ukb_populations/CEU_test.txt', delim = ' ', trim_ws = T)

combined_labels <- bind_rows(
  eur_train %>% mutate(pop = 'CEU_train'),
  eur_test %>% mutate(pop = 'CEU_test'),
  #AFR_df %>% mutate(pop = 'AFR'),
  #AMR_df %>% mutate(pop = 'AMR'),
  #EAS_df %>% mutate(pop = 'EAS'),
  #SAS_df %>% mutate(pop = 'SAS')
  CHB_df %>% mutate(pop = 'CHB'),
  JPT_df %>% mutate(pop = 'JPT'),
  CHS_df %>% mutate(pop = 'CHS'),
  CDX_df %>% mutate(pop = 'CDX'),
  KHV_df %>% mutate(pop = 'KHV'),
  #CHD_df %>% mutate(pop = 'CHD'),
  TSI_df %>% mutate(pop = 'TSI'),
  GBR_df %>% mutate(pop = 'GBR'),
  FIN_df %>% mutate(pop = 'FIN'),
  IBS_df %>% mutate(pop = 'IBS'),
  YRI_df %>% mutate(pop = 'YRI'),
  LWK_df %>% mutate(pop = 'LWK'),
  GWD_df %>% mutate(pop = 'GWD'),
  MSL_df %>% mutate(pop = 'MSL'),
  ESN_df %>% mutate(pop = 'ESN'),
  ASW_df %>% mutate(pop = 'ASW'),
  ACB_df %>% mutate(pop = 'ACB'),
  MXL_df %>% mutate(pop = 'MXL'),
  PUR_df %>% mutate(pop = 'PUR'),
  CLM_df %>% mutate(pop = 'CLM'),
  PEL_df %>% mutate(pop = 'PEL'),
  GIH_df %>% mutate(pop = 'GIH'),
  PJL_df %>% mutate(pop = 'PJL'),
  BEB_df %>% mutate(pop = 'BEB'),
  STU_df %>% mutate(pop = 'STU'),
  ITU_df %>% mutate(pop = 'ITU')
)
df_pop <- df %>%
  left_join(combined_labels, by = c('IID'))

df_pop <- cbind.data.frame(`#FID` = df_pop$`#FID`, df_pop %>% select(-`#FID`))

# Standard normalize with respect to training EUR population

covar <- df_pop %>% select(-`#FID`,-IID,-pop) %>% colnames()

train_eur_df <- df_pop %>% filter(pop=="CEU_train")

for (i in covar){
eur_mean <- train_eur_df %>% pull(i) %>% mean(na.rm=T)
eur_sd <- train_eur_df %>% pull(i) %>% sd(na.rm=T)
df_pop[,i] <- (df_pop[,i]-eur_mean)/eur_sd
}

# Export file

df_pop %>% write_tsv('data/ukb_merged/covar_all_samples.covar')

