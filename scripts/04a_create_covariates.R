
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
pc_df <- read_table('/work2/06568/joyce_w/stampede2/pgs_portability/data/extracted_phenotypes/ukb.pc.txt',
                     col_names=T)
cn = c()
for(i in 1:40){
  cn[i] = paste0("PC", i)
}
pc_df <- pc_df %>% separate(colnames(pc_df)[1], into = c("IID", cn), sep = "\t")
pc_df <- mutate_all(pc_df, function(x) as.numeric(x))

df <- sex_df %>% full_join(age_df,by="IID") %>% full_join(pc_df,by="IID")
df <- df %>%
	filter(!is.na(age)) %>%
# Code sex as 0 = missing, 1 = male, 2 = female, as in plink .sample files
	mutate(sex_covar=(2-is_male)) %>%
	mutate(sex_covar=replace(sex_covar,is.na(sex_covar),0)) %>%
	mutate(age_sq = age^2) %>%
	mutate(age_sex = age*sex_covar) %>%
	mutate(age_sq_sex = age_sq * sex_covar) %>%
	select(IID,sex_covar,age,age_sq,age_sex,age_sq_sex,contains("PC"))


# Import population files

NWB_df <- read_delim('data/ukb_populations/NWB_all.txt', delim = ' ', trim_ws = T)

eur_train <- read_delim('data/ukb_populations/WB_train.txt', delim = ' ', trim_ws = T)
eur_test <- read_delim('data/ukb_populations/WB_test.txt', delim = ' ', trim_ws = T)

combined_labels <- bind_rows(
  eur_train %>% mutate(pop = 'WB_train'),
  eur_test %>% mutate(pop = 'WB_test'),
  NWB_df %>% mutate(pop = 'NWB')
)
df_pop <- df %>%
  left_join(combined_labels, by = c('IID'))

df_pop <- cbind.data.frame(`#FID` = df_pop$`#FID`, df_pop %>% select(-`#FID`))

# Standard normalize with respect to training EUR population

covar <- df_pop %>% select(-`#FID`,-IID,-pop) %>% colnames()

train_eur_df <- df_pop %>% filter(pop=="WB_train")

for (i in covar){
eur_mean <- train_eur_df %>% pull(i) %>% mean(na.rm=T)
eur_sd <- train_eur_df %>% pull(i) %>% sd(na.rm=T)
df_pop[,i] <- (df_pop[,i]-eur_mean)/eur_sd
}

# Export file

df_pop %>% write_tsv('data/ukb_merged/covar_all_samples.covar')

