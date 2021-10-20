library(tidyverse)

# Table of each trait's name, code, field ID, and dtype (including eid)
traits_info <- read_delim('/work2/06568/joyce_w/stampede2/pgs_portability/data/martin_gwas_info.txt', delim = ' ',
                          col_types = cols_only(Trait = col_character(),
                                                UKBB_code = col_character())) %>%
    mutate(
        ukb_field = str_glue('{UKBB_code}-0.0') %>% as.character,
        dtype = 'd'
    ) %>%
    add_row(Trait = 'eid', UKBB_code = 'eid', ukb_field = 'eid', dtype = 'c')

# Load only the 18 columns of interest
field_to_dtype <- traits_info %>%
    select(ukb_field, dtype) %>%
    deframe %>%
    as.list

raw_phenotypes_df <- read_csv('/work2/06568/joyce_w/stampede2/pgs_portability/file-handlers/ukb45020.csv',
                              col_types = do.call(cols_only, field_to_dtype))

# Samples table (because not assured that FID == IID)
psam_df <- read.table(
  'data/ukb_merged/merged.fam'
) %>% mutate(V1=as.character(V1)) %>% mutate(V2=as.character(V2))
colnames(psam_df)[1:2] <- c("#FID","IID")

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
combined_labels$IID <- as.character(combined_labels$IID)
combined_labels$`#FID` <- as.character(combined_labels$`#FID`)

pheno_df <- raw_phenotypes_df %>%
    # Rename fields to their human-readable names
    pivot_longer(-eid, names_to = 'ukb_field') %>%
    inner_join(traits_info %>% select(-dtype, -UKBB_code), by = 'ukb_field') %>%
    pivot_wider(id_cols = eid, names_from = Trait, values_from = value) %>%
    rename(IID = eid) %>%
    left_join(psam_df, by = 'IID')  %>%
    left_join(combined_labels,by=c('IID','#FID'))

traits <- c("BMI","WBC","Height","RBC","MCV","MCH","Lymphocyte","Platelet","Monocyte","Eosinophil")

for (i in traits){
    eur_mean <- pheno_df[,c(i,"pop")] %>% filter(pop=="CEU_train") %>% pull(i) %>% mean(na.rm=T)
    eur_sd <- pheno_df[,c(i,"pop")] %>% filter(pop=="CEU_train") %>% pull(i) %>% sd(na.rm=T)

    pheno_df[,i] <- (pheno_df[,i]-eur_mean)/eur_sd
}

pheno_df %>% select('#FID', IID, all_of(traits)) %>%
	write_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', na = 'NA')
