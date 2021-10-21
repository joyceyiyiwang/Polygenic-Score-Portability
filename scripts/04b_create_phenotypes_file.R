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

NWB_df <- read_delim('data/ukb_populations/NWB_all.txt', delim = ' ', trim_ws = T)

eur_train <- read_delim('data/ukb_populations/WB_train.txt', delim = ' ', trim_ws = T)
eur_test <- read_delim('data/ukb_populations/WB_test.txt', delim = ' ', trim_ws = T)

combined_labels <- bind_rows(
  eur_train %>% mutate(pop = 'WB_train'),
  eur_test %>% mutate(pop = 'WB_test'),
  NWB_df %>% mutate(pop = 'NWB')
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
    eur_mean <- pheno_df[,c(i,"pop")] %>% filter(pop=="WB_train") %>% pull(i) %>% mean(na.rm=T)
    eur_sd <- pheno_df[,c(i,"pop")] %>% filter(pop=="WB_train") %>% pull(i) %>% sd(na.rm=T)

    pheno_df[,i] <- (pheno_df[,i]-eur_mean)/eur_sd
}

pheno_df %>% select('#FID', IID, all_of(traits)) %>%
	write_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', na = 'NA')
