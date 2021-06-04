library(tidyverse)

fam_file <- read_delim('data/ukb_merged/merged.fam',
                      col_names = c('#FID', 'IID', 'V1', 'V2', 'V3', 'V4'),
                       delim = ' ', trim_ws = T)

AFR_df <- read_delim('data/ukb_populations/AFR_all.txt', delim = ' ', trim_ws = T)
AMR_df <- read_delim('data/ukb_populations/AMR_all.txt', delim = ' ', trim_ws = T)
EAS_df <- read_delim('data/ukb_populations/EAS_all.txt', delim = ' ', trim_ws = T)
SAS_df <- read_delim('data/ukb_populations/SAS_all.txt', delim = ' ', trim_ws = T)

eur_train <- read_delim('data/ukb_populations/EUR_all.txt', delim = ' ', trim_ws = T)
eur_test <- read_delim('data/ukb_populations/EUR_test.txt', delim = ' ', trim_ws = T)
eur_train %>% nrow %>% print

combined_labels <- bind_rows(
  eur_train %>% mutate(pop = 'train'),
  eur_test %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  AFR_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  AMR_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  EAS_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  SAS_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
)

combined <- fam_file %>%
  left_join(combined_labels, by = c('#FID', 'IID'))

# Print a couple QC checks
#combined %>%
#  summarize_all(function(col) sum(is.na(col))) %>%
#  print


#combined %>%
#  group_by(pop) %>%
#  tally %>%
#  arrange(desc(n)) %>%
#  head %>%
#  print

combined %>%
  select(-'#FID') %>%
  select('#FID' = pop, 'IID', starts_with('V')) %>%
  write_tsv('data/ukb_merged/merged_fst.fam')
