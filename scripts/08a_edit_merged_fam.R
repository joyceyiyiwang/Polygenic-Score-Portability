library(tidyverse)

fam_file <- read_delim('data/ukb_merged/merged.fam',
                      col_names = c('#FID', 'IID', 'V1', 'V2', 'V3', 'V4'),
                       delim = ' ', trim_ws = T)

NWB_df <- read_delim('data/ukb_populations/NWB_all.txt', delim = ' ', trim_ws = T)

eur_train <- read_delim('data/ukb_populations/WB_train.txt', delim = ' ', trim_ws = T)
eur_test <- read_delim('data/ukb_populations/WB_test.txt', delim = ' ', trim_ws = T)
eur_train %>% nrow %>% print

combined_labels <- bind_rows(
  eur_train %>% mutate(pop = 'train'),
  eur_test %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  NWB_df %>% mutate(pop = str_glue('test{IID}') %>% as.character)
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
