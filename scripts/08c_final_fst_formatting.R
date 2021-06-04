library(tidyverse)

fst_df <- data.frame()

for (file in list.files(path = 'data/fst',
                          pattern = 'fst[0-9]+.est', full.names = T)) {
    
        fst_df_append <- read_tsv(file, col_names="Fst",col_types="c")
        fst_df_append$Fst <- fst_df_append$Fst %>% parse_number()
        fst_df_append$label <- rep(c("IID","Mean_Fst","Weighted_Fst"),nrow(fst_df_append)/3)
        fst_df <- rbind(fst_df,fst_df_append)
        rm(fst_df_append)
    
}                      
   
fst_df <- fst_df %>%
    pivot_wider(names_from=label,values_from=Fst)%>%
    arrange(IID)
                       
fst %>% write_tsv('data/fst/final_fst.tsv')
fst %>% write_tsv('data/phenotypes/final_fst_no_filter.tsv')


