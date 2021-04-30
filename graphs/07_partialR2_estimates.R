.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
#.libPaths("/moto/palab/users/jm4454/rpackages/")
library(asbio)
library(cowplot)
library(ggrepel)
library(ggplot2)
library(ggpmisc)
library(plyr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(bigsnpr)
library(Gmedian)
library(philentropy)

{

load_non_prs_df <- function(num_pc_groups,num_fst_groups, n_buckets=60) {
  # Loads all the covariates, phenotypes, and outcomes into a dataframe
  
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar_all_samples.covar')
  
  # Unweighted Medians
  medians <- covar_df %>% filter(pop=="EUR_train") %>% 
    select(starts_with("PC")) %>% Gmedian() %>% as.vector()
  covar_df <- covar_df %>% filter(pop!="EUR_train")
  
  dist_fun <- function(row){
    if (is.na(row)>0){
      return(NA)
    } 
    else{ 
      mat <- rbind(row,medians)
      return(distance(mat,method="euclidean"))
    }
  }
  
  eigenval <- read.table("data/kgp_merged/merged.eigenval",header=F)
  eigenval <- (eigenval/sum(eigenval)) %>% unlist() %>% as.vector()
  
  dist_fun_weighted <- function(row){
    if (is.na(row)>0){
      return(NA)
    } 
    else{ 
      mat <- sweep(rbind(row,medians),MARGIN=2,sqrt(eigenval), `*`)
      return(distance(mat,method="euclidean"))
      
    }
  }
  
  covar_df$PC_dist <- covar_df %>% select(starts_with("PC")) %>%
    as.matrix() %>%
    apply(1,dist_fun)
  covar_df$PC_dist_weighted <- covar_df %>% select(starts_with("PC")) %>%
    as.matrix() %>%
    apply(1,dist_fun_weighted)
  
  phenotypes_df <- read_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', trim_ws = T)
  
  # Load the individuals who are in the evaluation set for each population
  population_files <- c('data/ukb_populations/AFR_all.txt', 'data/ukb_populations/AMR_all.txt',
                        'data/ukb_populations/EAS_all.txt', 'data/ukb_populations/EUR_test.txt',
                        'data/ukb_populations/SAS_all.txt')
  populations_df <- data.frame()
  for (file in population_files) {
    pop_df <- read_delim(file, delim = ' ', trim_ws = T,
                         col_types = c('#FID' = col_integer(), 'IID' = col_integer())) %>%
      mutate(population = str_extract(file, '(?<=data/ukb_populations/)[A-Z]{3}'))
    populations_df <- bind_rows(populations_df, pop_df)
  }
  fst_values <- read_tsv('data/fst/final_fst.tsv')  
  #fst_values <- read_tsv('data/fst/final_fst.tsv')
  
  # Combine all the above tables into a table that will be joined with PRS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    left_join(fst_values,by=c("IID")) %>% #
    mutate(Weighted_Fst=ifelse(Weighted_Fst<=0,0,Weighted_Fst)) %>%
  #  mutate(Weighted_Fst=ifelse((is.na(Weighted_Fst)) & (population=="EUR_test"),0,Weighted_Fst)) %>%
    pivot_longer(BMI:Eosinophil, names_to = 'phenotype', values_to = 'phenotype_value')
  
  
  
 
  cutpoints_fst <<- seq(min(phenotypes_df$Weighted_Fst,na.rm=T),
                    max(phenotypes_df$Weighted_Fst,na.rm=T),
                    length=num_fst_groups+1)
  cutpoints_pcs <<- seq(min(phenotypes_df$PC_dist, na.rm=T),
                     max(phenotypes_df$PC_dist,na.rm=T),
                     length=num_pc_groups+1)
  cutpoints_wpcs <<- seq(min(phenotypes_df$PC_dist_weighted, na.rm=T),
                         max(phenotypes_df$PC_dist_weighted,na.rm=T),
                         length=num_pc_groups+1)
  
  phenotypes_df <- phenotypes_df %>%

    mutate(weighted_fst_groups = Weighted_Fst %>% ntile(n_buckets)) %>% #Divides into equally-sized components
    mutate(PC_groups = PC_dist %>% ntile(n_buckets)) %>% #Divides into equally-sized components
    mutate(WPC_groups = PC_dist_weighted %>% ntile(n_buckets)) %>% #Divides into equally-sized components
    add_count(PC_groups) %>%
    rename(PC_groups_count=n) %>%
    add_count(WPC_groups) %>%
    rename(WPC_groups_count=n) %>%
    add_count(weighted_fst_groups) %>%
    rename(fst_groups_count=n)
  
  return(phenotypes_df)
}

get_r2_values <- function(df,group_var) {
  # Computes the partial R^2 attributable to the PRS
  #    group_by(population, phenotype, threshold) %>%
  df$group <- df %>% select(any_of(group_var)) %>% as.data.frame()
  df1 <- df %>%
    group_by(group, phenotype, threshold) %>%
    do(
      # Regression on only covariates
      nested = lm(phenotype_value ~ sex_covar + age + age_sq + age_sex + age_sq_sex +
                    PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
                    PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + PC11_AVG + PC12_AVG +
                    PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                    PC19_AVG + PC20_AVG, data = .),
      
      # Regression now including the PRS + covariates
      full = lm(phenotype_value ~ prs + sex_covar + age + age_sq + age_sex + age_sq_sex +
                  PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
                  PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + PC11_AVG + PC12_AVG +
                  PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                  PC19_AVG + PC20_AVG,
                  data = .)
    ) %>%
    mutate(
      partial = partial.R2(nested, full),
      nested_r2 = summary(nested)$adj.r.squared,
      full_r2 = summary(full)$adj.r.squared,
      incremental_r2 = full_r2 - nested_r2
    ) %>%
    select(-nested, -full)
  df1$group <- 1:nrow(df1)
  colnames(df1)[1] <- group_var
  return(df1)
}

make_ldpred_evaluation_df <- function(non_prs_df,group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  prs_df <- data.frame()
  for (file in list.files(path = 'data/LDpred1/prs',
                          pattern = '[a-zA-Z]+_LDpred_scores_[0-9]{1}[0-9]?.sscore', full.names = T)) {
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', prs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/LDpred1/prs/)[A-Za-z]+(?=_)'),
        threshold = str_extract_all(string = file, pattern = '[0-9]+')[[1]][2] %>% as.integer
      ) %>%
      inner_join(non_prs_df, by =  c('#FID', 'IID','phenotype'))%>%
      get_r2_values(.,group_var=group_var_as_string)
    
    prs_df <- bind_rows(prs_df, this_df)
  }
  prs_df
}

make_prs_evaluation_df <- function(non_prs_df,group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  prs_df <- data.frame()
  for (file in list.files(path = 'data/prs1',
                          pattern = '[a-zA-Z]+_[0-9]_scores.sscore', full.names = T)) {
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', prs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/prs1/)[A-Za-z]+(?=_)'),
        threshold = str_extract_all(string = file, pattern = '[0-4]')[[1]][2] %>% as.integer
      ) %>%
      inner_join(non_prs_df, by =  c('#FID', 'IID','phenotype'))%>%
      get_r2_values(.,group_var=group_var_as_string)
    
    prs_df <- bind_rows(prs_df, this_df)
  }
  prs_df
}

# Comparison


####################
#Previously 60
non_prs_df <- load_non_prs_df(num_pc_groups = 8,num_fst_groups = 9, n_buckets=40)
ld_df_fst <- make_ldpred_evaluation_df(non_prs_df,"weighted_fst_groups")
#ld_df_PC <- make_ldpred_evaluation_df(non_prs_df,"PC_groups")
ld_df_WPC <- make_ldpred_evaluation_df(non_prs_df,"WPC_groups")
prs_df_fst  <- make_prs_evaluation_df(non_prs_df,"weighted_fst_groups")
#prs_df_PC <- make_prs_evaluation_df(non_prs_df,"PC_groups")
prs_df_WPC <- make_prs_evaluation_df(non_prs_df,"WPC_groups")

prs_df_pop <- make_prs_evaluation_df(non_prs_df,"population")  

ld_fst_df <- ld_df_fst %>% 
  arrange(phenotype,weighted_fst_groups,-threshold) %>%
  mutate(group_number = weighted_fst_groups) #%>%
#  filter(weighted_fst_groups!=10)

ld_pc_df <- ld_df_PC %>% 
  arrange(phenotype,PC_groups,-threshold) %>%
  mutate(group_number = PC_groups) #%>%
#  filter(PC_groups!=9)

ld_wpc_df <- ld_df_WPC %>% 
  arrange(phenotype,WPC_groups,-threshold) %>%
  mutate(group_number = WPC_groups)  #%>%
#  filter(WPC_groups!=9)  
                                   
prs_fst_df <- prs_df_fst %>% 
  arrange(phenotype,weighted_fst_groups,threshold) %>%
  mutate(group_number = weighted_fst_groups)  #%>%
#  filter(weighted_fst_groups!=10)

prs_pc_df <- prs_df_PC %>% 
  arrange(phenotype,PC_groups,threshold) %>%
  mutate(group_number = PC_groups) # %>%
#  filter(PC_groups!=9)

prs_wpc_df <- prs_df_WPC %>% 
  arrange(phenotype,WPC_groups,threshold) %>%
  mutate(group_number = WPC_groups) #%>%
#  filter(WPC_groups!=9)


prs_pop_df <- prs_df_pop %>%
  mutate(population = factor(population,levels=c(4,2,5,3,1),labels = c("EUR","AMR","SAS","EAS","AFR")))

non_prs_df %>% write_tsv('data/prs_comparisons2/non_prs_df.tsv')

prs_fst_df %>% write_tsv('data/prs_comparisons2/UKBB_geno_fst_CT.tsv')
ld_fst_df %>% write_tsv('data/prs_comparisons2/UKBB_geno_fst_LD.tsv')
#prs_pc_df %>% write_tsv('data/prs_comparisons2/UKBB_geno_pc_CT.tsv')
#ld_pc_df %>% write_tsv('data/prs_comparisons2/UKBB_geno_pc_LD.tsv')
prs_wpc_df %>% write_tsv('data/prs_comparisons2/UKBB_geno_wpc_CT.tsv')
ld_wpc_df %>% write_tsv('data/prs_comparisons2/UKBB_geno_wpc_LD.tsv')
prs_pop_df %>% write_tsv('data/prs_comparisons2/UKBB_geno_pop.tsv')

}

non_prs_df <- read_tsv('data/prs_comparisons2/non_prs_df.tsv')
prs_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_CT.tsv')
ld_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_LD.tsv')
#prs_pc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_pc_CT.tsv')
#ld_pc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_pc_LD.tsv')
prs_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_CT.tsv')
ld_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_LD.tsv')
prs_pop_df <- read_tsv('data/prs_comparisons2/UKBB_geno_pop.tsv')
#######################



groups <- unique(non_prs_df$weighted_fst_groups) %>% sort()
groups1 <- unique(non_prs_df$WPC_groups) %>% sort()

ids <- non_prs_df %>% 
  select("#FID","IID","weighted_fst_groups","WPC_groups") %>%
  distinct()

i <- 1

for (group in groups){
  sub <- ids %>% filter(weighted_fst_groups==group) %>% select(-weighted_fst_groups,-WPC_groups)
  if(nrow(sub) >250){
    sub_name <- paste0("data/theory/fst_",i,".txt")
    sub %>% write.table(sub_name, sep="\t",  col.names=FALSE, row.names = F)
    i <- i + 1
  }
}


i <- 1

for (group in groups1){
  sub <- ids %>% filter(WPC_groups==group) %>% select(-weighted_fst_groups,-WPC_groups)
 if(nrow(sub) >250){
    sub_name <- paste0("data/theory/wpc_",i,".txt")
    sub %>% write.table(sub_name, sep="\t",  col.names=FALSE, row.names = F)
    i <- i + 1
  }
}


#####################

mean_fst_values <- non_prs_df %>% 
  filter(phenotype=="Height") %>%
  select(Weighted_Fst,weighted_fst_groups,fst_groups_count) %>%
  group_by(weighted_fst_groups) %>%
  mutate(mean_fst = median(Weighted_Fst)) %>%
  ungroup() %>%
  distinct(weighted_fst_groups,mean_fst,.keep_all=T) %>%
  select(-Weighted_Fst) %>%
  na.omit() %>%
  arrange(mean_fst) %>%
  mutate(weighted_fst_groups = 1:n()) 

mean_pc_values <- non_prs_df %>% 
  filter(phenotype=="Height") %>%
  select(PC_dist,PC_groups,PC_groups_count) %>%
  group_by(PC_groups) %>%
  mutate(mean_pc = median(PC_dist)) %>%
  ungroup() %>%
  distinct(PC_groups,mean_pc,.keep_all=T) %>%
  select(-PC_dist) %>%
  na.omit() %>%
  mutate(mean_pc = as.numeric(mean_pc)) %>%
  arrange(mean_pc)%>%
  mutate(PC_groups = 1:n())

mean_wpc_values <- non_prs_df %>% 
  filter(phenotype=="Height") %>%
  select(PC_dist_weighted,WPC_groups,WPC_groups_count) %>%
  drop_na() %>%
  group_by(WPC_groups) %>%
  mutate(mean_wpc = median(PC_dist_weighted)) %>%
  ungroup() %>%
  distinct(WPC_groups,mean_wpc,.keep_all=T) %>%
  select(-PC_dist_weighted) %>%
  na.omit() %>%
  mutate(mean_wpc = as.numeric(mean_wpc)) %>%
  arrange(mean_wpc) %>%
  mutate(WPC_groups = 1:n()) 

##########################################


################### Initial sparsity values in LDpred
pheno <- function(phemeanno){
  return(phenotype)
}

df <- data.frame(p = seq_log(1e-4,0.9,20))
for (file in list.files(path = 'data/LDpred',
                        pattern = '[a-zA-Z]+_ldpred_auto_p_est.tsv', full.names = T)) {
  phenotype <- str_extract(string = file, pattern = '(?<=data/LDpred/)[A-Za-z]+(?=_)')
  add <- read_tsv(file,col_names=F,skip=1)
  add <- add %>% rename_with(pheno)
  df <- bind_cols(df,add)
  
}

df <- data.frame(h2 = 1:20)
for (file in list.files(path = 'data/LDpred',
                        pattern = '[a-zA-Z]+_ldpred_auto_h2_est.tsv', full.names = T)) {
  phenotype <- str_extract(string = file, pattern = '(?<=data/LDpred/)[A-Za-z]+(?=_)')
  add <- read_tsv(file,col_names=F,skip=1)
  add <- add %>% rename_with(pheno)
  df <- bind_cols(df,add)
  
}

median_p <- ldply(df,median) 
colnames(median_p) <- c("phenotype","prop")
median_p <- median_p[-1,]
median_p[,2] <- round(median_p[,2],3)




########### Remaining code sorted in other files
################################
all_list <- list(prs_fst_df,ld_fst_df,
                 #prs_pc_df,ld_pc_df,
                 prs_wpc_df,ld_wpc_df)

pheno_factor <- function(df){
  mutate(df,phenotype = factor(phenotype, levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                                   "Lymphocyte","WBC","Eosinophil")))
}
all_list <- lapply(all_list,pheno_factor) 


all_df <- ldply(all_list, rbind)
all_df$model <- c(rep("Fst",3200),rep("WPC",3200))
all_df %>% filter(threshold %in% c(1,2,3,4,29)) %>%
  group_by(model, phenotype) %>%
  mutate(max_h2 = max(partial)) %>%
  select(model,phenotype,threshold,max_h2) %>%
  distinct(model,phenotype,max_h2, .keep_all=T) %>%
  arrange(phenotype,model)


combined <- function(prs_df,ldpred_prs_df, type, line){
  

  ld_df <- ldpred_prs_df %>%
    filter(threshold==29)# %>%
 #   filter(phenotype!="Lymphocyte")
  
  plot_df <- rbind(prs_df,ld_df)
  
  plot_df$threshold <- factor(plot_df$threshold,levels=c(29,4,3,2,1,0),
                              labels=c("0-LDpred2","1e-2","1e-3","1e-4","1e-5","5e-8"),
                              ordered = T)

  plot_df <- plot_df %>%
    filter(threshold!="5e-8") %>%
    #filter(phenotype!="Lymphocyte") %>%
    group_by(phenotype,threshold) %>%
    mutate(
      eur_partial = max((partial*(group_number==1)),na.rm=T),
      relative_performance = partial / eur_partial) %>%
    ungroup() %>%
    arrange(threshold)
  #  filter(group_number!=5)
  
  
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(mean_fst_values,by="weighted_fst_groups") 
    plot_df <- plot_df %>% rename(mean_values = mean_fst)
    title1 = "Phenotype Variation explained by Polygenic Scores Across Fst Pools"
    xlabel1 = "Within-Group Median Fst"
    subtitle1 = "Fst was calculated between each individual against EUR training population."
    caption1 =  "Total of 40 groups each composed of around 750 observations."
    delta <- 0.05
    beta <- "\u03b2 * 0.05 (p)"
  }
  if (type=="PC") {
    plot_df <- plot_df %>% left_join(mean_pc_values,by="PC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_pc)
    title1 = "Phenotype Variation explained by Polygenic Scores PC Distance Pools"
    xlabel1 = "Within-Group Median PC Distance"
    subtitle1 = "PC distance was calculated between each individual's against the geometric median of the EUR training population."
    caption1 = "" #range_pc
  }
  if (type=="WPC") {
    plot_df <- plot_df %>% left_join(mean_wpc_values,by="WPC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_wpc)
    title1 = "Phenotype Variation explained by Polygenic Scores Across Weighted PC Distance Pools"
    xlabel1 = "Within-Group Median PC Distance"
    subtitle1 = "Weighted PC distance is Euclidean distance between each individual's 20 PCs against the geometric median of the EUR training population's PCs."
    caption1 = "Total of 40 groups each composed of around 750 observations."
    delta <- 50
    beta <- "\u03b2 * 50 (p)"
  }

  
  plot <-  ggplot(plot_df, aes(x=mean_values,y=partial, #or partial
                               color=phenotype))
  if(line ==1){
  plot <- plot +
    geom_smooth(size=1,method = lm,se=F,na.rm=T)
    
    tb <- plot_df %>% group_by(phenotype, threshold) %>%
      do(
        lm = lm(partial~mean_values, data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[2,c(1)]*delta,3),
                           " (",
                           round(summary(lm)$coef[2,c(4)],3),
                           ")")) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coef,`LM R2`) %>%
      rename("Phenotype"="phenotype") 
    
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
    h <- -.4
  }

  
  if(line ==5){
    plot <- plot +
      geom_smooth(size=1,method="loess",span=0.75,na.rm=T,se=F) #4
    
    
    tb <- plot_df %>%
      group_by(phenotype, threshold) %>%
      do(
        ls = loess(partial~mean_values, data=.)
      )  %>%
      ungroup() 
    tb$ls_fit <-  apply(sapply(tb$ls,fitted),2,list)
    tb <- tb %>%
      unnest(cols=ls_fit) %>%
      unnest(cols=ls_fit) %>%
      mutate(group_number = rep(1:40,10*5)) %>%
      select(group_number, phenotype, threshold, ls_fit) %>%
      right_join(plot_df , by=c("phenotype","threshold","group_number")) %>%
      group_by(phenotype, threshold) %>%
      do(
        lm = lm(ls_fit~mean_values,data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[2,c(1)]*delta,3),
                           " (",
                           round(summary(lm)$coef[2,c(4)],3),
                           ")")) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coef,`LM R2`) %>%
      rename("Phenotype"='phenotype')
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
    
    h <- -.4
  }
  
  tbs <- lapply(split(tb, tb$threshold), "[", -1)
  df <- tibble(x = rep(-Inf, length(tbs)), 
               y = rep(Inf, length(tbs)), 
               threshold = levels(as.factor(tb$threshold)), 
               tbl = tbs) %>%
    filter(threshold !="5e-8")
  
  prs_names <- c(
    `0-LDpred2` = "LDpred2",
    `1e-2` = "1e-2",
    `1e-3` = "1e-3",
    `1e-4` = "1e-4",
    `1e-5` = "1e-5"
  )
  
   plot <- plot + geom_point(size=1.5,color="white")+
      geom_point(size=1,alpha=0.75)+
      guides(color = guide_legend(title="Trait"),
             size=F)+
      theme_light() + 
     theme(legend.text=element_text(size=14),
           legend.title=element_text(size=16))+
 #     scale_color_brewer(palette="Set1")+
      xlab(xlabel1) +
      ylab(expression(Partial~R^2))+
      labs(
        #title=title1,
        #   subtitle=subtitle1,
           caption=caption1)+
     facet_grid(cols = vars(threshold), labeller = as_labeller(prs_names))+
     geom_table(data = df, aes(x = x, y = y, label = tbl),
                hjust = h, vjust = 1)
     


  return(plot)
}


for (i in c(1,5)){
comb_fst <- combined(as.data.frame(all_list[1]),as.data.frame(all_list[2]), "FST",i)
ggsave(paste0('img1/FINAL_COMB_FST',i,"_40.png"),comb_fst,width=20,height=10,dpi=1000)

comb_wpc <- combined(as.data.frame(all_list[3]),as.data.frame(all_list[4]), "WPC",i)
ggsave(paste0('img1/FINAL_COMB_WPC',i,"_40.png"),comb_wpc,width=20,height=10,dpi=1000)
}

combined_df <- function(prs_df,ldpred_prs_df, type, line){
  
  
  ld_df <- ldpred_prs_df %>%
    filter(threshold==29)
  
  plot_df <- rbind(prs_df,ld_df)
  
  plot_df$threshold <- factor(plot_df$threshold,levels=c(29,4,3,2,1,0),
                              labels=c("0-LDpred2","1e-2","1e-3","1e-4","1e-5","5e-8"),
                              ordered = T)
  
  plot_df <- plot_df %>%
    filter(threshold!="5e-8") %>%
    group_by(phenotype,threshold) %>%
    mutate(
      eur_partial = max((partial*(group_number==1)),na.rm=T),
      relative_performance = partial / eur_partial) %>%
    ungroup() %>%
    arrange(threshold)

  
  
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(mean_fst_values,by="weighted_fst_groups") 
    plot_df <- plot_df %>% rename(mean_values = mean_fst)
    delta <- 0.05
}
  if (type=="WPC") {
    plot_df <- plot_df %>% left_join(mean_wpc_values,by="WPC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_wpc)
    delta <- 50
  }
  
  if(line ==1){
    
    tb <- plot_df %>% group_by(phenotype, threshold) %>%
      do(
        lm = lm(partial~mean_values, data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coefficient = summary(lm)$coef[2,c(1)]*delta) %>%
      mutate(`p-val` = summary(lm)$coef[2,c(4)]) %>%
     mutate(Coef = paste0(round(summary(lm)$coef[2,c(1)],3),
                          " (",round(summary(lm)$coef[2,c(4)],3),
                          ")")) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coefficient,`p-val`,`LM R2`) %>%
      rename("Phenotype"="phenotype")
    h <- -.2
  }
  
  
  if(line ==5){
    
    tb <- plot_df %>%
      group_by(phenotype, threshold) %>%
      do(
        ls = loess(partial~mean_values, data=.)
      )  %>%
      ungroup() 
    tb$ls_fit <-  apply(sapply(tb$ls,fitted),2,list)
    tb <- tb %>%
      unnest(cols=ls_fit) %>%
      unnest(cols=ls_fit) %>%
      mutate(group_number = rep(1:40,10*5)) %>%
      select(group_number, phenotype, threshold, ls_fit) %>%
      right_join(plot_df , by=c("phenotype","threshold","group_number")) %>%
      group_by(phenotype, threshold) %>%
      do(
        lm = lm(ls_fit~mean_values,data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coefficient = summary(lm)$coef[2,1]*delta) %>%
      mutate(`p-val` = summary(lm)$coef[2,4]) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[2,c(1)]*delta,3),
                           " (",
                           round(summary(lm)$coef[2,c(4)],3),
                           ")")) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coefficient,`p-val`,`LM R2`) %>%
      rename("Phenotype"='phenotype')
    

  }
  return(tb)

}


comb <- data.frame()
for (i in c(1,5)){
  comb_fst <- combined_df(as.data.frame(all_list[1]),as.data.frame(all_list[2]), "FST",i)
  comb_wpc <- combined_df(as.data.frame(all_list[3]),as.data.frame(all_list[4]), "WPC",i)
  comb1 <- bind_rows(comb_fst,comb_wpc)
  comb <- rbind(comb,comb1)
  
}
comb$gd <- c(rep("FST",50), rep("WPC", 50),rep("FST", 50),rep("WPC",50))
comb$model <- c(rep("LM",100),rep("LOESS", 100))


combined_rc <- function(prs_df,ldpred_prs_df, type, line){
  
  
  ld_df <- ldpred_prs_df %>%
    filter(threshold==29) #%>%
   # filter(phenotype!="Lymphocyte")
  
  plot_df <- bind_rows(prs_df,ld_df)
  
  plot_df$threshold <- factor(plot_df$threshold,levels=c(29,4,3,2,1,0),
                              labels=c("0-LDpred2","1e-2","1e-3","1e-4","1e-5","5e-8"),
                              ordered = T)
  
  plot_df <- plot_df %>%
    filter(threshold!="5e-8") %>%
  #  filter(phenotype!="Lymphocyte") %>%
    group_by(phenotype,threshold) %>%
    mutate(
      eur_partial = max((partial*(group_number==1)),na.rm=T),
      relative_performance = partial / eur_partial) %>%
    ungroup() %>%
    arrange(phenotype,threshold)
  
  
  
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(mean_fst_values,by="weighted_fst_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_fst)
    title1 = "Relative Phenotype Variation Explained by Polygenic Scores Across Fst Pools"
    xlabel1 = "Within-Group Median Fst"
    subtitle1 = "Fst was calculated between each individual against EUR training population."
    caption1 ="Total of 40 groups each composed of around 750 observations."
    delta=0.05
    beta <- "\u03b2 * 0.05 (p)"
    
  }
  if (type=="PC") {
    plot_df <- plot_df %>% left_join(mean_pc_values,by="PC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_pc)
    title1 = "Relative Phenotype Variation Explained by Polygenic Scores PC Distance Pools"
    xlabel1 = "Within-Group Median PC Distance"
    subtitle1 = "PC distance was calculated between each individual's against the geometric median of the EUR training population."
    caption1 = "" # range_pc
  }
  if (type=="WPC") {
    plot_df <- plot_df %>% left_join(mean_wpc_values,by="WPC_groups") 
    plot_df <- plot_df %>% rename(mean_values = mean_wpc)
    title1 = "Relative Phenotype Variation Explained by Polygenic Scores Across Weighted PC Distance Pools"
    xlabel1 = "Within-Group Median PC Distance"
    subtitle1 = "Weighted PC distance is Euclidean distance between each individual's 20 PCs against the geometric median of the EUR training population's PCs."
    caption1 = "Total of 40 groups each composed of around 750 observations."
    delta=50
    beta <- "\u03b2 * 50 (p)"
  }
  
  min_value <- min(plot_df$mean_values)
  max_value <- max(plot_df$mean_values)
  
  plot_df$weight <- ifelse(plot_df$group_number==1,100000,1)

  plot <- ggplot(plot_df, aes(x=mean_values-min_value,y=relative_performance-1, #or partial
                              color=phenotype))
  if(line ==1){
    
    plot <- plot + geom_smooth(size=1,method = lm,formula=y ~ 0+x ,se=F,na.rm=T)
    
    tb <- plot_df %>% group_by(phenotype, threshold) %>%
      do(
        lm = lm(I(relative_performance-1)~0+I(mean_values-min_value), data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[1,c(1)]*delta,3),
                           " (",
                           round(summary(lm)$coef[1,c(4)],3),
                           ")")) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coef,`LM R2`) 
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
 #   colnames(tb)[colnames(tb) %in% "LM R2"] <- expression(LM~R2)
    h <- -.4
  }

  if(line ==5){
    
    plot <- plot + geom_smooth(aes(weight=weight),size=1,span=0.75,
                               method="loess",formula=y~0+x,na.rm=T,se=F) #4
    
    tb <- plot_df %>%
      group_by(phenotype, threshold) %>%
      do(
        ls = loess(I(relative_performance-1)~0+I(mean_values-min_value), data=.,
                   weights=weight)
      )  %>%
      ungroup() 
    tb$ls_fit <-  apply(sapply(tb$ls,fitted),2,list)
    tb <- tb %>%
      unnest(cols=ls_fit) %>%
      unnest(cols=ls_fit) %>%
      mutate(group_number = rep(1:40,10*5)) %>%
      select(group_number, phenotype, threshold, ls_fit) %>%
      right_join(plot_df , by=c("phenotype","threshold","group_number")) %>%
      group_by(phenotype, threshold) %>%
      do(
        lm = lm(I(ls_fit)~0+I(mean_values-min_value),data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[1,c(1)]*delta,3),
                           " (",
                           round(summary(lm)$coef[1,c(4)],3),
                           ")")) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coef,`LM R2`) %>%
      rename("Phenotype"='phenotype')
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
#    colnames(tb)[colnames(tb) %in% "LM R2"] <- expression(LM~R^2)
    h <- -.4

  }
  
    tbs <- lapply(split(tb, tb$threshold), "[", -1)
    df <- tibble(x = rep(-Inf, length(tbs)), 
               y = rep(Inf, length(tbs)), 
               threshold = levels(as.factor(tb$threshold)), 
               tbl = tbs) %>%
    filter(threshold !="5e-8")
  
    prs_names <- c(
    `0-LDpred2` = "LDpred2",
    `1e-2` = "1e-2",
    `1e-3` = "1e-3",
    `1e-4` = "1e-4",
    `1e-5` = "1e-5")
  
  
    plot <- plot + geom_point(size=1.5,
                              color="white")+
    geom_point(size=1,alpha=0.75)+
    guides(color = guide_legend(title="Trait"),
           size=F)+
    theme_light() +
      theme(legend.text=element_text(size=14),
            legend.title=element_text(size=16))+
    xlab(xlabel1) +
    geom_hline(yintercept=0, linetype="dashed") +
    ylab(expression(Partial~R^2~relative~to~Europeans))+
    labs(
      #title=title1,
      #   subtitle=subtitle1,
         caption=caption1)+
    coord_cartesian(ylim=c(-1, 1))+
    scale_y_continuous(breaks=c(-1,-.5,0,.5,1),labels = as.character(round(seq(0,2, by = 0.5),1))
                       )+
    scale_x_continuous(breaks=seq(0,max_value-min_value,length.out = 6),
                         labels = as.character(round(seq(min_value,max_value, length.out = 6),2)))+
    facet_grid(cols = vars(threshold), labeller = as_labeller(prs_names))

    plot <- plot + geom_table(data = df, aes(x = x, y = y, label = tbl),
                                hjust = h, vjust = 1) 
    
    return(plot)
}

for (i in c(1,5)){
  comb_fst <- combined_rc(all_list[[1]],all_list[[2]], "FST",i)
  ggsave(paste0('img1/FINAL_RC_COMB_FST',i,"_40.png"),comb_fst,width=20,height=10,dpi=1000)
  
  comb_wpc <- combined_rc(as.data.frame(all_list[3]),as.data.frame(all_list[4]), "WPC",i)
  ggsave(paste0('img1/FINAL_RC_COMB_WPC',i,"_40.png"),comb_wpc,width=20,height=10,dpi=1000)
}



combined_rc_df <- function(prs_df,ldpred_prs_df, type, line){
  
  
  ld_df <- ldpred_prs_df %>%
    filter(threshold==29) 
  
  plot_df <- rbind(prs_df,ld_df)
  
  plot_df$threshold <- factor(plot_df$threshold,levels=c(29,4,3,2,1,0),
                              labels=c("0-LDpred2","1e-2","1e-3","1e-4","1e-5","5e-8"),
                              ordered = T)
  
  plot_df <- plot_df %>%
    filter(threshold!="5e-8") %>%
    group_by(phenotype,threshold) %>%
    mutate(
      eur_partial = max((partial*(group_number==1)),na.rm=T),
      relative_performance = partial / eur_partial) %>%
    ungroup() %>%
    arrange(threshold)
  
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(mean_fst_values,by="weighted_fst_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_fst)
    delta=0.05
    
  }

  if (type=="WPC") {
    plot_df <- plot_df %>% left_join(mean_wpc_values,by="WPC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_wpc)
    delta=50
  }
  
  min_value <- min(plot_df$mean_values)
  max_value <- max(plot_df$mean_values)
  
  plot_df$weight <- ifelse(plot_df$group_number==1,100000,1)
  

  if(line ==1){
    
    tb <- plot_df %>% group_by(phenotype, threshold) %>%
      do(
        lm = lm(I(relative_performance-1)~0+I(mean_values-min_value), data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[1,c(1)]*delta,3),
                           " (",
                           round(summary(lm)$coef[1,c(4)],3),
                           ")")) %>%
      mutate(`Coefficient` = summary(lm)$coef[1,c(1)]*delta) %>%
      mutate(`p-val` = summary(lm)$coef[1,c(4)]) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coefficient,`p-val`,`LM R2`) %>%
      rename("Phenotype"=`phenotype`)
  }
  
  if(line ==5){
    tb <- plot_df %>%
      group_by(phenotype, threshold) %>%
      do(
        ls = loess(I(relative_performance-1)~0+I(mean_values-min_value), data=.,
                   weights=weight)
      )  %>%
      ungroup() 
    tb$ls_fit <-  apply(sapply(tb$ls,fitted),2,list)
    tb <- tb %>%
      unnest(cols=ls_fit) %>%
      unnest(cols=ls_fit) %>%
      mutate(group_number = rep(1:40,10*5)) %>%
      select(group_number, phenotype, threshold, ls_fit) %>%
      right_join(plot_df , by=c("phenotype","threshold","group_number")) %>%
      group_by(phenotype, threshold) %>%
      do(
        lm = lm(I(ls_fit)~0+I(mean_values-min_value),data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[1,c(1)]*delta,3),
                           " (",
                           round(summary(lm)$coef[1,c(4)],3),
                           ")")) %>%
      mutate(`Coefficient` = summary(lm)$coef[1,c(1)]*delta) %>%
      mutate(`p-val` = summary(lm)$coef[1,c(4)]) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coefficient,`p-val`,`LM R2`) %>%
      rename("Phenotype"='phenotype')
    
  }
  return(tb)
}

rc_comb <- data.frame()
for (i in c(1,5)){
  comb_fst <- combined_rc_df(as.data.frame(all_list[1]),as.data.frame(all_list[2]), "FST",i)
  comb_wpc <- combined_rc_df(as.data.frame(all_list[3]),as.data.frame(all_list[4]), "WPC",i)
  rc_comb1 <- bind_rows(comb_fst,comb_wpc)
  rc_comb <- rbind(rc_comb,rc_comb1)
  
}
rc_comb$gd <- c(rep("FST",50), rep("WPC", 50),rep("FST", 50),rep("WPC",50))
rc_comb$model <- c(rep("LM",100),rep("LOESS", 100))

all_comb <- rbind(comb,rc_comb) %>%
  mutate(type = c(rep("Raw",200),rep("Relative",200))) %>%
  group_by(Phenotype,gd,type) %>%
  mutate(max_beta = max(Coefficient)) %>%
  mutate(relative_beta = Coefficient/max_beta) %>%
  ungroup()


graph_comb <- all_comb %>%
  filter(type=="Raw") %>%
  ggplot(aes(y=relative_beta,x=`LM R2`,color=threshold, shape=model))+
  geom_point(size=4)+
  scale_color_brewer("Set1")+
  facet_grid(gd~Phenotype)
ggsave("img1/model_combinations.png",graph_comb,width=20,height=20,dpi=700)

###############################################

all_df <- prs_fst_df %>%
  bind_rows(ld_fst_df,
          #  prs_pc_df,ld_pc_df,
            prs_wpc_df,ld_wpc_df) %>%
  select(phenotype,threshold,partial,incremental_r2) %>%
  mutate(threshold=factor(threshold,labels=c("5e-8","1e-5","1e-4",
                                             "1e-3","1-e2","Naive","Unweighted LDpred",
                                             "Weighted LDpred")))# %>%
#  mutate(partial = replace_na(partial,0))
all_df$phenotype <- factor(all_df$phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
"Lymphocyte","WBC","Eosinophil"))

se_df <- all_df %>%
  group_by(phenotype) %>%
  mutate(mean_partial=mean(partial,na.rm=T),
         sem_partial=sd(partial,na.rm=T)/sqrt(length(unique(na.omit(all_df$threshold))))) %>%
  top_n(1,partial) %>%
  ungroup() 

traits <- all_df %>%
  ggplot(aes(x = phenotype, y = partial, color=phenotype, fill=phenotype)) +
  geom_violin(width = 1) +
  geom_jitter(aes(shape=as.factor(threshold)),color="gray44", size=1.5,height=0,width=0.3)+
  geom_errorbar(data=se_df, aes(x=phenotype,ymin=mean_partial - sem_partial, ymax=mean_partial + sem_partial),
                                            color='black', width=0.22,size=1.5)+
  theme_classic() +
  theme(legend.title=element_text(size=16), 
        legend.text=element_text(size=14),
        legend.position=c(0.9,0.8))+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  scale_shape_manual(values=1:8)+
  xlab('Traits') +
  ylab(expression(Partial~R^2))+
#  labs(title="Partial R2 across all PRS Models and groups",
#       caption="Traits are arranged in descending order of estimated heritability according to the Neale Lab.")+
  guides(fill = F,
         color= F,
         shape = guide_legend(title="Threshold"))

ggsave('img1/FINAL_traits_pr1.png',traits,width=15,height=10,dpi=1250)

############

df1 <- all_df[!is.na(all_df$partial),]
df1 <- df1[!is.na(df1$incremental_r2),]
cor(df1$partial,df1$incremental_r2)

##########################

p <- non_prs_df %>% mutate(population = factor(population,labels=c("EUR","AMR","SAS","EAS","AFR"))) %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  ggplot(aes(phenotype_value,color=population)) +
  geom_histogram(aes(y=..density..))+
  facet_grid(population~phenotype)+
  xlim(-5,5)+
  xlab("Standard Normalized Phenotype Value")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  ylab("Density")+
  theme_classic()+
  guides(color=F)+
  ggtitle("Distribution of Phenotypes across Major Population Groups")
ggsave("img/trait_distributions/traits.png",p,
       width=16,height=12,dpi=500)

q <- non_prs_df %>%
  mutate(`Fst Group`= factor(weighted_fst_groups)) %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  filter(phenotype!="Lymphocyte") %>%
  drop_na(`Fst Group`) %>%
  filter(weighted_fst_groups !="(0.107,0.121]") %>%
  ggplot(aes(phenotype_value,color=`Fst Group`)) +
#  geom_histogram() +
 geom_histogram(aes(y=..density..))+
  facet_grid(`Fst Group`~phenotype)+
  xlim(-5,5)+
  xlab("Standard Normalized Phenotype Value")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  ylab("Density")+
  theme_classic()+
  guides(color=F)+
  ggtitle("Distribution of Phenotypes across Fst Groups")
ggsave("img/traits_fst.png",q,
       width=16,height=12,dpi=500)

q <- non_prs_df %>%
  mutate(`WPC Groups`= factor(WPC_groups)) %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  filter(phenotype!="Lymphocyte") %>%
  drop_na(`WPC Groups`) %>%
#  filter(weighted_fst_groups !="(0.107,0.121]") %>%
  ggplot(aes(phenotype_value,color=`WPC Groups`)) +
  #  geom_histogram() +
  geom_histogram(aes(y=..density..))+
  facet_grid(`WPC Groups`~phenotype)+
  xlim(-5,5)+
  xlab("Standard Normalized Phenotype Value")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  ylab("Density")+
  theme_classic()+
  guides(color=F)+
  ggtitle("Distribution of Phenotypes across Weighted PC Groups")
ggsave("img/traits_wpc.png",q,
       width=16,height=12,dpi=500)


# Set up EUR train


non_prs_df_sub <- non_prs_df %>%
  left_join(mean_fst_values,by="weighted_fst_groups") %>%
  left_join(mean_wpc_values,by="WPC_groups")


r <- non_prs_df_sub %>%
  mutate(`Fst Group`= mean_fst) %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  drop_na(`Fst Group`) %>%
  group_by(`Fst Group`,phenotype) %>%
  mutate(pheno_mean = mean(phenotype_value,na.rm=T)) %>%
  ungroup() %>%
  distinct(`Fst Group`,phenotype,.keep_all=T) %>%
  ggplot(aes(x=mean_fst,y=pheno_mean,group=phenotype,color=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic()+
  ylab("Phenotype Mean")+
  xlab("Within-Group Median Fst")+
#  ggtitle("Phenotype Mean across Fst Pools")+
#  labs(caption="Phenotypes were standard-normalized with respect to the EUR train's mean and variance.")+
  scale_x_continuous(breaks = round(seq(0,
                                        0.13, by = 0.02),2)) +
  scale_y_continuous(breaks = round(seq(-1.5,
                                        1.5, by = .5),1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))+
  scale_color_discrete(name = "Phenotype")
ggsave("img1/FINAL_traits_mean_fst.png",r,
       width=16,height=12,dpi=1000)



covar_df  <- read_tsv('data/ukb_merged/covar_all_samples.covar') %>% select("#FID",IID,pop)
phenotypes_df <- read_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', trim_ws = T) %>%
  pivot_longer(BMI:Eosinophil,names_to="phenotype",values_to="phenotype_value") %>%
  left_join(covar_df,by=c("#FID","IID")) %>%
  filter(pop=="EUR_train") %>%
  group_by(phenotype) %>%
  summarize(EUR_train_var = var(phenotype_value,na.rm=T))


r <- non_prs_df_sub %>%
  mutate(`Fst Group`= mean_fst) %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  drop_na(`Fst Group`) %>%
  group_by(`Fst Group`,phenotype) %>%
  mutate(pheno_var = var(phenotype_value,na.rm=T)) %>%
  ungroup() %>%
  distinct(`Fst Group`,phenotype,.keep_all=T) %>%
  ggplot(aes(x=mean_fst,y=pheno_var,group=phenotype,color=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept=1, linetype="dashed") +
  theme_classic()+
  ylab("Phenotype Variance")+
  xlab("Within-Group Median Fst")+
#  ggtitle("Phenotype Variance across Fst Pools")+
#  labs(caption="Phenotypes were standard-normalized with respect to the EUR train's mean and variance.")+
  scale_x_continuous(breaks = round(seq(0,
                                        0.13, by = 0.02),2)) +
  scale_y_continuous(breaks = round(seq(0,
                                        3, by = .5),1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))+
  scale_color_discrete(name = "Phenotype")
ggsave("img1/FINAL_traits_var_fst.png",r,
       width=16,height=12,dpi=1000)



r <- non_prs_df_sub %>%
  rename(`WPC Group`= `WPC_groups`) %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  #  filter(phenotype!="Lymphocyte") %>%
  drop_na(`WPC Group`) %>%
  #  filter(`Fst Group` !="(0.107,0.121]") %>%
  group_by(`WPC Group`,phenotype) %>%
  mutate(pheno_mean = mean(phenotype_value,na.rm=T)) %>%
  ungroup() %>%
  distinct(`WPC Group`,phenotype,.keep_all=T) %>%
  ggplot(aes(x=mean_wpc,y=pheno_mean,group=phenotype,color=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic()+
  ylab("Phenotype Mean")+
  xlab("Within-Group Median Weighted PC Distance")+
  #ggtitle("Phenotype Mean across Weighted PC Pools")+
  #labs(caption="Phenotypes were standard-normalized with respect to the EUR train's mean and variance.")+
  scale_x_continuous(breaks = round(seq(0,
                                        140, by = 25),2)) +
  scale_y_continuous(breaks = round(seq(-1.5,
                                        1.5, by = .5),1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))+
  scale_color_discrete(name = "Phenotype")
ggsave("img1/FINAL_traits_mean_wpc.png",r,
       width=16,height=12,dpi=1000)



r <- non_prs_df_sub %>%
  rename(`WPC Group`= `WPC_groups`) %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  drop_na(`WPC Group`) %>%
  group_by(`WPC Group`,phenotype) %>%
  mutate(pheno_var = var(phenotype_value,na.rm=T)) %>%
  ungroup() %>%
  distinct(`WPC Group`,phenotype,.keep_all=T) %>%
  ggplot(aes(x=mean_wpc,y=pheno_var,group=phenotype,color=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept=1, linetype="dashed") +
  theme_classic()+
  ylab("Phenotype Variance")+
  xlab("Within-Group Median Weighted PC Distance")+
 # ggtitle("Phenotype Variance across Weighted PC Pools")+
#labs(caption="Phenotypes were standard-normalized with respect to the EUR train's mean and variance.")+
  scale_x_continuous(breaks = round(seq(0,
                                        140, by = 25),2)) +
  scale_y_continuous(breaks = round(seq(0,
                                        3, by = .5),1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))+
  scale_color_discrete(name = "Phenotype")
ggsave("img1/FINAL_traits_var_wpc.png",r,
       width=16,height=12,dpi=1000)


##########

covar_df  <- read_tsv('data/ukb_merged/covar_all_samples.covar') %>% select(`#FID`,IID,PC1_AVG,PC2_AVG,pop)
wb_df <- read.table('data/ukb_merged/wb.fam') %>% select(V1,V2)
colnames(wb_df) <- c("#FID","IID")
wb_df$`WB Status` <- "White British"

pc_df <- non_prs_df %>% 
  full_join(covar_df, by =c("#FID","IID")) %>%
  left_join(wb_df, by =c("#FID","IID")) %>%
  mutate(`WB Status` = ifelse(is.na(`WB Status`),"Not White British","White British")) %>%
  select(`#FID`,IID,PC1_AVG.y,PC2_AVG.y, pop.y,`WB Status`) %>%
  distinct(`#FID`,IID,.keep_all=T) %>%
  rename("Population" = "pop.y") %>%
  rename("PC1" = "PC1_AVG.y") %>%
  rename("PC2" = "PC2_AVG.y") %>%
  drop_na(Population) %>%
  mutate(Population=factor(Population,levels=c("EUR_test","AMR","SAS","EAS","AFR","EUR_train"))) %>%
  arrange(desc(Population))


medians <- data.frame(Population=unique(pc_df$Population),
                      PC1 = rep(NA,6),
                      PC2= rep(NA,6)) 
for(i in unique(pc_df$Population)){
  medians[medians$Population==i,c("PC1","PC2")] <- Gmedian(pc_df[pc_df$Population==i,c("PC1","PC2")])
}

sp <- pc_df %>%
  ggplot(aes(x=PC1,y=PC2,color=Population))+
  geom_point(size=0.6)+
  geom_point(data=medians,  mapping=aes(x = PC1, y = PC2,fill=Population),
             size=4,shape="square",color="black")+
  geom_point(data=medians,  mapping=aes(x = PC1, y = PC2,fill=Population),
             size=3.5,shape="square")+
  scale_colour_hue(l = 70, c = 150)+
  guides(color = guide_legend(override.aes = list(size = 3)),
         fill=F)+
  theme(legend.title=element_text(size=16), 
         legend.text=element_text(size=14),
        legend.position = c(0.1,0.8))

ggsave("img1/FINAL_PC1_PC2_all.png",sp,
       width=16,height=12,dpi=1000)

############################
scatter_plot <- non_prs_df %>%
  filter(phenotype=="BMI") %>%
  ggplot() +
  geom_point(aes(x=Weighted_Fst,y=PC_dist_weighted,color=population)) +
  theme_classic()+
  xlab("Weighted Fst") +
  ylab("Weighted PC Distance")+
#  ggtitle("Distribution of Weighted Fst and PC Distance in Test Individuals from Training Group")+
  theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16,face="bold"),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=14),
          legend.position=c(.9,.6))+
  scale_x_continuous(breaks = round(seq(min(non_prs_df$Weighted_Fst),
                                        max(non_prs_df$Weighted_Fst), by = 0.02),2)) +
  scale_y_continuous(breaks = round(seq(0,
                                        max(non_prs_df$PC_dist_weighted), by = 25),1))+
  scale_color_discrete(name="Population")

ggsave("img1/FINAL_fst_wpc.png",scatter_plot,
       width=16,height=12,dpi=1000)


smear <- non_prs_df %>%
  filter(phenotype=="BMI") %>%
  filter(((Weighted_Fst>=0.01) & (PC_dist_weighted<25) & (PC_dist_weighted >=20) & (population=="AMR"))|
           ((Weighted_Fst>=0.005) & (PC_dist_weighted<15) & (PC_dist_weighted >=10) & (population=="AMR"))|
           ((Weighted_Fst>=0.01) & (PC_dist_weighted<48) & (PC_dist_weighted >=40) & (population=="SAS"))|
           ((Weighted_Fst>=0.012) & (PC_dist_weighted<55) & (PC_dist_weighted >=50) & (population=="SAS"))|
           ((Weighted_Fst>=0.015) & (PC_dist_weighted<66) & (PC_dist_weighted >=60) & (population=="SAS"))|
           ((Weighted_Fst>=0.04) & (PC_dist_weighted<75) & (PC_dist_weighted >=70) & (population=="AFR"))
           
  ) %>%
  mutate(smears = 1) %>%
  select(`#FID`,IID, smears)

smear_non_prs_df <- non_prs_df %>%
  left_join(smear,by=c("#FID","IID"))  %>%
  filter(phenotype=="BMI")
smear_non_prs_df$smears[is.na(smear_non_prs_df$smears)] <- 0 

scatter_plot <-  ggplot() +
  geom_point(data=smear_non_prs_df %>% filter(smears==0),
             aes(x=PC1_AVG,y=PC2_AVG,color=population),
             size=2,alpha=0.2) +
  geom_point(data=smear_non_prs_df %>% filter(smears==1),
             aes(x=PC1_AVG,y=PC2_AVG),color="black",
             size=3.5,alpha=0.7) +
  geom_point(data=smear_non_prs_df %>% filter(smears==1),
             aes(x=PC1_AVG,y=PC2_AVG,color=population),
             size=3,alpha=0.9) +
  #  geom_vline(xintercept=unique(fst_values1), linetype="dashed") +
  #  geom_hline(yintercept=unique(WPC_values), linetype="dashed") +
  theme_classic()+
  xlab("PC1") +
  ylab("PC2")+
  ggtitle("PC1 vs PC2 in 'Smear' Test Individuals")+
  #  labs(subtitle="Lines refer to boundaries of genetic distance groups.")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))#+
 # scale_x_continuous(breaks = round(seq(min(non_prs_df$Weighted_Fst),
  #                                      max(non_prs_df$Weighted_Fst), by = 0.02),2)) +
  #scale_y_continuous(breaks = round(seq(0,
   #                                     max(non_prs_df$PC_dist_weighted), by = 25),1))


ggsave("img1/PC1_PC2_smears.png",scatter_plot,
       width=16,height=12,dpi=500)



scatter_plot <- ggplot() +
  geom_point(data=smear_non_prs_df %>% filter(smears==0),
             aes(x=Weighted_Fst,y=PC_dist_weighted,color=population),
             size=2,alpha=0.2) +
  geom_point(data=smear_non_prs_df %>% filter(smears==1),
             aes(x=Weighted_Fst,y=PC_dist_weighted),color="black",
             size=3.5,alpha=0.7) +
  geom_point(data=smear_non_prs_df %>% filter(smears==1),
             aes(x=Weighted_Fst,y=PC_dist_weighted,color=population),
             size=3,alpha=0.9) +
  theme_classic()+
  xlab("Weighted Fst") +
  ylab("Weighted PC Distance")+
  ggtitle("Distribution of Weighted Fst and PC Distance in Test Individuals from Training Group")+
  #  labs(subtitle="Lines refer to boundaries of genetic distance groups.")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))+
  scale_x_continuous(breaks = round(seq(min(non_prs_df$Weighted_Fst),
                                        max(non_prs_df$Weighted_Fst), by = 0.02),2)) +
  scale_y_continuous(breaks = round(seq(0,
                                        max(non_prs_df$PC_dist_weighted), by = 25),1))


ggsave("img1/fst_wpc_smeared.png",scatter_plot,
       width=16,height=12,dpi=500)


######

df <- non_prs_df[!is.na(non_prs_df$Weighted_Fst),]
df <- df[!is.na(df$PC_dist),]
cor(df$Weighted_Fst,df$PC_dist)

#######

