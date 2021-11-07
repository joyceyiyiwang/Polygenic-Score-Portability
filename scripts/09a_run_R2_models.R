.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
#.libPaths("/moto/palab/users/jm4454/rpackages/")
#library(asbio)
#library(cowplot)
#library(ggrepel)
library(ggplot2)
library(ggpmisc)
library(plyr)
library(tidyverse)
#library(gridExtra)
library(RColorBrewer)
library(bigsnpr)
library(Gmedian)
#library(philentropy)


load_non_prs_df <- function(fst_file_name = "final_fst",num_fst_groups, n_buckets=60, first = 1) {
  # Loads all the covariates, phenotypes, and outcomes into a dataframe
  
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar_all_samples.covar')
  covar_df <- covar_df %>% select(-c(28:47))
  
  # Unweighted Medians
  medians <- covar_df %>% filter(pop=="WB_train") %>% 
    select(starts_with("PC")) %>% Gmedian() %>% as.vector()
  covar_df <- covar_df %>% filter(pop!="WB_train")
  
  phenotypes_df <- read_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', trim_ws = T)
  
  # Load the individuals who are in the evaluation set for each population
  population_files <- c('data/ukb_populations/NWB_all.txt', 'data/ukb_populations/WB_train.txt', 
                        'data/ukb_populations/WB_test.txt')
  populations_df <- data.frame()
  for (file in population_files) {
    pop_df <- read_delim(file, delim = ' ', trim_ws = T,
                         col_types = c('#FID' = col_integer(), 'IID' = col_integer())) %>%
      mutate(population = str_extract(file, '(?<=data/ukb_populations/)[A-Z]{3}'))
    populations_df <- bind_rows(populations_df, pop_df)
  }
  fst_values <- read_tsv(paste0('data/fst/', fst_file_name, '.tsv'))
  
  # Combine all the above tables into a table that will be joined with PRS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    left_join(fst_values,by=c("IID")) %>% #
  #  mutate(Weighted_Fst=ifelse(Weighted_Fst<=0,0,Weighted_Fst)) %>%
  #  mutate(Weighted_Fst=ifelse((is.na(Weighted_Fst)) & (population=="EUR_test"),0,Weighted_Fst)) %>%
    pivot_longer(BMI:Eosinophil, names_to = 'phenotype', values_to = 'phenotype_value')
  
 
  #cutpoints_fst <<- seq(min(phenotypes_df$Weighted_Fst,na.rm=T),
  #                  max(phenotypes_df$Weighted_Fst,na.rm=T),
  #                  length=num_fst_groups+1)
  set.seed(3)
  # Randomize the IDs first because ntile() assign groups to the same values based on in order
  random_id <- sample(unique(phenotypes_df$IID))
  random_id <- rep(random_id, each = 10)
  random_id <- cbind.data.frame(IID = random_id,
                                phenotype = rep(c("BMI", "WBC", "Height", "RBC", "MCV", "MCH", 
                                                  "Lymphocyte", "Platelet", "Monocyte", 
                                                  "Eosinophil"), length(random_id) / 10))
  phenotypes_df <- random_id %>% left_join(phenotypes_df, by = c("IID", "phenotype"))
  phenotypes_df <- phenotypes_df[, c(3, 1:2, 4:ncol(phenotypes_df))]
  phenotypes_df <- phenotypes_df %>%
    mutate(weighted_fst_groups = Weighted_Fst %>% ntile(n_buckets + first - 1)) %>% #Divides into equally-sized components
    # Make first group 10x larger
    mutate(weighted_fst_groups = 
             ifelse(weighted_fst_groups <= first, 1, weighted_fst_groups - first + 1)) %>%
    arrange(IID) #%>%
    #add_count(weighted_fst_groups)# %>%
    #rename(fst_groups_count=n)
  
  return(phenotypes_df)
}

get_r2_values <- function(df,group_var) {
  # Computes the incremental R^2 attributable to the PRS
  #    group_by(population, phenotype, threshold) %>%
  df$group <- df %>% select(any_of(group_var)) %>% as.data.frame()
  df1 <- df %>%
    group_by(group, phenotype, threshold) %>%
    do(
      # Regression on only covariates
      nested = lm(phenotype_value ~ sex_covar + age + age_sq + age_sex + age_sq_sex +
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                    PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = .),
      
      # Regression now including the PRS + covariates
      full = lm(phenotype_value ~ prs + sex_covar + age + age_sq + age_sex + age_sq_sex +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
                  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = .)
    )
  nested_r2 = c()
  full_r2 = c()
  for(i in 1:nrow(df1)){
    nested_r2 = c(nested_r2, summary(df1[[4]][[i]])$r.squared)
    full_r2 = c(full_r2, summary(df1[[5]][[i]])$r.squared)
  }
  df1 <- df1 %>%
    cbind.data.frame(nested_r2 = nested_r2, full_r2 = full_r2) %>%
    mutate(
      #partial = partial.R2(nested, full),
      #nested_r2 = summary(nested)$adj.r.squared,
      #full_r2 = summary(full)$adj.r.squared,
      #nested_r2 = nested_r2,
      #full_r2 = full_r2,
      incremental_r2 = full_r2 - nested_r2
    ) %>%
    select(-nested, -full)
  df1$group <- 1:nrow(df1)
  colnames(df1)[1] <- group_var
  return(df1)
}

make_prs_evaluation_df <- function(non_prs_df,group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  prs_df <- data.frame()
  for (file in list.files(path = 'data/prs',
                          pattern = '[a-zA-Z]+_[0-9]_scores.sscore', full.names = T)) {
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', prs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/prs/)[A-Za-z]+(?=_)'),
        threshold = str_extract_all(string = file, pattern = '[0-4]')[[1]][1] %>% as.integer
      ) %>%
      inner_join(non_prs_df, by =  c('#FID', 'IID','phenotype'))%>%
      get_r2_values(.,group_var=group_var_as_string)
    
    prs_df <- bind_rows(prs_df, this_df)
  }
  prs_df
}
# Comparison


####################
#Previously 40
non_prs_df <- load_non_prs_df("final_fst",num_fst_groups = 100, n_buckets=100, first = 50)
# This is the file with FST calculated using all autosomes (but on only 10% individuals),
# Load it for comparison only
non_prs_full_df <- load_non_prs_df("final_fst_full",num_fst_groups = 100, n_buckets=100, first = 50)
prs_df_fst  <- make_prs_evaluation_df(non_prs_df,"weighted_fst_groups")

#prs_df_pop <- make_prs_evaluation_df(non_prs_df,"population")  

prs_fst_df <- prs_df_fst %>% 
  arrange(phenotype,weighted_fst_groups,threshold) %>%
  mutate(group_number = weighted_fst_groups)  #%>%
##  filter(weighted_fst_groups!=10)

#prs_pop_df <- prs_df_pop %>%
#  mutate(population = factor(population,levels=c(2,1),labels = c("WB","NWB")))

# Save the files
non_prs_df %>% write_tsv('data/prs_comparisons/non_prs_df.tsv')
non_prs_full_df %>% write_tsv('data/prs_comparisons/non_prs_full_df.tsv')
prs_fst_df %>% write_tsv('data/prs_comparisons/UKBB_geno_fst_CT.tsv')
#prs_pop_df %>% write_tsv('data/prs_comparisons/UKBB_geno_pop.tsv')


######### Separate FST groups for h2 calculations ######

#for(i in na.omit(unique(non_prs_df$weighted_fst_groups))){
#  group = non_prs_df %>% filter(weighted_fst_groups == i) %>%
#    select(`#FID`, IID) %>% distinct()
#  group %>% write.table(paste0('data/fst/group_', i, '_id.txt'), row.names=F, quote = FALSE)
#}





#########################################################

#non_prs_df <- read_tsv('data/prs_comparisons/non_prs_df.tsv')
#non_prs_full_df <- read_tsv('data/prs_comparisons/non_prs_full_df.tsv')
#prs_fst_df <- read_tsv('data/prs_comparisons/UKBB_geno_fst_CT.tsv')
#prs_pop_df <- read_tsv('data/prs_comparisons/UKBB_geno_pop.tsv')
#######################

# Save the IDs in each group
groups <- unique(non_prs_df$weighted_fst_groups) %>% sort()

ids <- non_prs_df %>% 
  select("#FID","IID","weighted_fst_groups") %>%
  distinct()

i <- 1

for (group in groups){
  sub <- ids %>% filter(weighted_fst_groups==group) %>% select(-weighted_fst_groups)
  if(nrow(sub) >250){
    sub_name <- paste0("data/theory/fst_",i,".txt")
    sub %>% write.table(sub_name, sep="\t",  col.names=FALSE, row.names = F)
    i <- i + 1
  }
}


#####################

# Calculate the median FST for each group
get_median_fst = function(file){
  median_fst_values <- file %>% 
    filter(phenotype =="Height") %>%
    select(weighted_fst_groups, Weighted_Fst) %>%
    group_by(weighted_fst_groups) %>%
    summarize(median_fst = median(Weighted_Fst),
              count = n()) %>%
    na.omit() %>%
    as.data.frame()
  return(median_fst_values)
}

median_fst = get_median_fst(non_prs_df)

# Create a scatterplot between FST values calculated using all autosomes and using just chromosome 1
non_prs_full_df <- non_prs_full_df[!is.na(non_prs_full_df$Weighted_Fst), ]
non_prs_full_df <- non_prs_full_df %>% 
  filter(phenotype == "Height") %>%
  select(IID, Weighted_Fst)
colnames(non_prs_full_df) <- c("IID", "Weighted_Fst_22")
chrom_compare = non_prs_df %>% 
  filter(phenotype == "Height") %>%
  select(IID, Weighted_Fst) %>% 
  inner_join(non_prs_full_df, by = "IID")
# Calculate and formate r^2 and p-value
r_2 = round(summary(lm(chrom_compare$Weighted_Fst_22~chrom_compare$Weighted_Fst))$r.squared, 3)
p = round(summary(lm(chrom_compare$Weighted_Fst_22~chrom_compare$Weighted_Fst))[[4]][2,4], 3)
p = sprintf(p, fmt = '%#.3f')
# Create and save the plot
corr = chrom_compare %>% ggplot(aes(x = Weighted_Fst, y = Weighted_Fst_22)) + 
  geom_point(size=4.5,color="white")+
  geom_point(size=4,alpha=0.75, color = "#B78AF3") +
  geom_smooth(method="lm",se=F,size=3, color = "white") +
  geom_smooth(method="lm",se=F,size=2.5, color = "black") +
  geom_smooth(method="lm",se=F, size = 2,color = "#B78AF3") +
  xlab("Chromosome 1") +
  ylab("All Autosomes") +
  ggtitle("Correlation between FST of 10% randomly chosen test individuals") +
  labs(subtitle = paste0("r \u00b2 = ", r_2, " (p = ", p, ")"),
        caption = paste("n =", nrow(chrom_compare)))+
  theme_light() + 
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        plot.title=element_text(size=28),
        plot.subtitle=element_text(size=20),
        plot.caption=element_text(size=16))
ggsave(paste0('img/correlation.png'),corr,width=20,height=10,dpi=1000)

##########################################



################################
all_list <- list(prs_fst_df)

pheno_factor <- function(df){
  mutate(df,phenotype = factor(phenotype, levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                                   "Lymphocyte","WBC","Eosinophil")))
}
all_list <- lapply(all_list,pheno_factor) 


all_df <- ldply(all_list, rbind)
all_df$model <- c(rep("Fst",5000))
all_df %>% filter(threshold %in% c(1,2,3,4,29)) %>%
  group_by(model, phenotype) %>%
  #mutate(max_h2 = max(partial)) %>%
  #select(model,phenotype,threshold,max_h2) %>%
  #distinct(model,phenotype,max_h2,.keep_all=T) %>%
  select(model,phenotype,threshold) %>%
  distinct(model,phenotype,.keep_all=T) %>%
  arrange(phenotype,model)

# Function for formatting decimals
keep_decimal_places <- function(x, num_digit){
  return(sprintf(paste0("%.", num_digit, "f"), x))
}

# Create the plots
# line: int, 1 for the regular linear regression; 
#            2 for 2 linear regressions separated at x = 0.01;
#            5 for splines
# c: str, hex code of a color
combined <- function(prs_df, type, line, trait, c, median_fst){

  plot_df <- prs_df
  
  plot_df$threshold <- factor(plot_df$threshold,levels=c(29,4,3,2,1,0),
                              labels=c("0-LDpred2","1e-2","1e-3","1e-4","1e-5","5e-8"),
                              ordered = T)

  denominator = plot_df %>%
    filter(threshold!="5e-8") %>%
    group_by(phenotype,threshold) %>% 
    filter(group_number == 1)
  denominator = denominator$incremental_r2
  denominator_new = c()
  for(i in 0:9){
    temp = denominator[(i * 4 + 1) : ((i + 1)*4)]
    denominator_new = c(denominator_new, rep(temp, 100))
  }
  
  
  plot_df <- plot_df %>%
    filter(threshold!="5e-8") %>%
    group_by(phenotype,threshold) %>%
    cbind.data.frame(eur_incremental_r2 = denominator_new) %>%
    mutate(
      #eur_incremental_r2 = max((incremental_r2*(group_number==1)),na.rm=T),
      #eur_incremental_r2 = denominator_new,
      relative_performance = incremental_r2 / eur_incremental_r2) %>%
    ungroup() %>%
    arrange(threshold)
  
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(median_fst,by="weighted_fst_groups") 
    #plot_df <- plot_df %>% rename(mean_values = mean_fst)
    plot_df$median_values = plot_df$median_fst
    plot_df <- plot_df %>% select(-median_fst)
    title1 = paste("Variation in", trait, "Explained by Polygenic Scores Across FST Bins")
    xlabel1 = "Within-Group Median FST"
    subtitle1 = "FST was calculated between each individual against WB training population."
    group_min = min(median_fst$count)
    group_max = max(median_fst$count[median_fst$weighted_fst_groups != 1])
    xtimes = round(median_fst$count[median_fst$weighted_fst_groups == 1] / group_min)
    caption1 =  paste0("Total of ", nrow(median_fst)," groups composed of ",
                       group_min, "-", group_max," each, except for the first group (", xtimes, 
                       "x larger).")
    #delta <- 0.05
    beta <- "\u03b2 (p)"
  }

  #plot <-  ggplot(plot_df, aes(x=mean_values,y=incremental_r2,#y=partial, #or partial
  #                             color=phenotype))
  plot <-  ggplot(plot_df, aes(x=median_values,y=relative_performance,#y=partial, #or partial
                               color=phenotype)) + 
    geom_hline(yintercept = 0, size = 0.5) +
    geom_vline(xintercept = 0, size = 0.5) +
    geom_point(data = subset(plot_df, phenotype == trait),size=4.5,color="white")+
    geom_point(data = subset(plot_df, phenotype == trait),size=4,alpha=0.75, color = c)
  if(line ==1){
  plot <- plot +
    geom_smooth(data = subset(plot_df, phenotype == trait),size=3,method = lm,se=F,na.rm=T,
                color = "white") +
    geom_smooth(data = subset(plot_df, phenotype == trait),size=2.5,method = lm,se=F,na.rm=T,
                color = "black") +
    geom_smooth(data = subset(plot_df, phenotype == trait),size=2,method = lm,se=F,na.rm=T,
                color = c)
    
    tb <- plot_df %>% group_by(phenotype, threshold) %>%
      # mutate(`LM R2` = 
      #          round(summary(lm(relative_performance~mean_values))$adj.r.squared,3)) %>%
      # mutate(Coef = paste0(round(summary(lm(relative_performance~mean_values))$coef[2,c(1)],3),
      #                      " (",
      #                      round(summary(lm(relative_performance~mean_values))$coef[2,c(4)],3),
      #                      ")")) %>%
       do(
         #lm = lm(incremental_r2~mean_values, data=.)
         lm = lm(relative_performance~median_values, data=.)
       )
      # mutate(`LM R2` = round(summary(lm[[1]])$adj.r.squared,3)) %>%
      # #mutate(Coef = paste0(round(summary(lm)$coef[2,c(1)]*delta,3),
      # mutate(Coef = paste0(round(summary(lm[[1]])$coef[2,c(1)],3),
      #                      " (",
      #                      round(summary(lm[[1]])$coef[2,c(4)],3),
      #                      ")")) %>%
    r = c()
    b = c()
    for(i in 1:nrow(tb)){
      r = c(r, keep_decimal_places(round(summary(tb[[3]][[i]])$r.squared,3), 3))
      b = c(b, paste0(keep_decimal_places(round(summary(tb[[3]][[i]])$coef[2,c(1)],3), 3),
                      " (",
                      keep_decimal_places(round(summary(tb[[3]][[i]])$coef[2,c(4)],3), 3),
                      ")"))
    }
    tb <- tb %>% ungroup() %>%
      mutate(r2 = r, Coef = b) %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coef,r2) #%>%
      #rename("Phenotype"="phenotype") 
    colnames(tb)[2] = "Phenotype"
    
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
    colnames(tb)[colnames(tb) %in% "r2"] <- "r \u00b2"
    h <- -.5
  }
  
  if(line ==2){
    plot <- plot +
      geom_smooth(data = subset(plot_df, phenotype == trait & median_values <= 0.01),size=3,method = lm,se=F,na.rm=T,
                  color = "white") +
      geom_smooth(data = subset(plot_df, phenotype == trait & median_values <= 0.01),size=2.5,method = lm,se=F,na.rm=T,
                  color = "black") +
      geom_smooth(data = subset(plot_df, phenotype == trait & median_values <= 0.01),size=2,method = lm,se=F,na.rm=T,
                  color = c) +
      geom_smooth(data = subset(plot_df, phenotype == trait & median_values > 0.01),size=3,method = lm,se=F,na.rm=T,
                  color = "white") +
      geom_smooth(data = subset(plot_df, phenotype == trait & median_values > 0.01),size=2.5,method = lm,se=F,na.rm=T,
                  color = "black") +
      geom_smooth(data = subset(plot_df, phenotype == trait & median_values > 0.01),size=2,method = lm,se=F,na.rm=T,
                  color = c)+
      geom_vline(xintercept = 0.01, size = 1, linetype = "dashed")
    
    tb <- plot_df %>% group_by(phenotype, threshold) %>%
      # mutate(`LM R2` = 
      #          round(summary(lm(relative_performance~mean_values))$adj.r.squared,3)) %>%
      # mutate(Coef = paste0(round(summary(lm(relative_performance~mean_values))$coef[2,c(1)],3),
      #                      " (",
      #                      round(summary(lm(relative_performance~mean_values))$coef[2,c(4)],3),
      #                      ")")) %>%
      do(
        #lm = lm(incremental_r2~mean_values, data=.)
        lm1 = lm(relative_performance~median_values, data=subset(., median_values <= 0.01)),
        lm2 = lm(relative_performance~median_values, data=subset(., median_values > 0.01))
      )
    # mutate(`LM R2` = round(summary(lm[[1]])$adj.r.squared,3)) %>%
    # #mutate(Coef = paste0(round(summary(lm)$coef[2,c(1)]*delta,3),
    # mutate(Coef = paste0(round(summary(lm[[1]])$coef[2,c(1)],3),
    #                      " (",
    #                      round(summary(lm[[1]])$coef[2,c(4)],3),
    #                      ")")) %>%
    r1 = c()
    b1 = c()
    r2 = c()
    b2 = c()
    for(i in 1:nrow(tb)){
      r1 = c(r1, keep_decimal_places(round(summary(tb[[3]][[i]])$r.squared,3), 3))
      b1 = c(b1, paste0(keep_decimal_places(round(summary(tb[[3]][[i]])$coef[2,c(1)],3), 3),
                      " (",
                      keep_decimal_places(round(summary(tb[[3]][[i]])$coef[2,c(4)],3), 3),
                      ")"))
      r2 = c(r2, keep_decimal_places(round(summary(tb[[4]][[i]])$r.squared,3), 3))
      b2 = c(b2, paste0(keep_decimal_places(round(summary(tb[[4]][[i]])$coef[2,c(1)],3), 3),
                        " (",
                        keep_decimal_places(round(summary(tb[[4]][[i]])$coef[2,c(4)],3), 3),
                        ")"))
    }
    tb <- tb %>% ungroup() %>%
      mutate(r21 = r1, Coef1 = b1, r22 = r2, Coef2 = b2) %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coef1,r21,Coef2,r22) #%>%
    #rename("Phenotype"="phenotype") 
    colnames(tb)[2] = "Phenotype"
    
    colnames(tb)[colnames(tb) %in% "Coef1"] <- "x \u2264 0.01: \u03b2 (p)"
    colnames(tb)[colnames(tb) %in% "r21"] <- "r \u00b2"
    colnames(tb)[colnames(tb) %in% "Coef2"] <- "x \u003e 0.01: \u03b2 (p)"
    colnames(tb)[colnames(tb) %in% "r22"] <- "r \u00b2"
    h <- -0.25
  }

  
  if(line ==5){
    plot <- plot +
      geom_smooth(data = subset(plot_df, phenotype == trait),size=3,
                  method="loess",span=0.75,na.rm=T,se=F, color = "white")+
      geom_smooth(data = subset(plot_df, phenotype == trait),size=2.5,
                  method="loess",span=0.75,na.rm=T,se=F, color = "black")+
      geom_smooth(data = subset(plot_df, phenotype == trait),size=2,
                  method="loess",span=0.75,na.rm=T,se=F, color = c) #4
    
    
    tb <- plot_df %>%
      group_by(phenotype, threshold) %>%
      do(
        #ls = loess(partial~mean_values, data=.)
        ls = loess(relative_performance~median_values, data=.)
      )
    tb$ls_fit <-  apply(sapply(tb$ls,fitted),2,list)
    tb <- tb %>%
      unnest(cols=ls_fit) %>%
      unnest(cols=ls_fit) %>%
      #mutate(group_number = rep(1:40,10*5)) %>%
      mutate(group_number = rep(1:100,10*4)) %>%
      select(group_number, phenotype, threshold, ls_fit) %>%
      right_join(plot_df , by=c("phenotype","threshold","group_number")) %>%
      group_by(phenotype, threshold) %>%
      do(
        lm = lm(ls_fit~median_values,data=.)
      )
    r = c()
    b = c()
    for(i in 1:nrow(tb)){
      r = c(r, keep_decimal_places(round(summary(tb[[3]][[i]])$r.squared,3), 3))
      b = c(b, paste0(keep_decimal_places(round(summary(tb[[3]][[i]])$coef[2,c(1)],3), 3),
                      " (",
                      keep_decimal_places(round(summary(tb[[3]][[i]])$coef[2,c(4)],3), 3),
                      ")"))
    }
    tb <- tb %>%
      #mutate(`LM R2` = r, Coef = b) %>%
      cbind.data.frame(r2 = r, Coef = b) %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coef,r2) #%>%
    #rename("Phenotype"="phenotype") 
    colnames(tb)[2] = "Phenotype"
    
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
    colnames(tb)[colnames(tb) %in% "r2"] <- "r \u00b2"
    h <- -.5
  }
  
  tbs <- lapply(split(tb, tb$threshold), "[", -1)
  tbs <- tbs[-1]
  df <- tibble(x = rep(-Inf, length(tbs)), 
               y = rep(Inf, length(tbs)), 
               threshold = levels(as.factor(tb$threshold))[-1], 
               tbl = tbs) %>%
    filter(threshold !="5e-8")
  
  prs_names <- c(
    #`0-LDpred2` = "LDpred2",
    `1e-2` = "1e-2",
    `1e-3` = "1e-3",
    `1e-4` = "1e-4",
    `1e-5` = "1e-5"
  )
  
  graph_labels = plot_df %>%
    filter(threshold!="5e-8" & phenotype == trait) %>%
    filter(group_number == 1) %>%
    group_by(threshold) %>%
    select(threshold, eur_incremental_r2) %>%
    mutate(eur_incremental_r2 = 
             paste("absolute inc. r \u00b2 in FST \n group 1 =", 
                   keep_decimal_places(round(eur_incremental_r2, 3), 3)))
  
   plot <- plot + 
      guides(color = guide_legend(title="Trait"))+#,
             #size="none")+
      theme_light() + 
     theme(axis.title=element_text(size=16),
           axis.text=element_text(size=12),
           plot.title=element_text(size=28),
           plot.caption=element_text(size=12))+
 #     scale_color_brewer(palette="Set1")+
      xlab(xlabel1) +
      #ylab(expression(Partial~R^2))+
      ylab(expression(Incremental~r^2~Relative~to~FST~Bin~1))+
      labs(
        title=title1,
        #   subtitle=subtitle1,
           caption=caption1)+
     facet_grid(cols = vars(threshold), labeller = as_labeller(prs_names))+
     geom_table(data = df, aes(x = x, y = y, label = tbl),
                hjust = h, vjust = 1) +
     geom_hline(yintercept = 1, size=1) +
     geom_text(data = graph_labels, aes(x = 0.095, y = 1, label = eur_incremental_r2), size=4,
               color = "black") + 
     scale_y_continuous(labels=function(x) sprintf("%.1f", x))
     
  return(plot)
}

colors = c("#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4",
           "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC")
pheno = c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
          "Lymphocyte","WBC","Eosinophil")
for(i in 1:10){
  for(j in c(1, 2, 5)){
    comb_fst <- combined(as.data.frame(all_list[[1]]),"FST",j, pheno[i], colors[i], median_fst)
    ggsave(paste0('img/FINAL_COMB_FST_',j, "_",pheno[i],"_", nrow(median_fst),"_",
                  round(median_fst$count[median_fst$weighted_fst_groups == 1] / 
                          min(median_fst$count[median_fst$weighted_fst_groups != 1])),
                  ".png"),comb_fst,width=20,height=10,dpi=1000)
  }
}


# Don't need the scripts below for now
##################################################################

# files <- lapply(Sys.glob("data/theoretical/fst_*_h2.tsv"), read.table, header=T)
# 
# pheno_fst <- function(fst, prs, i){
#   prs = prs %>% filter(threshold == 4 & weighted_fst_groups == i)
#   prs = left_join(prs, fst, by = "phenotype")
# }
# 
# plot_pheno_fst <- function(pheno_fst, i){
#   plot <- ggplot(pheno_fst, aes(x=h2,y=incremental_r2)) +
#     geom_smooth(size=1,method = lm,se=F,na.rm=T, color = "black")
#   
#   title1 = paste("Phenotype Variation explained by Heritability for Fst Group", i)
#   xlabel1 = "Additive Heritability"
#   #subtitle1 = "h-squared was calculated for each trait."
#   #caption1 =  "Total of 40 groups each composed of 92-93 observations."
#   delta <- 0.05
#   beta <- "\u03b2 * 0.05 (p)"
#   
#   tb <- pheno_fst %>% 
#     do(
#       #lm = lm(partial~mean_values, data=.)
#       lm = lm(incremental_r2~h2, data=.)
#     ) %>% 
#     mutate(`LM R2` = round(summary(lm[[1]])$adj.r.squared,3)) %>%
#     mutate(Coef = paste0(round(summary(lm[[1]])$coef[2,c(1)]*delta,3),
#                          " (",
#                          round(summary(lm[[1]])$coef[2,c(4)],3),
#                          ")")) %>%
#     select(Coef,`LM R2`)
#   
#   colnames(tb)[colnames(tb) %in% "Coef"] <- beta
#   h <- -.4
#   
#   tbs <- list(tb)
#   df <- tibble(x = rep(-Inf, length(tbs)), 
#                y = rep(Inf, length(tbs)), 
#                tbl = tbs)
#   
#   plot <- plot + geom_point(size=3,color="white")+
#     geom_point(aes(color = phenotype), size=2.5,alpha=0.75)+
#     guides(color = guide_legend(title="Trait"),
#            size="none")+
#     theme_light() + 
#     theme(legend.text=element_text(size=14),
#           legend.title=element_text(size=16))+
#     #     scale_color_brewer(palette="Set1")+
#     xlab(xlabel1) +
#     #ylab(expression(Partial~R^2))+
#     ylab(expression(Incremental~R^2))+
#     labs(
#       title=title1)+#,
#       #   subtitle=subtitle1,
#     #  caption=caption1)+
#     geom_table(data = df, aes(x = x, y = y, label = tbl),
#                hjust = h, vjust = 1.5) 
#   
#   return(plot)
# }
# 
# for(i in 1:length(files)){
#   phenotype_fst <- pheno_fst(files[[i]], prs_df_fst, i)
#   ggsave(paste0('img/R_h2_FST_',i,"_40.png"),plot_pheno_fst(phenotype_fst, i),width=10,height=5,dpi=1000)
# }


