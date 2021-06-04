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
library(philentropy)
library(khroma)
library(viridis)

non_prs_df <- read_tsv('data/prs_comparisons2/non_prs_df.tsv')
prs_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_CT.tsv')
ld_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_LD.tsv')
prs_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_CT.tsv')
ld_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_LD.tsv')

#############################################

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

#########################################

all_list <- list(prs_fst_df,ld_fst_df,
                 #prs_pc_df,ld_pc_df,
                 prs_wpc_df,ld_wpc_df)

pheno_factor <- function(df){
  mutate(df,phenotype = factor(phenotype, levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                                   "Lymphocyte","WBC","Eosinophil")))
  mutate(df, Heritability = as.numeric(as.character(factor(phenotype, levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                                       "Lymphocyte","WBC","Eosinophil"),
                                   labels=c(0.485,0.308,0.267,0.253,0.248,0.234,0.230,0.210,0.191,0.184)))))
}
all_list <- lapply(all_list,pheno_factor) 

sunset <- colour("sunset")
######################################
combined <- function(prs_df,ldpred_prs_df, type, line){
  
  
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
    title1 = "Phenotype Variation explained by Polygenic Scores Across Fst Pools"
    xlabel1 = "Within-Group Median Fst"
    subtitle1 = "Fst was calculated between each individual against EUR training population."
    caption1 =  "All models are associated with significant betas (p<0.05)."
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
    caption1 = "All models are associated with significant slopes (p<0.05)."
    delta <- 50
    beta <- "\u03b2 * 50 (p)"
  }
  
  tb <- plot_df %>% group_by(phenotype, threshold) %>%
    do(
      lm = lm(partial~mean_values, data=.)
    ) %>%
    mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
    mutate(Coef = paste0(round(summary(lm)$coef[2,c(1)]*delta,3),
                         " (",
                         round(summary(lm)$coef[2,c(4)],3),
                         ")")) %>%
    mutate(p = round(summary(lm)$coef[2,4],3)) %>%
    ungroup() %>%
    distinct(phenotype,threshold,.keep_all=T) %>%
    select(threshold,phenotype,p) #%>%
 #   rename("Phenotype"="phenotype") 
  
  plot_df <- plot_df %>%
    left_join(tb,by=c("phenotype","threshold")) %>%
    mutate(Significance = "N.S.") %>%
    mutate(Significance = ifelse(p<=0.05,"p<0.05",Significance)) %>%
    mutate(Heritability = ifelse(Heritability >= 0.30,0.30, Heritability)) %>%
    mutate(Heritability = ifelse(Heritability <= 0.20,0.20, Heritability))
  
  plot <-  ggplot(plot_df, aes(x=mean_values,y=partial, #or partial
                               group=phenotype,
                               color=Heritability))+
    geom_smooth(aes(linetype=Significance),size=1,method = lm,se=F,na.rm=T)
  
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
      mutate(p = round(summary(lm)$coef[2,4],3)) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,p) %>%
      rename("Phenotype"="phenotype") 
    
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
    h <- -.6
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
    guides(#color = guide_legend(title="Trait Heritability"),
           size=F,
           color=guide_colourbar(barwidth = 0.8, barheight = 7),
           linetype=F)+
    theme_bw() + 
    theme(legend.text=element_text(size=9),
          legend.title=element_text(size=14),
          strip.background=element_rect(fill="white"),
          strip.text = element_text(color="black",size=10))+ 
 #   scale_color_gradientn(colours = rev(sunset(9)))+
    scale_color_viridis(option = "D", breaks = c(0.20, 0.25, 0.30))+
    xlab(xlabel1) +
    ylab(expression(Partial~R^2))+
    coord_cartesian(ylim=c(0, 0.15))+
    theme(axis.text.x = element_text(angle = 45, vjust=.8, size=10),
          axis.text.y = element_text(size=10),
          axis.title=element_text(size=14))+
    labs(caption=caption1)+
    facet_grid(cols = vars(threshold), labeller = as_labeller(prs_names))
#+
#    geom_table(data = df, aes(x = x, y = y, label = tbl),
#               hjust = h, vjust = 1)
  
  
  
  return(plot)
}


for (i in c(2)){
  comb_fst <- combined(as.data.frame(all_list[1]),as.data.frame(all_list[2]), "FST",i)
  ggsave(paste0('img1/V4_COMB_FST',i,"_40.png"),comb_fst,width=20,height=10,dpi=1000)
  
  comb_wpc <- combined(as.data.frame(all_list[3]),as.data.frame(all_list[4]), "WPC",i)
  ggsave(paste0('img1/V4_COMB_WPC',i,"_40.png"),comb_wpc,width=20,height=10,dpi=1000)
}


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
    caption1 ="Dashed lines are associated with insignificant slopes (p>0.05)"
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
    caption1 = "Dashed lines are associated with insignificant slopes (p>0.05)"
    delta=50
    beta <- "\u03b2 * 50 (p)"
  }
  
  min_value <- min(plot_df$mean_values)
  max_value <- max(plot_df$mean_values)
  
  plot_df$weight <- ifelse(plot_df$group_number==1,100000,1)
  
  
  tb <- plot_df %>% group_by(phenotype, threshold) %>%
    do(
      lm = lm(I(relative_performance-1)~0+I(mean_values-min_value), data=.)
    ) %>%
    mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
    mutate(Coef = paste0(round(summary(lm)$coef[1,c(1)]*delta,3),
                         " (",
                         round(summary(lm)$coef[1,c(4)],3),
                         ")")) %>%
    mutate(p = round(summary(lm)$coef[1,4],3)) %>%
    ungroup() %>%
    distinct(phenotype,threshold,.keep_all=T) %>%
    select(threshold,phenotype,p)
  
  
  plot_df <- plot_df %>%
    left_join(tb,by=c("phenotype","threshold")) %>%
    mutate(Significance = "N.S.") %>%
    mutate(Significance = ifelse(p<=0.05,"p<0.05",Significance)) %>%
    mutate(Significance = factor(Significance,levels=c("p<0.05","N.S.")))%>%
    mutate(Heritability = ifelse(Heritability >= 0.30,0.30, Heritability)) %>%
    mutate(Heritability = ifelse(Heritability <= 0.20,0.20, Heritability))
  
  plot <- ggplot(plot_df, aes(x=mean_values-min_value,y=relative_performance-1,
                              color=Heritability,group=phenotype))+
    geom_hline(yintercept=0, size=0.5) +
    geom_smooth(aes(linetype=Significance),size=1,method = lm,formula=y ~ 0+x ,se=F,na.rm=T)
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
      mutate(p = round(summary(lm)$coef[1,4],3)) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,p) %>%
      rename("Trait"="phenotype") %>%
      mutate(Trait = factor(Trait,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                                   "Lymphocyte","WBC","Eosinophil"))) %>%
      arrange(threshold,Trait)
    colnames(tb)[colnames(tb) %in% "Coef"] <- beta
    #   colnames(tb)[colnames(tb) %in% "LM R2"] <- expression(LM~R2)
    h <- -.9
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
    h <- -.9
    
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
    theme_bw() +
    theme(legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          strip.background=element_rect(fill="white"),
          strip.text = element_text(color="black",size=10),
          axis.text.x = element_text(angle = 45, vjust=.8, size=10),
          axis.text.y = element_text(size=10),
          axis.title=element_text(size=14))+
    xlab(xlabel1) +
    scale_color_viridis(option = "D", breaks = c(0.20, 0.25, 0.30))+
    scale_linetype_manual(values=c("solid", "longdash"))+
    guides(
      size=F,
      color=guide_colourbar(barwidth = 0.8, barheight = 7),
      linetype=F)+
    ylab(expression(Relative~Partial~R^2))+
    labs(caption=caption1)+
    coord_cartesian(ylim=c(-1, 0.5))+
    scale_y_continuous(breaks=c(-1,-.5,0,.5,1),labels = as.character(round(seq(0,2, by = 0.5),1)))+
    scale_x_continuous(breaks=seq(0,max_value-min_value,length.out = 6),
                       labels = as.character(round(seq(min_value,max_value, length.out = 6),2)))+
    facet_grid(cols = vars(threshold), labeller = as_labeller(prs_names))
  
  return(plot)
}

for (i in c(2)){
  comb_fst <- combined_rc(all_list[[1]],all_list[[2]], "FST",i)
  ggsave(paste0('img1/V4_RC_COMB_FST',i,"_40.png"),comb_fst,width=20,height=10,dpi=1000)
  
  comb_wpc <- combined_rc(as.data.frame(all_list[3]),as.data.frame(all_list[4]), "WPC",i)
  ggsave(paste0('img1/V4_RC_COMB_WPC',i,"_40.png"),comb_wpc,width=20,height=10,dpi=1000)
}

#################################################


combined_rc_df <- function(prs_df,ldpred_prs_df, type, line){
 
  
  ld_df <- ldpred_prs_df %>%
    filter(threshold==29) 
  
  plot_df <- rbind(prs_df,ld_df)
  
  plot_df$threshold <- factor(plot_df$threshold,levels=c(29,4,3,2,1,0),
                              labels=c("LDpred2","1e-2","1e-3","1e-4","1e-5","5e-8"),
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
    
  }
  
  if (type=="WPC") {
    plot_df <- plot_df %>% left_join(mean_wpc_values,by="WPC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_wpc)
  }
  
  min_value <- min(plot_df$mean_values)
  max_value <- max(plot_df$mean_values)
  mid_value <- (max_value-min_value)/2
  plot_df$weight <- ifelse(plot_df$group_number==1,100000,1)
  
  
  if(line ==1){
    
    tb <- plot_df %>% group_by(phenotype, threshold) %>%
      do(
        lm = lm(I(relative_performance-1)~0+I(mean_values-min_value), data=.)
      ) %>%
      mutate(`LM R2` = round(summary(lm)$adj.r.squared,3)) %>%
      mutate(Coef = paste0(round(summary(lm)$coef[1,c(1)],3),
                           " (",
                           round(summary(lm)$coef[1,c(4)],3),
                           ")")) %>%
      mutate(`Coefficient` = summary(lm)$coef[1,c(1)]) %>%
      mutate(`p-val` = summary(lm)$coef[1,c(4)]) %>%
      ungroup() %>%
      distinct(phenotype,threshold,.keep_all=T) %>%
      select(threshold,phenotype,Coefficient,`p-val`,`LM R2`) %>%
      rename("Trait"=`phenotype`)
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
for (i in c(1)){
  comb_fst <- combined_rc_df(as.data.frame(all_list[1]),as.data.frame(all_list[2]), "FST",i)
  comb_wpc <- combined_rc_df(as.data.frame(all_list[3]),as.data.frame(all_list[4]), "WPC",i)
  rc_comb1 <- bind_rows(comb_fst,comb_wpc)
  rc_comb <- rbind(rc_comb,rc_comb1)
  
}
rc_comb$gd <- c(rep("FST",50), rep("WPC", 50))

rc_betas <- rc_comb %>%
  mutate(Heritability = as.numeric(as.character(factor(Trait, levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                                                               "Lymphocyte","WBC","Eosinophil"),
                                                           labels=c(0.485,0.308,0.267,0.253,0.248,0.234,0.230,0.210,0.191,0.184)))))

rc_betas_graph <- rc_betas %>%
  ggplot(aes(x=Heritability,y=Coefficient,size=-log10(`p-val`),color=Heritability)) +
  geom_point()+
  scale_color_gradientn(colours = rev(sunset(9)))+
  guides(color=F,size=F)+
  theme_bw() +
  theme(
        strip.background=element_rect(fill="white",size=1.5,color="black"),
        strip.text = element_text(color="black",size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  ylab("Slope")+
  xlab("Trait Heritability")+
  facet_wrap(~gd, scales = 'free_y')+
  scale_x_continuous(breaks = c(0.2,0.3,0.4,0.5), limits=c(0.15,0.5))
  #labs(caption="Slopes computed from a linear regression between relative partial R2 and within-group median genetic distance.")

ggsave("img1/rc_betas.png",rc_betas_graph,width=10,height=8,dpi=500)

