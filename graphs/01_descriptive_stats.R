.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
#.libPaths("/moto/palab/users/jm4454/rpackages/")
library(cowplot)
library(ggrepel)
library(ggplot2)
library(ggpmisc)
library(plyr)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(philentropy)
library(qqman)
library(boot)

non_prs_df <- read_tsv('data/prs_comparisons2/non_prs_df.tsv')
prs_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_CT.tsv')
ld_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_LD.tsv')
#prs_pc_df <- read_tsv('data/prs_comparisons1/UKBB_geno_pc_CT.tsv')
#ld_pc_df <- read_tsv('data/prs_comparisons1/UKBB_geno_pc_LD.tsv')
prs_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_CT.tsv')
ld_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_LD.tsv')
prs_pop_df <- read_tsv('data/prs_comparisons2/UKBB_geno_pop.tsv')

##################################################################


mean_fst_values <- non_prs_df %>% 
  filter(phenotype=="Height") %>%
  select(Weighted_Fst,weighted_fst_groups) %>%
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
  select(PC_dist_weighted,WPC_groups) %>%
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

###########################################################


non_prs_df_sub <- non_prs_df %>%
  left_join(mean_fst_values,by="weighted_fst_groups") %>%
  left_join(mean_wpc_values,by="WPC_groups") %>%
  mutate(`Fst Group`= mean_fst) %>%
  rename("WPC Group"= "WPC_groups") %>%
  mutate(phenotype=factor(phenotype,levels=c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                                             "Lymphocyte","WBC","Eosinophil"))) %>%
  select(`Fst Group`, mean_fst,
         `WPC Group`, mean_wpc,
         phenotype,phenotype_value)


mean_sample <- function(d,i){
  d2 <- d[i,] %>% pull(phenotype_value)
  return(mean(d2,na.rm=T))
}


var_sample <- function(d, i){
  d2 <- d[i,] %>% pull(phenotype_value)
  return(var(d2,na.rm=T))
}

extract_cis <- function(boot_obj){
  values <- boot_obj$percent[c(4,5)]
  return(paste0(values[1],",",values[2]))
}

set.seed(626)
df_boot_fst <- non_prs_df_sub %>%
  group_by(`Fst Group`, phenotype) %>%
  do(
    boot_mean = boot(.,mean_sample,R=500),
    boot_var = boot(.,var_sample,R=500)
  ) %>%
  ungroup() %>%
  mutate(mean_ci = unlist(lapply(lapply(boot_mean,boot.ci,type="perc"),extract_cis))) %>%
  mutate(var_ci = unlist(lapply(lapply(boot_var,boot.ci,type="perc"),extract_cis))) %>%
  select(-boot_mean,-boot_var) %>%
  pivot_longer(mean_ci:var_ci,names_to="stat",values_to="ci") %>%
  separate(ci, into = c("Min", "Max"), sep = ",") %>%
  mutate(stat = factor(stat,levels=c("mean_ci","var_ci"),labels=c("Mean","Variance"))) %>%
  mutate_if(is.character,as.numeric) %>%
  rename("Group"=`Fst Group`)

set.seed(626)
df_boot_wpc <- non_prs_df_sub %>%
  group_by(`WPC Group`, phenotype) %>%
  do(
    boot_mean = boot(.,mean_sample,R=500),
    boot_var = boot(.,var_sample,R=500)
  ) %>%
  ungroup() %>%
  mutate(mean_ci = unlist(lapply(lapply(boot_mean,boot.ci,type="perc"),extract_cis))) %>%
  mutate(var_ci = unlist(lapply(lapply(boot_var,boot.ci,type="perc"),extract_cis))) %>%
  select(-boot_mean,-boot_var) %>%
  pivot_longer(mean_ci:var_ci,names_to="stat",values_to="ci") %>%
  separate(ci, into = c("Min", "Max"), sep = ",") %>%
  mutate(stat = factor(stat,levels=c("mean_ci","var_ci"),labels=c("Mean","Variance"))) %>%
  mutate_if(is.character,as.numeric)%>%
  rename("Group"=`WPC Group`)


######################################################

non_prs_df_sub <- non_prs_df_sub %>%
  group_by(`Fst Group`,phenotype) %>%
  mutate(fst_pheno_mean = mean(phenotype_value,na.rm=T)) %>%
  mutate(fst_pheno_var = var(phenotype_value,na.rm=T)) %>%
  ungroup() %>%
  group_by(`WPC Group`,phenotype) %>%
  mutate(wpc_pheno_mean = mean(phenotype_value,na.rm=T)) %>%
  mutate(wpc_pheno_var = var(phenotype_value,na.rm=T)) %>%
  ungroup() %>%
  select(`Fst Group`,mean_fst, phenotype, fst_pheno_mean,fst_pheno_var,
         `WPC Group`,mean_wpc,phenotype,wpc_pheno_mean,wpc_pheno_var) 

top5 <- c("Height","Platelet","MCV","MCH","BMI")

fst_stats <- non_prs_df_sub %>%
  distinct(`Fst Group`,phenotype,.keep_all=T) %>%
  select(`Fst Group`,mean_fst,phenotype, fst_pheno_mean,fst_pheno_var) %>%
  pivot_longer(fst_pheno_mean:fst_pheno_var,names_to="stat",values_to="value") %>%
  mutate(stat = ifelse(stat=="fst_pheno_mean","Mean","Variance")) %>%
  mutate(stat = factor(stat,labels=c("Mean","Variance"))) %>%
  mutate(top5 = ifelse((phenotype %in% top5),"Top 5","Bottom 5")) %>%
  mutate(GD = "Fst") %>%
  rename("Group" = "Fst Group") %>%
  rename("Median" = "mean_fst") %>%
  left_join(df_boot_fst,by=c("Group","phenotype","stat"))

wpc_stats <- non_prs_df_sub %>%
  distinct(`WPC Group`,phenotype,.keep_all=T) %>%
  select(`WPC Group`,mean_wpc,phenotype, wpc_pheno_mean,wpc_pheno_var) %>%
  pivot_longer(wpc_pheno_mean:wpc_pheno_var,names_to="stat",values_to="value") %>%
  mutate(stat = ifelse(stat=="wpc_pheno_mean","Mean","Variance")) %>%
  mutate(stat = factor(stat,labels=c("Mean","Variance"))) %>%
  mutate(top5 = ifelse((phenotype %in% top5),"Top 5","Bottom 5")) %>%
  mutate(GD = "WPC") %>%
  rename("Group" = "WPC Group") %>%
  rename("Median" = "mean_wpc")%>%
  left_join(df_boot_wpc,by=c("Group","phenotype","stat"))

summary_stats <- fst_stats %>%
  bind_rows(wpc_stats)
summary_stats %>% write_tsv("data/prs_comparisons2/summary_stats.tsv")

for(i in c("Mean","Variance")){

  if(i=="Mean"){
    yint <- 0
    ylabel <- "Phenotype Mean"
    ymin1 <- -1.5
    ymax1 <- .75
    }
  else{
    yint <- 1
    ylabel <- "Phenotype Variance"
    ymin1 <- 0
    ymax1 <- 3
    }
  
pheno_graph <- fst_stats %>%
  bind_rows(wpc_stats) %>%
  mutate(top5 = factor(top5,levels=c("Top 5","Bottom 5"))) %>%
  filter(stat ==i) %>%
  ggplot(aes(x=Median,y=value,group=phenotype,color=phenotype)) +
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Min,ymax=Max),size=0.5, alpha=0.4)+
  geom_point(size=.5,color="black")+
#  geom_point(size=2, color="black")+
  geom_hline(yintercept=yint, linetype="dashed") +
  theme_classic()+
  ylab(ylabel)+
  coord_cartesian(ylim=c(ymin1, ymax1))+
  xlab("Within-Group Median Genetic Distance")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        legend.title=element_text(size=22), 
        legend.text=element_text(size=18),
        strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20)
        )+
  scale_color_discrete(name="Trait")+
#scale_color_brewer(palette="Paired")+
 # scale_color_discrete(name = "Phenotype")+
  facet_wrap(GD~top5, scales="free_x")+
  guides(color=guide_legend(override.aes = list(size=12)))

ggsave(paste0("img1/FINAL_combined_",i,"_errors1.png"),pheno_graph,
       width=20,height=16,dpi=1250)
}  

####################################################
r <- non_prs_df_sub %>%
  distinct(`Fst Group`,phenotype,.keep_all=T) %>%
  filter(phenotype %in% c("Height","Platelet","MCV","MCH","BMI"))
q <-  non_prs_df_sub %>%
  distinct(`Fst Group`,phenotype,.keep_all=T) %>%
  filter(!(phenotype %in% c("Height","Platelet","MCV","MCH","BMI")))

fst_mean <- function(df,i){
  
fst <- df %>%
  ggplot(aes(x=mean_fst,y=fst_pheno_mean,group=phenotype,color=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic()+
  ylab("Phenotype Mean")+
  xlab("Within-Group Median Fst")+
  scale_x_continuous(breaks = round(seq(0,
                                        0.13, by = 0.02),2)) +
  scale_y_continuous(breaks = round(seq(-1.5,
                                        1.5, by = .5),1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))+
  scale_color_discrete(name = "Phenotype")
ggsave(paste0("img1/FINAL_",i,"_mean_fst.png"),fst,
       width=16,height=12,dpi=1000)
}

fst_mean(r,"top5")
fst_mean(q,"bottom5")


fst_var <- function(df,i){
  
fst <- df %>%
  ggplot(aes(x=mean_fst,y=fst_pheno_var,group=phenotype,color=phenotype)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept=1, linetype="dashed") +
  theme_classic()+
  ylab("Phenotype Variance")+
  xlab("Within-Group Median Fst")+
  scale_x_continuous(breaks = round(seq(0,
                                        0.13, by = 0.02),2)) +
  scale_y_continuous(breaks = round(seq(0,
                                        3, by = .5),1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))+
  scale_color_discrete(name = "Phenotype")

ggsave(paste0("img1/FINAL_",i,"_var_fst.png"),fst,
       width=16,height=12,dpi=1000)
}
fst_var(r,"top5")
fst_var(q,"bottom5")

#########################################################

r <- non_prs_df_sub %>%
  distinct(`WPC Group`,phenotype,.keep_all=T) %>%
  filter(phenotype %in% c("Height","Platelet","MCV","MCH","BMI"))
q <-  non_prs_df_sub %>%
  distinct(`WPC Group`,phenotype,.keep_all=T) %>%
  filter(!(phenotype %in% c("Height","Platelet","MCV","MCH","BMI")))

wpc_mean <- function(df,i){
  
  fst <- df %>%
    ggplot(aes(x=mean_wpc,y=wpc_pheno_mean,group=phenotype,color=phenotype)) +
    geom_line()+
    geom_point()+
    geom_hline(yintercept=0, linetype="dashed") +
    theme_classic()+
    ylab("Phenotype Mean")+
    xlab("Within-Group Median WPC Distance")+
    scale_x_continuous(breaks = round(seq(0,
                                          140, by = 25),2)) +
    scale_y_continuous(breaks = round(seq(-1.5,
                                          1.5, by = .5),1))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16,face="bold"),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=14))+
    scale_color_discrete(name = "Phenotype")
  ggsave(paste0("img1/FINAL_",i,"_mean_wpc.png"),fst,
         width=16,height=12,dpi=1000)
}

wpc_mean(r,"top5")
wpc_mean(q,"bottom5")


wpc_var <- function(df,i){
  
  fst <- df %>%
    ggplot(aes(x=mean_wpc,y=wpc_pheno_var,group=phenotype,color=phenotype)) +
    geom_line()+
    geom_point()+
    geom_hline(yintercept=1, linetype="dashed") +
    theme_classic()+
    ylab("Phenotype Variance")+
    xlab("Within-Group Median WPC Distance")+
    scale_x_continuous(breaks = round(seq(0,
                                          140, by = 25),2)) +
    scale_y_continuous(breaks = round(seq(0,
                                          3, by = .5),1))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16,face="bold"),
          legend.title=element_text(size=16), 
          legend.text=element_text(size=14))+
    scale_color_discrete(name = "Phenotype")
  
  ggsave(paste0("img1/FINAL_",i,"_var_wpc.png"),fst,
         width=16,height=12,dpi=1000)
}
wpc_var(r,"top5")
wpc_var(q,"bottom5")

##########################################################

scatter_plot <- non_prs_df %>%
  filter(phenotype=="BMI") %>%
  ggplot() +
  geom_point(aes(x=Weighted_Fst,y=PC_dist_weighted,color=population)) +
  theme_classic()+
  xlab("Weighted Fst") +
  ylab("Weighted PC Distance")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22,face="bold"),
        legend.title=element_text(size=22), 
        legend.text=element_text(size=18),
        legend.position=c(.825,.25))+
  guides(color=guide_legend(override.aes = list(size=14,shape="circle")))+
  scale_x_continuous(breaks = round(seq(min(non_prs_df$Weighted_Fst),
                                        max(non_prs_df$Weighted_Fst), by = 0.02),2)) +
  scale_y_continuous(breaks = round(seq(0,
                                        max(non_prs_df$PC_dist_weighted), by = 25),1))+
  scale_color_discrete(name="Population")

ggsave("img1/FINAL_fst_wpc.png",scatter_plot,
       width=16,height=12,dpi=1000)

#####################################################

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
  mutate(Population=factor(Population,levels=c("EUR_train","EUR_test","AMR","SAS","EAS","AFR"),
                           labels=c("EUR (train)","EUR (test)","AMR","SAS","EAS","AFR"))) %>%
  arrange(Population)


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
             size=4,shape="square",color="black",stroke=2)+
  geom_point(data=medians,  mapping=aes(x = PC1, y = PC2,fill=Population),
             size=3.5,shape="square")+
  scale_colour_hue(l = 70, c = 150)+
  guides(color=guide_legend(override.aes = list(size=8,shape="circle")),
         fill=F)+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=22,face="bold"),
        legend.title=element_text(size=22), 
        legend.text=element_text(size=18),
        legend.position = c(0.1,0.8))

ggsave("img1/FINAL_PC1_PC2_all.png",sp,
       width=16,height=12,dpi=1000)

######################################################

height_gwas <- tibble()
for (file in list.files(path = 'data/gwas_results1',
                        pattern = 'Height*.glm.linear', full.names = T)) {
  
  add <- read_tsv(file)
  height_gwas <- height_gwas %>% bind_rows(add)
  
}

height_gwas <- height_gwas %>%
  select(ID,`#CHROM`,POS,P) %>%
  rename("SNP"="ID") %>%
  rename("CHR"="#CHROM") %>%
  rename("BP"="POS")

png("img1/Height_manhattan.png",width=1020,height=800)
manhattan(height_gwas, ylim = c(0, 15), cex = 0.6, cex.axis = 0.9, 
          col = c("purple", "orange"))
dev.off()


#######################################

pop_df <- prs_pop_df %>%
  filter(threshold!=0) %>%
  mutate(threshold=factor(threshold,levels=4:1,labels=c("1e-2","1e-3","1e-4","1e-5"))) %>%
  group_by(population,phenotype) %>%
  mutate(max_partial=max(partial)) %>% 
  ungroup() %>%
  group_by(phenotype) %>%
  mutate(eur_partial = max(partial*(population=="EUR"))) %>%
  mutate(rel_performance = max_partial/eur_partial) %>%
  ungroup() %>%
  group_by(population,phenotype) %>%
  top_n(1,partial)

pop_df <- prs_pop_df %>%
  filter(threshold!=0) %>%
  mutate(threshold=factor(threshold,levels=4:1,labels=c("1e-2","1e-3","1e-4","1e-5"))) %>%
  group_by(phenotype,threshold) %>%
  mutate(eur_partial = max(partial*(population=="EUR"))) %>%
  mutate(rel_performance = partial/eur_partial) %>%
  ungroup() %>%
  select(-partial,-nested_r2,-full_r2,-incremental_r2)


pop_df_se <- pop_df %>%
  ungroup() %>%
  group_by(population) %>%
  dplyr::mutate(mean_rel_eur=mean(rel_performance), sem_rel_eur=sd(rel_performance)/sqrt(10)) %>%
  subset(phenotype=='BMI') %>%
  arrange(desc(mean_rel_eur)) %>%
  select(population, mean_rel_eur,sem_rel_eur)

hex_codes1 <- c("black",hue_pal()(4))
hex_codes1   

pop_violin <- pop_df %>% 
  mutate(population=factor(population,levels=c("EUR","AMR","SAS","EAS","AFR"))) %>%
  ggplot(aes(x=population, fill=population, color=population)) +
  geom_violin(aes(y=rel_performance))+
  geom_point(aes(y=rel_performance,shape=threshold),color='black',size=4) +
  geom_errorbar(data=pop_df_se %>% filter(population!="EUR"), aes(ymin=mean_rel_eur - sem_rel_eur,
                                    ymax=mean_rel_eur + sem_rel_eur),
                color='black', width=0.15,size=1.1) +
  labs(x='Population', y='Prediction accuracy\n(relative to test Europeans)') +
  theme_bw() +
  scale_fill_manual(values=hex_codes1, name='Population') +
  scale_color_manual(values=hex_codes1, name='Population') +
  guides(fill=F, color=F,shape=guide_legend(title="C+T Cutoff")) +
  theme(axis.text = element_text(color='black', size=16),
        text = element_text(size=18),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom")
ggsave("img1/FINAL_population_violin.png",pop_violin,width=15,height=12,dpi=1000)

############################################


pop_df <- prs_pop_df %>%
  filter(threshold!=0) %>%
  mutate(threshold=factor(threshold,levels=4:1,labels=c("1e-2","1e-3","1e-4","1e-5"))) %>%
  group_by(phenotype,threshold) %>%
  mutate(eur_partial = max(partial*(population=="EUR"))) %>%
  mutate(rel_performance = partial/eur_partial) %>%
  ungroup() %>%
  select(-partial,-nested_r2,-full_r2,-incremental_r2)


pop_df_se <- pop_df %>%
  ungroup() %>%
  group_by(population,threshold) %>%
  dplyr::mutate(mean_rel_eur=mean(rel_performance), sem_rel_eur=sd(rel_performance)/sqrt(10)) %>%
  subset(phenotype=='BMI') %>%
  arrange(desc(mean_rel_eur)) %>%
  select(population, threshold,mean_rel_eur,sem_rel_eur)

hex_codes1 <- c("black",hue_pal()(4))
hex_codes1   

pop_violin <- pop_df %>% 
  mutate(population=factor(population,levels=c("EUR","AMR","SAS","EAS","AFR"))) %>%
  ggplot(aes(x=population, fill=population, color=population)) +
  geom_violin(aes(y=rel_performance))+
  geom_point(aes(y=rel_performance),color='black',size=2) +
  geom_errorbar(data=pop_df_se %>% filter(population!="EUR"), aes(ymin=mean_rel_eur - sem_rel_eur,
                                                                  ymax=mean_rel_eur + sem_rel_eur),
                color='black', width=0.15,size=1.1) +
  labs(x='Population', y='Prediction accuracy\n(relative to test Europeans)') +
  theme_bw() +
  scale_fill_manual(values=hex_codes1, name='Population') +
  scale_color_manual(values=hex_codes1, name='Population') +
  guides(fill=F, color=F,shape=guide_legend(title="C+T Cutoff")) +
  theme(axis.text = element_text(color='black', size=16),
        text = element_text(size=18),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom")+
  facet_wrap(~threshold,ncol=2)
ggsave("img1/FINAL_population_violin_panel.png",pop_violin,width=15,height=12,dpi=1000)


