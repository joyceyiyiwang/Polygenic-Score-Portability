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
library(bigsnpr)
library(philentropy)


non_prs_df <- read_tsv('data/prs_comparisons2/non_prs_df.tsv')

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
  mutate(weighted_fst_groups = 1:n()) %>%
  mutate(group_type = "Fst") %>%
  rename("group_number" = "weighted_fst_groups") %>%
  rename("median" = "mean_fst") %>%
  select(group_type, group_number, median)


median_values <- non_prs_df %>% 
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
  mutate(WPC_groups = 1:n()) %>%
  mutate(group_type = "WPC") %>%
  rename("group_number" = "WPC_groups") %>%
  rename("median" = "mean_wpc") %>%
  select(group_type, group_number, median) %>%
  bind_rows(mean_fst_values)

#################################################################

prs_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_CT.tsv') %>%
  mutate(group_type="Fst") %>%
  filter(threshold!=0) %>%
mutate(threshold = as.factor(threshold)) %>%
select(-weighted_fst_groups)

ld_fst_df <- read_tsv('data/prs_comparisons2/UKBB_geno_fst_LD.tsv')%>%
  mutate(group_type="Fst") %>%
  filter(threshold==29) %>%
  mutate(threshold = factor(threshold,label="LDpred")) %>%
  select(-weighted_fst_groups)

prs_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_CT.tsv') %>%
  mutate(group_type="WPC") %>%
  filter(threshold!=0) %>%
  mutate(threshold = as.factor(threshold))%>%
  select(-WPC_groups)

ld_wpc_df <- read_tsv('data/prs_comparisons2/UKBB_geno_wpc_LD.tsv') %>%
  mutate(group_type="WPC") %>%
  filter(threshold==29) %>%
  mutate(threshold = factor(threshold,label="LDpred"))%>%
  select(-WPC_groups)

all_dfs <- prs_fst_df %>%
  bind_rows(ld_fst_df) %>%
  bind_rows(prs_wpc_df) %>%
  bind_rows(ld_wpc_df) %>%
  filter(phenotype %in% c("Height","Platelet","MCV")) %>%
  group_by(group_type,phenotype,threshold) %>%
  mutate(eur_partial = max((partial*(group_number==1)),na.rm=T),
         relative_performance = partial / eur_partial) %>%
  ungroup() %>%
  select(group_type,group_number,phenotype,threshold, relative_performance)


#################################################################

pop_med <- non_prs_df %>%
  group_by(population) %>%
  mutate(med_fst = median(Weighted_Fst,na.rm=T)) %>%
  mutate(med_wpc = median(PC_dist_weighted,na.rm=T)) %>%
  ungroup() %>%
  distinct(population, med_fst,med_wpc) %>%
  filter(population!="EUR") %>%
  pivot_longer(med_fst:med_wpc,names_to="median",values_to="pop_median") %>%
  mutate(median = factor(median,levels=c("med_fst","med_wpc"),
                         labels=c("Fst","WPC")))


######################

sumstats_ct <- read_csv("theor_sumstats_ct/sumstats_ct.txt",
                        col_names=c("GD","group_number","phenotype","threshold","sumstat")) %>%
  mutate(sumstat = parse_number(sumstat)) %>%
  mutate(threshold = as.character(threshold)) %>%
  mutate(GD=str_extract(GD,"[a-zA-Z]{3}"))

sumstats_height <- read_csv("theor_sumstats/height_ldpred_sumstats.txt",
                     col_names=c("GD","group_number","sumstat"))%>%
  mutate(sumstat=parse_number(sumstat)) %>%
  mutate(phenotype="Height") %>%
  mutate(GD=str_extract(GD,"[a-zA-Z]{3}"))

sumstats <- read_csv("theor_sumstats/platelet_mcv_ldpred_sumstats.txt",
                     col_names=c("GD","group_number","sumstat","phenotype")) %>%
  mutate(phenotype=str_extract(phenotype,"[a-zA-Z]+")) %>%
  mutate(GD = str_extract(GD,"[a-zA-Z]{3}")) %>%
  bind_rows(sumstats_height) %>%
  mutate(threshold="LDpred") %>%
  bind_rows(sumstats_ct) %>%
  mutate(GD = factor(GD,labels=c("Fst","WPC"))) %>%
  rename("group_type"="GD")

####################

ld <- read_csv("theor_ld/ld.txt",
               col_names=c("GD","group_number","LD"))%>%
  mutate(LD=parse_number(LD)) %>%
  mutate(GD = factor(GD,labels=c("Fst","WPC"))) %>%
  arrange(GD,group_number) %>%
  rename("group_type"="GD")


###################

herit_tibble <- tibble()

for (file in list.files(path = 'theor_herit',
                        pattern = '[a-z]{3}_[0-9]+_h2.tsv', full.names = T)) {
  
  group_type <- str_extract(file,"fst|wpc")
  add <- read_tsv(file) %>%
    mutate(group_type = group_type) %>%
    rename("group_number"="group") %>%
    select(group_type,group_number,phenotype,h2) 
  herit_tibble <- bind_rows(herit_tibble,add)
  
}
herit <- herit_tibble %>%
  mutate(group_type=ifelse(group_type=="fst","Fst","WPC")) %>%
  mutate(group_type = factor(group_type)) %>%
  group_by(group_type,phenotype) %>%
  mutate(h2_1 = max(h2*(group_number==1))) %>%
  mutate(rel_h2 = h2/h2_1) %>%
  ungroup() %>%
  mutate(group_number = as.numeric(group_number)) %>%
  filter(group_number!=1) %>%
  arrange(group_type,group_number,phenotype)

###################

pred <- sumstats %>%
  left_join(ld, by=c("group_type","group_number")) %>%
  left_join(herit, by=c("group_type","group_number","phenotype")) %>%
  left_join(median_values,by=c("group_type","group_number")) %>%
  mutate(eq_herit = LD*sumstat) %>%
  mutate(uneq_herit = rel_h2*LD*sumstat) %>%
  select(-LD,-sumstat,-h2,-h2_1,-rel_h2) %>%
  mutate(weight=1)


for(i in c("Fst","WPC")){
  for(j in c("Height","Platelet","MCV")){
    for(k in c("LDpred",as.character(1:4))){
      pred <- pred %>%
        add_row(group_type=i,group_number=1,phenotype=j,threshold=k,median=0,eq_herit=1,
                uneq_herit=1,weight=1000)
    }
  }
}

##########################


pred2 <- pred %>%
  left_join(all_dfs, by=c("group_type","group_number","phenotype","threshold")) %>%
  pivot_longer(cols=c("eq_herit","uneq_herit","relative_performance"),names_to="model",
               values_to="rel_pred") %>%
  mutate(model = factor(model,levels=c("eq_herit","uneq_herit","relative_performance"),
                        labels=c("Theor (eq h2)",
                                 "Theor (uneq h2)",
                                 "Observed"))) %>%
  mutate(threshold = factor(threshold,levels=c("LDpred",as.character(4:1)),
                            labels=c("LDpred","1e-2","1e-3","1e-4","1e-5"))) %>%
  mutate(phenotype = factor(phenotype, levels=c("Height","Platelet","MCV")))

#####################


eq <- pred2 %>%
  group_by(group_type,phenotype,threshold,model) %>%
  mutate(min_value = min(median)) %>%
  do(
    lm = lm(I(rel_pred-1)~0+I(median-min_value),data=.)
  ) %>%
  mutate(beta = coefficients((lm))[1] %>% round(5)) %>%
  mutate(p = summary(lm)$coefficients[1,4] %>% round(5)) %>%
  mutate(r2 = summary(lm)$adj.r.squared %>% round(3)) %>%
  ungroup() %>%
  select(-lm,-r2)

pred3 <- pred2 %>% 
  left_join(eq, by =c("group_type","phenotype","threshold","model")) %>%
  mutate(significance = ifelse(p<=0.05,"0.05","N.S.")) %>%
 # mutate(significance = ifelse(p<0.01,"0.01",significance)) %>%
  mutate(significance = factor(significance,levels=c("0.05","N.S."),labels=c("p<0.05","N.S."))) %>%
  rename("Model" = "model")

for (i in c("Fst","WPC")){

  if(i =="Fst") {
    x_lab <- "Within-Group Median Fst"
    ylim2 <- 1.75
    min_median <- 0
    }
  if(i =="WPC") {
    x_lab <- "Within-Group Median WPC Distance"
    ylim2 <- 0.25
    min_median <- 2.15
    }
  
graph_df <- pred3 %>%
  left_join(pop_med ,by=c("group_type"="median")) %>%
  filter(group_type==i)

graph <- graph_df %>%
  ggplot(aes(x=median,y=rel_pred-1,group=Model, color=Model))+
  geom_hline(yintercept=0,size=0.5)+
  geom_point(size=0.5)+
  geom_smooth(aes(linetype=significance),se=T,method="lm", level=0.95,
              formula=y~0+I(x-min_median))+
  scale_color_brewer(palette="Set1")+
  coord_cartesian(ylim=c(-1, ylim2))+
  scale_y_continuous(breaks=c(-1,-.5,0,.5,1,1.5,2),
                     labels = as.character(round(seq(0,3, by = 0.5),1)))+
#  scale_linetype_manual(values=c("solid", "longdash"))+
  facet_grid(phenotype~threshold,scales="free_x")+
  ylab("Relative Prediction Accuracy")+
  xlab(x_lab)+
  theme_bw()+
  guides(linetype=guide_legend(title=NULL,override.aes = list(color="black")))+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=16),
#        strip.background=element_rect(fill="white"),
        strip.text = element_text(color="black",size=10),
        axis.text.x = element_text(angle = 45, vjust=0.75,size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=14))+
  labs(caption="The gray regions represent the 95% confidence interval of the fitted line from the theoretical predictions assuming unequal heritability.")

name <- paste0("img1/V2_theor_",i,".png")
ggsave(name,graph,dpi=1000,width=12,height=10)

}





