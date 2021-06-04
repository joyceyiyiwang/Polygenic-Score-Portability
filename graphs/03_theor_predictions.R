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

sumstats_ct <- read_csv("data/theoretical/ct_sumstats.txt",
                        col_names=c("GD","group_number","phenotype","threshold","sumstat")) %>%
  mutate(sumstat = parse_number(sumstat)) %>%
  mutate(threshold = as.character(threshold)) %>%
  mutate(GD=str_extract(GD,"[a-zA-Z]{3}"))

#sumstats_height <- read_csv("data/theoretical/height_ldpred_sumstats.txt",
#                     col_names=c("GD","group_number","sumstat"))%>%
#  mutate(sumstat=parse_number(sumstat)) %>%
#  mutate(phenotype="Height") %>%
#  mutate(GD=str_extract(GD,"[a-zA-Z]{3}"))

sumstats <- read_csv("data/theoretical/ldpred_sumstats.txt",
                     col_names=c("GD","group_number","sumstat","phenotype")) %>%
  mutate(phenotype=str_extract(phenotype,"[a-zA-Z]+")) %>%
  mutate(GD = str_extract(GD,"[a-zA-Z]{3}")) %>%
 # bind_rows(sumstats_height) %>%
  mutate(threshold="LDpred") %>%
  bind_rows(sumstats_ct) %>%
  mutate(GD = factor(GD,labels=c("Fst","WPC"))) %>%
  rename("group_type"="GD")

####################

ld <- read_csv("data/theoretical/ld.txt",
               col_names=c("GD","group_number","LD"))%>%
  mutate(LD=parse_number(LD)) %>%
  mutate(GD = factor(GD,labels=c("Fst","WPC"))) %>%
  arrange(GD,group_number) %>%
  rename("group_type"="GD")

components <- function(){
#Need to extract num and den from output files
ld_components <- read_csv("theor_ld/ld_components.txt",
               col_names=c("GD","group_number","LD","num","den"))%>%
  mutate(GD = factor(GD,labels=c("Fst","WPC"))) %>%
  mutate(den = parse_number(den)) %>%
  arrange(GD,group_number) %>%
  rename("group_type"="GD") %>%
  left_join(median_values, by=c("group_type","group_number"))

ld_comp <- ld_components %>%
  mutate(clr = ifelse(group_type=="Fst","red","blue")) %>%
  ggplot(aes(x=median,y=num,color=group_type)) +
  geom_line(size=1)+
  geom_point(size=1.5,color="black")+
  geom_hline(aes(yintercept=den), linetype="dashed")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"),
        strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)
  )+
  xlab("Within-Group Median Genetic Distance")+
  ylab("Weighted LD Correlation")+
  facet_wrap(~group_type,scales="free_x")+
  scale_color_brewer(palette="Set1")+
  guides(color=F)
  
ggsave("img1/ld_comp_numerator.png",ld_comp, width=8,height=6)
}

###################

herit_tibble <- tibble()

for (file in list.files(path = 'data/theoretical/',
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

sumstats_fst <- sumstats %>% filter(group_type=="Fst") %>%left_join(mean_fst_values,by="group_number")
s_fst <- sumstats_fst %>%
  ggplot(aes(x=median,y=sumstat,color=threshold))+
  geom_point(size=1)+
  xlab("Within-Group Median Fst")+
  ylab("Term iv")+
  geom_vline(xintercept=quantile(sumstats_fst$median,probs=c(0.25,0.5,0.75)), linetype="dashed")+
  facet_wrap(~phenotype)

ggsave("img1/theor_term_iv.png",s_fst)

ld_fst <- ld %>% filter(group_type=="Fst") %>% select(-group_type) %>% left_join(mean_fst_values,by="group_number")
ld_fst1 <- ld_fst %>%
  ggplot(aes(x=median,y=LD))+
  geom_point(size=1)+
  xlab("Within-Group Median Fst")+
  ylab("Term iii")+
  geom_vline(xintercept=quantile(ld_fst$median,probs=c(0.25,0.5,0.75)), linetype="dashed")

ggsave("img1/theor_term_iii.png",ld_fst1)

herit_fst <- herit %>% filter(group_type=="Fst") %>% select(-group_type) %>% left_join(mean_fst_values,by="group_number")
h_fst1 <- herit_fst %>%
  ggplot(aes(x=median,y=rel_h2))+
  geom_point(size=1)+
  xlab("Within-Group Median Fst")+
  ylab("Term ii (Relative Heritability")+
  facet_wrap(~phenotype)+
  geom_vline(xintercept=quantile(ld_fst$median,probs=c(0.25,0.5,0.75)), linetype="dashed")

ggsave("img1/theor_term_ii.png",h_fst1)
#######################

het <- tibble()

for (file in list.files(path = 'data/theory',
                        pattern = '[a-z]{3}_[0-9]+.het', full.names = T)) {
  het_file <- read.table(file, sep = "" , header = T,stringsAsFactors = F) %>%
  mutate(group_number = parse_number(file)) %>%
  mutate(group_type = str_extract(file,"fst|wpc"))
  het <- het %>% bind_rows(het_file)
  
}
  
het1 <- het %>% group_by(group_type,group_number) %>%
  summarize(mean_O = mean(O.HOM.,na.rm=T),
            mean_E = mean(E.HOM.,na.rm=T)) %>%
  ungroup() %>%
  mutate(mean_O_1 = max(mean_O*(group_number==1))) %>%
  mutate(rel_mean_O = mean_O / mean_O_1) %>%
  mutate(group_type = factor(group_type,labels=c("Fst","WPC"), levels=c("fst","wpc"))) %>%
  left_join(median_values,by=c("group_type","group_number"))

homo_graph <- het1 %>%
#  pivot_longer(mean_O:mean_E,names_to="Homozygosity",values_to="Values") %>%
#  mutate(Homozygosity = factor(Homozygosity,levels=c("mean_E","mean_O"),labels=c("Expected","Observed")))%>%
#  filter(Homozygosity =="Observed") %>%
  ggplot(aes(x=median,y=rel_mean_O,color=group_type))+
  geom_line(size=1)+
  geom_point(size=1.5,color="black")+
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  xlab("Within-Group Median Genetic Distance")+
  ylab("Relative Average of Observed # of Homozyg. Genotypes per Sample")+
  facet_wrap(~group_type,scales="free_x")
ggsave("img1/homozygosity.png",homo_graph)

####################


##################

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





