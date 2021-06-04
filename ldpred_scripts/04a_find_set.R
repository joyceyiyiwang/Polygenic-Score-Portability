library(tidyverse)
.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(bigsnpr)


phenotypes <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
                "Lymphocyte","WBC","Eosinophil","DBP","SBP","Hb","Ht","MCHC","Neutrophil","Basophil" )

for (pheno in phenotypes){

  fam <- read.table('data/LDpred1/LD_EUR_train_1.fam') %>% as_tibble() %>% select(V1,V2)
  colnames(fam) <- c("#FID","IID")
  
  for (file in list.files(path = 'data/LDpred1/val_prs',
                          pattern = '[a-zA-Z]+_LDpred_scores_[0-9]+.sscore', full.names = T)) {
  
    phenotype <- str_extract(string = file, pattern = '(?<=data/LDpred1/val_prs/)[A-Za-z]+(?=_)')
  
    if(phenotype==pheno){
      score <- str_extract_all(string = file, pattern = '[0-9]+')[[1]][2] #Latter part may be affected by #s in directory
      if (nchar(score)==1) {score <- paste("0",score,sep="")}
      add <- read_tsv(file,col_names=T,col_types = cols()) %>% select(`#FID`,IID,SCORE1_AVG)
      colnames(add)[3] <- paste("SCORE",score,"_AVG",sep = "")
  
      fam <- fam %>% full_join(add,by=c("#FID","IID"))
    }
  }
 fam <- fam %>% 
   select(sort(names(.)))
 
 pred_auto <- fam %>% select(-"#FID",-"IID"
 #,-"SCORE06_AVG",-"SCORE28_AVG"
 )
  
  #Filtering conditions of coefficients done similarly by Prive et al
  sc <- apply(pred_auto, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  
  # Select sets of betas that are approved
  sumstats_file <- paste0('data/LDpred1/prs/',pheno,'_LDpred_sumstats.tsv')
  sumstats <- read_tsv(sumstats_file, col_types=cols())
    
  final_beta_auto <- sumstats %>% select(starts_with("V"))
  final_beta_auto <- rowMeans(final_beta_auto[, keep])
  
  sumstats$LDpred <- final_beta_auto
  
  final_sumstats <- paste0('data/LDpred1/',pheno,"_final_LDpred.tsv")
  sumstats %>% 
	# select(chr,rsid,pos,beta,LDpred) %>%
	 write_tsv(final_sumstats)
}

