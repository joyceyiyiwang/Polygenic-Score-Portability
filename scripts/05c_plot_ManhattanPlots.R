#install.packages("qqman", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("dplyr", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("rlang", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("vctrs", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("glue", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("magrittr", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("ellipsis", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("fansi", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("tibble", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("purrr", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")
#install.packages("Rcpp", repos=c("https://cloud.r-project.org"), lib="/home1/06568/joyce_w/R_libs/")

.libPaths("/home1/06568/joyce_w/R_libs/")
suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
options(bitmapType='cairo')

library("qqman",lib="/home1/06568/joyce_w/R_libs/")
library("dplyr",lib="/home1/06568/joyce_w/R_libs/")

phenotypes <- c("BMI", "Eosinophil", "Height", "Lymphocyte",
  "MCH", "MCV", "Monocyte", "Platelet", "RBC", "WBC")

for (i in 1:length(phenotypes)){
	chr1 <- paste0("data/gwas_results/",phenotypes[i],".chr1.",phenotypes[i],".glm.linear")
	results_as <- read.csv(chr1,sep="\t")
	
	#Bind linear regression output together of all autosomes
	for (j in 2:22){
        	chr_bind <- paste0("data/gwas_results/",phenotypes[i],".chr",j,".",phenotypes[i],".glm.linear")
        	results_bind <- read.csv(chr_bind,sep="\t")
        	results_as <- bind_rows(results_as,results_bind)
	}
	results_as <- results_as[!is.na(results_as$P),]
	results_as <- results_as[results_as$P != 0,]

	# Manhattan Plot
	png(paste0("img/",phenotypes[i],"_linear_manhattan.png"))
	manhattan(results_as,chr="X.CHROM",bp="POS",p="P",snp="ID",
		 main = paste0("Manhattan plot: ",phenotypes[i]),
		 col = c("blue4", "orange3"), suggestiveline=T, genomewideline=T, cex=0.4)
	dev.off()
	
	png(paste0("img/",phenotypes[i],"_linear_manhattan_zoomed.png"))
	manhattan(results_as,chr="X.CHROM",bp="POS",p="P",snp="ID",
	          main = paste0("Manhattan plot: ",phenotypes[i]),
	          col = c("blue4", "orange3"), suggestiveline=T, genomewideline=T, cex=0.4,
	          ylim = c(0, 20))
	dev.off()
	
	#QQ Plot of P-values
        png(paste0("img/",phenotypes[i],"_linear_qqplot.png"))
	qq(results_as$P, main = paste0("Q-Q plot of GWAS p-values for ",phenotypes[i]),
		 xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)
	dev.off()
}

