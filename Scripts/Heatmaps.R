##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Create heatmaps from samples 
## Date: 02.07.19
## Author: Mahmoud Hallal
##################################################

## Load libraries
library(Biobase)

#install.packages("ggplot2", repos="http://cran.rstudio.com/")#
library(ggplot2)

library(pheatmap)

library(yaml)

library(matrixStats)

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation
conds <- unlist(strsplit(params$Conditions,","))
control <- params$control

## Load ExpressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))
data <- exprs(test_eSet)

## Choose topVariances
topVarianceGenes <- head(order(rowVars(data), decreasing=T),1000)
matrix <- data[ topVarianceGenes, ]

## plot output pdf
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Heatmap_DEsites_",params$cell_line,".pdf"), width = 5)
  
pheatmap(data, 
           show_rownames = F, 
           fontsize_col = 8,
           color=colorRampPalette(c("blue", "white", "red"))(nrow(data)),
           border_color = NA,
           scale='row',
           na_col="grey",
         #silac fribourg 05.04.20
         cluster_rows = T
         )
dev.off()
