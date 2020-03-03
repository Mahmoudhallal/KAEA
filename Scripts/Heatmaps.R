##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Create heatmaps from samples 
## Date: 02.07.19
## Author: Mahmoud Hallal
##################################################
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

## Load expressionSet
# load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))
# load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Volcano_plot_material_",params$cell_line,".Rda"))
# 
# ##select DE phosphosites
# for (x in 1:length(conds)){
#   samp <- conds[x]
#   DE_sites <- volcano_material[[samp]]
#   keep <- rownames(DE_sites)[which(DE_sites$padj<0.05)]
#   
#   ###
#   matrix <- exprs(test_eSet)
#   matrix2 <- matrix[rownames(matrix) %in% keep,]
#   #matrix2[matrix2 == 0] <- NA
#   
#   #plot
#   pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Heatmap_DEsites_",samp,'_',params$cell_line,".pdf"), width = 8)
#   
#   pheatmap(matrix2, 
#            show_rownames = F, 
#            fontsize_row = 8,
#            main=paste0("K562 - top ",nrow(matrix2)), 
#            color=colorRampPalette(c("blue", "white", "red"))(nrow(matrix2)),
#            border_color = NA,
#            scale='row',
#            na_col="grey")
#   dev.off()
# }


###
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))
data <- exprs(test_eSet)
#nn <- c("sMOLM13_MIDO1", "sMOLM13_MIDO2" ,"sMOLM13_MIDO3", "sMOLM13_CTRL1", "sMOLM13_CTRL2", "sMOLM13_CTRL3", "rMOLM13_CTRL1", "rMOLM13_CTRL2",
#        "rMOLM13_CTRL3", "rMOLM13_MIDO1", "rMOLM13_MIDO2", "rMOLM13_MIDO3")
  
#colnames(data) <- nn
#heatmap(data,scale='row')

topVarianceGenes <- head(order(rowVars(data), decreasing=T),1000)
matrix <- data[ topVarianceGenes, ]

#plot
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Heatmap_DEsites_",params$cell_line,".pdf"), width = 5)
  
pheatmap(data, 
           show_rownames = F, 
           fontsize_col = 8,
         #angle_col = 315,
           #main=paste0("K562 - top ",nrow(matrix2)), 
           color=colorRampPalette(c("blue", "white", "red"))(nrow(data)),
           border_color = NA,
           scale='row',
           na_col="grey")
dev.off()

