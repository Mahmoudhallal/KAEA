##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Produce heat map
## Date: 09.03.17
## Author: Cedric Simillion
##################################################

## Load libraries
#install.packages("gplots", repos="http://cran.rstudio.com/")
library(gplots)
#install.packages("pheatmap", repos="http://cran.rstudio.com/")
library(pheatmap)   
library(ggplot2)
library(RColorBrewer)
library(yaml)
library(Rcpp)

#install.packages("dplyr", repos="http://cran.rstudio.com/")
library(dplyr)

## Load parameters
params <- read_yaml("./config.yaml")

## Load input topTables of under-expressed phosphosites
inputFile_less = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_less/pathways.txt")
outputFile2 = paste0(params$CWD,"/results/",params$cell_line,"_",params$pvalue_cutoff,"P_",params$fdr_cutoff,"FDR/","Heatmap_",params$cell_line,"_" ,params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR.pdf',sep='')

results_less = read.table(inputFile_less, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
pathwayMatrix_less = as.matrix(results_less[,1:(ncol(results_less)-3)])
rownames(pathwayMatrix_less) <- results_less$description

## Load input topTables of over-expressed phosphosites
inputFile1_greater = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_greater/pathways.txt")

results_greater = read.table(inputFile1_greater, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
pathwayMatrix_greater = as.matrix(results_greater[,1:(ncol(results_greater)-3)])
rownames(pathwayMatrix_greater) <- results_greater$description


pathwayMatrix_less <- -pathwayMatrix_less
pathwayMatrix <- rbind(pathwayMatrix_greater, pathwayMatrix_less)

## Sum kinase values if they are present in over- and under- enriched analysis 
dd1 <- as.data.frame(pathwayMatrix_less)
dd2 <- as.data.frame(pathwayMatrix_greater)

ll <- bind_rows(dd1 %>% add_rownames(), 
                dd2 %>% add_rownames()) %>% 
  # evaluate following calls for each value in the rowname column
  group_by(rowname) %>% 
  # add all non-grouping variables
  summarise_all(sum)
cc <- as.data.frame(ll)
pathwayMatrix <- cc
rownames(pathwayMatrix) <- cc$rowname

contrasts <- unlist(strsplit(params$Conditions,","))
control <- params$control
if (params$SILAC == "T"){
  samples = contrasts
} else {
  samples <- contrasts[!(contrasts == control)]
}

colnames(pathwayMatrix) <- c('rowname',samples)
## Create output file for dataframe to be used for heatmap
write.table(pathwayMatrix, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Material_for_heatmap_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".csv"))

## Plot Heatmap
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
pathwayMatrix1 <- as.data.frame(pathwayMatrix[,-1])
rownames(pathwayMatrix1) <- rownames(pathwayMatrix)
colnames(pathwayMatrix1) <- samples

pdf(outputFile2, width=4, height=6)
pheatmap(pathwayMatrix1,color = my_palette, cluster_rows = T, cluster_cols = F, fontsize_row = 9, fontsize_col = 12, onefile=TRUE)
dev.off()



