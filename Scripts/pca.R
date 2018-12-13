##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: PCA analysis 
## Date: 05.03.17
## Author: Mahmoud Hallal/ Cedric Simillion
##################################################
## Load libraries
#source("https://bioconductor.org/biocLite.R")

library(ggplot2)
library(Biobase)

library(yaml)

## Load parameters file
params <- read_yaml("./config.yaml")

#########################################################################################
#perform a MAD ranking and cutoff to take the proteins with the hight variance only
MAD_ranking <- function(e, z_cutoff, mad_cutoff){
  exprMatrix1 = data.frame(exprs(e))
  exprMatrix1_cleaned <- exprMatrix1[rowSums(exprMatrix1 > 0)>z_cutoff,]
  
  exprMatrix1_cleaned <- exprMatrix1
  mads = apply(exprMatrix1_cleaned, 1, mad)
  cutoff = sort(mads, decreasing = TRUE)[mad_cutoff]
  new_exprMatrix = exprMatrix1_cleaned[mads >= cutoff,]
  return(new_exprMatrix)
}


## PCA function
pcaPlot <- function(e, title) {
  #Define input df
  e_expr <- as.data.frame(exprs(e))
  #apply pca function 
  pca = prcomp(t(e_expr))
  plotData = cbind(pData(e), pca$x[,1:3])
  plotData$group = as.factor(plotData$group)
  plot(pca, main="cutoff = 100")
  
  summary <- summary(pca)

  #plot PC1 and PC2
  p = ggplot(data=plotData, 
             aes(x=PC1, y=PC2, color=cell_line)) + 
    geom_point()+ ggtitle("PCA of cell lines data") +
    geom_text(aes(label=plotData$group), size=3) +
    labs(x = paste0("PC1(", round(as.data.frame(summary$importance)$PC1[2]*100,3),"%)"),
         y = paste0("PC2(", round(as.data.frame(summary$importance)$PC2[2]*100,3),"%)"))
  plot(p)
  
  #plot PC1 and PC3
  p2 = ggplot(data=plotData, 
             aes(x=PC1, y=PC3, color=cell_line)) + 
    geom_point()+ ggtitle("PCA cell lines data") +
    geom_text(aes(label=plotData$group), size=3) +
    labs(x = paste0("PC1(", round(as.data.frame(summary$importance)$PC1[2]*100,3),"%)"),
         y = paste0("PC3(", round(as.data.frame(summary$importance)$PC3[2]*100,3),"%)"))
   plot(p2)
}

## Load expressionSet 
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"test_eSet1_",params$cell_line,".Rda"))

## Produce output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"PCA_plot_filtered_",params$cell_line,".pdf"), onefile = TRUE)
pcaPlot(test_eSet, "cell_lines_PCA")
dev.off()


##
# rr <- exprs(test_eSet)
# rr <- rr[-duplicated(rr),]
# cc <- rbind(rr[,c(1,2,3)],rr[,c(4,5,6)])
# tsne_out <- Rtsne(cc,check_duplicates=FALSE) # Run TSNE
# #dnames <- c(rep("DasR_Par", nrow(rr)), rep("PazR_Par",nrow(rr)))
# dnames <- c(rep("Molm13_CTRL", nrow(rr)), rep("Molm13_DRG",nrow(rr)))
# 
# tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = dnames)
# ggplot(tsne_plot) + 
#   geom_point(aes(x=x, y=y, color=col))
# 
# #pca
# e_expr <- as.data.frame(exprs(test_eSet))
# #apply pca function 
# pca = prcomp(t(e_expr))
# ll <- as.data.frame(pca$rotation)
# #plot PC1 and PC2
# p = ggplot(data=ll, 
#            aes(x=PC1, y=PC2)) + 
#   geom_point()+ ggtitle("PCA of cell lines data") +
#   #geom_text(aes(label=plotData$group), size=3) +
#   labs(x = paste0("PC1(", round(as.data.frame(summary$importance)$PC1[2]*100,3),"%)"),
#        y = paste0("PC2(", round(as.data.frame(summary$importance)$PC2[2]*100,3),"%)"))
# plot(p)
