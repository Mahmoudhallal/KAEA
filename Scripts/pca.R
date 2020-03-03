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
#BiocManager::install("SummarizedExperiment", version = "3.8")
#library(SummarizedExperiment)

## Load parameters file
params <- read_yaml("./config.yaml")
imp <- params$Imputation

#########################################################################################
#perform a MAD ranking and cutoff to take the proteins with the hight variance only
MAD_ranking <- function(e, z_cutoff, mad_cutoff){
  exprMatrix1 = data.frame(exprs(e))
  exprMatrix1_cleaned <- exprMatrix1[rowSums(exprMatrix1 > 0,na.rm=T)>z_cutoff,]
  exprMatrix1_cleaned[is.na(exprMatrix1_cleaned)] <- 0
  exprMatrix1_cleaned <- exprMatrix1
  mads = apply(exprMatrix1_cleaned, 1, mad, na.rm=TRUE)
  cutoff = sort(mads, decreasing = TRUE)[mad_cutoff]
  new_exprMatrix = exprMatrix1_cleaned[mads >= cutoff,]
  return(new_exprMatrix)
}


## PCA function
pcaPlot <- function(e, title) {
  #old usage
  #Define input df
  # e_expr <- as.data.frame(exprs(e))
  # #apply pca function 
  # pca = prcomp(e_expr)
  # plotData = cbind(pData(e), pca$rotation[,1:3])
  # plotData$group = as.factor(plotData$group)
  # plot(pca, main="cutoff = 100")
  # 
  # summary <- summary(pca)
  
  ## new 20.03.2019
  #mad_cutoff <- 500
  #e_expr <- as.data.frame(exprs(e))
  e_expr2 <- MAD_ranking(e, 0, nrow(exprs(e)))
  e_expr2[is.na(e_expr2)] <- 0
  #rv <- rowVars(as.matrix(e_expr))
  #select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  #pca <- prcomp(t(e_expr[select, ]))
  pca = prcomp(t(e_expr2))
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], name = rownames(pca$x))
  d$group <- gsub('(X\\d+)_(.*)','\\2',d$name)
  #Proportion of Variance is nothing else than normalized standard deviations.
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  attr(d, "percentVar") <- percentVar[1:2]

  #plot PC histogram
  plot(pca)
  
  #plot pca
  p <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] *100), "% variance")) + 
    #coord_fixed() +
    geom_text(aes(label=gsub('(X)(.*)','\\2',d$name)), size=3, hjust=0.5, vjust=-0.3)
  plot(p)


#__impF_26.02.2019
  #plot PC1 and PC2
  # p = ggplot(data=plotData,
  #            aes(x=PC1, y=PC2, color=cell_line)) +
  #   geom_point()+ ggtitle("PCA of cell lines data") +
  #   geom_text(aes(label=plotData$group), size=3) +
  #   labs(x = paste0("PC1(", round(as.data.frame(summary$importance)$PC1[2]*100,3),"%)"),
  #        y = paste0("PC2(", round(as.data.frame(summary$importance)$PC2[2]*100,3),"%)"))
  # plot(p)
  
  #plot PC1 and PC3
  # p2 = ggplot(data=plotData, 
  #            aes(x=PC1, y=PC3, color=cell_line)) + 
  #   geom_point()+ ggtitle("PCA cell lines data") +
  #   geom_text(aes(label=plotData$group), size=3) +
  #   labs(x = paste0("PC1(", round(as.data.frame(summary$importance)$PC1[2]*100,3),"%)"),
  #        y = paste0("PC3(", round(as.data.frame(summary$importance)$PC3[2]*100,3),"%)"))
  #  plot(p2)
}

## Load expressionSet 
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

## Produce output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/PCA_plot_filtered_",params$cell_line,".pdf"), onefile = TRUE)
pcaPlot(test_eSet, "cell_lines_PCA")
dev.off()

## tSNE
# load('../results/K562_0.01P_0.05FDR/test_eSet1_K562.Rda')
# load('../results/K562_0.01P_0.05FDR/topTables_K562_less.Rda')
# rr <- exprs(test_eSet)
# rr <- as.data.frame(topTables_less$K562_DRG$topTable)
# rr <- rr[-duplicated(rr),]
# cc <- as.matrix(rr)
# cc <- t(cc)
# cc <- rbind(rr[,c(1,2,3)],rr[,c(4,5,6)])
# #cc <- as.matrix(cc$value)
# tsne_out <- Rtsne(cc,check_duplicates=FALSE,perplexity = 500) # Run TSNE
# tsne_out <- Rtsne(cc, dims = 2, perplexity=5000, verbose=TRUE, check_duplicates=FALSE, max_iter = 500)
# 
# #dnames <- c(rep("DasR_Par", nrow(rr)), rep("PazR_Par",nrow(rr)))
# dnames <- c(rep("Molm13_CTRL", nrow(rr)), rep("Molm13_DRG",nrow(rr)))
# dnames <- rownames(cc)
# #dnames <- ll$Var2
# tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],col = dnames)#, 
# ggplot(tsne_plot) +
#   geom_point(aes(x=x, y=y,col=col))
# 
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
