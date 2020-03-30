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
imp <- params$Imputation

#########################################################################################
## Perform a MAD ranking and cutoff to take the proteins with the hight variance only
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
  
  e_expr2 <- MAD_ranking(e, 0, nrow(exprs(e)))
  e_expr2[is.na(e_expr2)] <- 0
  
  #perform pca
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
    geom_text(aes(label=gsub('(X)(.*)','\\2',d$name)), size=3, hjust=0.5, vjust=-0.3)
  
  plot(p)
}

## Load expressionSet 
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

## Produce output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/PCA_plot_filtered_",params$cell_line,".pdf"), onefile = TRUE)
pcaPlot(test_eSet, "cell_lines_PCA")
dev.off()
