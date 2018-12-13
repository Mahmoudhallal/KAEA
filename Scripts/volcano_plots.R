##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Plot volcano plots
## Date: 28.08.2018
## Author: Mahmoud Hallal
##################################################
## Load libraries
source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library(Biobase)

library(ggplot2)

library(yaml)


## Load parameters file
params <- read_yaml("./config.yaml")

## Define logFC function
getLogFC <- function(row, phenoTable, groupA, groupB, groupVar) {
  indicesA = phenoTable[[groupVar]] == groupA
  indicesB = phenoTable[[groupVar]] != groupA
  mean(row[indicesA]) - mean(row[indicesB])
}

## Define PValue function
getPValue <- function(row, phenoTable, groupA, groupB, groupVar, ttest) {
  row = row + runif(length(row), -1e-6, 1e-6)
  indicesA = phenoTable[[groupVar]] == groupA
  indicesB = phenoTable[[groupVar]] != groupA
  #cz we are only looking at overrepresentation
  t.test(row[indicesA], row[indicesB], alternative = ttest)$p.value
}


## Define a function to create a table of p-value, p-adj and logFC
getTopTable <- function(contrasts, eSet, groupVar, x, ttest) {
  groupA = contrasts[x]
  groupB = contrasts[-x]
  topTable = fData(eSet)
  eSet_modified <- exprs(eSet)
  #keep complete cases only
  eSet_modified[eSet_modified==0] <- NA
  eSet_modified <- eSet_modified[complete.cases(eSet_modified),]
  #calculate FC
  topTable = apply(eSet_modified, 1, getLogFC, pData(eSet), groupA, groupB, groupVar)
  topTable <- as.data.frame(topTable)
  #calculate pvalue
  topTable$p = apply(eSet_modified, 1, getPValue, pData(eSet), groupA, groupB, groupVar,ttest)
  topTable$padj = p.adjust(topTable$p, method="fdr")
  topTable = topTable[order(topTable$padj, topTable$p),]
  topTable
}


## Load input expressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"test_eSet1_",params$cell_line,".Rda"))

## Define the conditions to compare
contrasts <- unlist(strsplit(params$Conditions,","))

## Calculate two-tailed t-test
topTables_ts <- lapply(1:length(contrasts), function(x) getTopTable(contrasts, test_eSet, "cell_line", x, "two.sided"))
names(topTables_ts) <- contrasts

## Plot volcano plots for every condition
for (x in 1:length(topTables_ts)){
  data <- topTables_ts
  #Mark threshold on 0.05pvalue
  data[[x]]$threshold = as.factor(data[[x]]$padj < 0.05)
  data[[x]]$prots <- rownames(data[[x]])
  
  #plot the volcano plot
  g <- ggplot(data=data[[x]],
              aes(x=topTable, y =-log10(p),colour=threshold, label=rownames(data[[x]]))) +
    geom_point(alpha=0.4, size=1.75) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle("Volcano plot",names(data)[[x]]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    theme(axis.title=element_text(size=24))

  #create pdf output
  pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Volcano_plot_",params$cell_line,"_x_",names(data)[[x]],".pdf"), width = 8)
  plot(g)
  dev.off()
}

## Save output for Shiny
data_volcano <- topTables_ts
volcano_material <- lapply(1:length(topTables_ts), function(x){
  #Mark threshold on 0.05pvalue
  data_volcano[[x]]$threshold <- as.factor(data_volcano[[x]]$p < 0.05)
  data_volcano[[x]]$prots <- rownames(data_volcano[[x]])
  data_volcano[[x]]
})
names(volcano_material) <- names(topTables_ts)
save(volcano_material, file = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Volcano_plot_material_",params$cell_line,".Rda"))
