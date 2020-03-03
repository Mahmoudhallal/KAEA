##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Plot volcano plots
## Date: 28.08.2018
## Author: Mahmoud Hallal
##################################################
## Load libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library(Biobase)

library(ggplot2)

library(yaml)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
#library(devtools)
#devtools::install_github("r-lib/xml2")

#library(xml2)

## Load parameters file
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Define logFC function
getLogFC <- function(row, phenoTable, groupA, groupB, groupVar) {
  indicesA = phenoTable[[groupVar]] == groupA
  indicesB = phenoTable[[groupVar]] != groupA
  mean(row[indicesA],na.rm=T) - mean(row[indicesB], na.rm=T)
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
  print(x)
  groupA = contrasts[x]
  groupB = contrasts[-x]
  topTable = fData(eSet)
  eSet_modified <- exprs(eSet)
  #keep complete cases only
  #eSet_modified[eSet_modified==0] <- NA
  #eSet_modified <- eSet_modified[complete.cases(eSet_modified),]
  #calculate FC
  topTable = apply(eSet_modified, 1, getLogFC, pData(eSet), groupA, groupB, groupVar)
  topTable <- as.data.frame(topTable)
  print(topTable)
  
  #calculate pvalue
  topTable$p = apply(eSet_modified, 1, getPValue, pData(eSet), groupA, groupB, groupVar,ttest)
  topTable$padj = p.adjust(topTable$p, method="fdr")
  topTable = topTable[order(topTable$padj, topTable$p),]
  topTable
}


## Load input expressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

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
              aes(x=topTable, y =-log10(padj),colour=threshold, label=rownames(data[[x]]))) +
    geom_point(alpha=0.4, size=1.75) +
    xlab("log2 fold change") + ylab("-log10 p-value") +
    theme_bw() +
    theme(legend.position="none") +
    ggtitle("Volcano plot",names(data)[[x]]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    theme(axis.title=element_text(size=24))

  #create pdf output
  pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Volcano_plot_",params$cell_line,"_x_",names(data)[[x]],".pdf"), width = 8)
  plot(g)
  dev.off()
}

## Map uniprot to gene names 
## Transform uniprot to protein names
library(biomaRt)
rnms <- rownames(test_eSet)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "uswest.ensembl.org")
prots <- gsub('(^.*)(\\-[d])?_([STY]\\d+)','\\1',rnms)
unique_prots <- unique(prots)

rr <- getBM(
  attributes= c("uniprotswissprot","hgnc_symbol"),
  filters = 'uniprotswissprot',
  values= unique_prots,
  mart= mart)

## Save output for Shiny
data_volcano <- topTables_ts
volcano_material <- lapply(1:length(topTables_ts), function(x){
  #Mark threshold on 0.05pvalue
  data_volcano[[x]]$threshold <- as.factor(data_volcano[[x]]$padj < 0.05)
  data_volcano[[x]]$prots <- rownames(data_volcano[[x]])
  data_volcano[[x]]$only_prots <- gsub('(^.*)_([STY]\\d+)','\\1',data_volcano[[x]]$prots)
  
  data_volcano[[x]]$gene <- rr$hgnc_symbol[match(data_volcano[[x]]$only_prots, rr$uniprotswissprot)]
  data_volcano[[x]]
})
names(volcano_material) <- names(topTables_ts)
save(volcano_material, file = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Volcano_plot_material_",params$cell_line,".Rda"))


############################################################################
## Create topTables for specific volcano plot of barplot in shiny
## Define a function to create a table of p-value, p-adj and logFC

## IF phosphosite has one or more kinases
all_dbs <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/all_dbs.csv"))
one_kinase <- names(table(all_dbs$geneID)[table(all_dbs$geneID) == 1])
five_kinase <-  names(table(all_dbs$geneID)[which(table(all_dbs$geneID) > 1 & table(all_dbs$geneID) <= 5)])

#more_kinase <- all_dbs$geneID[table(all_dbs$geneID) != 1]

#'Q9UQ35_T983' %in% one_kinase
## out of order 11.04.2019 because we use the normal one up
# getTopTable_special <- function(contrasts, eSet, groupVar, x, ttest) {
#   groupA = contrasts[x]
#   groupB = contrasts[-x]
#   topTable = fData(eSet)
#   eSet_modified <- exprs(eSet)
#   #keep complete cases only
#   eSet_modified[eSet_modified==0] <- 1
#   eSet_modified <- eSet_modified[complete.cases(eSet_modified),]
#   #calculate FC
#   topTable = apply(eSet_modified, 1, getLogFC, pData(eSet), groupA, groupB, groupVar)
#   topTable <- as.data.frame(topTable)
#   #calculate pvalue
#   topTable$p = apply(eSet_modified, 1, getPValue, pData(eSet), groupA, groupB, groupVar,ttest)
#   topTable$padj = p.adjust(topTable$p, method="fdr")
#   topTable = topTable[order(topTable$padj, topTable$p),]
#   topTable
# }

#topTables_ts_special <- lapply(1:length(contrasts), function(x) getTopTable_special(contrasts, test_eSet, "cell_line", x, "two.sided"))
topTables_ts_special <- topTables_ts
names(topTables_ts_special) <- contrasts

data_volcano_special <- topTables_ts_special
volcano_material_special <- lapply(1:length(topTables_ts_special), function(x){
  #Mark threshold on 0.05pvalue
  data_volcano_special[[x]]$threshold <- as.factor(data_volcano_special[[x]]$padj < 0.05)
  data_volcano_special[[x]]$prots <- rownames(data_volcano_special[[x]])
  data_volcano_special[[x]]$only_prots <- gsub('(^\\w+)_([STY]\\d+)','\\1',data_volcano_special[[x]]$prots)
  
  data_volcano_special[[x]]$gene <- rr$hgnc_symbol[match(data_volcano_special[[x]]$only_prots, rr$uniprotswissprot)]
  
  data_volcano_special[[x]]$share <- data_volcano_special[[x]]$prots %in% one_kinase
  data_volcano_special[[x]]$shape <- 'x'
  data_volcano_special[[x]][data_volcano_special[[x]]$share == "TRUE",]$share <- 'Unique'
  data_volcano_special[[x]][data_volcano_special[[x]]$share == "Unique",]$shape <- 'circle'
  
  data_volcano_special[[x]][which(data_volcano_special[[x]]$prots %in% five_kinase),]$share <- 'Shared (2-5)'
  data_volcano_special[[x]][which(data_volcano_special[[x]]$prots %in% five_kinase),]$shape <- 'o'
  
  data_volcano_special[[x]][which(data_volcano_special[[x]]$share == "FALSE"),]$share <- 'Not specific (>5)'
  #data_volcano_special[[x]][which(data_volcano_special[[x]]$share == "FALSE"),]$shape <- 'x'
  
  #
  data_volcano_special[[x]]$shape <- as.factor(data_volcano_special[[x]]$shape)
  data_volcano_special[[x]]$share <- as.factor(data_volcano_special[[x]]$share)
  
  data_volcano_special[[x]]
})
names(volcano_material_special) <- names(topTables_ts_special)
save(volcano_material_special, file = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Volcano_plot_material_",params$cell_line,"_special.Rda"))
