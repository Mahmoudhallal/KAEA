##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Perform setRank enrichment and export network
## Date: 09.03.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
#install.packages("SetRank", repos="http://cran.rstudio.com/")
library(SetRank)

#install.packages("yaml", repos="http://cran.rstudio.com/")
library(yaml)

## Load parameters 
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Load collection
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/collection_",params$cell_line,".Rda"))

## Load topTables files
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_less.Rda"))
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_greater.Rda"))

## Define the number of cores for setRank 
options(mc.cores=1)

## Underactive kinases
## Define enrichment analysis
enrichment_analysis <- 
  lapply(1:length(topTables_less), function(x) {
    setRankAnalysis(rownames(topTables_less[[x]]), collection, use.ranks= TRUE, setPCutoff = params$pvalue_cutoff , fdrCutoff = params$fdr_cutoff, delete = TRUE)
  })

gene_lists <- lapply(1:length(topTables_less), function(x) rownames(topTables_less[[x]]))
names(gene_lists) <- names(topTables_less)

names(enrichment_analysis) <- names(gene_lists)

all_nets <- enrichment_analysis
all_gene_lists <- gene_lists

## Export output file
exportMultipleResults(all_nets, all_gene_lists, collection ,IDConverter = NULL, paste(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_less/",sep=""))

print("DONEEEEEE")
## Overactive kianses
enrichment_analysis1 <- 
  lapply(1:length(topTables_greater), function(x) {
    setRankAnalysis(rownames(topTables_greater[[x]]), collection, use.ranks= TRUE, setPCutoff = params$pvalue_cutoff , fdrCutoff = params$fdr_cutoff, delete = TRUE)
  })

gene_lists <- lapply(1:length(topTables_greater), function(x) rownames(topTables_greater[[x]]))
names(gene_lists) <- names(topTables_greater)

names(enrichment_analysis1) <- names(gene_lists)

all_nets <- enrichment_analysis1#
all_gene_lists <- gene_lists

## Export output file
exportMultipleResults(all_nets, all_gene_lists, collection ,IDConverter = NULL, paste(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_greater/",sep=""))

