##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: create collection for SetRank analysis
## Date: 09.03.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library(Biobase)

#install.packages("SetRank", repos="http://cran.rstudio.com/")
library(SetRank)

library(yaml)

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation
maxSetSize <- 500

## Define the fucntion to create a "collection"
create_collection <- function(list_prots, all_dbs, maxSetSize){
  #Define reference: union of all phos proteins of all conditions
  reference_union_all_phos_prots <- unique(list_prots)
  
  #create the collection
  options(mc.cores=2)
  collection_union_all_phos_prots500 <- buildSetCollection(all_dbs, referenceSet = reference_union_all_phos_prots, maxSetSize = maxSetSize)
  
  #return the collection object
  return(collection_union_all_phos_prots500)
}

## Load input expressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

## Load database 
all_dbs <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/all_dbs.csv"))

## Create the collection
collection <- create_collection(rownames(exprs(test_eSet)), all_dbs, maxSetSize = maxSetSize)

## Write output file
save(collection, file = paste(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/collection_",params$cell_line,".Rda",sep=""))
