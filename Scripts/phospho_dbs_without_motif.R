##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Construct the final databse when not using motif analysis
## Date: 09.03.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library(Biobase)

library(yaml)

## Load parameters
params <- read_yaml("./config.yaml")

## Print that no MOTIF analysis will be used
print("PHOSPHO_DBS_WITHOUT_MOTIF_WILL_BE_USED!!")

## Load DB without motif analysis
all_dbs2 <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"all_dbs1.csv"))

## Redefine the name
all_dbs <- all_dbs2

## Remove duplicated cases
all_dbs <- all_dbs[!duplicated(all_dbs),]

## Write output file
write.csv(all_dbs, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"all_dbs.csv"))
#