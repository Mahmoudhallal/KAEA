##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Filter out proteins based on the row Medians, replace by zero all the values of a specific protein if the median of the 3 replicates is 0
## Date: 13.06.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
#source("https://bioconductor.org/biocLite.R")
set.seed(22)

#install.packages("norm", repos = "https://cran.rstudio.com")
library(norm)

#install.packages("RNetCDF", repos = "https://cran.rstudio.com")
#library(RNetCDF)

#biocLite("Biobase")
library(Biobase)

#install.packages("Rcpp", repos = "https://cran.rstudio.com")
library(Rcpp)

#biocLite("Rhdf5lib")
#library(Rhdf5lib)

#biocLite("mzR")
library(mzR)

#biocLite("MSnbase")
library(MSnbase)

#install.packages("/Users/Mahmoud.Hallal/Desktop/PhD/new_scripts/results_test_sm/SetTools_1.01.tar.gz", repos = NULL, type = "source")
library(SetTools)

#install.packages("/Users/Mahmoud.Hallal/Desktop/PhD/new_scripts/results_test_sm/OmicsTools_0.0.22.tar.gz", repos = NULL, type = "source")
library(OmicsTools)

#install.packages("mvtnorm", repos = "https://cran.rstudio.com")
#library(mvtnorm)

library(yaml)

## Load parameters file
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Define annotated dataframe function
aDataFrame <- function(input) {
  metaData = data.frame(labelDescription=colnames(input),
                        row.names=colnames(input))
  new("AnnotatedDataFrame", data=input, varMetadata=metaData)
}

## Normalisation
## Prepare MSnSet object
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet_notfiltered_",params$cell_line,".Rda"))
cell_lines <- strsplit(params$Conditions,",")[[1]]
exprs_all <- exprs(test_eSet)

#Get phosphosites with only 1 quantification among all replicates of all conditions
indices_to_remove_1_entry_only <- (apply(exprs_all != 0, 1, sum) == 1)
exprs(test_eSet) <- exprs_all

## Make 0s NAs before normalisation
exprs(test_eSet)[is.nan(exprs(test_eSet))] <- NA
exprs(test_eSet)[exprs(test_eSet)== 0] <- NA

## MSnSet ExpressionSet
yy <- as.MSnSet.ExpressionSet(test_eSet)

## Normalise by "quantiles"
msnset.nrm <- normalise(yy, "quantiles")
exprs_all <- exprs(msnset.nrm)
exprs_all[is.na(exprs_all)] <- 0
exprs(test_eSet) <- exprs_all


# Imputation based on the condition
## Impute based on cell line separately
## Replace 0 by NA if 1 or 2 values are missing for a triplicate
if (params$Imputation == "T"){
  print("## IMPUTATION START ##")
  options("expressions"=500000)
  for (x in 1:length(cell_lines)) {
    indices <- which(test_eSet$cell_line == cell_lines[x])
    #assign(paste0('indices_0_',cell_lines[x]), (apply(exprs_all[, indices] == 0, 1, sum) > 1))
    
    #select phosphosites with 1/3 missing values
    indices_NA <-  (apply(exprs_all[, indices] == 0, 1, sum) == 1)
    #indices_NA <-  (apply(exprs_all[, indices] == 0, 1, sum) < (2/3)*ncol(exprs_all[, indices]))
    #exprs_all[indices_NA,indices][exprs_all[indices_NA,indices] == 0]  <- NA
    
    #replace all 0s with NAs
    exprs_all[, indices][exprs_all[, indices]==0] <- NA
    
    assign(paste0('indices_0_',cell_lines[x]), indices_NA)
  }
  
  exprs(test_eSet) <- exprs_all
  exprs_all2 <- exprs_all
  
  #create MSNset
  yy <- as.MSnSet.ExpressionSet(test_eSet)
  
  ## impute separately for every cell line and them add all to test_eSet
  ## mixed model
  for (x in 1:length(cell_lines)) {
    indices <- which(yy$cell_line == cell_lines[x])
    indices0 <- get(paste0('indices_0_',cell_lines[x]))
    
    exprs_all2[,indices] <- exprs(impute(yy[,indices],method = "mixed",
                                         randna = indices0,
                                         mar = "MLE", mnar = "MinProb"))
    
  }
  
  #remove phosphosites with 1 entry in all replicates of all conditions (selected above)
  exprs_all4 <- exprs_all2[!indices_to_remove_1_entry_only,]
  exprs_all4
  print("## IMPUTATION FINISH ##")
  
} else {
  
  ## Filter: Keep if rowMedian of biological replicates is >0
  for (x in 1:length(cell_lines)) {
    indices <- which(test_eSet$cell_line == cell_lines[x])
    indices_tokeep <- rowMedians(abs(as.matrix((exprs_all[, indices])))) != 0
    exprs_all[!indices_tokeep,indices] <- 0
  }
  exprs(test_eSet) <- exprs_all
  
  #remove rows with only 0s
  exprs_all <- exprs_all[rowSums(exprs_all) != 0,]
  
  #introduce NAs
  for (x in 1:length(cell_lines)) {
    indices <- which(test_eSet$cell_line == cell_lines[x])
    assign(paste0('indices_0_',cell_lines[x]), (apply(as.data.frame(exprs_all[, indices]) == 0, 1, sum) == ncol(as.data.frame(exprs_all[, indices]))))
    
    indices_NA <-  (apply(as.data.frame(exprs_all[, indices]) == 0, 1, sum) < ncol(as.data.frame(exprs_all[, indices])))
    exprs_all[indices_NA,indices][exprs_all[indices_NA,indices] == 0]  <- NA
    
  }
  exprs_all4 <- exprs_all
  exprs_all4
}


#########################################################
## Recreate the ExpressionSet
test_eSet <- ExpressionSet(assayData = exprs_all4, aDataFrame(pData(test_eSet)))

## Save expressionSet
save(test_eSet, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

