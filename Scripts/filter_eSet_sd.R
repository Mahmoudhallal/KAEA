##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Filter out proteins based on the row Medians, replace by zero all the values of a specific protein if the median of the 3 replicates is 0
## Date: 13.06.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
source("https://bioconductor.org/biocLite.R")

install.packages("norm", repos = "https://cran.rstudio.com")
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

install.packages("/home/user/KAEA/Packages/SetTools_1.01.tar.gz", repos = NULL, type = "source")
library(SetTools)

install.packages("/home/user/KAEA/Packages/OmicsTools_0.0.22.tar.gz", repos = NULL, type = "source")
library(OmicsTools)

library(yaml)

## Load parameters file
params <- read_yaml("./config.yaml")

### FOR normalisation
## Prepare MSnSet object
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"test_eSet_notfiltered_",params$cell_line,".Rda"))
cell_lines <- strsplit(params$Conditions,",")[[1]]
exprs_all <- exprs(test_eSet)


## Keep if rowMedian of biological replicates is >0
# for (x in 1:length(cell_lines)) {
#   indices <- which(test_eSet$cell_line == cell_lines[x])
#   indices_tokeep <- rowMedians(as.matrix((exprs_all[, indices]))) =! 0
#   exprs_all[!indices_tokeep,indices] <- 0
# }
# exprs(test_eSet) <- exprs_all
#
##Make 0s NAs before normalisation
exprs(test_eSet)[exprs(test_eSet)== 0] <- NA
yy <- as.MSnSet.ExpressionSet(test_eSet)

## Normalise by diff.mean
msnset.nrm <- normalise(yy, "quantiles")
exprs_all <- exprs(msnset.nrm)
exprs_all[is.na(exprs_all)] <- 0
## Remove lines with 0s only
exprs_all2 <- exprs_all[rowSums(exprs_all) != 0,]

##### end
## Load expressionSet
# load(paste0(params$CWD,"/results/test_eSet_notfiltered_",params$cell_line,".Rda"))
# cell_lines <- strsplit(params$Conditions,",")[[1]]
# exprs_all <- exprs(test_eSet)
# 
# ## Keep if rowMedian of biological replicates is >0
# for (x in 1:length(cell_lines)) {
#   indices <- which(test_eSet$cell_line == cell_lines[x])
#   indices_tokeep <- rowMedians(as.matrix((exprs_all[, indices]))) > 0
#   exprs_all[!indices_tokeep,indices] <- 0
#  }
# 
# ## Remove lines with 0s only
# exprs_all2 <- exprs_all[rowSums(exprs_all) != 0,]

#exprs(test_eSet) <- exprs_all2

## Imputation based on the condition
#Impute based on cell line separately
#Replace 0 by NA if 1 or 2 values are missing for a triplicate
if (params$Imputation == "T"){
  print("## IMPUTATION START ##")
  options("expressions"=500000)
  for (x in 1:length(cell_lines)) {
    indices <- which(test_eSet$cell_line == cell_lines[x])
    #indices_NA <-  (apply(exprs_all[, indices] == 0, 1, sum) == 1) | (apply(exprs_all[, indices] == 0, 1, sum) == 2)
    indices_NA <-  (apply(exprs_all[, indices] == 0, 1, sum) < ncol(exprs_all[, indices]))
    exprs_all[indices_NA,indices][exprs_all[indices_NA,indices] == 0]  <- NA
  }
  exprs(test_eSet) <- exprs_all
  exprs_all2 <- exprs_all
  yy <- as.MSnSet.ExpressionSet(test_eSet)
  
  # impute separately for every cell line and them add all to test_eSet
  # #MLE
  for (x in 1:length(cell_lines)) {
    indices <- which(yy$cell_line == cell_lines[x])
    exprs_all2[,indices] <- exprs(impute(yy[,indices], method = "MLE"))
  }
  exprs_all2 <- exprs_all2[!rowSums(exprs_all2) ==0,]
  print("## IMPUTATION FINISH ##")
  
}
#exprs(test_eSet) <- exprs_all2

##KNN
# data <- data.frame(exprs(yy))
# row_data <-seq(1:nrow(data))
# ind <- split(row_data, ceiling(seq_along(row_data)/2000))
# 
# #ind<-as.factor(c(gl(round(row_data/clus)-1,clus),rep(round(row_data/clus)-1+1, row_data-length(gl(round(row_data/clus)-1,clus)))))
# #newMat <- split(data, ind)
# #newMat2 <- newMat
# for(y in 1:length(ind)){
#   for (x in 1:length(cell_lines)) {
#     indices <- which(yy$cell_line == cell_lines[x])
#     exprs_all2[unlist(ind[y], use.names=FALSE),indices] <- exprs(impute(yy[unlist(ind[y], use.names=FALSE),indices], method = "knn"))
#   }
# }
# exprs(test_eSet) <- exprs_all2


# myimpute <- function(data,clus = 2000) # clus = Number or row in splited matrix
# {
#   library(impute)
#   data <- data.frame(exprs(yy))
#   row_data <-nrow(data)
#   ind<-as.factor(c(gl(round(row_data/clus)-1,clus),rep(round(row_data/clus)-1+1, row_data-length(gl(round(row_data/clus)-1,clus)))))
#   newMat <- split(data, ind)
#   res <- lapply(newMat,function(x)impute.knn(as.matrix(x)))
#   res <- lapply(res,"[[","data")
#   res <- do.call(rbind, res)
#   res
# }

####
#Recreate the expression set
aDataFrame <- function(input) {
  metaData = data.frame(labelDescription=colnames(input),
                        row.names=colnames(input))
  new("AnnotatedDataFrame", data=input, varMetadata=metaData)
}

test_eSet <- ExpressionSet(assayData = exprs_all2, aDataFrame(pData(test_eSet)))


## Save expressionSet
save(test_eSet, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"test_eSet1_",params$cell_line,".Rda"))

