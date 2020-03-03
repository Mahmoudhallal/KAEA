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
#install.packages("imputeLCMD", repos = "https://cran.rstudio.com")

library(yaml)

## Load parameters file
params <- read_yaml("./config.yaml")
imp <- params$Imputation

aDataFrame <- function(input) {
  metaData = data.frame(labelDescription=colnames(input),
                        row.names=colnames(input))
  new("AnnotatedDataFrame", data=input, varMetadata=metaData)
}


### FOR normalisation
## Prepare MSnSet object
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet_notfiltered_",params$cell_line,".Rda"))
cell_lines <- strsplit(params$Conditions,",")[[1]]
exprs_all <- exprs(test_eSet)


## Keep if rowMedian of biological replicates is >0
#
##Make 0s NAs before normalisation
exprs(test_eSet)[is.nan(exprs(test_eSet))] <- NA
exprs(test_eSet)[exprs(test_eSet)== 0] <- NA

########
#MSnSet ExpressionSet
yy <- as.MSnSet.ExpressionSet(test_eSet)

## Normalise by diff.mean
msnset.nrm <- normalise(yy, "quantiles")
exprs_all <- exprs(msnset.nrm)
exprs_all[is.na(exprs_all)] <- 0

##filter
for (x in 1:length(cell_lines)) {
  indices <- which(test_eSet$cell_line == cell_lines[x])
  indices_tokeep <- rowMedians(abs(as.matrix((exprs_all[, indices])))) != 0
  exprs_all[!indices_tokeep,indices] <- 0
}
exprs(test_eSet) <- exprs_all



### Imputation based on the condition
#Impute based on cell line separately
#Replace 0 by NA if 1 or 2 values are missing for a triplicate
if (params$Imputation == "T"){
  print("## IMPUTATION START ##")
  options("expressions"=500000)
  for (x in 1:length(cell_lines)) {
    indices <- which(test_eSet$cell_line == cell_lines[x])
    #indices_NA <-  (apply(exprs_all[, indices] == 0, 1, sum) == 1) | (apply(exprs_all[, indices] == 0, 1, sum) == 2)
    assign(paste0('indices_0_',cell_lines[x]), (apply(exprs_all[, indices] == 0, 1, sum) == ncol(exprs_all[, indices])))
    
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
    indices0 <- get(paste0('indices_0_',cell_lines[x]))
    exprs_all2[!indices0,indices] <- exprs(impute(yy[!indices0,indices], method = "MLE"))
  }
  exprs_all2 <- exprs_all2[!rowSums(exprs_all2) ==0,]
  
  
  #imputation2
  library(imputeLCMD)
  exprs_all2[exprs_all2 == 0] <- NA
  exprs_all3 <- exprs_all2
  test_eSet <- ExpressionSet(assayData = exprs_all3, aDataFrame(pData(test_eSet)))
  #exprs(test_eSet) <- exprs_all3
  yy2 <- as.MSnSet.ExpressionSet(test_eSet)
  for (x in 1:length(cell_lines)) {
    indices <- which(yy2$cell_line == cell_lines[x])
    exprs_all3[,indices] <- exprs(impute(yy2[,indices], method = "MinDet"))
  }
  
  #exprs_all2[is.na(exprs_all2)] <- min(exprs_all2, na.rm=TRUE)
  exprs_all3
  print("## IMPUTATION FINISH ##")
  
} else {
  exprs_all <- exprs_all[rowSums(exprs_all) != 0,]
  for (x in 1:length(cell_lines)) {
    indices <- which(test_eSet$cell_line == cell_lines[x])
    #indices_NA <-  (apply(exprs_all[, indices] == 0, 1, sum) == 1) | (apply(exprs_all[, indices] == 0, 1, sum) == 2)
    assign(paste0('indices_0_',cell_lines[x]), (apply(as.data.frame(exprs_all[, indices]) == 0, 1, sum) == ncol(as.data.frame(exprs_all[, indices]))))
    
    indices_NA <-  (apply(as.data.frame(exprs_all[, indices]) == 0, 1, sum) < ncol(as.data.frame(exprs_all[, indices])))
    exprs_all[indices_NA,indices][exprs_all[indices_NA,indices] == 0]  <- NA
    
    
  }
  exprs_all3 <- exprs_all
  exprs_all3
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

# ## Normalise according to SN proteome
# protein_groups <- read.delim('/Users/Mahmoud.Hallal/Desktop/PhD/cell_lines_data/Stimulated_data/Molm13_Rh_Sn_05.07.2019_all_samples/20190709_MQres_MOLM13_SN_MH/proteinGroups.txt')
# protein_groups1 <- protein_groups[protein_groups$Reverse != "+",]
# protein_groups1 <- protein_groups1[protein_groups1$Potential.contaminant != "+",]
# protein_groups1$cleaned_prots <- gsub("[|].*","", protein_groups1$Protein)
# 
# keep_cols <- c('cleaned_prots',"Intensity.R_C1","Intensity.R_C2","Intensity.R_C3","Intensity.R_D1","Intensity.R_D2","Intensity.R_D3","Intensity.S_C1","Intensity.S_C2","Intensity.S_C3","Intensity.S_D1","Intensity.S_D2","Intensity.S_D3")
# protein_groups2 <- protein_groups1[,keep_cols]
# protein_groups2[protein_groups2 == 0] <- 1
# protein_groups2[,c(2,3,4,5,6,7,8,9,10,11,12,13)] <- log2(protein_groups2[,c(2,3,4,5,6,7,8,9,10,11,12,13)])
# colnames(protein_groups2) <- c("cleaned_prots","X1_MOLM13H_RC", "X2_MOLM13H_RC","X3_MOLM13H_RC",
#                                "X1_MOLM13H_RD", "X2_MOLM13H_RD" ,"X3_MOLM13H_RD",
#                                "X1_MOLM13H_SC", "X2_MOLM13H_SC", "X3_MOLM13H_SC",
#                                "X1_MOLM13H_SD", "X2_MOLM13H_SD", "X3_MOLM13H_SD")
# 
# protein_groups3 <- protein_groups2[,-1]
# 
# #phosphoproteome
# phos_data <- as.data.frame(exprs_all3)
# phos_data$protein <- gsub("(.*)_.*","\\1",rownames(phos_data))
# 
# phos_data2 <- phos_data[,-ncol(phos_data)]
# phos_data22 <- phos_data[,-ncol(phos_data)]
# 
# 
# #
# phosphoproteome <- unique(phos_data$protein)
# #3162
# 
# proteome <- unique(protein_groups2$cleaned_prots)
# #4435
# 
# intersection <- intersect(proteome,phosphoproteome)
# #1359
# 
# pos <- lapply(1:nrow(phos_data), function(x){
#   prot <- phos_data$protein[x]
#   prot2 <- paste0('^',prot,'$')
#   grep(prot2, protein_groups2$cleaned_prots)
# })
# 
# for(y in 1:length(pos)){
#   print(y)
#   if (isEmpty(pos[[y]])){
#   }else{
#     phos_data2[y,] <- phos_data2[y,]/protein_groups3[pos[[y]],]
#   }
# }
# #phos_data3 <- phos_data2[!is.na(phos_data2)]
# 
# phos_data3 <- as.matrix(phos_data2)
# 






#########################################################
#Recreate the expression set
aDataFrame <- function(input) {
  metaData = data.frame(labelDescription=colnames(input),
                        row.names=colnames(input))
  new("AnnotatedDataFrame", data=input, varMetadata=metaData)
}

test_eSet <- ExpressionSet(assayData = exprs_all3, aDataFrame(pData(test_eSet)))


## Save expressionSet
save(test_eSet, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

