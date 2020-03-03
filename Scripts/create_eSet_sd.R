###########################################################-
# Objective: create ExpressionSet object from dataframe with Intensity, experiment and proteins columns
# Author: Mahmoud Hallal
# Date created: 03.03.17
###########################################################-

## Load libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rsamtools")
library(Rsamtools)
#biocLite("AnnotationDbi")
library(AnnotationDbi)

#biocLite("Biobase")
library(Biobase)
 
#biocLite("GenomicFeatures")
library(GenomicFeatures)
 
#biocLite("affy")
#library(affy)
 
#biocLite("limma")
#library(limma)
 
#biocLite("BSgenome")
library('BSgenome')

#biocLite("BSgenome.Hsapiens.UCSC.hg19")
library('BSgenome.Hsapiens.UCSC.hg19')

#biocLite("Homo.sapiens")
library(Homo.sapiens)

install.packages("/Users/Mahmoud.Hallal/Desktop/PhD/new_scripts/results_test_sm/SetTools_1.01.tar.gz", repos = NULL, type = "source")
library(SetTools)

install.packages("/Users/Mahmoud.Hallal/Desktop/PhD/new_scripts/results_test_sm/OmicsTools_0.0.22.tar.gz", repos = NULL, type = "source")
library(OmicsTools)

library(yaml)

## Load parameters file
params <- read_yaml("./config.yaml")

## cell line
cell_line <- params$cell_line
imp <- params$Imputation

## Define input file
inputFile <- paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/eSet_data_test1_",params$cell_line,".csv")

## Annotation function
aDataFrame <- function(input) {
  metaData = data.frame(labelDescription=colnames(input), 
                        row.names=colnames(input))
  new("AnnotatedDataFrame", data=input, varMetadata=metaData)
}

## Read the table content
rawData2 = read.table(inputFile, sep="\t", header=TRUE, 
                      stringsAsFactors = FALSE, quote="", comment.char = "", check.names = FALSE)

## Define the expression set (intensities) here they all start with X (remove proteins)
expression2 = as.matrix(rawData2[,grep(cell_line, colnames(rawData2))])#21.01.2019: remove X from grep
#expression2[expression2 == 0] <- 1
#expression2 <- log2(expression2)

## reference#
#tt2 <- expression2[order(rowVars(expression2), decreasing = F),]
#tt2 <- as.data.frame(tt2)

#select top1
#tt2 <- tt2[1,]
#expression3 <- as.data.frame(expression2)/as.data.frame(tt2)
#expression3 <- sweep(expression2, 2, tt2, '-')
#expression2 <- expression3
#with no X
#expression2 = as.matrix(rawData2[,-17])

## Define rowNames as phosphorylation positions
rownames(expression2) <- rawData2$Proteins

## PhenoTable create
phenoTable2 = data.frame(
  group = sub("X(\\w+)", "\\1", colnames(expression2)),
  replicates = sub("X(\\d+)_(\\w+)", "\\1_\\2", colnames(expression2)),
  cell_line = sub("X(\\d+)_(\\w+)", "\\2", colnames(expression2)),
  batch = sub("X(\\d+)_(.+)", "\\1", colnames(expression2)),#changed 21.01.2019
  row.names = colnames(expression2))

phenoTable2$cell_line = as.character(phenoTable2$cell_line)
phenoTable2$group = as.character(phenoTable2$group)
phenoTable2$replicates = as.character(phenoTable2$replicates)
phenoTable2$batch = as.character(phenoTable2$batch)

## Create the ExpresionSet
test_eSet <- ExpressionSet(assayData = expression2, aDataFrame(phenoTable2))

## Save the expressionSet
save(test_eSet, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet_notfiltered_",params$cell_line,".Rda"))

