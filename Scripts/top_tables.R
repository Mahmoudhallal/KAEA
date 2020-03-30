##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Create table with p-values from one tail (greater) t-test (cell line x vs all other cell lines)
## Date: 09.07.17
## Author: Mahmoud Hallal
##################################################

##Load libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library(Biobase)

#install.packages("ggplot2", repos="http://cran.rstudio.com/")#
library(ggplot2)

library(yaml)

## Load parameters
params <- read_yaml("./config.yaml")
biological_reps <- strsplit(params$Biological_replicates,"-")[[1]]
imp <- params$Imputation

## TopTables
if (params$SILAC == "T"){
  
  ## calculate logFC
  getLogFC <- function(row, phenoTable, groupA, groupVar) {
    indicesA = phenoTable[[groupVar]] == groupA
    mean(row[indicesA]) 
    
  }
  
  #calculate p-value
  getPValue <- function(row, phenoTable, groupA, groupVar, ttest) {
    row = row + runif(length(row), -1e-10, 1e-10)
    indicesA = phenoTable[[groupVar]] == groupA
    
    print(ttest)
    t.test(row[indicesA],alternative = ttest)$p.value
  }
  
  
  #function to create a table of p-value, p-adj and logFC
  #group A is the desired cell line and groupB is the control
  getTopTable <- function(samples, control, eSet, groupVar, x, ttest) {
    groupA = samples[x]
    
    #when only selecting the proteins that are present in every cell line
    topTable = data.frame()
    
    eSet2 <- eSet[rowSums(exprs(eSet)[,pData(eSet)[[groupVar]] == groupA]) != 0,]
    #eSet2 <- eSet
    topTable = apply(exprs(eSet2), 1, getLogFC, pData(eSet2), groupA, groupVar)
    topTable <- as.data.frame(topTable)
    topTable$p = apply(exprs(eSet2), 1, getPValue, pData(eSet2), groupA, groupVar, ttest)
    topTable$padj = p.adjust(topTable$p, method="fdr")
    topTable = topTable[order(topTable$padj, topTable$p),]
    topTable
    }
  } else {
    
    ## calculate logFC
    getLogFC <- function(row, phenoTable, groupA, groupB, groupVar) {
      indicesA = phenoTable[[groupVar]] == groupA
      indicesB = phenoTable[[groupVar]] == groupB
      mean(row[indicesA], na.rm=T) - mean(row[indicesB], na.rm=T)
      }

  ## calculate p-value 
  getPValue <- function(row, phenoTable, groupA, groupB, groupVar, ttest) {
      row = row + runif(length(row), -1e-10, 1e-10)#
      indicesA = phenoTable[[groupVar]] == groupA
      indicesB = phenoTable[[groupVar]] == groupB
      
      t.test(row[indicesA], row[indicesB], alternative = ttest)$p.value
      }


  ## function to create a table of p-value, p-adj and logFC
  ## group A is the desired cell line and groupB is the control
  getTopTable <- function(samples, control, eSet, groupVar, x, ttest, parm) {
    
    #Define conditions  
    groupA = samples[x]
    if (control != 'ALL'){
      groupB = control
    } else {
      groupB = samples[x] 
    } 
    
    #Define fData
    topTable = fData(eSet)
    
    #Define dataset
    eSet_modified <- exprs(eSet)
    
    #Calculate FC
    if (biological_reps[x] > 1){
      topTable = apply(eSet_modified, 1, getLogFC, pData(eSet), groupA, groupB, groupVar)
      topTable <- as.data.frame(topTable)
      
      #Calculate p-value
      topTable$p = apply(eSet_modified, 1, getPValue, pData(eSet), groupA, groupB, groupVar, ttest)
      topTable$padj = p.adjust(topTable$p, method="fdr")
      topTable = topTable[order(topTable$padj, topTable$p),]
      #topTable
      if (parm == "less"){
        topTable_up <- topTable[topTable$topTable > 0,]
        topTable_up = topTable_up[order(topTable_up$padj, topTable_up$p,decreasing = TRUE),]
        topTable_down <- topTable[topTable$topTable < 0,]
        topTable_down = topTable_down[order(topTable_down$padj, topTable_down$p,decreasing = FALSE),]
        topTable <- rbind(topTable_down,topTable_up)
      } else {
        topTable_up <- topTable[topTable$topTable > 0,]
        topTable_up = topTable_up[order(topTable_up$padj, topTable_up$p,decreasing = FALSE),]
        topTable_down <- topTable[topTable$topTable < 0,]
        topTable_down = topTable_down[order(topTable_down$padj, topTable_down$p,decreasing = TRUE),]
        topTable <- rbind(topTable_up,topTable_down)
        }
      topTable
    } else {
      if (ttest == "less"){
        topTable = apply(exprs(eSet), 1, getLogFC, pData(eSet), groupA, groupB, groupVar)
        topTable <- as.data.frame(topTable)
        topTable = topTable[order(topTable$topTable,decreasing = FALSE),,drop=FALSE]
        topTable
      } else {
        topTable = apply(exprs(eSet), 1, getLogFC, pData(eSet), groupA, groupB, groupVar)
        topTable <- as.data.frame(topTable)
        topTable = topTable[order(topTable$topTable,decreasing = TRUE),,drop=FALSE]
        topTable
      }
    }
    topTable
    }
}

## Load expressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

## Define contrasts
contrasts <- unlist(strsplit(params$Conditions,","))
control <- params$control
if (params$SILAC == "T"){
  samples = contrasts
} else {
  samples <- contrasts[!(contrasts == control)]
}

## Create topTables for over- and unver-expressed phosphosites
topTables_less <- lapply(1:length(samples), function(x) getTopTable(samples, control, test_eSet, "cell_line", x, "two.sided","less"))
topTables_greater <- lapply(1:length(samples), function(x) getTopTable(samples, control,test_eSet, "cell_line", x, "two.sided","greater"))

## Name the topTables according to conditions
names(topTables_less) <- samples
names(topTables_greater) <- samples

## Save output files separately
save(topTables_less, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_less.Rda"))
save(topTables_greater, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_greater.Rda"))
