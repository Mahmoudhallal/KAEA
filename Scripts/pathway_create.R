##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Produce heat map
## Date: 09.03.17
## Author: Cedric Simillion
##################################################

## Load libraries
#install.packages("gplots", repos="http://cran.rstudio.com/")
library(gplots)
#install.packages("pheatmap", repos="http://cran.rstudio.com/")
library(pheatmap)   
library(ggplot2)
library(RColorBrewer)
library(yaml)
library(Rcpp)

#install.packages("dplyr", repos="http://cran.rstudio.com/")
library(dplyr)

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation

#######################################################
#######################################################
#######################################################
### We changed the heatmap so that it is prduced using p-values of the waterfall and not the scores of setRank
### The reason is that there was some differences in the order of kinases according to scores or corrected p values
### this is changed on 22.02.2019

## Load input topTables of under-expressed phosphosites
# inputFile_less = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_less/pathways.txt")
# outputFile2 = paste0(params$CWD,"/results/",params$cell_line,"_",params$pvalue_cutoff,"P_",params$fdr_cutoff,"FDR/","Heatmap_",params$cell_line,"_" ,params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR.pdf',sep='')
# 
# results_less = read.table(inputFile_less, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
# pathwayMatrix_less = as.matrix(results_less[,1:(ncol(results_less)-3)])
# rownames(pathwayMatrix_less) <- results_less$description
# 
# ## Load input topTables of over-expressed phosphosites
# inputFile1_greater = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_greater/pathways.txt")
# 
# results_greater = read.table(inputFile1_greater, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
# pathwayMatrix_greater = as.matrix(results_greater[,1:(ncol(results_greater)-3)])
# rownames(pathwayMatrix_greater) <- results_greater$description
# 
# 
# pathwayMatrix_less <- -pathwayMatrix_less
# pathwayMatrix <- rbind(pathwayMatrix_greater, pathwayMatrix_less)
# 
# ## Sum kinase values if they are present in over- and under- enriched analysis 
# dd1 <- as.data.frame(pathwayMatrix_less)
# dd2 <- as.data.frame(pathwayMatrix_greater)
# 
# ll <- bind_rows(dd1 %>% add_rownames(), 
#                 dd2 %>% add_rownames()) %>% 
#   # evaluate following calls for each value in the rowname column
#   group_by(rowname) %>% 
#   # add all non-grouping variables
#   summarise_all(sum)
# cc <- as.data.frame(ll)
# pathwayMatrix <- cc
# rownames(pathwayMatrix) <- cc$rowname
# 
# contrasts <- unlist(strsplit(params$Conditions,","))
# control <- params$control
# if (params$SILAC == "T"){
#   samples = contrasts
# } else {
#   samples <- contrasts[!(contrasts == control)]
# }
# 
# colnames(pathwayMatrix) <- c('rowname',samples)

##################################################
##################################################
##################################################
## Load parameters
params <- read_yaml("./config.yaml")
outputFile2 = paste0(params$CWD,"/results/",params$cell_line,"_",params$pvalue_cutoff,"P_",params$fdr_cutoff,"FDR_imp",imp,"/Heatmap_",params$cell_line,"_" ,params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR.pdf',sep='')

contrasts <- unlist(strsplit(params$Conditions,","))
control <- params$control
if (params$SILAC == "T"){
  samples = contrasts
} else {
  samples <- contrasts[!(contrasts == control)]
}

#conditions <- strsplit(params$Conditions, split=",")[[1]]
conditions <- samples

## Load input file from underactive kianses
pathwayMatrix_less <- lapply(1:length(conditions), function(x){
  nme <- paste0('inputFile_',strsplit(params$Conditions, split=",")[[1]][x])
  assign(nme, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_less/",conditions[x] ,"_pathways.txt"))
  results = read.table(get(nme), sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
  if (nrow(results) != 0){
    pathwayMatrix = as.matrix(results[,c("correctedPValue")])
    rownames(pathwayMatrix) <- results$description
    colnames(pathwayMatrix) <- conditions[x]
    pathwayMatrix <- as.data.frame(pathwayMatrix)
    pathwayMatrix$color <- "blue"
    pathwayMatrix$kinase <- rownames(pathwayMatrix)
    pathwayMatrix <- pathwayMatrix[order(pathwayMatrix[[conditions[x]]], decreasing = T),]
    pathwayMatrix[[conditions[x]]]<- -log10(pathwayMatrix[[conditions[x]]])
    pathwayMatrix[[conditions[x]]] <- round(pathwayMatrix[[conditions[x]]], digits = 5)
    pathwayMatrix[[conditions[x]]] <- -pathwayMatrix[[conditions[x]]]
    pathwayMatrix
  } else{
    pathwayMatrix <- data.frame()
    pathwayMatrix
  }
})

## Load the input file for overactive kinases
pathwayMatrix_greater <- lapply(1:length(conditions), function(x){
  nme <- paste0('inputFile_',strsplit(params$Conditions, split=",")[[1]][x])
  assign(nme, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_greater/",conditions[x] ,"_pathways.txt"))
  results = read.table(get(nme), sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
  
  if (nrow(results) != 0){
    pathwayMatrix = as.matrix(results[,c("correctedPValue")])
    rownames(pathwayMatrix) <- results$description
    colnames(pathwayMatrix) <- conditions[x]
    pathwayMatrix <- as.data.frame(pathwayMatrix)
    pathwayMatrix$color <- "red"
    pathwayMatrix$kinase <- rownames(pathwayMatrix)
    pathwayMatrix <- pathwayMatrix[order(pathwayMatrix[[conditions[x]]], decreasing = T),]
    pathwayMatrix[[conditions[x]]]<- -log10(pathwayMatrix[[conditions[x]]])
    pathwayMatrix[[conditions[x]]] <- round(pathwayMatrix[[conditions[x]]], digits = 5)
    
    pathwayMatrix
  } else {
    pathwayMatrix <- data.frame()
    pathwayMatrix
  }
})

### Heatmap
# dd1 <- lapply(1:length(conditions), function(x){
#   as.data.frame(pathwayMatrix_less[[x]][,grep(params$cell_line,colnames(pathwayMatrix_less[[x]])),drop=FALSE])
# })
# #dd1 <- as.data.frame(dd1$MOLM13_DRG,drop=FALSE)
# dd2 <- lapply(1:length(conditions), function(x){
#   as.data.frame(pathwayMatrix_greater[[x]][,grep(params$cell_line,colnames(pathwayMatrix_greater[[x]])),drop=FALSE])
# })
# 
# ll <- lapply(1:length(conditions), function(x){
#   bind_rows(dd1[[x]] %>% add_rownames(), 
#             dd2[[x]] %>% add_rownames()) %>% 
#     # evaluate following calls for each value in the rowname column
#     group_by(rowname) %>% 
#     # add all non-grouping variables
#     summarise_all(sum)
# })
#   
# cc <- lapply(1:length(conditions), function(x){
#   as.data.frame(ll[[x]])
# })
# pathwayMatrix <- cc
# rownames(pathwayMatrix) <- cc$rowname

########
########
########
both <- lapply(1:length(pathwayMatrix_greater), function(x){
  dd1 <- as.data.frame(pathwayMatrix_less[[x]][,grep(params$cell_line,colnames(pathwayMatrix_less[[x]])),drop=FALSE])
  #dd1 <- as.data.frame(dd1$MOLM13_DRG,drop=FALSE)
  dd2 <- as.data.frame(pathwayMatrix_greater[[x]][,grep(params$cell_line,colnames(pathwayMatrix_greater[[x]])),drop=FALSE])
  ll <- bind_rows(dd1 %>% add_rownames(), 
                  dd2 %>% add_rownames()) %>% 
    # evaluate following calls for each value in the rowname column
    group_by(rowname) %>% 
    # add all non-grouping variables
    summarise_all(sum)
  cc <- as.data.frame(ll)
  material_for_waterfall <- cc
  rownames(material_for_waterfall) <- cc$rowname
  for (y in 1:nrow(material_for_waterfall)){if (material_for_waterfall[[conditions[[x]]]][y] > 0){material_for_waterfall$Color[y] <- "red"} else {material_for_waterfall$Color[y] <- "blue"}}
  colnames(material_for_waterfall) <- c("Kinase","Enrichment.score","Color")
  material_for_waterfall <- material_for_waterfall[order(material_for_waterfall$Enrichment.score,decreasing=F),]
  material_for_waterfall
})

names(both) <- conditions
material_for_heatmap <- both


## Create output file for dataframe to be used for heatmap
save(material_for_heatmap, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Material_for_heatmap_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".Rda"))

## Plot Heatmap
for (x in names(both)){
  outputFile = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Heatmap_", params$pvalue_cutoff,"P_",params$fdr_cutoff,'FDR_', x,'.pdf')
  pathwayMatrix1 <- as.data.frame(both[[x]][,c('Enrichment.score'),drop=F])
  
  
  ### New 21.01.2020 from 
  moitie = nrow(pathwayMatrix1)/2
  Min = min(pathwayMatrix1[,1])
  Max = max(pathwayMatrix1[,1])
  Center = 0
  
  ## Negative values
  neg = colorRampPalette(colors = c("blue", "white"), space="Lab")(moitie)    
  ## POsitive values
  pos = colorRampPalette(colors = c("white", "red"), space="Lab")(moitie)
  rampcols = c(neg, pos)
  ## In your example, this line sets the color for values between 49 and 51. 
  rampcols[c(moitie, moitie+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
  
  rb1 = seq(Min, Center, length.out=moitie+1)
  rb2 = seq(Center, Max, length.out=moitie+1)[-1]
  rampbreaks = c(rb1, rb2)
  colnames(pathwayMatrix1) <- x
  
  ###
  pdf(outputFile, width=4, height=8)
  pheatmap(pathwayMatrix1, breaks=rampbreaks,color = rampcols, cluster_rows = T, cluster_cols = F, fontsize_row = 9, fontsize_col = 12, onefile=TRUE)
  dev.off()
}


