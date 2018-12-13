##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: PRoduce heat map
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

## Load parameters
params <- read_yaml("./config.yaml")

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
  assign(nme, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_less/",conditions[x] ,"_pathways.txt"))
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
  assign(nme, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"output_networks_with_positions_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P_",params$cell_line,"_greater/",conditions[x] ,"_pathways.txt"))
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


## Merge over- and under- active kinases
both <- lapply(1:length(pathwayMatrix_greater), function(x){
  df <- rbind(pathwayMatrix_less[[x]],pathwayMatrix_greater[[x]])
  colnames(df) <- c("Enrichment.score","Color","Kinase")
  df
})

names(both) <- conditions
material_for_waterfall <- both

## Save the merged dataframe as dataframe for Waterfall
save(material_for_waterfall, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Material_for_waterfall_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".Rda"))


## Plot waterfall (barplot)
for (x in names(both)){
  outputFile = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Waterfall_",params$pvalue_cutoff,"P_",params$fdr_cutoff,'FDR_', x,'.pdf')
  
  pdf(outputFile, width=5, height=8)
  
  water_fall <- ggplot(data=both[[x]], aes(x=reorder(Kinase,Enrichment.score), y=Enrichment.score)) +
    geom_bar(stat="identity", fill=both[[x]]$Color) +
    coord_flip() +
    ggtitle("Kinase enrichments") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 18)) +
    theme(axis.text.y = element_text(size = 12)) +
    xlab("Kinase") +
    theme(#axis.text=element_text(size=14),
      axis.title=element_text(size=20))#
  
  plot(water_fall)
  
  dev.off()
}
