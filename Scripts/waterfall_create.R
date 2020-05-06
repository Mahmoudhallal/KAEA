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

library(dplyr)

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Define contrasts and control
contrasts <- unlist(strsplit(params$Conditions,","))
control <- params$control
if (params$SILAC == "T"){
  samples = contrasts
} else {
  samples <- contrasts[!(contrasts == control)]
}

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

## Get Oncogene/tumor suppressor CancerMine database
cancerMine <- read.delim('../Phospho_DBs/cancermine_collated.tsv')
leukemia_genes <- grep('leukemia',cancerMine$cancer_normalized)
cancerMine_leukemia <- cancerMine[leukemia_genes,]
TS_leukemia <- cancerMine_leukemia[cancerMine_leukemia$role == "Tumor_Suppressor",]
ONC_leukemia <- cancerMine_leukemia[cancerMine_leukemia$role == "Oncogene",]
DR_leukemia <- cancerMine_leukemia[cancerMine_leukemia$role == "Driver",]
#Intersection between TS and ONC
Intersection <- intersect(ONC_leukemia$gene_normalized,TS_leukemia$gene_normalized)

## NO Overlapping kkinases
### Heatmap
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
  #direcion
  for (y in 1:nrow(material_for_waterfall)){if (material_for_waterfall[[conditions[[x]]]][y] > 0){material_for_waterfall$Color[y] <- "red"} else {material_for_waterfall$Color[y] <- "blue"}}
  
  #hjust
  for (y in 1:nrow(material_for_waterfall)){if (material_for_waterfall[[conditions[[x]]]][y] > 0){material_for_waterfall$hjust[y] <- "1"} else {material_for_waterfall$hjust[y] <- "0"}}
  
  #Category: TS, ONC or Driver
  for (y in 1:nrow(material_for_waterfall)){
    if (material_for_waterfall[['rowname']][y] %in% TS_leukemia$gene_normalized){
      material_for_waterfall$Category[y] <- "TS"
      } else if (material_for_waterfall[['rowname']][y] %in% ONC_leukemia$gene_normalized) {
        material_for_waterfall$Category[y] <- "ONC"
        } else if (material_for_waterfall[['rowname']][y] %in% DR_leukemia$gene_normalized){
          material_for_waterfall$Category[y] <- "DR"
          } else if (material_for_waterfall[['rowname']][y] %in% Intersection){
            material_for_waterfall$Category[y] <- "POTSF"
          } else {
            material_for_waterfall$Category[y] <- ""
          }
    }
  
  colnames(material_for_waterfall) <- c("Kinase","Enrichment.score","Color","hjust","category")
  material_for_waterfall <- material_for_waterfall[order(material_for_waterfall$Enrichment.score,decreasing=F),]
  material_for_waterfall
})

names(both) <- conditions
material_for_waterfall <- both

## Save the merged dataframe as dataframe for Waterfall
save(material_for_waterfall, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Material_for_waterfall_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".Rda"))


## Plot waterfall (barplot)
for (x in names(both)){
  outputFile = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Waterfall_",params$pvalue_cutoff,"P_",params$fdr_cutoff,'FDR_', x,'.pdf')
  
  pdf(outputFile, width=5, height=8)
  
  water_fall <- ggplot(data=both[[x]], aes(x=reorder(Kinase,Enrichment.score), y=Enrichment.score)) +
    geom_bar(stat="identity", fill=both[[x]]$Color) +
    coord_flip() +
    ggtitle("Kinase enrichments") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24)) +
    theme(axis.text.y = element_text(size = 14),axis.title=element_text(size=22)) +
    xlab("Kinase") +
    ylab("-log10 p-value")+
    geom_text(aes(label=both[[x]]$category), hjust = as.numeric(both[[x]]$hjust), color="white", position = position_dodge(1), size=3.5)

  
  plot(water_fall)
  
  dev.off()
}

