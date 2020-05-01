##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Prepare output file for Shiny app
## Date: 21.09.2018
## Author: Mahmoud Hallal
##################################################

## Load libraries
#source("http://bioconductor.org/biocLite.R")
#biocLite("pathview")
library(pathview)

#install.packages("filesstrings", repos="http://cran.rstudio.com/")
library(filesstrings)
data("gene.idtype.bods")

library(yaml)

#library(dplyr)

#library(plotly)

#detach("package:dplyr", unload=TRUE)

#detach("package:plotly", unload=TRUE)

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Load material for heatmap
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Material_for_heatmap_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".Rda"))

## Load material for waterfall plots
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Material_for_waterfall_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".Rda"))

## Create table of pathways using waterfall input
to_sbmt <- lapply(1:length(material_for_waterfall), function(x){
  df1 <- material_for_waterfall[[x]][!duplicated(material_for_waterfall[[x]]$Kinase),]
  
  #make dataframe
  df <- as.data.frame(df1$Enrichment.score)
  rownames(df) <- df1$Kinase
  colnames(df) <- "value"
  dd1 <- df
  dd2 <- as.data.frame(dd1)
  
  dd2
})

names(to_sbmt) <- names(material_for_waterfall)
files <- names(to_sbmt)


## Define directory and KEGG pathways
directory <- paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/pathview")
unlink(directory,recursive = TRUE)
dir.create(directory)
if (params$Species == "Human"){
  pathways = c("hsa04140","hsa04150","hsa05221","hsa05200","hsa04014","hsa04010","hsa04012","hsa04310","hsa04350","hsa04630","hsa04151","hsa05220")
  pathways_names <- c('Autophagy animal','mTOR Signaling Pathway','Acute Myeloid Leukemia','Pathways in Cancer','RAS Signaling Pathway',
                      'MAPK Signaling Pathway','ERBB Signaling Pathway','WNT Signaling Pathway','TGF-BETA Signaling Pathway','JAK-STAT Signaling Pathway',
                      'PI3K-AKT Signaling Pathway','Chronic Myeloid Leukemia')
  } else {
  pathways = c("sce04138","sce04011","sce04392")
  pathways_names <- c('Autophagy animal','MAPK Signaling Pathway','Hippo signaling pathway')
}
  
## Create all pathview.png for all pathways and conditions 
all_files <- lapply(1:length(files), function(y){
  sp_directory <- paste0(directory,'/',files[y],'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR')
  dir.create(sp_directory)
  ## assign the species according to input, this is for KEGG mapping
  if (params$Species == "Human"){
    sp = "hsa"
    type = "SYMBOL"
    
  } else {
    sp = "sce"
    type = "GENENAME"
    
  }
  table_con <-  lapply(1:length(pathways), function(x) {
    pp <- pathview(gene.data = to_sbmt[[files[y]]], pathway.id = pathways[x],
                   #we changed gene.idtype from SYMBOL to GENENAME 22.01.2019
             species = sp, gene.idtype = type, 
             limit = list(gene = c(min(to_sbmt[[files[y]]]), max(to_sbmt[[files[y]]])), cpd = 1),
             kegg.dir = sp_directory)
    
    file.move(paste(pathways[x], ".pathview.png", sep = ""), sp_directory, overwrite = TRUE)
    file.copy(paste0(sp_directory,'/',pathways[x],'.pathview.png'), paste0(sp_directory,'/',pathways_names[x],'.pathview.png'))
    pp
    })
})

## Give the file names
names(all_files) <- files

## Count the number of kinases present in every pathway
enriched_kin <- lapply(1:length(files), function(y){ lapply(1:length(all_files[[files[y]]]), function(x){sum(!is.na(all_files[[files[y]]][[x]]$plot.data.gene$ge1))})})

## Give pathway names
for(x in 1:length(enriched_kin)){names(enriched_kin[[x]]) <- pathways}
names(enriched_kin) <- files

## Add pathway names and reurun final table for every condition
enriched_kin1 <- lapply(1:length(enriched_kin), function(x){
  df1 <- t(as.data.frame(enriched_kin[[x]]))
  df2 <- as.data.frame(cbind(rownames(df1),df1))
  if (params$Species == "Human"){
  df2$Pathway <- c('Autophagy animal','mTOR Signaling Pathway','Acute Myeloid Leukemia','Pathways in Cancer','RAS Signaling Pathway',
                                  'MAPK Signaling Pathway','ERBB Signaling Pathway','WNT Signaling Pathway','TGF-BETA Signaling Pathway','JAK-STAT Signaling Pathway',
                                  'PI3K-AKT Signaling Pathway','Chronic Myeloid Leukemia')
  } else {
    df2$Pathway <- c('Autophagy animal','MAPK Signaling Pathway','Hippo signaling pathway')
  }
  df2 <- df2[,c(3,1,2)]
  colnames(df2) <- c('Pathway','KEGG pathway','#Enriched kinases')
  df2$`#Enriched kinases` <- as.numeric(as.character(df2$`#Enriched kinases`))
  df2 <- df2[order(df2$'#Enriched kinases',decreasing=TRUE),]
  df2
  })


names(enriched_kin1) <- files
material_for_table <- enriched_kin1


## LoadExpressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

## Load Quality report material
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Quality_Report_Shiny_",params$cell_line,".Rda"))

## Load Volcano plot material
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Volcano_plot_material_",params$cell_line,".Rda"))

## Load Coefficient of variation
CVs <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/CV_",params$cell_line,".csv"))

## Load special topTables for volcano plot of barplot
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Volcano_plot_material_",params$cell_line,"_special.Rda"))

## Get the kinase-substrate database
db <- read.delim(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,'/all_dbs.csv'),sep=',')

## Get the enriched pathways
#load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,'/Enriched_pathways_',params$cell_line,'.Rda'))


## Save the heatmap, waterfall and table material in a list format
shiny_results <- list(Quality_rep$hist_sites, Quality_rep$hist_peps, Quality_rep$hist_prots, Quality_rep$venn_sites, Quality_rep$venn_peps, Quality_rep$venn_prots,
                      Quality_rep$table_counts, material_for_heatmap, material_for_waterfall, material_for_table, volcano_material, volcano_material_special, exprs(test_eSet),CVs,db)
names(shiny_results) <- c(names(Quality_rep),'Heatmap_rep', 'waterfall_rep','table_rep','Volcano','Volcano_special','ExpressSet','CV','DB')
save(shiny_results, file = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/results_shiny_",params$cell_line,"_",params$pvalue_cutoff,"P_",params$fdr_cutoff,"FDR.Rda"))
