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
# library(dplyr)
# library(plotly)
# detach("package:dplyr", unload=TRUE)
# detach("package:plotly", unload=TRUE)
## Load parameters
params <- read_yaml("./config.yaml")

## Load material for heatmap
Heatmap_rep <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Material_for_heatmap_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".csv"), row.names = 1,sep=" ")

## Load material for waterfall plots
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Material_for_waterfall_", params$fdr_cutoff,"FDR_", params$pvalue_cutoff,"P_",params$cell_line,".Rda"))

## Create table of pathways using waterfall input
to_sbmt <- lapply(1:length(material_for_waterfall), function(x){
  df1 <- material_for_waterfall[[x]][!duplicated(material_for_waterfall[[x]]$Kinase),]
  #make dataframe
  df <- as.data.frame(df1$Enrichment.score)
  rownames(df) <- df1$Kinase
  colnames(df) <- "value"
  #dd <- df$value
  #names(dd) <- gsub('(\\w+)[+-]','\\1',rownames(df))
  #dd1 <- dd[!duplicated(names(dd))]
  dd1 <- df
  dd2 <- as.data.frame(dd1)
  #rownames(dd2) <- names(dd1)
  #colnames(dd2) <- "value"
  dd2
})
names(to_sbmt) <- names(material_for_waterfall)
files <- names(to_sbmt)


## Define directory and KEGG pathways
directory <- paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"pathview")
unlink(directory,recursive = TRUE)
dir.create(directory)
pathways = c("hsa04140","hsa04150","hsa05221","hsa05200","hsa04014","hsa04010","hsa04012","hsa04310","hsa04350","hsa04630","hsa04151","hsa05220")
pathways_names <- c('Autophagy animal','mTOR Signaling Pathway','Acute Myeloid Leukemia','Pathways in Cancer','RAS Signaling Pathway',
                    'MAPK Signaling Pathway','ERBB Signaling Pathway','WNT Signaling Pathway','TGF-BETA Signaling Pathway','JAK-STAT Signaling Pathway',
                    'PI3K-AKT Signaling Pathway','Chronic Myeloid Leukemia')
## Create all pathview.png for all pathways and conditions 
all_files <- lapply(1:length(files), function(y){
  sp_directory <- paste0(directory,'/',files[y],'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR')
  dir.create(sp_directory)
  table_con <-  lapply(1:length(pathways), function(x) {
    pp <- pathview(gene.data = to_sbmt[[files[y]]], pathway.id = pathways[x],
             species = "hsa", gene.idtype = "SYMBOL", 
             limit = list(gene = 1, cpd = 1),
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

#Give pathway names
for(x in 1:length(enriched_kin)){names(enriched_kin[[x]]) <- pathways}
names(enriched_kin) <- files

## Add pathway names and reurun final table for every condition
enriched_kin1 <- lapply(1:length(enriched_kin), function(x){
  df1 <- t(as.data.frame(enriched_kin[[x]]))
  df2 <- as.data.frame(cbind(rownames(df1),df1))
  df2$Pathway <- c('Autophagy animal','mTOR Signaling Pathway','Acute Myeloid Leukemia','Pathways in Cancer','RAS Signaling Pathway',
                                  'MAPK Signaling Pathway','ERBB Signaling Pathway','WNT Signaling Pathway','TGF-BETA Signaling Pathway','JAK-STAT Signaling Pathway',
                                  'PI3K-AKT Signaling Pathway','Chronic Myeloid Leukemia')
  df2 <- df2[,c(3,1,2)]
  colnames(df2) <- c('Pathway','KEGG pathway','#Enriched kinases')
  df2$`#Enriched kinases` <- as.numeric(as.character(df2$`#Enriched kinases`))
  df2 <- df2[order(df2$'#Enriched kinases',decreasing=TRUE),]
  df2
  })

names(enriched_kin1) <- files
material_for_table <- enriched_kin1

## Load Quality report material
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Quality_Report_Shiny_",params$cell_line,".Rda"))

## Load Volcano plot material
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Volcano_plot_material_",params$cell_line,".Rda"))

## LoadExpressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"test_eSet1_",params$cell_line,".Rda"))

## Load Coefficient of variation
CVs <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"CV_",params$cell_line,".csv"))

## Save the heatmap, waterfall and table material in a list format
shiny_results <- list(Quality_rep$hist_sites, Quality_rep$hist_peps, Quality_rep$hist_prots, Quality_rep$venn_sites, Quality_rep$venn_peps, Quality_rep$venn_prots,
                      Quality_rep$table_counts, Heatmap_rep, material_for_waterfall, material_for_table, volcano_material, exprs(test_eSet),CVs)
names(shiny_results) <- c(names(Quality_rep),'Heatmap_rep', 'waterfall_rep','table_rep','Volcano','ExpressSet','CV')
save(shiny_results, file = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$pvalue_cutoff,"P_",params$fdr_cutoff,"FDR.Rda"))

#