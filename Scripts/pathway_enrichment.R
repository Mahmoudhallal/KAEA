#################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Apply pathway enrichment on identified proteins
## Date: 08.02.2019
## Author: Mahmoud Hallal
##################################################
#GSEA
#install.packages("devtools", repos="http://cran.rstudio.com/")
#library(devtools)
#devtools::install_github('C3c6e6/SetRank',force=TRUE)
#devtools::install_github('cran/SetRank')

#install.packages("SetRank", repos="http://cran.rstudio.com/")
library(SetRank)
library(igraph)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Homo.sapiens")
library(Homo.sapiens)

#utils::install.packages('/Users/Mahmoud.Hallal/Desktop/PhD/new_scripts/GeneSets/GeneSets.Homo.sapiens_16.5.3.tar.gz')
library(GeneSets.Homo.sapiens)

## Load libraries
library("reshape2")
#install.packages("yaml", repos="http://cran.rstudio.com/")
library(yaml)
#library(plotly)
#library(dplyr)
#detach("package:reshape2", unload=TRUE)
#detach("package:dplyr", unload=TRUE)

## Load parameters file
params <- read_yaml("./config.yaml")
samples <- unlist(strsplit(params$Samples,","))
imp <- params$Imputation

## Load expressionSet 
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))
#load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Volcano_plot_material_",params$cell_line,".Rda"))
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_less.Rda"))
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_greater.Rda"))

#define reference
union_all_phos_prots <- unique(rownames(exprs(test_eSet)))
union_all_phos_prots <- gsub("(\\w+)(_[STY]\\d+)", "\\1" ,union_all_phos_prots)
union_all_phos_prots <- unique(union_all_phos_prots)

## Define uniprot to entrezid converter
uniprot2EntrezID <- createIDConverter("Homo.sapiens", "UNIPROT","ENTREZID")

## COnvert reference
reference_union_all_phos_prots <- uniprot2EntrezID(union_all_phos_prots)

## COnvert desired dataset
front_proteins <- lapply(1:length(samples), function(x){
  #less
  ttl <- topTables_less[[samples[[x]]]][topTables_less[[samples[[x]]]]$topTable < 0,]
  proteins_test_less <- unique(gsub("(\\w+)(_[STY]\\d+)", "\\1" ,rownames(ttl)))
  proteins_test_less <- uniprot2EntrezID(proteins_test_less)
  #greater
  ttl <- topTables_greater[[samples[[x]]]][topTables_greater[[samples[[x]]]]$topTable > 0,]
  proteins_test_greater <- unique(gsub("(\\w+)(_[STY]\\d+)", "\\1" ,rownames(topTables_greater[[samples[[x]]]])))
  proteins_test_greater <- uniprot2EntrezID(proteins_test_greater)
  
  list(proteins_test_less, proteins_test_greater)
})

# test_MOLM13_DRG <- unique(gsub("(\\w+)(_[STY]\\d+)", "\\1" ,rownames(volcano_material[samples])))
# test_MOLM13_DRG <- uniprot2EntrezID(test_MOLM13_DRG)

#create the collection
collection_union_all_phos_prots <- buildSetCollection(KEGG , referenceSet = reference_union_all_phos_prots, maxSetSize = 500)
collection_pathways <- collection_union_all_phos_prots
#saveRDS(collection_union_all_phos_prots, "collection_allDBs.rds")
#collection_pathways <- readRDS('~/Desktop/PhD/new_scripts/collection_union_all_phos_prots.rds')

#PATHWAYS
#cell_lines <- params$cell_line
options(mc.cores=2)
enriched_pathways <- lapply(1:length(samples), function(x){
  less_greater <- lapply(1:length(front_proteins[[x]]), function(y){
    res <- setRankAnalysis(front_proteins[[x]][[y]], collection_pathways, use.ranks= TRUE, setPCutoff = 0.05, fdrCutoff = 0.05)
    enriched_pathways <- as_data_frame(res, what="vertices")
    enriched_dbs <- unique(enriched_pathways$database)
    separate_dbs <- lapply(1:length(enriched_dbs), function(x) {enriched_pathways[enriched_pathways$database == enriched_dbs[x],]})
    #View(enriched_pathways)
    #exportSingleResult(res, get(paste0('test_',cell_lines[x])), collection_pathways, paste0("results_",cell_lines[x]), IDConverter, paste0("./", cell_lines[x], "_enrichments_pathways_0.01"))
    enriched_pathways
    })
  names(less_greater) <- c("less","greater")
  less_greater
  })

names(enriched_pathways) <- samples

save(enriched_pathways,  file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Enriched_pathways_", params$cell_line,".Rda"))


load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Enriched_pathways_", params$cell_line,".Rda"))
enriched_pathways$MOLM13H_RD$less
enriched_pathways$MOLM13H_RD$greater
