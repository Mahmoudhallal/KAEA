#################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Apply pathway enrichment using topGO package on identified phosphosites
## Date: 24.07.2019
## Author: Mahmoud Hallal
##################################################

#load packages
install.packages("SetRank", repos="http://cran.rstudio.com/")
library(topGO)
library(plyr)
library(dplyr)
#library(SetRank)
#library(DESeq2) # 1.22.2
library(topGO) # 2.34.0
library(yaml)
library(biomaRt)
library("org.Hs.eg.db")

library(yaml)

## Load parameters file
params <- read_yaml("./config.yaml")
samples <- unlist(strsplit(params$Samples,","))
imp <- params$Imputation
#samples <- 'MOLM13H_RD'




RunTopGO <- function(fdr.threshold=0.05, database="org.Hs.eg.db", direction, sample = samples){
  
  
  #Load topTables
  #load('../results/MOLM13H_0.01P_0.05FDR_impF_bgRC_nomatch/topTables_MOLM13H_greater.Rda')
  #load('../results/MOLM13H_0.01P_0.05FDR_impF_bgSC_nomatch/topTables_MOLM13H_greater.Rda')
  
  load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_less.Rda"))
  load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/topTables_",params$cell_line,"_greater.Rda"))
  
  if (direction == "greater"){
    topTables <- topTables_greater
    } else {
      topTables <- topTables_less
      }
  
  pp <- lapply(1:length(samples), function(t){
  sample <- as.character(samples[t])
  #DEResults_less<-topTables_less$MOLM13H_RD
  DEResults_greater<-topTables[[sample]]

  #Define background
  ttl <- DEResults_greater#[DEResults_greater$topTable < 0,]
  proteins_test_greater <- unique(gsub("(\\w+)(-\\d+)?(_[STY]\\d+)", "\\1" ,rownames(ttl)))
  
  #use mart to convert uniprot to Ensembl ID
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "uswest.ensembl.org",ensemblRedirect = FALSE)
  #prots <- unique(HPRD_df$substrate_gene_symbol)
  rr2 <- getBM(
    attributes= c("ensembl_gene_id","uniprotswissprot"),
    filters = 'uniprotswissprot',
    values= proteins_test_greater,
    mart= mart)
  
  
  rr3 <- rr2[!rr2$uniprotswissprot=="",]
  gene.universe <- unique(rr3$ensembl_gene_id)
  #write.csv(gene.universe,'../gene.universe.csv')
  # Define foreground
  if (direction == "less"){
    ttl2 <- ttl[ttl$topTable<0 & ttl$padj <0.05,]
  } else {
    ttl2 <- ttl[ttl$topTable>0 & ttl$padj <0.05,]
  }
  
  proteins_test_greater <- unique(gsub("(\\w+)(-\\d+)?(_[STY]\\d+)", "\\1" ,rownames(ttl2)))
  
  proteins_test_greater <- rr3$ensembl_gene_id[match(proteins_test_greater,rr2$uniprotswissprot)]
  print(proteins_test_greater)
  DE_genes <- proteins_test_greater[!is.na(proteins_test_greater)]
  geneList <- factor(as.integer(gene.universe %in% DE_genes))
  names(geneList) <- gene.universe
  #print(geneList)
  # create topGOdata object
  onts<-c("BP","MF")
  tab<-as.list(onts)
  names(tab)<-onts
  
  ann.genes<-list(MF=list(), BP=list())
  names(ann.genes)<-onts
  
  for(i in 1:2){
  GOdata <- new("topGOdata", 
                ontology = onts[i], 
                allGenes = geneList,
                annot = annFUN.org,
                mapping=database,
                ID="ensembl")
  
 
  #perform enrichment
  resultTopGo.weight01.Fisher<-runTest(GOdata, algorithm="weight01", statistic="Fisher")
  resultTopGo.classic.Fisher<-runTest(GOdata, algorithm="classic", statistic="Fisher")
    
  #create table
  tab[[i]]<-GenTable(GOdata, 
                     weight01.Fisher=resultTopGo.weight01.Fisher, 
                     classic.Fisher=resultTopGo.classic.Fisher, 
                     orderBy="weight01.Fisher", 
                     ranksOf="classic.Fisher", 
                     topNodes=length(resultTopGo.weight01.Fisher@score))
  tab[[i]]$ontology<-onts[i]
  }

  #tt <- pp[pp$weight01.Fisher<0.05,]
  #View(tt)
  #return(pp)
  tab2 <- rbind_all(tab)
  tab2
  })
  }
  
all_samples_enr <- lapply(1:length(samples), function(x){
  enr <- RunTopGO(fdr.threshold=0.05, database="org.Hs.eg.db", direction="less", sample = samples[x])
  names(enr) <- samples[x]
  all_less <- rbind_all(enr)
  all_less <- all_less[all_less$classic.Fisher<0.05,]
  #write.csv(all_less, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/TopGO_less_cutoff0.05_background_selected_", params$cell_line,".csv"))
  
  enr2 <- RunTopGO(fdr.threshold=0.05, database="org.Hs.eg.db", direction="greater", sample = samples[x])
  names(enr2) <- samples[x]
  all_greater <- rbind_all(enr2)
  all_greater <- all_greater[all_greater$classic.Fisher<0.05,]
  #write.csv(all_greater, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/TopGO_greater_cutoff0.05_background_selected_", params$cell_line,".csv"))
  
  less_greater <- list(all_less, all_greater)
  names(less_greater) <- c("less","greater")
  less_greater
  
})
#pvalue_greater <- rbind_all(enr)
#pp <- all_samples_enr
#topGO_results_greater_all_genes <- enr
#topGO_results_greater
#save(topGO_results_less,file='./topGO_results_less.Rda')
names(all_samples_enr) <- samples
enriched_pathways <- all_samples_enr

save(enriched_pathways,  file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Enriched_pathways_", params$cell_line,".Rda"))
