##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Convert protein Ensembl to uniprot ID
## Date: 07.11.2018
## Author: Mahmoud Hallal
##################################################
library(SetRank)
## Load the predictions table from NetworKIN
db <- read.delim('../old_results/networkin_human_predictions_3.1.tsv')
db2 <- db[,c("X.substrate","position","id","sequence","networkin_score")]

db2$sub_gene <- gsub('(\\w+)(\\s)(\\(ENSP\\d+\\))','\\1',db2$X.substrate)
db2$sub_ensp <- gsub('(\\w+)(\\s)\\((ENSP\\d+)\\)','\\3',db2$X.substrate)

## ID converter method
IDConverter = createIDConverter("Homo.sapiens", "SYMBOL", "UNIPROT")
IDConverter2 = createIDConverter("Homo.sapiens", "ENSEMBLPROT", "UNIPROT")

db2$sub_uniprot1 <- IDConverter(as.character(db2$sub_gene))
db2$sub_uniprot1 <- IDConverter2(as.character(db2$sub_ensp))

## Biomart method

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
prots <- unique(db2$sub_ensp)
rr <- getBM(
  attributes= c("ensembl_peptide_id","uniprotswissprot"),
  filters = 'ensembl_peptide_id',
  values= prots,
  mart= mart)
db2$sub_uninprot <- rr$uniprotswissprot[match(db2$sub_ensp, rr$ensembl_peptide_id)]

# ll <- sapply(1:nrow(db2), function(x)
#   getBM(attributes= c("ensembl_peptide_id","uniprotswissprot"),
#         filters = 'ensembl_peptide_id', 
#         values= as.character(db2$sub_ensp[x]),
#         mart= mart))


## Get peptides
ll <- sapply(1:nrow(db2), function(x)
  toupper(strsplit(as.character(db2$sequence[x]),'')[[1]][grep('[sty]', strsplit(as.character(db2$sequence[x]),'')[[1]])])
)

db2$aa <- ll
db2$merge <- paste0(db2$sub_uninprot,'_',db2$aa,db2$position)
db2 <- db2[!is.na(db2$sub_uninprot),]
db2 <- db2[!db2$sub_uninprot=="",]

unique_peps <- unique(db2$merge)

##Take the highest score
pp <- lapply(1:length(unique_peps), function(x){
  print(x)
  dd <- db2[db2$merge == unique_peps[x],]
  rr <- dd[which(dd$networkin_score == max(dd$networkin_score)),]
  print(rr)
  rr
})

library(dplyr)
cc <- rbind_all(pp)

write.csv(cc, '~/Desktop/PhD/Generalized_pipeline/results/Updated_NetworKIN_predictions.csv',row.names = F)
