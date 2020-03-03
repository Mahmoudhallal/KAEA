##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Build database of kinase;substrate relationship from multiple databases
## Date: 09.03.17
## Author: Mahmoud Hallal
##################################################

## Load libraries##
#source("https://bioconductor.org/biocLite.R")

#biocLite("Rsamtools")
library(Rsamtools)
#biocLite("AnnotationDbi")
library(AnnotationDbi)

#biocLite("Biobase")
library(Biobase)

#biocLite("GenomicFeatures")
library(GenomicFeatures)
#biocLite("BSgenome")
library('BSgenome')
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
library('BSgenome.Hsapiens.UCSC.hg19')
#biocLite("Homo.sapiens")
library(Homo.sapiens)
#detach("package:dplyr", unload=TRUE)

library(yaml)
library(biomaRt)
#install.packages("igraph", repos="http://cran.rstudio.com/")

library(igraph)
install.packages("SetRank", repos="http://cran.rstudio.com/")
library(SetRank)#

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Load kinase-substrate databases
if (params$Species == "Human"){
  phospho_database <- read.delim(paste0(params$CWD,"/Phospho_DBs/Human_dbs.csv"), sep=',')
} else if (params$Species == "Mouse"){
  phospho_database <- read.delim(paste0(params$CWD,"/Phospho_DBs/Mouse_dbs.csv"), sep=',')
}

column_names <- c("geneID","termID","termName","dbName","description")


####################################################################################
## NetworKIN#
nwkin_input <- params$NWKIN_Input
nwkin <- read.delim(nwkin_input,sep=",")
nwkin$residue <- substr(nwkin$sequence,6,6)
nwkin <- nwkin[!nwkin$residue == "?",]
nwkin$residue <- toupper(nwkin$residue)
nwkin$protein_with_pos <- paste0(nwkin[,1],'_', nwkin$residue, nwkin$position)
nwkin <- nwkin[!duplicated(nwkin),]


nwkin <- nwkin[,c('protein_with_pos','id')]

motif_df <- data.frame(unique(nwkin$id),paste("MTF_",seq(1:length(unique(nwkin$id))), sep=""))
nwkin$substrateID <- motif_df$paste..MTF_...seq.1.length.unique.nwkin.id.....sep......[match(nwkin$id, motif_df$unique.nwkin.id.)]
nwkin$dbName <- "MTF"

# Rearrange columns#
motif_new <- nwkin[,c(1,3,2,4)]
motif_new$domain <- 'domain'
colnames(motif_new) <- column_names
motif_new <- motif_new[complete.cases(motif_new),]
motif_new <- motif_new[!duplicated(motif_new),]

#MusiteDeep#
MD_dbs <- read.delim('/Users/Mahmoud.Hallal/Desktop/PhD/MusiteDeep/DBs_MusiteDeep_17.07.2019/Molm13H_allsample_MusiteDeep_DBs_above10_17.07.2019_95cutoff.csv',sep=',')
MD_dbs1 <- MD_dbs
subs_MD_df <- data.frame(unique(MD_dbs1$Kinase_only),paste("MD_",seq(1:length(unique(MD_dbs1$Kinase_only))), sep=""))
MD_dbs1$substrateID <- subs_MD_df$paste..MD_...seq.1.length.unique.MD_dbs1.Kinase_only.....sep......[match(MD_dbs1$Kinase_only, subs_MD_df$unique.MD_dbs1.Kinase_only.)]
MD_dbs1$dbName <- 'MusiteDeep'
MD_dbs1$description <- 'domain'

#Rearrange columns
MD_dbs2 <- MD_dbs1[,c(8,10,9,11,12)]
colnames(MD_dbs2) <- column_names

####################################################################################
#Merge the all databases
all_dbs <- rbind(phospho_database, motif_new)
all_dbs$termName <- as.character(all_dbs$termName)
    
## remove PTPNs
#all_dbs <- all_dbs[-grep('PTPN',all_dbs$termName),]
all_dbs$termName[grep('^PTPN',all_dbs$termName)] <- ""

## Unify names given by different databases
unique(all_dbs$termName[grep('^CK',all_dbs$termName)])
unique(all_dbs$termName[grep('^CSNK',all_dbs$termName)])

all_dbs$termName[grep('^CK2_group$',all_dbs$termName)] <- "CSNK2"

all_dbs$termName[grep('^CK2_?alpha$',all_dbs$termName)] <- "CSNK2A"
all_dbs$termName[grep('^CK2a1$',all_dbs$termName)] <- "CSNK2A1"
all_dbs$termName[grep('^CK2a2$',all_dbs$termName)] <- "CSNK2A2"
all_dbs$termName[grep('^CK2a$',all_dbs$termName)] <- "CSNK2A"

all_dbs$termName[grep('^CSNK2A1$',all_dbs$termName)] <- "CSNK2A1"
all_dbs$termName[grep('^CSNK2A2$',all_dbs$termName)] <- "CSNK2A2"

all_dbs$termName[grep('^CK2B$',all_dbs$termName)] <- "CSNK2B"
all_dbs$termName[grep('^CK2_beta$',all_dbs$termName)]<- "CSNK2B"

#CK1
unique(all_dbs$termName[grep('^CSNK$',all_dbs$termName)])
unique(all_dbs$termName[grep('^CK1$',all_dbs$termName)])

all_dbs$termName[grep('^CK1alpha$',all_dbs$termName)] <- "CSNK1A"
all_dbs$termName[grep('^CK1a$',all_dbs$termName)] <- "CSNK1A"
all_dbs$termName[grep('^CSNK1A1$',all_dbs$termName)] <- "CSNK1A1"
all_dbs$termName[grep('^CK1_delta$',all_dbs$termName)] <- "CSNK1D"
all_dbs$termName[grep('^CK1d$',all_dbs$termName)] <- "CSNK1D"
all_dbs$termName[grep('^CK1delta$',all_dbs$termName)] <- "CSNK1D"
all_dbs$termName[grep('^CK1e$',all_dbs$termName)] <- "CSNK1E"
all_dbs$termName[grep('^CK1_epsilon$',all_dbs$termName)] <- "CSNK1E"
all_dbs$termName[grep('^CK1gamma2$',all_dbs$termName)] <- "CSNK1G2"

##DOUBLE CHECK THIS #06.11.2018
all_dbs$termName[grep('^CK$',all_dbs$termName)] <- "CSNK2"

#MRKCb
unique(all_dbs$termName[grep('^CDC42',all_dbs$termName)])
all_dbs$termName[grep('^MRCKb$',all_dbs$termName)] <- "CDC42BPB"

#SRC
unique(all_dbs$termName[grep('^SRC',all_dbs$termName)])
all_dbs$termName[grep('^SRC_group$',all_dbs$termName)] <- "SRC"
all_dbs$termName[grep('^SRC-type Tyr-kinases$',all_dbs$termName)] <- "SRC"


#GSK3
unique(all_dbs$termName[grep('^GSK',all_dbs$termName)])

all_dbs$termName[grep('^GSK-3_alpha$',all_dbs$termName)] <- "GSK3A"
all_dbs$termName[grep('^GSK3alpha$',all_dbs$termName)] <- "GSK3A"
#
all_dbs$termName[grep('^GSK-3_beta$',all_dbs$termName)] <- "GSK3B"
all_dbs$termName[grep('^GSK3beta$',all_dbs$termName)] <- "GSK3B"
#
all_dbs$termName[grep('^GSK-3_group$',all_dbs$termName)] <- "GSK3"
all_dbs$termName[grep('^GSK3_group$',all_dbs$termName)] <- "GSK3"

#Aurora
#A
unique(all_dbs$termName[grep('^Aur',all_dbs$termName)])
unique(all_dbs$termName[grep('^AUR',all_dbs$termName)])

all_dbs$termName[grep('Aurora A$',all_dbs$termName)] <- "AURKA"
all_dbs$termName[grep('AurA$',all_dbs$termName)] <- "AURKA"
all_dbs$termName[grep('AuroraA$',all_dbs$termName)] <- "AURKA"

#B
all_dbs$termName[grep('Aurora B$',all_dbs$termName)] <- "AURKB"
all_dbs$termName[grep('AurB$',all_dbs$termName)] <- "AURKB"

#PKB
unique(all_dbs$termName[grep('^PKB',all_dbs$termName)])

all_dbs$termName[grep('PKB_group$',all_dbs$termName)] <- "AKT"
all_dbs$termName[grep('PKBalpha$',all_dbs$termName)] <- "AKT1"
all_dbs$termName[grep('PKB_beta$',all_dbs$termName)] <- "AKT2"
all_dbs$termName[grep('PKBbeta$',all_dbs$termName)] <- "AKT2"

#PKC
unique(all_dbs$termName[grep('PKC',all_dbs$termName)])
all_dbs$termName[grep('PKCalpha$',all_dbs$termName)] <- "PRKCA"
all_dbs$termName[grep('PKCa$',all_dbs$termName)] <- "PRKCA"
all_dbs$termName[grep('PKC_alpha$',all_dbs$termName)] <- "PRKCA"

all_dbs$termName[grep('PKCbeta$',all_dbs$termName)] <- "PRKCB"
all_dbs$termName[grep('PKCb$',all_dbs$termName)] <- "PRKCB"
all_dbs$termName[grep('PKC_beta$',all_dbs$termName)] <- "PRKCB"

all_dbs$termName[grep('PKCdelta$',all_dbs$termName)] <- "PRKCD"
all_dbs$termName[grep('PKCd$',all_dbs$termName)] <- "PRKCD"
all_dbs$termName[grep('PKC_delta$',all_dbs$termName)] <- "PRKCD"

all_dbs$termName[grep('PKCepsilon$',all_dbs$termName)] <- "PRKCE"
all_dbs$termName[grep('PKCe$',all_dbs$termName)] <- "PRKCE"
all_dbs$termName[grep('PKC_epsilon$',all_dbs$termName)] <- "PRKCE"

all_dbs$termName[grep('PKCtheta$',all_dbs$termName)] <- "PRKCT"
all_dbs$termName[grep('PKCt',all_dbs$termName)] <- "PRKCT"
all_dbs$termName[grep('PKC_theta',all_dbs$termName)] <- "PRKCT"

all_dbs$termName[grep('PKCgamma',all_dbs$termName)] <- "PRKCG"
all_dbs$termName[grep('PKCg$',all_dbs$termName)] <- "PRKCG"
all_dbs$termName[grep('PKC_gamma$',all_dbs$termName)] <- "PRKCG"

all_dbs$termName[grep('PKCz$',all_dbs$termName)] <- "PRKCZ"
all_dbs$termName[grep('PKC_zeta$',all_dbs$termName)] <- "PRKCZ"

all_dbs$termName[grep('PKCiota$',all_dbs$termName)] <- "PRKCI"
all_dbs$termName[grep('PKCi$',all_dbs$termName)] <- "PRKCI"
all_dbs$termName[grep('PKC_iota$',all_dbs$termName)] <- "PRKCI"

all_dbs$termName[grep('PKC_eta$',all_dbs$termName)] <- "PRKCH"
all_dbs$termName[grep('PKCh$',all_dbs$termName)] <- "PRKCH"

all_dbs$termName[grep('PKC_group$',all_dbs$termName)] <- "PRKC"

##IKK
unique(all_dbs$termName[grep('IKK',all_dbs$termName)])

all_dbs$termName[grep('IKKa$',all_dbs$termName)] <- "IKKA"
all_dbs$termName[grep('CHUK$',all_dbs$termName)] <- "IKKA"
all_dbs$termName[grep('IKK_alpha$',all_dbs$termName)] <- "IKKA"

all_dbs$termName[grep('IKKb$',all_dbs$termName)] <- "IKKB"
all_dbs$termName[grep('IKBKB$',all_dbs$termName)] <- "IKKB"
all_dbs$termName[grep('IKK_beta$',all_dbs$termName)] <- "IKKB"

all_dbs$termName[grep('IKKe$',all_dbs$termName)] <- "IKKE"
all_dbs$termName[grep('IKBKE$',all_dbs$termName)] <- "IKKE"
all_dbs$termName[grep('IKK_epsilon$',all_dbs$termName)] <- "IKKE"

all_dbs$termName[grep('IKK-complex$',all_dbs$termName)] <- "IKK"
all_dbs$termName[grep('IKK_group$',all_dbs$termName)] <- "IKK"


##BCR/ABL
unique(all_dbs[grep('Abl',all_dbs$termName),])
all_dbs$termName[grep('BCR/ABL$',all_dbs$termName)] <- "BCR/ABL"
all_dbs$termName[grep('BCR-ABL$',all_dbs$termName)] <- "BCR/ABL"

all_dbs$termName[grep('^ABL$',all_dbs$termName)] <- "ABL"
all_dbs$termName[grep('^Abl$',all_dbs$termName)] <- "ABL"

#PKA
all_dbs$termName[grep('^PKAalpha$',all_dbs$termName)] <- "PRKACA"
all_dbs$termName[grep('^PKACa$',all_dbs$termName)] <- "PRKACA"
all_dbs$termName[grep('^PKA_alpha$',all_dbs$termName)] <- "PRKACA"

#CAMK2
unique(all_dbs$termName[grep('CaM-K$',all_dbs$termName)])
unique(all_dbs$termName[grep('CaMK$',all_dbs$termName)])

all_dbs$termName[grep('^CaMKIIdelta$',all_dbs$termName)] <- "CaMK2D"

all_dbs$termName[grep('^CaMKIIbeta$',all_dbs$termName)] <- "CaMK2B"
all_dbs$termName[grep('^CaMK2b$',all_dbs$termName)] <- "CaMK2B"

all_dbs$termName[grep('^CaMKIIalpha$',all_dbs$termName)] <- "CAMK2A"
all_dbs$termName[grep('^CaMK2a$',all_dbs$termName)] <- "CAMK2A"
all_dbs$termName[grep('^CaMKIIa$',all_dbs$termName)] <- "CAMK2A"

all_dbs$termName[grep('^CaM-KII_alpha$',all_dbs$termName)] <- "CAMK2A"
all_dbs$termName[grep('^CaM-KII_group$',all_dbs$termName)] <- "CAMK2"
all_dbs$termName[grep('^CaM-CaMK2_group$',all_dbs$termName)] <- "CAMK2"

all_dbs$termName[grep('^CaM-KI_alpha$',all_dbs$termName)] <- "CaMK1a"
all_dbs$termName[grep('^CaMK1a$',all_dbs$termName)] <- "CAMK1A"

all_dbs$termName[grep('^CaM-KIV',all_dbs$termName)] <- "CAMK4"
all_dbs$termName[grep('^CaMK4',all_dbs$termName)] <- "CAMK4"

all_dbs$termName[grep('^CaM-KK_alpha$',all_dbs$termName)] <- "CAMKK1"
all_dbs$termName[grep('^CaMKK1$',all_dbs$termName)] <- "CAMKK1"

#PDK (PDK-1 and PDK-2 belong to ELM and they represent PDPK1 and PDPK2)
unique(all_dbs$termName[grep('PDK',all_dbs$termName)])
all_dbs$termName[grep('^PDK-1$',all_dbs$termName)] <- "PDPK1"
all_dbs$termName[grep('^PDK-2$',all_dbs$termName)] <- "PDPK2"

#RSK
unique(all_dbs$termName[grep('RSK',all_dbs$termName)])
all_dbs$termName[grep('^RSK-1$',all_dbs$termName)] <- "RSK1"
all_dbs$termName[grep('^RSK-2$',all_dbs$termName)] <- "RSK2"
all_dbs$termName[grep('^RSK-3$',all_dbs$termName)] <- "RSK3"
all_dbs$termName[grep('^RSK-5$',all_dbs$termName)] <- "RSK5"

#GRK
unique(all_dbs$termName[grep('PIM',all_dbs$termName)])
all_dbs$termName[grep('^GRK-1$',all_dbs$termName)] <- "GRK1"
all_dbs$termName[grep('^GRK-2$',all_dbs$termName)] <- "GRK2"
all_dbs$termName[grep('^GRK-3$',all_dbs$termName)] <- "GRK3"
all_dbs$termName[grep('^GRK-4$',all_dbs$termName)] <- "GRK4"
all_dbs$termName[grep('^GRK-5$',all_dbs$termName)] <- "GRK5"
all_dbs$termName[grep('^GRK-6$',all_dbs$termName)] <- "GRK6"

#PIM
unique(all_dbs$termName[grep('PIM',all_dbs$termName)])
all_dbs$termName[grep('^PIM-1$',all_dbs$termName)] <- "PIM1"

#PDGFR
unique(all_dbs$termName[grep('PDGFR',all_dbs$termName)])
all_dbs$termName[grep('^PDGFR_alpha$',all_dbs$termName)] <- "PDGFRA"
all_dbs$termName[grep('^PDGFR_beta$',all_dbs$termName)] <- "PDGFRB"
all_dbs$termName[grep('^PDGFRa$',all_dbs$termName)] <- "PDGFRA"
all_dbs$termName[grep('^PDGFRb$',all_dbs$termName)] <- "PDGFRB"
all_dbs$termName[grep('^PDGFR_group$',all_dbs$termName)] <- "PDGFR"


#FRAP and mTOR
unique(all_dbs$termName[grep('RAFT',all_dbs$termName)])
all_dbs$termName[grep('^FRAP$',all_dbs$termName)] <- "MTOR"
all_dbs$termName[grep('^mTOR$',all_dbs$termName)] <- "MTOR"

#FYN
unique(all_dbs$termName[grep('Fyn',all_dbs$termName)])
all_dbs$termName[grep('^Fyn$',all_dbs$termName)] <- "FYN"

#CDK
unique(all_dbs$termName[grep('Cyclin',all_dbs$termName)])
unique(all_dbs$termName[grep('CDK',all_dbs$termName)])

all_dbs$termName[grep('^CyclinD/CDK4$',all_dbs$termName)] <- "CDK4"
all_dbs$termName[grep('^CyclinE/CDK2$',all_dbs$termName)] <- "CDK2"
all_dbs$termName[grep('^CyclinB/CDK1$',all_dbs$termName)] <- "CDK1"
all_dbs$termName[grep('^CDC2$',all_dbs$termName)] <- "CDK1"
all_dbs$termName[grep('^CDC42BPA$',all_dbs$termName)] <- "MRCKa"


all_dbs$termName[grep('^CyclinA2/CDK2$',all_dbs$termName)] <- "CDK2"

#ADRBK1, BARK, GRK2
all_dbs$termName[grep('^BARK1$',all_dbs$termName)] <- "ADRBK1"
all_dbs$termName[grep('^GRK2$',all_dbs$termName)] <- "ADRBK1"

#WNK4
#all_dbs$termName[grep('^Wnk4pp$',all_dbs)] <- "WNK4"

#DNA-PK
all_dbs$termName[grep('^DNA-PK$',all_dbs$termName)] <- "DNAPK"


#Remove autocatalysis
all_dbs <- all_dbs[!all_dbs$termName=="autocatalysis",]
all_dbs <- all_dbs[!all_dbs$termName=="inisoform HMG-I",]

# Remove phosphatases
#all_dbs[grep('^PTPRE$',all_dbs$termName),]$termName <- ""
all_dbs$termName[grep('^PTPR',all_dbs$termName)] <- ""
all_dbs$termName[grep('^CDC14B$',all_dbs$termName)] <- ""
all_dbs$termName[grep('^CDC14A$',all_dbs$termName)] <- ""

all_dbs$termName[grep('^PPP1CA$',all_dbs$termName)] <- ""
all_dbs$termName[grep('^PPP2CA$',all_dbs$termName)] <- ""
all_dbs$termName[grep('^PPP3CA$',all_dbs$termName)] <- ""
all_dbs$termName[grep('^PPP2CB$',all_dbs$termName)] <- ""
all_dbs$termName[grep('^PPP3CC$',all_dbs$termName)] <- ""
all_dbs$termName[grep('^PPM1D$',all_dbs$termName)] <- ""


all_dbs$termName[grep('^DUSP4$',all_dbs$termName)] <- ""
all_dbs <- all_dbs[!all_dbs$termName=="",]

## Remove kinases with one entry
table_count <- table(all_dbs$termName)
table_count_1 <- table_count[table_count == 1]
all_dbs <- all_dbs[!all_dbs$termName %in% names(table_count_1),]

## Lyn
all_dbs$termName[grep('^Lyn$',all_dbs$termName)] <- "LYN"
all_dbs$termName[grep('^Lck$',all_dbs$termName)] <- "LCK"
all_dbs$termName[grep('^Met$',all_dbs$termName)] <- "MET"
all_dbs$termName[grep('^TGFbR2$',all_dbs$termName)] <- "TGFBR2"
all_dbs$termName[grep('^PKD1$',all_dbs$termName)] <- "PRKD1"

# Make groups look same
groups <- grep('[A-Za-z0-9]group',all_dbs$termName)
#above100
all_dbs$termName[groups] <- gsub('(.*)(group)','\\1_\\2',all_dbs$termName[groups])


## Unified database#
all_dbs$dbName <- "DBS"
all_dbs <- subset(all_dbs, select = -termID)
all_dbs <- all_dbs[!duplicated(all_dbs),]
dd <- data.frame(unique(all_dbs$termName),paste("DBS_",seq(1:length(unique(all_dbs$termName))), sep=""))
all_dbs$termID <- dd$paste..DBS_...seq.1.length.unique.all_dbs.termName.....sep......[match(all_dbs$termName, dd$unique.all_dbs.termName.)]
all_dbs <- all_dbs[, c(1,5,2,3,4)]
all_dbs <- all_dbs[!duplicated(all_dbs),]
all_dbs$termName <- as.factor(all_dbs$termName)
#write.csv(all_dbs,'../../MusiteDeep/all_dbs_02.07.2019.csv' , row.names = FALSE)
## Write output DB#
#write.table(all_dbs, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/all_dbs.csv"), row.names = FALSE, sep=',')
write.csv(all_dbs, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/all_dbs.csv"), row.names = FALSE)
#write.csv(all_dbs,'~/Desktop/PhD/MusiteDeep/all_dbs_02.07.2019_corrected_05.08.2019.csv',row.names=FALSE)
