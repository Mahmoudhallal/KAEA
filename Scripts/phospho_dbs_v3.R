##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Build database of kinase;substrate relationship from multiple databases
## Date: 09.03.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
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
library(SetRank)

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Load kinase-substrate databases
if (params$Species == "Human"){
  phospho_database <- read.delim(paste0(params$CWD,"/Phospho_DBs/Human_dbs.csv"), sep=',')
} else if (params$Species == "Mouse"){
  phospho_database <- read.delim(paste0(params$CWD,"/Phospho_DBs/Mouse_dbs.csv"), sep=',')
}

## column names
column_names <- c("geneID","termID","termName","dbName","description")

####################################################################################
## NetworKIN
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

## Rearrange columns
motif_new <- nwkin[,c(1,3,2,4)]
motif_new$domain <- 'domain'
colnames(motif_new) <- column_names
motif_new <- motif_new[complete.cases(motif_new),]
motif_new <- motif_new[!duplicated(motif_new),]

#####################################################################################
## Merge all databases
all_dbs <- rbind(phospho_database, motif_new)
all_dbs$termName <- as.character(all_dbs$termName)
    
## remove PTPNs
all_dbs$termName[grep('^PTPN',all_dbs$termName)] <- ""

## Unify names given by different databases
unique(all_dbs$termName[grep('^CK',all_dbs$termName)])
unique(all_dbs$termName[grep('^CSNK',all_dbs$termName)])

all_dbs$termName[grep('^CK2_group$',all_dbs$termName)] <- "CSNK2"

all_dbs$termName[grep('^CK2_?alpha$',all_dbs$termName)] <- "CSNK2A"
all_dbs$termName[grep('^CK2a1$',all_dbs$termName)] <- "CSNK2A1"
all_dbs$termName[grep('^CK2a2$',all_dbs$termName)] <- "CSNK2A2"
all_dbs$termName[grep('^CK2a$',all_dbs$termName)] <- "CSNK2A"

all_dbs$termName[grep('^CK2B$',all_dbs$termName)] <- "CSNK2B"
all_dbs$termName[grep('^CK2_beta$',all_dbs$termName)]<- "CSNK2B"

#CK1
unique(all_dbs$termName[grep('^CSNK$',all_dbs$termName)])
unique(all_dbs$termName[grep('^CK1$',all_dbs$termName)])

all_dbs$termName[grep('^CK1_group$',all_dbs$termName)] <- "CSNK1"
all_dbs$termName[grep('^CK1alpha$',all_dbs$termName)] <- "CSNK1A"
all_dbs$termName[grep('^CK1_alpha$',all_dbs$termName)] <- "CSNK1A"
all_dbs$termName[grep('^CK1a$',all_dbs$termName)] <- "CSNK1A"
all_dbs$termName[grep('^CK1_delta$',all_dbs$termName)] <- "CSNK1D"
all_dbs$termName[grep('^CK1d$',all_dbs$termName)] <- "CSNK1D"
all_dbs$termName[grep('^CK1delta$',all_dbs$termName)] <- "CSNK1D"
all_dbs$termName[grep('^CK1e$',all_dbs$termName)] <- "CSNK1E"
all_dbs$termName[grep('^CK1_epsilon$',all_dbs$termName)] <- "CSNK1E"
all_dbs$termName[grep('^CK1gamma2$',all_dbs$termName)] <- "CSNK1G2"
all_dbs$termName[grep('^CK1gamma3$',all_dbs$termName)] <- "CSNK1G3"
all_dbs$termName[grep('^CK1epsilon$',all_dbs$termName)] <- "CSNK1E"

##
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

all_dbs$termName[grep('^Aurora A$',all_dbs$termName)] <- "AURKA"
all_dbs$termName[grep('^AurA$',all_dbs$termName)] <- "AURKA"
all_dbs$termName[grep('^AuroraA$',all_dbs$termName)] <- "AURKA"

#B
all_dbs$termName[grep('^Aurora B$',all_dbs$termName)] <- "AURKB"
all_dbs$termName[grep('^AurB$',all_dbs$termName)] <- "AURKB"

#C
all_dbs$termName[grep('^AurC$',all_dbs$termName)] <- "AURKC"



#PKB
unique(all_dbs$termName[grep('^PKB',all_dbs$termName)])

all_dbs$termName[grep('^PKB_group$',all_dbs$termName)] <- "AKT"
all_dbs$termName[grep('^PKBalpha$',all_dbs$termName)] <- "AKT1"
all_dbs$termName[grep('^PKB_beta$',all_dbs$termName)] <- "AKT2"
all_dbs$termName[grep('^PKBbeta$',all_dbs$termName)] <- "AKT2"
all_dbs$termName[grep('^PKBgamma$',all_dbs$termName)] <- "AKT3"

#PKC
unique(all_dbs$termName[grep('PKC',all_dbs$termName)])
all_dbs$termName[grep('^PKCalpha$',all_dbs$termName)] <- "PRKCA"
all_dbs$termName[grep('^PKCa$',all_dbs$termName)] <- "PRKCA"
all_dbs$termName[grep('^PKC_alpha$',all_dbs$termName)] <- "PRKCA"

all_dbs$termName[grep('^PKCbeta$',all_dbs$termName)] <- "PRKCB"
all_dbs$termName[grep('^PKCb$',all_dbs$termName)] <- "PRKCB"
all_dbs$termName[grep('^PKC_beta$',all_dbs$termName)] <- "PRKCB"

all_dbs$termName[grep('^PKCdelta$',all_dbs$termName)] <- "PRKCD"
all_dbs$termName[grep('^PKCd$',all_dbs$termName)] <- "PRKCD"
all_dbs$termName[grep('^PKC_delta$',all_dbs$termName)] <- "PRKCD"

all_dbs$termName[grep('^PKCepsilon$',all_dbs$termName)] <- "PRKCE"
all_dbs$termName[grep('^PKCe$',all_dbs$termName)] <- "PRKCE"
all_dbs$termName[grep('^PKC_epsilon$',all_dbs$termName)] <- "PRKCE"

all_dbs$termName[grep('^PKCtheta$',all_dbs$termName)] <- "PRKCQ"
all_dbs$termName[grep('^PKCt',all_dbs$termName)] <- "PRKCQ"
all_dbs$termName[grep('^PKC_theta',all_dbs$termName)] <- "PRKCQ"

all_dbs$termName[grep('^PKCgamma',all_dbs$termName)] <- "PRKCG"
all_dbs$termName[grep('^PKCg$',all_dbs$termName)] <- "PRKCG"
all_dbs$termName[grep('^PKC_gamma$',all_dbs$termName)] <- "PRKCG"

all_dbs$termName[grep('^PKCz$',all_dbs$termName)] <- "PRKCZ"
all_dbs$termName[grep('^PKC_zeta$',all_dbs$termName)] <- "PRKCZ"

all_dbs$termName[grep('^PKCiota$',all_dbs$termName)] <- "PRKCI"
all_dbs$termName[grep('^PKCi$',all_dbs$termName)] <- "PRKCI"
all_dbs$termName[grep('^PKC_iota$',all_dbs$termName)] <- "PRKCI"

all_dbs$termName[grep('^PKC_eta$',all_dbs$termName)] <- "PRKCH"
all_dbs$termName[grep('^PKCh$',all_dbs$termName)] <- "PRKCH"

all_dbs$termName[grep('^PKC_group$',all_dbs$termName)] <- "PRKC"

##IKK
unique(all_dbs$termName[grep('IKK',all_dbs$termName)])

all_dbs$termName[grep('IKKa$',all_dbs$termName)] <- "CHUK"
all_dbs$termName[grep('IKKA$',all_dbs$termName)] <- "CHUK"
all_dbs$termName[grep('IKK_alpha$',all_dbs$termName)] <- "CHUK"
all_dbs$termName[grep('IKKalpha$',all_dbs$termName)] <- "CHUK"

all_dbs$termName[grep('IKKb$',all_dbs$termName)] <- "IKBKB"
all_dbs$termName[grep('IKKB$',all_dbs$termName)] <- "IKBKB"
all_dbs$termName[grep('IKK_beta$',all_dbs$termName)] <- "IKBKB"
all_dbs$termName[grep('IKKbeta$',all_dbs$termName)] <- "IKBKB"

all_dbs$termName[grep('IKKe$',all_dbs$termName)] <- "IKBKE"
all_dbs$termName[grep('IKKE$',all_dbs$termName)] <- "IKBKE"
all_dbs$termName[grep('IKK_epsilon$',all_dbs$termName)] <- "IKBKE"

all_dbs$termName[grep('IKK-complex$',all_dbs$termName)] <- "IKK"
all_dbs$termName[grep('IKK_group$',all_dbs$termName)] <- "IKK"


##BCR/ABL
unique(all_dbs[grep('Abl',all_dbs$termName),])
all_dbs$termName[grep('^BCR-ABL$',all_dbs$termName)] <- "BCR/ABL"

all_dbs$termName[grep('^ABL$',all_dbs$termName)] <- "ABL"
all_dbs$termName[grep('^Abl$',all_dbs$termName)] <- "ABL1"

#DNA-PK
all_dbs$termName[grep('^DNA-PK$',all_dbs$termName)] <- "DNAPK"


#CAMK2
unique(all_dbs$termName[grep('CaM-K$',all_dbs$termName)])
unique(all_dbs$termName[grep('CaMK$',all_dbs$termName)])

all_dbs$termName[grep('^CaMKIIdelta$',all_dbs$termName)] <- "CAMK2D"

all_dbs$termName[grep('^CaMKIIbeta$',all_dbs$termName)] <- "CAMK2B"
all_dbs$termName[grep('^CaMK2b$',all_dbs$termName)] <- "CAMK2B"

all_dbs$termName[grep('^CaMKIIalpha$',all_dbs$termName)] <- "CAMK2A"
all_dbs$termName[grep('^CaMK2a$',all_dbs$termName)] <- "CAMK2A"
all_dbs$termName[grep('^CaMKIIa$',all_dbs$termName)] <- "CAMK2A"

all_dbs$termName[grep('^CaM-KII_alpha$',all_dbs$termName)] <- "CAMK2A"
all_dbs$termName[grep('^CaM-KII_group$',all_dbs$termName)] <- "CAMK2"

all_dbs$termName[grep('^CaMKIV$',all_dbs$termName)] <- "CAMK4"

all_dbs$termName[grep('^CaMKIIgamma$',all_dbs$termName)] <- "CAMK2G"

#CAMK1
all_dbs$termName[grep('^CaM-KI_alpha$',all_dbs$termName)] <- "CaMK1a"
all_dbs$termName[grep('^CaMK1a$',all_dbs$termName)] <- "CAMK1A"

all_dbs$termName[grep('^CaM-KIV',all_dbs$termName)] <- "CAMK4"
all_dbs$termName[grep('^CaMK4',all_dbs$termName)] <- "CAMK4"

all_dbs$termName[grep('^CaM-KK_alpha$',all_dbs$termName)] <- "CAMKK1"
all_dbs$termName[grep('^CaMKK1$',all_dbs$termName)] <- "CAMKK1"

#PDK (PDK-1 and PDK-2 belong to ELM and they represent PDPK1 and PDPK2)
#unique(all_dbs$termName[grep('PDK',all_dbs$termName)])
#all_dbs$termName[grep('^PDK-1$',all_dbs$termName)] <- "PDPK1"
#all_dbs$termName[grep('^PDK-2$',all_dbs$termName)] <- "PDPK2"

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
all_dbs$termName[grep('^MRCKa$',all_dbs$termName)] <- "CDC42BPA"


all_dbs$termName[grep('^CyclinA2/CDK2$',all_dbs$termName)] <- "CDK2"

#ADRBK1, BARK, GRK2
all_dbs$termName[grep('^BARK1$',all_dbs$termName)] <- "GRK2"
all_dbs$termName[grep('^ADRBK1$',all_dbs$termName)] <- "GRK2"
all_dbs$termName[grep('^ADRBK2$',all_dbs$termName)] <- "GRK3"
all_dbs$termName[grep('^BARK2$',all_dbs$termName)] <- "GRK3"

#WNK4
#all_dbs$termName[grep('^Wnk4pp$',all_dbs)] <- "WNK4"

#DNA-PK
all_dbs$termName[grep('^DNA-PK$',all_dbs$termName)] <- "DNAPK"

#Lyn
all_dbs$termName[grep('^Lyn$',all_dbs$termName)] <- "LYN"
all_dbs$termName[grep('^Lck$',all_dbs$termName)] <- "LCK"
all_dbs$termName[grep('^Met$',all_dbs$termName)] <- "MET"
all_dbs$termName[grep('^TGFbR2$',all_dbs$termName)] <- "TGFBR2"
all_dbs$termName[grep('^Tyk2$',all_dbs$termName)] <- "TYK2"

all_dbs$termName[grep('^PKD1$',all_dbs$termName)] <- "PRKD1"
all_dbs$termName[grep('^PKD2$',all_dbs$termName)] <- "PRKD2"
all_dbs$termName[grep('^PKD3$',all_dbs$termName)] <- "PRKD3"

## Make groups look same
groups <- grep('[A-Za-z0-9]group',all_dbs$termName)
all_dbs$termName[groups] <- gsub('(.*)(group)','\\1_\\2',all_dbs$termName[groups])

#PKA
all_dbs$termName[grep('^PKAalpha$',all_dbs$termName)] <- "PRKACA"
all_dbs$termName[grep('^PKACa$',all_dbs$termName)] <- "PRKACA"
all_dbs$termName[grep('^PKA_alpha$',all_dbs$termName)] <- "PRKACA"
all_dbs$termName[grep('^PKAgamma$',all_dbs$termName)] <- "PRKACG"
all_dbs$termName[grep('^PKAbeta$',all_dbs$termName)] <- "PRKACB"

## Add p42229_Y694 of FLT3
#addition <- t(data.frame(c("P42229_Y694","OWN","FLT3","OWN2","domain")))
#colnames(addition) <- colnames(all_dbs)
#rownames(addition) <- NULL
#all_dbs <- rbind(all_dbs, addition)

#Remove autocatalysis
all_dbs <- all_dbs[!all_dbs$termName=="autocatalysis",]
all_dbs <- all_dbs[!all_dbs$termName=="inisoform HMG-I",]

# Remove phosphatases
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

#14.05.20
all_dbs$termName[grep('^Abl2$',all_dbs$termName)] <- "ABL2"
all_dbs$termName[grep('^eEF2K$',all_dbs$termName)] <- "EEF2K"
all_dbs$termName[grep('^Brk$',all_dbs$termName)] <- "PTK6"
all_dbs$termName[grep('^BRK$',all_dbs$termName)] <- "PTK6"
all_dbs$termName[grep('^Yes$',all_dbs$termName)] <- "YES1"
all_dbs$termName[grep('^YES$',all_dbs$termName)] <- "YES1"
all_dbs$termName[grep('^Csk$',all_dbs$termName)] <- "CSK"
all_dbs$termName[grep('^Fes$',all_dbs$termName)] <- "FES"
all_dbs$termName[grep('^Mer$',all_dbs$termName)] <- "MERTK"
all_dbs$termName[grep('^MER$',all_dbs$termName)] <- "MERTK"
all_dbs$termName[grep('^Axl$',all_dbs$termName)] <- "AXL"
all_dbs$termName[grep('^Fgr$',all_dbs$termName)] <- "FGR"
all_dbs$termName[grep('^Kit$',all_dbs$termName)] <- "KIT"
all_dbs$termName[grep('^Wnk1$',all_dbs$termName)] <- "WNK1"
all_dbs$termName[grep('^NuaK1$',all_dbs$termName)] <- "NUAK1"
all_dbs$termName[grep('^ErbB2$',all_dbs$termName)] <- "ERBB2"
all_dbs$termName[grep('^TGFbR1$',all_dbs$termName)] <- "TGFBR1"
all_dbs$termName[grep('^HRI$',all_dbs$termName)] <- "EIF2AK1"
all_dbs$termName[grep('^AMPKa1$',all_dbs$termName)] <- "PRKAA1"
all_dbs$termName[grep('^AMPKa2$',all_dbs$termName)] <- "PRKAA2"
all_dbs$termName[grep('^CHK1$',all_dbs$termName)] <- "CHEK1"
all_dbs$termName[grep('^OSR1$',all_dbs$termName)] <- "OXSR1"
all_dbs$termName[grep('^TRKB$',all_dbs$termName)] <- "NTRK2"
all_dbs$termName[grep('^MSK2$',all_dbs$termName)] <- "RPS6KA4"
all_dbs$termName[grep('^PYK2$',all_dbs$termName)] <- "PTK2B"
all_dbs$termName[grep('^P70S6KB$',all_dbs$termName)] <- "RPS6KB2"
all_dbs$termName[grep('^Brk$',all_dbs$termName)] <- "PTK6"
all_dbs$termName[grep('^BRK$',all_dbs$termName)] <- "PTK6"
all_dbs$termName[grep('^CHK2$',all_dbs$termName)] <- "CHEK2"
all_dbs$termName[grep('^DNAPK$',all_dbs$termName)] <- "PRKDC"
all_dbs$termName[grep('^EphA2$',all_dbs$termName)] <- "EPHA2"
all_dbs$termName[grep('^LKB1$',all_dbs$termName)] <- "STK11"
all_dbs$termName[grep('^RON$',all_dbs$termName)] <- "MST1R"
all_dbs$termName[grep('^KIS$',all_dbs$termName)] <- "UHMK1"
all_dbs$termName[grep('^PKG1$',all_dbs$termName)] <- "PRKG1"
all_dbs$termName[grep('^FAK$',all_dbs$termName)] <- "PTK2"
all_dbs$termName[grep('^PRP4$',all_dbs$termName)] <- "PRPF4B"
all_dbs$termName[grep('^MST2$',all_dbs$termName)] <- "STK3"
all_dbs$termName[grep('^MLK3$',all_dbs$termName)] <- "MAP3K11"
all_dbs$termName[grep('^NIK$',all_dbs$termName)] <- "MAP3K14"
all_dbs$termName[grep('^p70S6K$',all_dbs$termName)] <- "RPS6KB1"
all_dbs$termName[grep('^ChaK1$',all_dbs$termName)] <- "TRPM7"
all_dbs$termName[grep('^Chak1$',all_dbs$termName)] <- "TRPM7"
all_dbs$termName[grep('^Mnk1$',all_dbs$termName)] <- "MKNK1"
all_dbs$termName[grep('^MNK1$',all_dbs$termName)] <- "MKNK1"
all_dbs$termName[grep('^MST1$',all_dbs$termName)] <- "STK4"
all_dbs$termName[grep('^EphA4$',all_dbs$termName)] <- "EPHA4"
all_dbs$termName[grep('^smMLCK$',all_dbs$termName)] <- "MYLK"
all_dbs$termName[grep('^ACK$',all_dbs$termName)] <- "TNK2"
all_dbs$termName[grep('^MSK1$',all_dbs$termName)] <- "RPS6KA5"
all_dbs$termName[grep('^RSK2$',all_dbs$termName)] <- "RPS6KA3"
all_dbs$termName[grep('^EphA3$',all_dbs$termName)] <- "EPHA3"
all_dbs$termName[grep('^EphA8$',all_dbs$termName)] <- "EPHA8"
all_dbs$termName[grep('^EphB1$',all_dbs$termName)] <- "EPHB1"
all_dbs$termName[grep('^EphB2$',all_dbs$termName)] <- "EPHB2"
all_dbs$termName[grep('^EphB3$',all_dbs$termName)] <- "EPHB3"
all_dbs$termName[grep('^EphB4$',all_dbs$termName)] <- "EPHB4"
all_dbs$termName[grep('^ALK4$',all_dbs$termName)] <- "ACVR1B"
all_dbs$termName[grep('^TRKB$',all_dbs$termName)] <- "NTRK2"
all_dbs$termName[grep('^RON$',all_dbs$termName)] <- "MST1R"
all_dbs$termName[grep('^PHKg1$',all_dbs$termName)] <- "PHKG1"
all_dbs$termName[grep('^NDR2$',all_dbs$termName)] <- "STK38L"
all_dbs$termName[grep('^PYK2$',all_dbs$termName)] <- "PTK2B"
all_dbs$termName[grep('^BRK$',all_dbs$termName)] <- "PTK6"
all_dbs$termName[grep('^Brk$',all_dbs$termName)] <- "PTK6"
all_dbs$termName[grep('^TRKA$',all_dbs$termName)] <- "NTRK1"
all_dbs$termName[grep('^GPRK6$',all_dbs$termName)] <- "GRK6"
all_dbs$termName[grep('^GPRK5$',all_dbs$termName)] <- "GRK5"
all_dbs$termName[grep('^GPRK4$',all_dbs$termName)] <- "GRK4"
all_dbs$termName[grep('^Eg3 kinase$',all_dbs$termName)] <- "MELK"
all_dbs$termName[grep('^TIE2$',all_dbs$termName)] <- "TEK"
all_dbs$termName[grep('^YSK1$',all_dbs$termName)] <- "STK25"
all_dbs$termName[grep('^RSK3$',all_dbs$termName)] <- "RPS6KA2"
all_dbs$termName[grep('^CDC2$',all_dbs$termName)] <- "CDK1"
all_dbs$termName[grep('^QIK$',all_dbs$termName)] <- "SIK2"
all_dbs$termName[grep('^RHOK$',all_dbs$termName)] <- "GRK1"
all_dbs$termName[grep('^LOK$',all_dbs$termName)] <- "STK10"
all_dbs$termName[grep('^ACK$',all_dbs$termName)] <- "TNK2"
all_dbs$termName[grep('^QSK$',all_dbs$termName)] <- "SIK3"
all_dbs$termName[grep('^NDR1$',all_dbs$termName)] <- "STK38"
all_dbs$termName[grep('^PKR$',all_dbs$termName)] <- "EIF2AK2"
all_dbs$termName[grep('^MST3$',all_dbs$termName)] <- "STK24"
all_dbs$termName[grep('^CRIK$',all_dbs$termName)] <- "CIT"
all_dbs$termName[grep('^Wnk4$',all_dbs$termName)] <- "WNK4"
all_dbs$termName[grep('^p70S6Kb$',all_dbs$termName)] <- "RPS6KB2"
all_dbs$termName[grep('^PKG2$',all_dbs$termName)] <- "PRKG2"
all_dbs$termName[grep('^PKG1$',all_dbs$termName)] <- "PRKG1"

all_dbs$termName[grep('^MRCKb$',all_dbs$termName)] <- "CDC42BPB"
all_dbs$termName[grep('^PKG1cGKI$',all_dbs$termName)] <- "PKG1/cGK-I"
all_dbs$termName[grep('^PKG2cGKII$',all_dbs$termName)] <- "PKG2/cGK-II"

## Unified database
all_dbs$dbName <- "DBS"
all_dbs <- subset(all_dbs, select = -termID)
all_dbs <- all_dbs[!duplicated(all_dbs),]
dd <- data.frame(unique(all_dbs$termName),paste("DBS_",seq(1:length(unique(all_dbs$termName))), sep=""))
all_dbs$termID <- dd$paste..DBS_...seq.1.length.unique.all_dbs.termName.....sep......[match(all_dbs$termName, dd$unique.all_dbs.termName.)]
all_dbs <- all_dbs[, c(1,5,2,3,4)]
all_dbs <- all_dbs[!duplicated(all_dbs),]
all_dbs$termName <- as.factor(all_dbs$termName)

## Remove kinases with one entry
table_count <- table(all_dbs$termName)
table_count_1 <- table_count[table_count == 1]
all_dbs <- all_dbs[!all_dbs$termName %in% names(table_count_1),]

## Write output DB
write.csv(all_dbs, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/all_dbs.csv"), row.names = FALSE)
