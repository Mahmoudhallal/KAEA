##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Build database of kinase;substrate relationship from multiple databases
## Date: 09.03.17
## Author: Mahmoud Hallal
##################################################

## Load libraries##
#install.packages("SetRank", repos="http://cran.rstudio.com/")
library(SetRank)#

source("https://bioconductor.org/biocLite.R")

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

## Load parameters
params <- read_yaml("./config.yaml")

## Load kinase-substrate databases
phosphoELM_vertebrate <- read.delim("/home/user/KAEA/Phospho_database/phosphoELM_db/phosphoELM_vertebrate_2015-04.dump")
Kinase_Substrate_Dataset <- read.delim("/home/user/KAEA/Phospho_database/phosphoSitePlus/Kinase_Substrate_Dataset_v2")
HPRD <- read.delim('/home/user/KAEA/Phospho_database/HPRD/FLAT_FILES_072010/POST_TRANSLATIONAL_MODIFICATIONS.txt', header = FALSE)
RegPhos <- read.delim('/home/user/KAEA/Phospho_database/RegPhos/RegPhos_Phos_human.txt')
Signor <- read.delim('/home/user/KAEA/Phospho_database/SIGNOR/human_phosphorylations_12_07_18.csv', sep = ",")

## Create unified column names
column_names <- c("geneID","termID","termName","dbName","description")

####################################################################################
## ELM
#only homo sapiens
phosphoELM_vertebrate <- phosphoELM_vertebrate[phosphoELM_vertebrate$species == "Homo sapiens",]
#remove empty cases and filter
phosphoELM_vertebrate <- phosphoELM_vertebrate[!phosphoELM_vertebrate$kinases == "",]
phosphoELM_vertebrate$description <- paste(phosphoELM_vertebrate$code,phosphoELM_vertebrate$position,sep="_")
ELM <- phosphoELM_vertebrate[,c("acc","kinases","position","code")]
ELM$acc <- sub("(\\w+)(-\\d?)", "\\1", ELM$acc)
ELM <- ELM[!duplicated(ELM),]
ELM$domain <- "domain"
ELM$dbName <- "ELM"
ELM$acc <- paste(ELM$acc, '_', ELM$code, ELM$position, sep="")

#Add substrateID table
subs_df <- data.frame(unique(ELM$kinases),paste("ELM_",seq(1:length(unique(ELM$kinases))), sep=""))
ELM$substrateID <- subs_df$paste..ELM_...seq.1.length.unique.ELM.kinases....[match(ELM$kinases, subs_df$unique.ELM.kinases.)]

#Rearrange columns
ELM <- ELM[,c(1,7,2,6,5)]
colnames(ELM) <- column_names

####################################################################################
## PhosphoSitePlus (PSP)
#In vitro and In vivo were taken into consideration
Kinase_Substrate_Dataset <- Kinase_Substrate_Dataset[Kinase_Substrate_Dataset$SUB_ORGANISM == "human",]
Kinase_Substrate_Dataset <- Kinase_Substrate_Dataset[Kinase_Substrate_Dataset$KIN_ORGANISM == "human",]

# Arrange and filter data
PSP <- Kinase_Substrate_Dataset[,c('SUB_ACC_ID',"GENE","SUB_MOD_RSD")]
PSP <- PSP[!duplicated(PSP),]
PSP$domain <- "domain"
PSP$dbName <- "PhosphoSitePlus"
#isoforms
PSP$SUB_ACC_ID <- sub("(\\w+)(-\\d?)", "\\1", PSP$SUB_ACC_ID)
subs_PSP_df <- data.frame(unique(PSP$GENE),paste("PSP_",seq(1:length(unique(PSP$GENE))), sep=""))
PSP$substrateID <- subs_PSP_df$paste..PSP_...seq.1.length.unique.PSP.GENE....[match(PSP$GENE, subs_PSP_df$unique.PSP.GENE.)]

PSP$SUB_ACC_ID <- paste(PSP$SUB_ACC_ID, PSP$SUB_MOD_RSD, sep = "_")

#Rearrange columns
PSP <- PSP[,c(1,6,2,5,4)]
colnames(PSP) <- column_names



####################################################################################
## HPRD

#Define converter from Symbol to Uniprot
IDConverter = createIDConverter("Homo.sapiens", "SYMBOL", "UNIPROT")

# Select columns and filter
hprd_colnames <- c("substrate_hprd_id","substrate_gene_symbol","substrate_isoform_id","substrate_refseq_id","site","residue","enzyme_name","enzyme_hprd_id","modification_type","experiment_type","reference_id")
colnames(HPRD) <- hprd_colnames
keep_hprd_colnames <- c("substrate_gene_symbol","site","residue","enzyme_name","modification_type","site","residue")
HPRD_df <- HPRD[,keep_hprd_colnames]
HPRD_df <- HPRD_df[HPRD_df$modification_type == "Phosphorylation",]
HPRD_df <- HPRD_df[!HPRD_df$enzyme_name == "-",]


ll <- sapply(1:nrow(HPRD_df), function(x) IDConverter(as.character(HPRD_df$substrate_gene_symbol[x])))
ll_length <- sapply(1:length(ll), function(x) length(unlist(ll[x])))
kinases <- rep(HPRD_df$enzyme_name ,ll_length)
phos_pos <- paste(HPRD_df$residue, HPRD_df$site, sep="")
phos_pos <- rep(phos_pos, ll_length)
HPRD_df_new <- data.frame(kinases, as.data.frame(unlist(ll)), phos_pos)
HPRD_df_new <- HPRD_df_new[!duplicated(HPRD_df_new),]
HPRD_df_new$unlist.ll. <- paste(HPRD_df_new$unlist.ll., HPRD_df_new$phos_pos, sep="_")


HPRD_df2 <- data.frame(unique(HPRD_df_new$kinases),paste("HPRD_",seq(1:length(unique(HPRD_df_new$kinases))), sep=""))
HPRD_df_new$substrateID <- HPRD_df2$paste..HPRD_...seq.1.length.unique.HPRD_df_new.kinases....[match(HPRD_df_new$kinases, HPRD_df2$unique.HPRD_df_new.kinases.)]
HPRD_df_new$domain <- "domain"
HPRD_df_new$dbName <- "HPRD"

#remove ;- from some proteins
HPRD_df_new$unlist.ll. <- gsub(";-","", HPRD_df_new$unlist.ll.)

# Rearrange columns
HPRD_df_new <- HPRD_df_new[,c(2,4,1,6,5)]
colnames(HPRD_df_new) <- column_names

####################################################################################
## RegPhos
RegPhos_new <- RegPhos[,c("ID", "AC", "position", "catalytic.kinase", "code")]
RegPhos_new$domain <- paste(RegPhos_new$code, RegPhos_new$position, sep="")
RegPhos_new$ID <- sub("(\\w+)_HUMAN", '\\1', RegPhos_new$ID)
RegPhos_new <- RegPhos_new[!RegPhos_new$catalytic.kinase=="",]

#Change MAPK format
only_MAPKs <- RegPhos_new$catalytic.kinase[grep("\\(MAPK", RegPhos_new$catalytic.kinase)]
RegPhos_new$catalytic.kinase <- as.character(RegPhos_new$catalytic.kinase)
RegPhos_new$catalytic.kinase[grep("\\(MAPK", RegPhos_new$catalytic.kinase)] <-
  sub("(\\w+)\\((\\w*)\\)", "\\2", RegPhos_new$catalytic.kinase[grep("\\(MAPK", RegPhos_new$catalytic.kinase)])

#Change ABL format
only_ABL <- RegPhos_new$catalytic.kinase[grep("ABL[0-9]", RegPhos_new$catalytic.kinase)]
RegPhos_new$catalytic.kinase[grep("ABL[0-9]", RegPhos_new$catalytic.kinase)] <-
  sub("(\\w+)\\((\\w*)\\)", "\\1", RegPhos_new$catalytic.kinase[grep("ABL[0-9]", RegPhos_new$catalytic.kinase)])


#
RegPhos_df <- data.frame(unique(RegPhos_new$catalytic.kinase),paste("REGPHOS_",seq(1:length(unique(RegPhos_new$catalytic.kinase))), sep=""))
RegPhos_new$substrateID <- RegPhos_df$paste..REGPHOS_...seq.1.length.unique.RegPhos_new.catalytic.kinase....[match(RegPhos_new$catalytic.kinase, RegPhos_df$unique.RegPhos_new.catalytic.kinase.)]
RegPhos_new$dbName <- "RegPhos"
RegPhos_new$AC <- paste(RegPhos_new$AC, RegPhos_new$domain, sep = "_")
#
RegPhos_new$domain <- "domain"
RegPhos_new2 <- RegPhos_new[,c(2,7,4,8,6)]
colnames(RegPhos_new2) <- column_names
RegPhos_new2 <- RegPhos_new2[!duplicated(RegPhos_new2),]

####################################################################################
## Signor
Signor_phos <- Signor[Signor$MECHANISM == "phosphorylation",]
Signor_human <- Signor[Signor_phos$TAX_ID == "9606",]
Signor_select_cols <- Signor_human[,c("ENTITYA","IDA","ENTITYB","IDB","RESIDUE")]
colnames(Signor_select_cols) <- c("Kinase_name","kinase_ID","Substrate_name","Substrate_ID","Residue")
Signor_select_cols$Residue <- as.character(Signor_select_cols$Residue)
Signor_select_cols$Residue[grep("Tyr", Signor_select_cols$Residue)] <- sub("Tyr","Y", Signor_select_cols$Residue[grep("Tyr", Signor_select_cols$Residue)])
Signor_select_cols$Residue[grep("Ser", Signor_select_cols$Residue)] <- sub("Ser","S", Signor_select_cols$Residue[grep("Ser", Signor_select_cols$Residue)])
Signor_select_cols$Residue[grep("Thr", Signor_select_cols$Residue)] <- sub("Thr","T", Signor_select_cols$Residue[grep("Thr", Signor_select_cols$Residue)])
Signor_select_cols$Substrate_ID <- paste(Signor_select_cols$Substrate_ID, Signor_select_cols$Residue, sep = "_")
Signor_reduced <- Signor_select_cols[,c(4,1,5)]

Signor_all <- Signor_reduced[!Signor_reduced$Residue == "",]

#
Signor_df <- data.frame(unique(Signor_all$Kinase_name),paste("SIGNOR_",seq(1:length(unique(Signor_all$Kinase_name))), sep=""))
Signor_all$substrateID <- Signor_df$paste..SIGNOR_...seq.1.length.unique.Signor_all.Kinase_name.....[match(Signor_all$Kinase_name, Signor_df$unique.Signor_all.Kinase_name.)]
Signor_all$dbName <- "SIGNOR"

# Rearrange columns
Signor_new <- Signor_all[,c(1,4,2,5)]
Signor_new$domain <- 'domain'
colnames(Signor_new) <- column_names
Signor_new <- Signor_new[complete.cases(Signor_new),]
Signor_new <- Signor_new[-grep("SIGNOR", Signor_new$geneID),]
Signor_new <- Signor_new[!duplicated(Signor_new),]

####################################################################################
## NetworKIN#
nwkin_input <- params$NWKIN_Input
nwkin <- read.delim(nwkin_input,sep=",")
nwkin$residue <- substr(nwkin$sequence,6,6)
nwkin <- nwkin[!nwkin$residue == "?",]
nwkin$residue <- toupper(nwkin$residue)
nwkin$X.substrate <- gsub("(\\w+)(-\\d+)?","\\1",nwkin$X....substrate)
nwkin$protein_with_pos <- paste0(nwkin$X.substrate,'_', nwkin$residue, nwkin$position)
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

# nwkin_input <- params$NWKIN_Input
# nwkin <- read.delim(nwkin_input,sep=",")
# nwkin <- nwkin[!nwkin$sub_uninprot == "",]
# nwkin$residue <- toupper(as.character(nwkin$aa))
# nwkin$protein_with_pos <- paste0(nwkin$sub_uninprot,'_', nwkin$residue, nwkin$position)
# nwkin <- nwkin[!duplicated(nwkin),]
# 
# nwkin <- nwkin[,c('protein_with_pos','id')]
# 
# motif_df <- data.frame(unique(nwkin$id),paste("MTF_",seq(1:length(unique(nwkin$id))), sep=""))
# nwkin$substrateID <- motif_df$paste..MTF_...seq.1.length.unique.nwkin.id.....sep......[match(nwkin$id, motif_df$unique.nwkin.id.)]
# nwkin$dbName <- "MTF"
# 
# # Rearrange columns
# motif_new <- nwkin[,c(1,3,2,4)]
# motif_new$domain <- 'domain'
# colnames(motif_new) <- column_names
# motif_new <- motif_new[complete.cases(motif_new),]
# motif_new <- motif_new[!duplicated(motif_new),]

####################################################################################
#Merge the all databases
all_dbs <- rbind(PSP, ELM, HPRD_df_new, RegPhos_new2, Signor_new,motif_new)
all_dbs$termName <- as.character(all_dbs$termName)

## Unify names given by different databases
unique(all_dbs[grep('^CK',all_dbs$termName),]$termName)
unique(all_dbs[grep('^CSNK',all_dbs$termName),]$termName)

all_dbs[grep('^CK2_group',all_dbs$termName),]$termName <- "CSNK2"

all_dbs[grep('^CK2_?alpha',all_dbs$termName),]$termName <- "CSNK2A"
all_dbs[grep('^CK2a',all_dbs$termName),]$termName <- "CSNK2A"
all_dbs[grep('^CSNK2A1',all_dbs$termName),]$termName <- "CSNK2A"
all_dbs[grep('^CSNK2A2',all_dbs$termName),]$termName <- "CSNK2A"

all_dbs[grep('^CK2B',all_dbs$termName),]$termName <- "CSNK2B"
all_dbs[grep('^CK2_beta',all_dbs$termName),]$termName <- "CSNK2B"

#CK1
unique(all_dbs[grep('^CSNK',all_dbs$termName),]$termName)
unique(all_dbs[grep('^CK1',all_dbs$termName),]$termName)

all_dbs[grep('^CK1alpha',all_dbs$termName),]$termName <- "CSNK1A"
all_dbs[grep('^CK1a',all_dbs$termName),]$termName <- "CSNK1A"
all_dbs[grep('^CSNK1A1',all_dbs$termName),]$termName <- "CSNK1A"
all_dbs[grep('^CK1_delta',all_dbs$termName),]$termName <- "CSNK1D"
all_dbs[grep('^CK1d',all_dbs$termName),]$termName <- "CSNK1D"
#all_dbs[grep('^CK1delta',all_dbs$termName),]$termName <- "CSNK1D"
all_dbs[grep('^CK1e',all_dbs$termName),]$termName <- "CSNK1E"
all_dbs[grep('^CK1_epsilon',all_dbs$termName),]$termName <- "CSNK1E"

#all_dbs[grep('^CK1gamma2',all_dbs$termName),]$termName <- "CSNK1G2"
##DOUBLE CHECK THIS #06.11.2018
all_dbs[grep('^CK$',all_dbs$termName),]$termName <- "CSNK2"

#MRKCb
unique(all_dbs[grep('^CDC42',all_dbs$termName),]$termName)
all_dbs[grep('^MRCKb$',all_dbs$termName),]$termName <- "CDC42BPB"

#SRC
unique(all_dbs[grep('^SRC',all_dbs$termName),]$termName)
all_dbs[grep('^SRC_group$',all_dbs$termName),]$termName <- "SRC"
all_dbs[grep('^SRC-type Tyr-kinases$',all_dbs$termName),]$termName <- "SRC"


#GSK3
unique(all_dbs[grep('^GSK',all_dbs$termName),]$termName)

all_dbs[grep('^GSK-3_alpha',all_dbs$termName),]$termName <- "GSK3A"
all_dbs[grep('^GSK3alpha',all_dbs$termName),]$termName <- "GSK3A"
#
all_dbs[grep('^GSK-3_beta',all_dbs$termName),]$termName <- "GSK3B"
all_dbs[grep('^GSK3beta',all_dbs$termName),]$termName <- "GSK3B"
#
all_dbs[grep('^GSK-3_group',all_dbs$termName),]$termName <- "GSK3"
all_dbs[grep('^GSK3_group',all_dbs$termName),]$termName <- "GSK3"

#Aurora
#A
unique(all_dbs[grep('^Aur',all_dbs$termName),]$termName)
unique(all_dbs[grep('^AUR',all_dbs$termName),]$termName)

all_dbs[grep('Aurora A',all_dbs$termName),]$termName <- "AURKA"
all_dbs[grep('AurA',all_dbs$termName),]$termName <- "AURKA"
all_dbs[grep('AuroraA',all_dbs$termName),]$termName <- "AURKA"

#B
all_dbs[grep('Aurora B',all_dbs$termName),]$termName <- "AURKB"
all_dbs[grep('AurB',all_dbs$termName),]$termName <- "AURKB"

#PKB
unique(all_dbs[grep('^PKB',all_dbs$termName),]$termName)

all_dbs[grep('PKB_group',all_dbs$termName),]$termName <- "AKT"
all_dbs[grep('PKBalpha',all_dbs$termName),]$termName <- "AKT1"
#all_dbs[grep('PKB_beta',all_dbs$termName),]$termName <- "AKT2"
#all_dbs[grep('PKBbeta',all_dbs$termName),]$termName <- "AKT2"

#PKC
unique(all_dbs[grep('PKC',all_dbs$termName),]$termName)
all_dbs[grep('PKCalpha',all_dbs$termName),]$termName <- "PRKCA"
all_dbs[grep('PKCa',all_dbs$termName),]$termName <- "PRKCA"
all_dbs[grep('PKC_alpha',all_dbs$termName),]$termName <- "PRKCA"

all_dbs[grep('PKCbeta',all_dbs$termName),]$termName <- "PRKCB"
all_dbs[grep('PKCb',all_dbs$termName),]$termName <- "PRKCB"
all_dbs[grep('PKC_beta',all_dbs$termName),]$termName <- "PRKCB"

all_dbs[grep('PKCdelta',all_dbs$termName),]$termName <- "PRKCD"
all_dbs[grep('PKCd',all_dbs$termName),]$termName <- "PRKCD"
all_dbs[grep('PKC_delta',all_dbs$termName),]$termName <- "PRKCD"

all_dbs[grep('PKCepsilon',all_dbs$termName),]$termName <- "PRKCE"
all_dbs[grep('PKCe',all_dbs$termName),]$termName <- "PRKCE"
all_dbs[grep('PKC_epsilon',all_dbs$termName),]$termName <- "PRKCE"

all_dbs[grep('PKCtheta',all_dbs$termName),]$termName <- "PRKCT"
all_dbs[grep('PKCt',all_dbs$termName),]$termName <- "PRKCT"
all_dbs[grep('PKC_theta',all_dbs$termName),]$termName <- "PRKCT"

all_dbs[grep('PKCgamma',all_dbs$termName),]$termName <- "PRKCG"
all_dbs[grep('PKCg',all_dbs$termName),]$termName <- "PRKCG"
all_dbs[grep('PKC_gamma',all_dbs$termName),]$termName <- "PRKCG"

all_dbs[grep('PKCz',all_dbs$termName),]$termName <- "PRKCZ"
all_dbs[grep('PKC_zeta',all_dbs$termName),]$termName <- "PRKCZ"

all_dbs[grep('PKCiota',all_dbs$termName),]$termName <- "PRKCI"
all_dbs[grep('PKCi',all_dbs$termName),]$termName <- "PRKCI"
all_dbs[grep('PKC_iota',all_dbs$termName),]$termName <- "PRKCI"

all_dbs[grep('PKC_eta',all_dbs$termName),]$termName <- "PRKCH"
all_dbs[grep('PKCh',all_dbs$termName),]$termName <- "PRKCH"

all_dbs[grep('PKC_group',all_dbs$termName),]$termName <- "PRKC"

##IKK
unique(all_dbs[grep('IKK',all_dbs$termName),]$termName)

all_dbs[grep('IKKa',all_dbs$termName),]$termName <- "IKKA"
all_dbs[grep('CHUK',all_dbs$termName),]$termName <- "IKKA"
all_dbs[grep('IKK_alpha',all_dbs$termName),]$termName <- "IKKA"

all_dbs[grep('IKKb',all_dbs$termName),]$termName <- "IKKB"
all_dbs[grep('IKBKB',all_dbs$termName),]$termName <- "IKKB"
all_dbs[grep('IKK_beta',all_dbs$termName),]$termName <- "IKKB"

all_dbs[grep('IKKe',all_dbs$termName),]$termName <- "IKKE"
all_dbs[grep('IKBKE',all_dbs$termName),]$termName <- "IKKE"
all_dbs[grep('IKK_epsilon',all_dbs$termName),]$termName <- "IKKE"

all_dbs[grep('IKK-complex',all_dbs$termName),]$termName <- "IKK"
all_dbs[grep('IKK_group',all_dbs$termName),]$termName <- "IKK"


##BCR/ABL
unique(all_dbs[grep('ABL',all_dbs$termName),]$termName)
all_dbs[grep('BCR/ABL',all_dbs$termName),]$termName <- "ABL1"
all_dbs[grep('BCR-ABL',all_dbs$termName),]$termName <- "ABL1"

all_dbs[grep('^ABL$',all_dbs$termName),]$termName <- "ABL1"
all_dbs[grep('^Abl$',all_dbs$termName),]$termName <- "ABL1"

#PKA
all_dbs[grep('^PKAalpha',all_dbs$termName),]$termName <- "PRKACA"
all_dbs[grep('^PKACa',all_dbs$termName),]$termName <- "PRKACA"
all_dbs[grep('^PKA_alpha',all_dbs$termName),]$termName <- "PRKACA"

#CAMK2
unique(all_dbs[grep('CaM-K',all_dbs$termName),]$termName)
unique(all_dbs[grep('CaMK',all_dbs$termName),]$termName)

all_dbs[grep('^CaMKIIdelta',all_dbs$termName),]$termName <- "CaMK2D"

all_dbs[grep('^CaMKIIbeta',all_dbs$termName),]$termName <- "CaMK2B"
all_dbs[grep('^CaMK2b',all_dbs$termName),]$termName <- "CaMK2B"

all_dbs[grep('^CaMKIIalpha',all_dbs$termName),]$termName <- "CAMK2A"
all_dbs[grep('^CaMK2a',all_dbs$termName),]$termName <- "CAMK2A"

#all_dbs[grep('^CaMKIIa',all_dbs$termName),]$termName <- "CaMK2a"

all_dbs[grep('^CaM-KII_alpha',all_dbs$termName),]$termName <- "CAMK2A"
all_dbs[grep('^CaM-KII_group',all_dbs$termName),]$termName <- "CAMK2"
#all_dbs[grep('^CaM-CaMK2_group',all_dbs$termName),]$termName <- "CAMK2"

all_dbs[grep('^CaM-KI_alpha',all_dbs$termName),]$termName <- "CaMK1a"
all_dbs[grep('^CaMK1a',all_dbs$termName),]$termName <- "CAMK1A"

all_dbs[grep('^CaM-KIV',all_dbs$termName),]$termName <- "CAMK4"
all_dbs[grep('^CaMK4',all_dbs$termName),]$termName <- "CAMK4"

all_dbs[grep('^CaM-KK_alpha',all_dbs$termName),]$termName <- "CAMKK1"
all_dbs[grep('^CaMKK1',all_dbs$termName),]$termName <- "CAMKK1"

#PDK
unique(all_dbs[grep('PDK',all_dbs$termName),]$termName)
all_dbs[grep('^PDK-1',all_dbs$termName),]$termName <- "PDK1"
all_dbs[grep('^PDK-2',all_dbs$termName),]$termName <- "PDK2"

#RSK
unique(all_dbs[grep('RSK',all_dbs$termName),]$termName)
all_dbs[grep('^RSK-1',all_dbs$termName),]$termName <- "RSK1"
all_dbs[grep('^RSK-2',all_dbs$termName),]$termName <- "RSK2"
all_dbs[grep('^RSK-3',all_dbs$termName),]$termName <- "RSK3"
all_dbs[grep('^RSK-5',all_dbs$termName),]$termName <- "RSK5"

#GRK
unique(all_dbs[grep('PIM',all_dbs$termName),]$termName)
all_dbs[grep('^GRK-1',all_dbs$termName),]$termName <- "GRK1"
all_dbs[grep('^GRK-2',all_dbs$termName),]$termName <- "GRK2"
all_dbs[grep('^GRK-3',all_dbs$termName),]$termName <- "GRK3"
all_dbs[grep('^GRK-4',all_dbs$termName),]$termName <- "GRK4"
all_dbs[grep('^GRK-5',all_dbs$termName),]$termName <- "GRK5"
all_dbs[grep('^GRK-6',all_dbs$termName),]$termName <- "GRK6"

#PIM
unique(all_dbs[grep('PIM',all_dbs$termName),]$termName)
all_dbs[grep('^PIM-1',all_dbs$termName),]$termName <- "PIM1"

#PDGFR
unique(all_dbs[grep('PDGFR',all_dbs$termName),]$termName)
all_dbs[grep('^PDGFR_alpha',all_dbs$termName),]$termName <- "PDGFRA"
all_dbs[grep('^PDGFR_beta',all_dbs$termName),]$termName <- "PDGFRB"
all_dbs[grep('^PDGFRa',all_dbs$termName),]$termName <- "PDGFRA"
all_dbs[grep('^PDGFRb',all_dbs$termName),]$termName <- "PDGFRB"
all_dbs[grep('^PDGFR_group',all_dbs$termName),]$termName <- "PDGFR"


#FRAP and mTOR
unique(all_dbs[grep('RAFT',all_dbs$termName),]$termName)
all_dbs[grep('^FRAP$',all_dbs$termName),]$termName <- "MTOR"
all_dbs[grep('^mTOR$',all_dbs$termName),]$termName <- "MTOR"

#FYN
unique(all_dbs[grep('Fyn',all_dbs$termName),]$termName)
all_dbs[grep('^Fyn$',all_dbs$termName),]$termName <- "FYN"

#CDK
unique(all_dbs[grep('Cyclin',all_dbs$termName),]$termName)
unique(all_dbs[grep('CDK',all_dbs$termName),]$termName)

all_dbs[grep('^CyclinD/CDK4$',all_dbs$termName),]$termName <- "CDK4"
all_dbs[grep('^CyclinE/CDK2$',all_dbs$termName),]$termName <- "CDK2"
all_dbs[grep('^CyclinB/CDK1$',all_dbs$termName),]$termName <- "CDK1"
all_dbs[grep('^CyclinA2/CDK2$',all_dbs$termName),]$termName <- "CDK2"

#Remove autocatalysis
all_dbs <- all_dbs[!all_dbs$termName=="autocatalysis",]

# Remove phosphatases
all_dbs[grep('^PTPRE$',all_dbs$termName),]$termName <- ""
all_dbs[grep('^PTPN11$',all_dbs$termName),]$termName <- ""
all_dbs[grep('^DUSP4$',all_dbs$termName),]$termName <- ""
all_dbs <- all_dbs[!all_dbs$termName=="",]

## Unified database#
all_dbs$dbName <- "DBS"
all_dbs <- subset(all_dbs, select = -termID)
all_dbs <- all_dbs[!duplicated(all_dbs),]
dd <- data.frame(unique(all_dbs$termName),paste("DBS_",seq(1:length(unique(all_dbs$termName))), sep=""))
all_dbs$termID <- dd$paste..DBS_...seq.1.length.unique.all_dbs.termName.....sep......[match(all_dbs$termName, dd$unique.all_dbs.termName.)]
all_dbs <- all_dbs[, c(1,5,2,3,4)]
all_dbs <- all_dbs[!duplicated(all_dbs),]
all_dbs$termName <- as.factor(all_dbs$termName)

## Write output DB
write.csv(all_dbs, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"all_dbs1.csv"), row.names = FALSE)


