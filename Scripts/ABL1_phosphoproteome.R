#################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Test the levels of ABL1 in the SN nonphosphorylated part
## Date: 22.08.17
## Author: Mahmoud Hallal
##################################################
#the protein of ABL1 is P00519
#With no phosphorylation sites
params <- read_yaml("./config.yaml")
biological_reps <- strsplit(params$Biological_replicates,"-")[[1]]
imp <- params$Imputation

## Load input expressionSet
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))

## Load database 
all_dbs <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/all_dbs.csv"))


##
rr <- read.delim('../../cell_lines_data/Stimulated_data/K562_Nilotinib_t2/SN/proteinGroups.txt')
View(rr)
rr2 <- read.delim('../../cell_lines_data/Stimulated_data/K562_Nilotinib_t2/SN/evidence.txt')

##########################################################################################
#WIth the phosphorylation sites 
cc <- as.data.frame(exprs(test_eSet)[grep("P51692", rownames(exprs(test_eSet))),,drop=F])

##
cc <- as.data.frame(t(cc))
cc$Replicate <- gsub("X","",rownames(cc))
cc$Replicate <- factor(cc$Replicate, levels = cc$Replicate)

cc <- melt(cc)
#cc <- cc[-13,]
cc$value <- 2^cc$value

pp <- ggplot(cc, aes(x = Replicate, y = value, fill = variable )) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Intensity levels of SIRT1 (\"Q96EB6\") phosphorylation sites in the PP") +
  theme(plot.title = element_text(size = 12)) 
  #scale_x_discrete(limits=c("1_K562","2_K562","3_K562","1_NB4","2_NB4","3_NB4",
  #                          "1_THP1","2_THP1","3_THP1","1_OCI","2_OCI","3_OCI",
  #                          "1_MOLM13","2_MOLM13","3_MOLM13"))
  
  #scale_x_discrete(limits=c("1_MOLM13_CTRL","2_MOLM13_CTRL","3_MOLM13_CTRL", "1_MOLM13_DRG","2_MOLM13_DRG","3_MOLM13_DRG"))
pp

#BL-BQ
#parallel coordinates plot
library(GGally)
cc <- as.data.frame(exprs(test_eSet)[grep("FLT3", rownames(exprs(test_eSet))),])
cc$Replicate <- gsub("X","",rownames(cc))
#cc <- cc[,c("X1_K562_CTRL","X2_K562_CTRL","X3_K562_CTRL", "X1_K562_DRG","X2_K562_DRG","X3_K562_DRG","Replicate")]
#colnames(cc) <- c("K562_CTRL1","K562_CTRL2","K562_CTRL3", "K562_DRG1","K562_DRG2","K562_DRG3","Replicate")
ggparcoord(cc, columns = 1:6, groupColumn = "Replicate", scale = "globalminmax") +theme(axis.text.x = element_text(angle = 90, hjust = 1))

rr <- list(FYN_subs,LYN_subs,HCK_subs)
names(rr) <- c('FYN','LYN','HCK')
ll <- venn.diagram(rr,
             fill = 2:4,
             filename = NULL,
             alpha = rep(0.5,3),
             cat.fontface = 4,
             lty =2,
             cex = 1.1,
             scaled = TRUE,
             euler.d = TRUE ,
             main = "Phosphorylated Proteins_phosphorylation site of K562 Ctrl and Drug",
             main.cex = 1.2)

grid.draw(ll)

#Get all the substrates of ABL1
abl1_subs <- unique(all_dbs[grep('^FLT3$', all_dbs$termName),]$geneID)
cc <- as.data.frame(exprs(test_eSet))[rownames(exprs(test_eSet)) %in% abl1_subs,]
cc <- as.data.frame(t(cc))
cc$Replicate <- gsub("X","",rownames(cc))
cc$Replicate <- factor(cc$Replicate, levels = cc$Replicate)


cc <- melt(cc)
cc <- cc[!is.na(cc$value),]
cc$value <- 2^cc$value
#cc <- cc[cc$variable == "P42229_Y694",]
#colnames(cc) <- c("Replicate","Phosphosite","Quantification")



#FYN_subs <- unique(colnames(cc))
#LYN_subs <- unique(colnames(cc))
#HCK_subs <- unique(colnames(cc))
#load('../results/K562_0.025P_0.025FDR/topTables_K562_less.Rda')
#pp <- topTables_less$K562_DRG
#pp2 <- pp[LYN_subs,]
# plot_ly(x = cc$Replicate, y=cc$value,name=cc$variable, type = "bar") %>%
#   layout(title = paste0("Intensities of substrates of CDK1"), barmode = 'stack')

pp <- ggplot(cc, aes(x = Replicate, y = value, fill = variable )) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #ggtitle("Intensity levels of ABL1 (\"P00519\") substrate phosphorylation sites in the PP") +
  theme(plot.title = element_text(size = 12))  #+ guides(fill=FALSE)+
  #scale_x_discrete(limits=c("1_K562","2_K562","3_K562","1_NB4","2_NB4","3_NB4",
  #                          "1_THP1","2_THP1","3_THP1","1_OCI","2_OCI","3_OCI",
  #                          "1_MOLM13","2_MOLM13","3_MOLM13"))
  #scale_x_discrete(limits=c("1_K562_CTRL","2_K562_CTRL","3_K562_CTRL", "1_K562_DRG","2_K562_DRG","3_K562_DRG"))
  #scale_x_discrete(limits=c("1_MOLM13_CTRL","2_MOLM13_CTRL","3_MOLM13_CTRL", "1_MOLM13_DRG","2_MOLM13_DRG","3_MOLM13_DRG"))
  #scale_x_discrete(limits=c("1_MOLM13H_SC","2_MOLM13H_SC","3_MOLM13H_SC", "1_MOLM13H_SD","2_MOLM13H_SD","3_MOLM13H_SD",
  #                          "1_MOLM13H_RC","2_MOLM13H_RC","3_MOLM13H_RC","1_MOLM13H_RD","2_MOLM13H_RD","3_MOLM13H_RD"))+
  #theme(legend.position = "none")
  
pp

#silac#
library(ggbeeswarm)
ggplot(cc) + 
  geom_quasirandom(aes(x = Replicate, y = value),size=0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
#


colnames(cc) <- c("Intensity", "Replicate")
pp <- ggplot(cc, aes(x = Replicate, y = Intensity)) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Intensity levels of FLT3 (\"P36888\") in the different replicates in the Phosphoproteome") +
  theme(plot.title = element_text(size = 12))
pp

#
#parallel coordinates plot
library(GGally)
cc <- as.data.frame(exprs(test_eSet)[rownames(exprs(test_eSet)) %in% abl1_subs,])
cc$Replicate <- gsub("X","",rownames(cc))
ggparcoord(cc, columns = 1:6, groupColumn = "Replicate", scale = "globalminmax") +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

