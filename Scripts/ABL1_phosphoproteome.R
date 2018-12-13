#################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Test the levels of ABL1 in the SN nonphosphorylated part
## Date: 22.08.17
## Author: Mahmoud Hallal
##################################################
#the protein of ABL1 is P00519
#With no phosphorylation sites
load("/Users/Mahmoud.Hallal/Desktop/PhD/Generalized_pipeline/results/CML_0.01P_0.05FDR/test_eSet1_CML.Rda")

##########################################################################################
#WIth the phosphorylation sites 
cc <- as.data.frame(exprs(test_eSet)[grep("P42229", rownames(exprs(test_eSet))),])

##
cc <- as.data.frame(t(cc))
cc$Replicate <- gsub("X","",rownames(cc))
cc <- melt(cc)
#cc <- cc[-13,]
cc$value <- 10^cc$value
pp <- ggplot(cc, aes(x = Replicate, y = value, fill = variable )) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Intensity levels of ABL1 (\"P00519\") phosphorylation sites in the PP") +
  theme(plot.title = element_text(size = 12)) +
  #scale_x_discrete(limits=c("1_K562","2_K562","3_K562","1_NB4","2_NB4","3_NB4",
  #                          "1_THP1","2_THP1","3_THP1","1_OCI","2_OCI","3_OCI",
  #                          "1_MOLM13","2_MOLM13","3_MOLM13"))
  
  scale_x_discrete(limits=c("1_MOLM13_CTRL","2_MOLM13_CTRL","3_MOLM13_CTRL", "1_MOLM13_DRG","2_MOLM13_DRG","3_MOLM13_DRG"))
pp


#parallel coordinates plot
library(GGally)
cc <- as.data.frame(exprs(test_eSet)[grep("P00519", rownames(exprs(test_eSet))),])
cc$Replicate <- gsub("X","",rownames(cc))
cc <- cc[,c("X1_K562_CTRL","X2_K562_CTRL","X3_K562_CTRL", "X1_K562_DRG","X2_K562_DRG","X3_K562_DRG","Replicate")]
colnames(cc) <- c("K562_CTRL1","K562_CTRL2","K562_CTRL3", "K562_DRG1","K562_DRG2","K562_DRG3","Replicate")
ggparcoord(cc, columns = 1:6, groupColumn = "Replicate", scale = "globalminmax") +theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Get all the substrates of ABL1
abl1_subs <- unique(all_dbs[grep('^ABL1$', all_dbs$termName),]$geneID)

cc <- as.data.frame(exprs(test_eSet))[rownames(exprs(test_eSet)) %in% abl1_subs,]
cc <- as.data.frame(t(cc))
cc$Replicate <- gsub("X","",rownames(cc))
cc <- melt(cc)
cc$value <- 2^cc$value

plot_ly(x = cc$Replicate, y=cc$value,name=cc$variable, type = "bar") %>%
  layout(title = paste0("Intensities of substrates of CDK1"), barmode = 'stack')

pp <- ggplot(cc, aes(x = Replicate, y = value, fill = variable )) + 
  geom_bar(stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Intensity levels of ABL1 (\"P00519\") substrate phosphorylation sites in the PP") +
  theme(plot.title = element_text(size = 12)) 
  #scale_x_discrete(limits=c("1_K562","2_K562","3_K562","1_NB4","2_NB4","3_NB4",
  #                          "1_THP1","2_THP1","3_THP1","1_OCI","2_OCI","3_OCI",
  #                          "1_MOLM13","2_MOLM13","3_MOLM13"))
  scale_x_discrete(limits=c("1_MOLM13_CTRL","2_MOLM13_CTRL","3_MOLM13_CTRL", "1_MOLM13_DRG","2_MOLM13_DRG","3_MOLM13_DRG"))
pp



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

