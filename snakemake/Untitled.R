phos_STY_INC1 <- read.delim('/Users/Mahmoud.Hallal/Desktop/PhD/cell_lines_data/Stimulated_data/Fribourg/Phospho (STY)Sites.txt')

#Choose the column names to use
to_keep_STY <- c("Protein", "Amino.acid", "Position",  "Gene.names" , "Phospho..STY..Probabilities", "Reverse", "Contaminant", "Localization.prob" ,
                 "Ratio.H.L.normalized.Exp1" , "Ratio.H.L.normalized.Exp2" ,"Ratio.H.L.normalized.Exp3" ,"Ratio.M.L.normalized.Exp4","Ratio.H.L.normalized.Exp4",
                 "Ratio.H.M.normalized.Exp4","Ratio.M.L.normalized.Exp5","Ratio.H.L.normalized.Exp5","Ratio.H.M.normalized.Exp5","Sequence.window")

#define cell lines
cell_lines <- c("U2OS_Rap_WT","U2OS_RapConA_WT","U2OS_RapConA_Rap")

dd <- phos_STY_INC1[,c("Localization.prob.Exp1","Localization.prob.Exp2","Localization.prob.Exp3","Localization.prob.Exp4","Localization.prob.Exp5","Localization.prob")]
rr <- phos_STY_INC1[,c("Localization.prob.Exp1","Ratio.H.L.normalized.Exp1")]
#Select the columns
#colnames(phos)[grep("Intensity", colnames(phos))]
files_to_use_STY <- c("phos_STY_INC1")
for (x in 1:length(files_to_use_STY)) {assign(paste0(files_to_use_STY[x],"_new"), get(files_to_use_STY[x])[,to_keep_STY])}

#add INC columns#
phos_STY_INC1_new$incubation <- "INC1"
phos_STY_INC1_new$Ratio.H.L.normalized.Exp3 <- 1/phos_STY_INC1_new$Ratio.H.L.normalized.Exp3
phos_STY_INC1_new$Ratio.M.L.normalized.Exp5 <- 1/phos_STY_INC1_new$Ratio.M.L.normalized.Exp5
phos_STY_INC1_new$Ratio.H.L.normalized.Exp5 <- 1/phos_STY_INC1_new$Ratio.H.L.normalized.Exp5
#
phos_STY_INC1_new[is.na(phos_STY_INC1_new)] <- 0
#Merge the 2 dfs since they have the same order and column names and
#then divide them according to 3 batches
all_phos_STY <- rbind(phos_STY_INC1_new)
#
##
#Duplicate the P entries to have duplicates, add +1 to each to have difference
extra <- subset(all_phos_STY, select=-c(Ratio.H.L.normalized.Exp1 , Ratio.H.L.normalized.Exp2 ,Ratio.H.L.normalized.Exp3 ,Ratio.M.L.normalized.Exp4, Ratio.H.L.normalized.Exp4, Ratio.H.M.normalized.Exp5, Ratio.M.L.normalized.Exp5,Ratio.H.L.normalized.Exp5,Ratio.H.M.normalized.Exp4))
# extra$Intensity.P9_rep2 <- all_phos_STY$Intensity.P9 + runif(1, 1e-6, 2e-6)
# extra$Intensity.P10_rep2 <- all_phos_STY$Intensity.P10 + runif(1, 1e-6, 2e-6)
# extra$Intensity.P11_rep2 <- all_phos_STY$Intensity.P11 + runif(1, 1e-6, 2e-6)
# extra$Intensity.P12_rep2 <- all_phos_STY$Intensity.P12 + runif(1, 1e-6, 2e-6)

#
extra$"U2OS_Rap_WT_rep1" <- all_phos_STY$Ratio.H.L.normalized.Exp1 #+runif(length(extra$Intensity.P9), -1e-6, 1e-6)
extra$"U2OS_Rap_WT_rep2" <- all_phos_STY$Ratio.H.L.normalized.Exp2 #+runif(length(extra$Intensity.P10), -1e-6, 1e-6)
extra$"U2OS_Rap_WT_rep3" <- all_phos_STY$Ratio.H.L.normalized.Exp3 #+runif(length(extra$Intensity.P11), -1e-6, 1e-6)
extra$"U2OS_Rap_WT_rep4" <- all_phos_STY$Ratio.M.L.normalized.Exp4 #+runif(length(extra$Intensity.P12), -1e-6, 1e-6)
extra$"U2OS_Rap_WT_rep5" <- all_phos_STY$Ratio.H.M.normalized.Exp5 #+runif(length(extra$Intensity.P12), -1e-6, 1e-6)
extra$"U2OS_RapConA_WT_rep1" <- all_phos_STY$Ratio.H.L.normalized.Exp4 #+runif(length(extra$Intensity.P12), -1e-6, 1e-6)
extra$"U2OS_RapConA_WT_rep2" <- all_phos_STY$Ratio.M.L.normalized.Exp5 #+runif(length(extra$Intensity.P12), -1e-6, 1e-6)
extra$"U2OS_RapConA_Rap_rep1" <- all_phos_STY$Ratio.H.M.normalized.Exp4 #+runif(length(extra$Intensity.P12), -1e-6, 1e-6)
extra$"U2OS_RapConA_Rap_rep2" <- all_phos_STY$Ratio.H.L.normalized.Exp5 #+runif(length(extra$Intensity.P12), -1e-6, 1e-6)
extra$"U2OS_RapConA_WT_rep3" <- 0
extra$"U2OS_RapConA_WT_rep4" <- 0
extra$"U2OS_RapConA_WT_rep5" <- 0
extra$"U2OS_RapConA_Rap_rep3" <- 0
extra$"U2OS_RapConA_Rap_rep4" <- 0
extra$"U2OS_RapConA_Rap_rep5" <- 0


#
all_phos_STY <- extra
tt <- all_phos_STY
cc <- c("U2OS_Rap_WT_rep1","U2OS_Rap_WT_rep2","U2OS_Rap_WT_rep3","U2OS_Rap_WT_rep4","U2OS_Rap_WT_rep5","U2OS_RapConA_WT_rep1","U2OS_RapConA_WT_rep2","U2OS_RapConA_WT_rep3","U2OS_RapConA_WT_rep4","U2OS_RapConA_WT_rep5","U2OS_RapConA_Rap_rep1","U2OS_RapConA_Rap_rep2","U2OS_RapConA_Rap_rep3","U2OS_RapConA_Rap_rep4","U2OS_RapConA_Rap_rep5")

gsub("(U2OS)_(\\w+)_(\\w+)_rep(\\d{1})","Intensity.\\4_\\1_\\2_\\3" ,colnames(tt[,cc]))

colnames(tt[,cc])
colnames(tt)[colnames(tt) %in% cc] <- gsub("(U2OS)_(\\w+)_(\\w+)_rep(\\d{1})","Intensity.\\4_\\1_\\2_\\3" ,colnames(tt)[colnames(tt) %in% cc])
write.table(tt,'../Fribourg_data_pipeline.txt', quote=FALSE, sep='\t', row.names = FALSE)
