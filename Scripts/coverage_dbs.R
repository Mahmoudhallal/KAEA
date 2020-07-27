##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Prepare output file for Shiny app
## Date: 21.09.2018
## Author: Mahmoud Hallal
##################################################

## Load libraries
library(reshape2)

library(ggplot2)

library(yaml)

library(Biobase)

library(dplyr)

## Load parameters
params <- read_yaml("./config.yaml")
imp <- params$Imputation

## Load data
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))
all_dbs <- read.csv(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/all_dbs.csv"))

## Assign proteins to samples
sample_sites <- lapply(1:length(colnames(exprs(test_eSet))), function(x){
  data <- as.data.frame(exprs(test_eSet))
  data_sample <- as.data.frame(data[,colnames(data)[x], drop=F])
  #turn 0s into NA and remove them
  data_sample[data_sample == 0] <- NA
  proteins <- rownames(data_sample[!is.na(data_sample),,drop=F])
  unique(proteins)
})

## yes = covered, no = not covered
total_cov_prots <- lapply(1:length(sample_sites), function(x) length(sample_sites[[x]]))
coverage_yes <- lapply(1:length(sample_sites), function(x) sum(sample_sites[[x]] %in% unique(all_dbs$geneID)))
coverage_no <- lapply(1:length(sample_sites), function(x) sum(!sample_sites[[x]] %in% unique(all_dbs$geneID)))

coverage <- lapply(1:length(sample_sites), function(x){
  cc <- c(coverage_no[[x]],coverage_yes[[x]])
  names(cc) <- c('no','yes')
  cc
  })

names(coverage) <- colnames(exprs(test_eSet))
#remove t() only for docker version
coverage2 <- bind_rows(coverage)
colnames(coverage2) <- c('no','yes')
#added only for docker version
rownames(coverage2) <- names(coverage)
coverage2 <- as.data.frame(coverage2)
coverage2$sample <- rownames(coverage2)


colnames(coverage2) <- c("Absent", "Present","cell_line")
coverage2$cell_line <- as.factor(coverage2$cell_line)
Total <- as.numeric(coverage2$Absent) + as.numeric(coverage2$Present)
(as.numeric(coverage2$Present)*100)/Total

coverage2$Present <- as.numeric(coverage2$Present)
coverage2$Absent <- as.numeric(coverage2$Absent)

coverage_df_new <- melt(coverage2)

colnames(coverage_df_new) <- c("cell_line", "Status", "value")

## Create barplot
dd <- ggplot(data = coverage_df_new, aes(x = cell_line, y = value, fill= Status))+
  geom_bar(stat="identity")+
  scale_x_discrete(name="Samples") +
  ggtitle("Number of present/absent phosphorylation sites in the reference DB") +
  geom_text(aes(label=as.numeric(value)),vjust=2, size = 3) +
  ylab("#Phosphorylated sites") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) 
  
## Create output pdf
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Coverage_Database_",params$cell_line,".pdf"), width = 8)
plot(dd)
dev.off()

