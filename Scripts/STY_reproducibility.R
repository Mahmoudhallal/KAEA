##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Plot violin plots and venn diagrams
## Date: 28.08.2018
## Author: Mahmoud Hallal
##################################################

## Load libraries
library(yaml)

#install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)

#install.packages("VennDiagram", repos="http://cran.rstudio.com/")
library(VennDiagram)

#install.packages("reshape2", repos="http://cran.rstudio.com/")
library("reshape2")

library(Biobase)

library(plyr)

flog.threshold(ERROR)

## load parameters file
params <- read_yaml("./config.yaml")
biological_reps <- strsplit(params$Biological_replicates,"-")[[1]]
imp <- params$Imputation

## Define conditions used
Conditions <- unlist(strsplit(params$Conditions,","))
biological_reps <- strsplit(params$Biological_replicates,"-")[[1]]

###################################################################################################
## Coeffiecient of variation
## load expressionSet 
if (biological_reps[1] > 1){
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))
all_quants <- as.data.frame(exprs(test_eSet))

## Calculate CV between replicates of every conditions
## CV rule: (SD/mean)*100
cc <- lapply (1:length(Conditions), function(x){
  all_quants2 <- all_quants[,grep(Conditions[x],colnames(all_quants),fixed=TRUE)]
  all_quants2$Proteins <- rownames(all_quants)
  all_quants2 <- melt(all_quants2)
  all_quants2$cell_line <- gsub("\\d_(\\w+)",'\\1',all_quants2$variable)
  all_quants2$factor <- paste(all_quants2$Proteins, all_quants2$cell_line, sep='_') 
  
  all_quants2 <- all_quants2[,c('factor','value')]
  all_quants2 <- all_quants2[!all_quants2$value == 0,]
  mean_df <- ddply(all_quants2, .(factor), colwise(mean))
  sd_df <- ddply(all_quants2, .(factor), colwise(sd))
  final_mean <- mean_df[,c('factor','value')]
  colnames(final_mean) <- c('factor','mean')
  final_sd <- sd_df[,c('factor','value')]
  colnames(final_sd) <- c('factor','sd')
  final_cv <- cbind(final_mean, final_sd)
  final_cv$CV <- (final_cv$sd/final_cv$mean)*100
  final_cv <- final_cv[!is.na(final_cv$CV),]
  final_cv$cell_line <- Conditions[x]
  final_cv
  
})

## Plot the violin plot of every condition
cc2 <- do.call(rbind, cc)
cc2$cell_line <- as.factor(cc2$cell_line)
cc2$CV <- as.numeric(cc2$CV)
cc2 <- cc2[,c('cell_line','CV')]

#write as csv file
write.csv(cc2, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/CV_",params$cell_line,".csv"))

## ggplot for CV
g <- ggplot(cc2, aes(cell_line, CV, fill = cell_line))
g <- g + geom_violin(alpha = 0.5, draw_quantiles = c(0.5)) + 
  labs(title="CV between triplicates", 
       x="Condition",
       y="CV %") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title=element_text(size=24))

# Output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Phos_proteins_CV_violin_plot_",params$cell_line,".pdf"), width = 8)
plot(g)
dev.off()
}
