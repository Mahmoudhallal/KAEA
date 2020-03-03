##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Plot violin plots and venn diagrams
## Date: 28.08.2018
## Author: Mahmoud Hallal
##################################################

## Load Packages
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

## Take every batch separately and make some statistics of proteome and phosphorpteoome
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/batches_all_data_",params$cell_line,".Rda"))
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/batches_all_phospho_data_",params$cell_line,".Rda"))

## Define conditions used
Conditions <- unlist(strsplit(params$Conditions,","))
biological_reps <- strsplit(params$Biological_replicates,"-")[[1]]

## GIve statistics on phosphopeptides and proteins
peptide_counts <- lapply(1:length(Conditions), function(x) {
  lapply(1:biological_reps[x], function(y) {
    a = length(unique(batches_all_data_STY[[x]][[y]]$Modified.sequence.position))
    b = length(unique(batches_all_phospho_data_STY[[x]][[y]]$Modified.sequence.position))
    diff = a-b
    c = batches_all_phospho_data_STY[[x]][[y]]$Modified.sequence.position
    d = batches_all_phospho_data_STY[[x]][[y]]$Protein
    return(list(all = a, phos = b, diff = diff, c = c, d = d))
  })
})
#a = numberof all peptides
#b = number of phosphorylated peptides
#diff = difference between the 2 values above
#c = the modified sequences of this replicate cell line
#d = the proteins of this replicate cell line

## Proteins
proteins_counts <- lapply(1:length(Conditions), function(x) {
  lapply(1:biological_reps[x], function(y) {
    a = length(unique(batches_all_data_STY[[x]][[y]]$Protein))
    b = length(unique(batches_all_phospho_data_STY[[x]][[y]]$Protein))
    diff = a-b
    c = batches_all_phospho_data_STY[[x]][[y]]$Modified.sequence.position
    d = batches_all_phospho_data_STY[[x]][[y]]$Protein
    e = batches_all_phospho_data_STY[[x]][[y]]$Protein_with_position
    
    return(list(all = a, phos = b, diff = diff, c = c, d = d, e = e))
  })
})


###################################################################################################
#Draw venn diagrams of overlaps between replicates of every condition
# for(cond in 1:length(Conditions)){
#   pr <- c(biological_reps)
#   rr <- lapply(1:pr[cond],function(x) proteins_counts[[cond]][[x]]$e)
#   names(rr) <- lapply(1:pr[cond], function(x) paste0(Conditions[cond],'_',x))
#   #K562
#   venn_diag <- venn.diagram(rr,
#     fill = 2:(as.numeric(pr[cond])+1),filename = NULL, alpha = rep(0.4,pr[cond]),
#     cat.fontface = 4,lty =2, cat.cex = 1.6, cex = 2.2, scaled = TRUE, euler.d = TRUE,
#     main = "Phosphorylated sites between K562 CTRL replicates", main.cex = 2)
#   
#   pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Phos_proteins_",params$cell_line,"_venn_diagram_", Conditions[cond],".pdf"), width = 8)
#   grid.draw(venn_diag) 
#   dev.off()
#   
# }

###################################################################################################
## Coeffiecient of variation
## load expressionSet 
if (biological_reps[1] > 1){
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/test_eSet1_",params$cell_line,".Rda"))
all_quants <- as.data.frame(exprs(test_eSet))

## Calculate CV between replicates of every conditions
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
  #plot(density(final_cv$CV))
  
})

## Plot the violin plot of every condition
cc2 <- do.call(rbind, cc)
cc2$cell_line <- as.factor(cc2$cell_line)
cc2$CV <- as.numeric(cc2$CV)
cc2 <- cc2[,c('cell_line','CV')]

write.csv(cc2, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/CV_",params$cell_line,".csv"))

g <- ggplot(cc2, aes(cell_line, CV, fill = cell_line))
g <- g + geom_violin(alpha = 0.5,draw_quantiles = c(0.5)) + 
  labs(title="CV between triplicates", 
       x="Condition",
       y="CV %") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(#axis.text=element_text(size=14),
    axis.title=element_text(size=24))


pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/Phos_proteins_CV_violin_plot_",params$cell_line,".pdf"), width = 8)
plot(g)
dev.off()
}
