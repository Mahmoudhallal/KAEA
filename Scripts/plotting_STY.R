##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Create stacked barplots and venn diagrams of peptides and proteins
## Date: 28.03.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
#install.packages("ggplot2", repos="http://cran.rstudio.com/")
library(ggplot2)
#install.packages("VennDiagram", repos="http://cran.rstudio.com/")
library(VennDiagram)
#install.packages("reshape2", repos="http://cran.rstudio.com/")
library("reshape2")
flog.threshold(ERROR)
library(yaml)

## Load parameters file
params <- read_yaml("./config.yaml")


#Load phospho- and proteome data
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"batches_all_data_",params$cell_line,".Rda"))
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"batches_all_phospho_data_",params$cell_line,".Rda"))

## Define conditions
Conditions <- strsplit(params$Conditions,",")[[1]]
biological_reps <- strsplit(params$Biological_replicates,"-")[[1]]

## Define statistics
proteins_counts <- lapply(1:length(Conditions), function(x) {
  lapply(1:biological_reps[x], function(y) {
    a = length(unique(batches_all_phospho_data_STY[[x]][[y]]$Protein))
    b = length(unique(batches_all_phospho_data_STY[[x]][[y]]$Protein_with_position))
    c = batches_all_phospho_data_STY[[x]][[y]]$Protein_with_position
    d = batches_all_phospho_data_STY[[x]][[y]]$Protein
    e = batches_all_phospho_data_STY[[x]][[y]]$Modified.sequence
    return(list(all = a, phos = b, c = c, d = d, e = e))
  })
})


###############################################################################################
##Stacked barplot for all batches of number of phosphorylated and non-phosphorylated peptides and proteins
## Phosphorylated proteins with positions
## Define labels of replicates
all_prots_labels <- lapply(1:length(Conditions), function(x){
  lapply(1:biological_reps[x], function(y){paste0(Conditions[[x]],'_',y) })
})
all_prots_labels <- unlist(all_prots_labels)

## Get phosphorylation sites numbers
phos_sites <- lapply(1:length(Conditions), function(x){
  lapply(1:biological_reps[x], function(y){proteins_counts[[x]][[y]]$phos })
})

## Prepare data to be plotted
phos_sites <- as.data.frame(unlist(phos_sites))
colnames(phos_sites) <- "prots"
phos_sites$row <- c("Phos")
phos_sites$rep <- all_prots_labels
phos_sites$rep <- factor(phos_sites$rep, levels=unique(phos_sites$rep))
phos_sites$Cell_line <- gsub("(.*)_.*","\\1",phos_sites$rep)

dd <- ggplot(phos_sites, aes(x=rep, y=prots, fill = Cell_line)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=prots),vjust=-1, size=3) +
  xlab("Replicates") +
  ylab("# sites") +
  labs(title="Phosphorylated sites") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title=element_text(size=24))

## Produce output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Phos_prots_with_position_histogram_",params$cell_line,".pdf"),height=7.5)
plot(dd)
dev.off()

## Get phosphorylated peptides numbers
phos_peps <- lapply(1:length(Conditions), function(x){
  lapply(1:biological_reps[x], function(y){length(unique(proteins_counts[[x]][[y]]$e))})
})

phos_peps <- as.data.frame(unlist(phos_peps))
colnames(phos_peps) <- "prots"
phos_peps$row <- c("Phos")
phos_peps$rep <- all_prots_labels
phos_peps$rep <- factor(phos_peps$rep, levels=unique(phos_peps$rep))
phos_peps$Cell_line <- gsub("(.*)_.*","\\1",phos_peps$rep)

dd <- ggplot(phos_peps, aes(x=rep, y=prots, fill = Cell_line)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=prots),vjust=-1, size=3) +
  xlab("Replicates") +
  ylab("# proteins") +
  labs(title="Phosphorylated proteins") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title=element_text(size=24))

## Produce output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Phos_peps_histogram_",params$cell_line,".pdf"))
plot(dd)
dev.off()

## Get phosphorylated proteins numbers
phos_prots <- lapply(1:length(Conditions), function(x){
  lapply(1:biological_reps[x], function(y){proteins_counts[[x]][[y]]$all })
})

phos_prots <- as.data.frame(unlist(phos_prots))
colnames(phos_prots) <- "prots"
phos_prots$row <- c("Phos")
phos_prots$rep <- all_prots_labels
phos_prots$rep <- factor(phos_prots$rep, levels=unique(phos_prots$rep))
phos_prots$Cell_line <- gsub("(.*)_.*","\\1",phos_prots$rep)

dd <- ggplot(phos_prots, aes(x=rep, y=prots, fill = Cell_line)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=prots),vjust=-1, size=3) +
  xlab("Replicates") +
  ylab("# proteins") +
  labs(title="Phosphorylated proteins") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
  theme(axis.text.y = element_text(size = 22)) +
  theme(axis.title=element_text(size=24))

## Produce output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Phos_prots_histogram_",params$cell_line,".pdf"))
plot(dd)
dev.off()


###################################################################################################
## Compare overlap between conditions
## Only phosphorylated peptides 
all_peps <- lapply(1:length(Conditions), function(x){
  unique(unlist(lapply(1:biological_reps[x], function(y){
    proteins_counts[[x]][[y]]$e })))
})
names(all_peps) <- Conditions

## Produce venn diagram
dd <- venn.diagram(all_peps,
                   fill = 2:(1+length(all_peps)),
                   filename = NULL, 
                   alpha = rep(0.5,length(all_peps)),
                   cat.fontface = 4,
                   lty =2, cex = 1.1, 
                   scaled = TRUE, euler.d = TRUE,
                   main = "Phosphorylated Peptides of K562 Ctrl and Drug", 
                   main.cex = 1.2)

## Produce output PDF
pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Phos_peptides_venn_diagram_",params$cell_line,".pdf"), width = 8)
grid.draw(dd) 
dev.off()

## Proteins with positions
all_sites <- lapply(1:length(Conditions), function(x){
  unique(unlist(lapply(1:biological_reps[x], function(y){
    proteins_counts[[x]][[y]]$c })))
})
names(all_sites) <- Conditions

union_sites <- unique(unlist(all_sites))

dd <- venn.diagram(all_sites,
                   fill = 2:(1+length(all_sites)),
                   filename = NULL, 
                   alpha = rep(0.5,length(all_sites)),
                   cat.fontface = 4,
                   lty =2, 
                   cex = 1.1, 
                   scaled = TRUE, 
                   euler.d = TRUE ,
                   main = "Phosphorylated Proteins_phosphorylation site of K562 Ctrl and Drug", 
                   main.cex = 1.2)

pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Phos_prots_with_position_venn_diagram_",params$cell_line,".pdf"), width = 8)
grid.draw(dd) 
dev.off()

## Proteins with positions
all_prots <- lapply(1:length(Conditions), function(x){
  unique(unlist(lapply(1:biological_reps[x], function(y){
    proteins_counts[[x]][[y]]$d})))
})
names(all_prots) <- Conditions

union_prots <- unique(unlist(all_prots))

dd <- venn.diagram(all_prots,
                   fill = 2:(1+length(all_peps)),
                   filename = NULL, 
                   alpha = rep(0.5,length(all_peps)),
                   cat.fontface = 4,
                   lty =2, 
                   cex = 1.1, 
                   scaled = TRUE, 
                   euler.d = TRUE ,
                   main = "Phosphorylated Proteins_phosphorylation site of K562 Ctrl and Drug", 
                   main.cex = 1.2)

pdf(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Phos_prots_venn_diagram_",params$cell_line,".pdf"), width = 8)
grid.draw(dd) 
dev.off()

#########################################################
## Make a table of counts for peptides and proteins
table_counts_all <- 
  data.frame(cbind(
    Conditions,
    unlist(lapply(1:length(all_peps),function(x){length(unique(all_peps[[x]]))})),#peps
    unlist(lapply(1:length(all_prots),function(x){length(unique(all_prots[[x]]))})),#prots with phos sites
    unlist(lapply(1:length(all_prots),function(x){length(unique(gsub('_.*','',all_prots[[x]])))}))
  ))
      
colnames(table_counts_all) <- c("Cell line","# Peps", "# Prots_with_position", "# Prots")
write.table(table_counts_all, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"table_counts_all_",params$cell_line,".csv"), sep=",")

## Save for Shiny
Quality_rep <- list(phos_sites, phos_peps, phos_prots, all_sites, all_peps, all_prots, table_counts_all)
names(Quality_rep) <- c('hist_sites','hist_peps','hist_prots','venn_sites','venn_peps','venn_prots','table_counts')
save(Quality_rep, file = paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"Quality_Report_Shiny_",params$cell_line,".Rda"))
