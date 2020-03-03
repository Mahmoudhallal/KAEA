
library(pathview)
library(filesstrings)

##
to_sbmt <- lapply(1:length(topTables_greater), function(x){
  df1 <- topTables_greater[[x]]
  #make dataframe
  df <- as.data.frame(df1$log2FoldChange,drop=FALSE)
  rownames(df) <- rownames(df1)
  df <- df[df$`df1$log2FoldChange` >1,,drop=FALSE]
  colnames(df) <- "value"
  #dd <- df$value
  #names(dd) <- gsub('(\\w+)[+-]','\\1',rownames(df))
  #dd1 <- dd[!duplicated(names(dd))]
  dd1 <- df
  dd2 <- as.data.frame(dd1)
  #rownames(dd2) <- names(dd1)

  
  dd2
  
})
write.csv(to_sbmt[[1]],'../../submit_data.csv')
## Define directory and KEGG pathways
directory <- '../../Generalized_pipeline/result_pathview/'
unlink(directory,recursive = TRUE)
dir.create(directory)

pathways = c("01100")
pathways_names <- c('Metabolic pathways')

## Create all pathview.png for all pathways and conditions 
files <- c('Molm13R')
all_files <- lapply(1:length(files), function(y){
  sp_directory <- paste0(directory,'/','test1')
  dir.create(sp_directory)
  ## assign the species according to input, this is for KEGG mapping
  sp = "hsa"
  type = "SYMBOL"
    
  
  table_con <-  lapply(1:length(pathways), function(x) {
    pp <- pathview(gene.data = to_sbmt[[1]], pathway.id = pathways[x],
                   #we changed gene.idtype from SYMBOL to GENENAME 22.01.2019
                   species = sp, gene.idtype = type, 
                   limit = list(gene = 0.05, cpd = 0.05),
                   kegg.dir = sp_directory)
    file.move(paste('hsa',pathways[x], ".pathview.png", sep = ""), sp_directory, overwrite = TRUE)
    #file.copy(paste0(sp_directory,'/',pathways[x],'.pathview.png'), paste0(sp_directory,'/',pathways_names[x],'.pathview.png'))
    pp
  })
})

## Give the file names
names(all_files) <- files
