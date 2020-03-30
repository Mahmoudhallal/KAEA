##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Use MaxQuant STY files to arrange phosphorylated proteins and peptides
## Date: 13.06.17
## Author: Mahmoud Hallal
##################################################

## Load libraries
library("reshape2")
install.packages("yaml", repos="http://cran.rstudio.com/")
library(yaml)
library(S4Vectors)

## Load parameters file
params <- read_yaml("./config.yaml")
biological_reps <- unlist(strsplit(params$Biological_replicates, "-"))
imp <- params$Imputation

## Load the input file
phos_STY_INC1 <- read.delim(params$path_Input)#,sep = ','

## Change contaminant name
if ("Contaminant" %in% colnames(phos_STY_INC1)){
  colnames(phos_STY_INC1)[colnames(phos_STY_INC1) == 'Contaminant'] <- 'Potential.contaminant'
}

## Choose the column names to use#
Intensity_columns <- params$Intensity_columns
to_keep_STY <- c("Protein", "Amino.acid", "Position","Sequence.window" ,"Phospho..STY..Probabilities", "Reverse", "Potential.contaminant", "Localization.prob",Intensity_columns )

## define cell lines
cell_lines <- params$cell_line
Conditions <- strsplit(params$Conditions,",")[[1]]

## Select the columns#
files_to_use_STY <- c("phos_STY_INC1")
for (x in 1:length(files_to_use_STY)) {assign(paste0(files_to_use_STY[x],"_new"), get(files_to_use_STY[x])[,to_keep_STY])}

## Merge inputs if available and filter <0.75, contaminants and Reverse sites
all_phos_STY <- rbind(phos_STY_INC1_new)
all_phos_STY <- all_phos_STY[all_phos_STY$Localization.prob >= 0.75,]
all_phos_STY <- all_phos_STY[!all_phos_STY$Reverse == "+",]
all_phos_STY <- all_phos_STY[!all_phos_STY$Potential.contaminant == "+",]
all_phos_STY <- all_phos_STY[!is.na(all_phos_STY$Position),]

## Sum the injection intensity values for every replicate separately
if (params$Sum_conditions == "T"){
  for (x in 1:length(params$Conditions_to_sum)){
    new_col <- params$Conditions_to_sum[[x]][1]
    all_phos_STY[[new_col]] <- rowSums(all_phos_STY[,c(params$Conditions_to_sum[[x]][-1])])
  }
}

## clear all columns that are from the injections
if (isEmpty(grep('_INC', colnames(all_phos_STY))) == F){
  all_phos_STY <- all_phos_STY[,-grep('_INC', colnames(all_phos_STY))]
}

## Separate the X batches of biological replicates
## format of Intensity column name should be "Intensity.#replicate_Conditionname"
all_batches_STY <- lapply(1:length(Conditions), function(x){
  nested_ll <- lapply(1:biological_reps[x], function(y){
    cbind(all_phos_STY[,-grep(paste(cell_lines,collapse="|"), colnames(all_phos_STY))],
          setNames(data.frame(all_phos_STY[[paste0("Intensity.",y,'_',Conditions[x])]], 
                              rep(Conditions[x], nrow(all_phos_STY))), c("Intensity", "cell_line")),
          setNames(data.frame(paste0(y,'_',Conditions[x])),c("Exp_short")))
    
  })
  new_file <- paste0("batch_new_STY_",x)
  assign(paste0("batch_new_STY_",x),do.call(rbind,nested_ll))
  get(paste0("batch_new_STY_",x))
})


################################################################################################
## Clean the Protein and sequences accordint to what we need
new_all_batches_STY <- lapply(1:length(all_batches_STY), function(x){
  
  #get gene name
  all_batches_STY[[x]]$Gene <- gsub("(^sp\\|)(.*)\\|(\\w+)(_HUMAN)","\\3", all_batches_STY[[x]]$Protein)
  
  #trim accession number
  all_batches_STY[[x]]$Protein <- gsub("(^sp\\|)?","", all_batches_STY[[x]]$Protein)
  all_batches_STY[[x]]$Protein <- gsub("[|].*","", all_batches_STY[[x]]$Protein)
  
   
  #with or without isoforms
  #all_batches_STY[[x]]$Protein <- gsub("-[0-9]{1,2}","",all_batches_STY[[x]]$Protein)
  
  #add phos site and position
  all_batches_STY[[x]]$Protein_with_position <- paste0(all_batches_STY[[x]]$Protein, "_",all_batches_STY[[x]]$Amino.acid, all_batches_STY[[x]]$Position)
  
  #remove PTMs except phosphorylations
  all_batches_STY[[x]]$Modified.sequence <- gsub("\\(\\d+(\\.\\d{1,3})?\\)","", all_batches_STY[[x]]$Phospho..STY..Probabilities)
  all_batches_STY[[x]]$Modified.sequence.position <- paste0(all_batches_STY[[x]]$Modified.sequence, "_",all_batches_STY[[x]]$Amino.acid, all_batches_STY[[x]]$Position)
  
  #all_batches_STY[[x]]$Protein_with_position <- all_batches_STY[[x]]$Sequence.window 
  
  all_batches_STY[[x]] <- all_batches_STY[[x]][!is.na(all_batches_STY[[x]]$Intensity),]
  
  all_batches_STY[[x]]
})

## Divide the batchesin a separate variable each
for(x in 1:length(all_batches_STY)) {assign(paste0('batch_STY_',x), new_all_batches_STY[[x]])}

## Export only phosphoproteome for qualitative analysis
batches_all_phospho_data_STY <- lapply(1:length(Conditions), function(x){
  cell_lines_all_data_STY <- lapply(1:biological_reps[x], function(y) {
    new <- new_all_batches_STY[[x]][new_all_batches_STY[[x]]$Exp_short == paste0(y,'_',Conditions[x]),]
    new <- new[new$Intensity != 0,]
    
  })
})
save(batches_all_phospho_data_STY, file=paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/batches_all_phospho_data_",params$cell_line,".Rda"))


###############################################################################################
## Prepare to create eSet object
## Merge X batches
all_batches_eSet_only_phos_STY <- do.call(rbind, lapply(1:length(all_batches_STY), function(x){get(paste0('batch_STY_',x))}))

## Choose columns of interest
eSet_data1_STY <- all_batches_eSet_only_phos_STY[,c("Protein_with_position","Exp_short", "Intensity")]

## remove zero intensity values
eSet_data1_STY <- eSet_data1_STY[!eSet_data1_STY$Intensity == 0,]

## Sum the intensities according to protein_position
eSet_data_STY <- dcast(eSet_data1_STY, Exp_short ~ Protein_with_position, fun.aggregate = sum)
rownames(eSet_data_STY) <- eSet_data_STY[,1]

## Inverse the data frame 
eSet_data_STY <- t(eSet_data_STY)
eSet_data_STY <- data.frame(eSet_data_STY[-1,])
eSet_data_STY$Proteins <- rownames(eSet_data_STY)

## Write the table in csv format
write.table(eSet_data_STY, paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR_imp',imp,"/eSet_data_test1_",params$cell_line,".csv"),sep="\t",row.names = FALSE,quote = FALSE)
