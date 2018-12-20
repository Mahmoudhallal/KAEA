#source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library(Biobase)

print("PHOSPHO_DBS_WITH_MOTIF_WILL_BE_USED!!")

all_dbs2 <- read.csv("/Users/Mahmoud.Hallal/Desktop/PhD/Generalized_pipeline/results/all_dbs1.csv")
#motif_db <- read.csv('/Users/Mahmoud.Hallal/Desktop/PhD/Motif_analysis/all_dbs.csv')
motif_db <- data.frame()

all_dbs2 <- all_dbs2[,-1]
#motif_db <- motif_db[,-1]

all_dbs <- rbind(all_dbs2, motif_db)
#all_dbs <- all_dbs2
all_dbs <- all_dbs[!duplicated(all_dbs),]

write.csv(all_dbs, "/Users/Mahmoud.Hallal/Desktop/PhD/Generalized_pipeline/results/all_dbs.csv")
#
