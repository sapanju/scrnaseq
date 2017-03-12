
setwd("~/Documents/e09")

library(readr)
cellnames <- read_csv("~/Documents/e09/cellnames.csv", col_names = FALSE);

raw_counts = read.table("genes.rsem_counts.txt", row.names = 1);

mydata = raw_counts;
colnames(mydata) = t(cellnames)

revisedHg19annotationFileSummary <- read.table("~/Documents/E12/Analysis/revisedHg19annotationFileSummary.txt", quote="\"")
names(revisedHg19annotationFileSummary)[1] = "EnsembleId"
names(revisedHg19annotationFileSummary)[2] = "GeneType"
names(revisedHg19annotationFileSummary)[3] = "GeneName"

# Add column for gene_type and gene name, rename accordingly
mydata = data.frame(mydata[,1], mydata[,1], mydata[,1], mydata)
names(mydata)[1] = "EnsembleId"
names(mydata)[2] = "GeneName"
names(mydata)[3] = "GeneType"
mydata$EnsembleId = rownames(mydata)
mydata$GeneName = revisedHg19annotationFileSummary$GeneName[match(mydata$EnsembleId, revisedHg19annotationFileSummary$EnsembleId)]
mydata$GeneType = revisedHg19annotationFileSummary$GeneType[match(mydata$EnsembleId, revisedHg19annotationFileSummary$EnsembleId)]

# number of duplicates: length(which(duplicated(mydata_coding_only$GeneName)))
# = 103
# Nov 10 (e12), Mar 07 (e09)

# Remove all but protein coding genes
mydata_coding_only = mydata[which(mydata$GeneType == "protein_coding"),]

# Add column for entrezid
mydata_coding_only = data.frame(mydata_coding_only[,1], mydata_coding_only)
names(mydata_coding_only)[1] = "EnsembleId"
names(mydata_coding_only)[2] = "EntrezId"

source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

# call remove duplicates script



# # Match entrez id by gene name (1466 NA's)
# xx <- as.list(org.Hs.egSYMBOL2EG)
# xx2 <- cbind("symbol"=names(xx), "entrez"=sapply(xx, function (x) { return (x[1])} ))
# mydata_coding_only$EntrezId = xx2[,2][match(mydata_coding_only$GeneName, xx2[,1])]
# # number of duplicated entrez id's: 101
# # number of duplicated gene names: 101
# 
# # Match entrez id by ensemble id (1431 NA's)
# # Remove decimal point values from Ensemble ids first
# ensembleids <- read.table("~/Documents/E12/Analysis/ensembleids.txt", quote="\"", comment.char="")
# View(ensembleids)
# mydata_coding_only[,1] = as.list(ensembleids)
# xx <- as.list(org.Hs.egENSEMBL2EG)
# xx2 <- cbind("symbol"=names(xx), "entrez"=sapply(xx, function (x) { return (x[1])} ))
# mydata_coding_only$EntrezId = xx2[,2][match(mydata_coding_only$EnsembleId, xx2[,1])]
# 
# # Remove NA's
# mydata_coding_only_noNA = mydata_coding_only[!is.na(mydata_coding_only$EntrezId),]
# 
# # number of duplicated entrez id's: 64
# # number of duplicated gene names: 11