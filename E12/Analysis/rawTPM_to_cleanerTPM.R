setwd("~/Documents/E12/Analysis")

# Load original readcounts matrix (tpm) and the annotation file
E12_original <- read.delim("~/Documents/E12/Analysis/E12_original.txt")
revisedHg19annotationFileSummary <- read.table("~/Documents/E12/Analysis/revisedHg19annotationFileSummary.txt", quote="\"")
names(revisedHg19annotationFileSummary)[1] = "EnsembleId"
names(revisedHg19annotationFileSummary)[2] = "GeneType"
names(revisedHg19annotationFileSummary)[3] = "GeneName"

# Add column for gene_type and gene name, rename accordingly
E12 = data.frame(E12_original[,1], E12_original[,1], E12_original)
names(E12)[1] = "EnsembleId"
names(E12)[2] = "GeneName"
names(E12)[3] = "GeneType"
E12$GeneName = revisedHg19annotationFileSummary$GeneName[match(E12$EnsembleId, revisedHg19annotationFileSummary$EnsembleId)]
E12$GeneType = revisedHg19annotationFileSummary$GeneType[match(E12$EnsembleId, revisedHg19annotationFileSummary$EnsembleId)]

# number of duplicates: length(which(duplicated(E12_coding_only$GeneName)))
# = 103
# Nov 10

# Remove all but protein coding genes
E12_coding_only = E12[which(E12$GeneType == "protein_coding"),]

# Add column for entrezid
E12_coding_only = data.frame(E12_coding_only[,1], E12_coding_only)
names(E12_coding_only)[1] = "EnsembleId"
names(E12_coding_only)[2] = "EntrezId"

source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

# Match entrez id by gene name (1466 NA's)
xx <- as.list(org.Hs.egSYMBOL2EG)
xx2 <- cbind("symbol"=names(xx), "entrez"=sapply(xx, function (x) { return (x[1])} ))
E12_coding_only$EntrezId = xx2[,2][match(E12_coding_only$GeneName, xx2[,1])]
# number of duplicated entrez id's: 101
# number of duplicated gene names: 101

# Match entrez id by ensemble id (1431 NA's)
# Remove decimal point values from Ensemble ids first
ensembleids <- read.table("~/Documents/E12/Analysis/ensembleids.txt", quote="\"", comment.char="")
View(ensembleids)
E12_coding_only[,1] = as.list(ensembleids)
xx <- as.list(org.Hs.egENSEMBL2EG)
xx2 <- cbind("symbol"=names(xx), "entrez"=sapply(xx, function (x) { return (x[1])} ))
E12_coding_only$EntrezId = xx2[,2][match(E12_coding_only$EnsembleId, xx2[,1])]

# Remove NA's
E12_coding_only_noNA = E12_coding_only[!is.na(E12_coding_only$EntrezId),]

# number of duplicated entrez id's: 64
# number of duplicated gene names: 11