source("https://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")

load("~/PSets/GDSC1000.RData")
setwd("~/Documents/E12/drugprediction")

# load("GDSCexpression.RData")
# load("GDSCexpression_breast.RData")
# load("GDSC1000_auc.RData")
# load("GDSC1000_breast_auc.RData")

top_1000_genes_eset <- read.delim("~/Documents/E12/QC/top_1000_genes_eset.txt")

tissueDescriptor <- cellInfo(GDSC1000)
breastCells = tissueDescriptor[which(tissueDescriptor$GDSC.Tissue.descriptor.1 == 'breast'),"Sample.Name"]

GDSC1000.auc <- summarizeSensitivityProfiles(
  pSet=GDSC1000,
  sensitivity.measure='auc_recomputed',
  cell.lines = cellNames(GDSC1000),
  summary.stat="median",
  verbose=TRUE)

GDSC1000_breast.auc <- summarizeSensitivityProfiles(
  pSet=GDSC1000,
  sensitivity.measure='auc_recomputed',
  cell.lines = breastCells,
  summary.stat="median",
  verbose=TRUE)

# GDSC1000.ic50 <- summarizeSensitivityProfiles(
#   pSet=common$GDSC,
#   sensitivity.measure='ic50_recomputed',
#   summary.stat="median",
#   verbose=FALSE)

# Convert back to ensembleid .... grrr.... and add _at to match with the format of ensembleid names in GDSC1000
annotation = read.table("~/Documents/E12/Analysis/revisedHg19annotationFileSummary.txt", quote="\"")
annotation$V1 = strtrim(annotation$V1, 15) # trim off the decimal places in the ensemble id name
names(annotation)[1] = "EnsembleId"
names(annotation)[2] = "GeneType"
names(annotation)[3] = "GeneName"
gene_names = rownames(top_1000_genes_eset)
feature_names = annotation$EnsembleId[match(gene_names, annotation$GeneName)]
feature_names = paste(feature_names, "_at", sep = "")

GDSCexpression <- summarizeMolecularProfiles(GDSC1000,
                                             cellNames(GDSC1000),
                                             mDataType="rna",
                                             features = feature_names,
                                             verbose=TRUE)
# Only 898/1000 features can be found
GDSCexpression_breast <- summarizeMolecularProfiles(GDSC1000,
                                             breastCells,
                                             mDataType="rna",
                                             features = feature_names,
                                             verbose=TRUE)

