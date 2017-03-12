#QC

setwd("~/Documents/e09/qc")
load("~/Documents/e09/rawdata/e09_raw_counts_eset.RData")

mydata = e09_raw_counts_eset
mydata = apply(mydata, 2, as.numeric)
rownames(mydata) = rownames(e09_raw_counts_eset)

sce <- newSCESet(countData=as.matrix(mydata))
dim(sce)

is.spike <- grepl("^ERCC", rownames(sce))
is.mito <- grepl("^MT-", rownames(sce))

sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
head(colnames(pData(sce)))

## [1] "total_counts"             "log10_total_counts"       "filter_on_total_counts"  
## [4] "total_features"           "log10_total_features"     "filter_on_total_features"
library(scran)
isSpike(sce) <- is.spike

par(mfrow=c(1,2))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)

par(mfrow=c(1,2))
hist(sce$pct_exprs_feature_controls_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_exprs_feature_controls_ERCC, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

mito.drop <- isOutlier(sce$pct_exprs_feature_controls_Mt, nmads=3, type="higher")
spike.drop <- isOutlier(sce$pct_exprs_feature_controls_ERCC, nmads=3, type="higher")

sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce))

# ByLibSize ByFeature ByMito BySpike Remaining
# Samples         1         3      4       1        32


e09_clean_counts_eset = counts(sce)
genesums = rowSums(e09_clean_counts_eset)
e09_clean_counts_eset = e09_clean_counts_eset[-(which(genesums < 50)),]



# USING TPM DATA:
load(e09_raw_tpm_eset)
exprs = log2(e09_raw_tpm_eset + 1)
sce = newSCESet(exprsData = as.matrix(exprs), tpmData= as.matrix(e09_raw_tpm_eset))
is_exprs(sce) <- exprs > 0.1
## ByLibSize ByFeature ByMito BySpike Remaining
## Samples         0         0      5       5        70

ave.tpm <- rowMeans(tpm(sce))
keep <- ave.tpm >= 1
sum(keep)

e09_clean_tpm_eset = tpm(sce)
e09_clean_tpm_eset = e09_clean_tpm_eset[keep,]

## [1] 9398
# hist(log10(ave.tpm), breaks=100, main="", col="grey80",
#      xlab=expression(Log[10]~"average count"))
# abline(v=log10(1), col="blue", lwd=2, lty=2)

# REALLY SHOULD HAVE DONE THRESHOLD OF 0.5, NOT 1
#ave.counts <- rowMeans(counts(sce))
#keep <- ave.counts >= 0.5
#sum(keep)
#hist(log10(ave.counts), breaks=100, main="", col="grey80",
#     xlab=expression(Log[10]~"average count"))
#abline(v=log10(0.5), col="blue", lwd=2, lty=2)


#TOP 1000 most variable genes based on IQR
matrix = data.frame(e09_clean_counts_eset)
matrix[,33] = rownames(e09_clean_counts_eset)
rownames(matrix) = 1:10406
iqr_list = apply(matrix[,1:32], 1, IQR)
iqr_list = sort(iqr_list, TRUE)
indices_to_remove = as.numeric(names(iqr_list[1001:length(iqr_list)]))
matrix = matrix[-indices_to_remove,]
top_1000_genes_eset = as.data.frame(matrix)
rownames(top_1000_genes_eset) = top_1000_genes_eset[,33]
top_1000_genes_eset = top_1000_genes_eset[,1:32]

# ## repeated for the tpm eset
matrix = e09_clean_tpm_eset
matrix[,33] = rownames(e09_clean_tpm_eset)
rownames(matrix) = 1:9478
iqr_list = apply(matrix[,1:32], 1, IQR)
iqr_list = sort(iqr_list, TRUE)
indices_to_remove = as.numeric(names(iqr_list[1001:length(iqr_list)]))
matrix = matrix[-indices_to_remove,]
top_1000_genes_eset = as.data.frame(matrix)
rownames(top_1000_genes_eset) = top_1000_genes_eset[,33]
top_1000_genes_eset = top_1000_genes_eset[,1:32]

library(pvclust)
pvcl <- pvclust(data=as.matrix(top_1000_genes_eset), method.hclust="complete", method.dist="correlation",use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)
library(gplots)

heatmap.2(as.matrix(top_1000_genes_eset),Colv = as.dendrogram(pvcl$hclust),dendrogram = "column", scale = "row", trace = "none", cexRow = 0.4, cexCol = 0.4)
