#QC
#load E12_eset_counts
E12_eset_counts = read.table("~/Documents/E12/QC/E12_eset_counts.txt", sep = "\t", header = TRUE)

library(scater)
sce <- newSCESet(countData=E12_eset_counts)
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

## ByLibSize ByFeature ByMito BySpike Remaining
## Samples         1         1      5       5         74

# USING TPM DATA:
# Load tpm_eset.txt
exprs = log2(real_tpm_eset + 1)
sce = newSCESet(exprsData = as.matrix(exprs), tpmData= as.matrix(real_tpm_eset))
is_exprs(sce) <- exprs > 0.1
## ByLibSize ByFeature ByMito BySpike Remaining
## Samples         0         0      5       5        70

ave.tpm <- rowMeans(tpm(sce))
keep <- ave.tpm >= 1
sum(keep)
## [1] 11494
hist(log10(ave.tpm), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)

# REALLY SHOULD HAVE DONE THRESHOLD OF 0.5, NOT 1
#ave.counts <- rowMeans(counts(sce))
#keep <- ave.counts >= 0.5
#sum(keep)
#hist(log10(ave.counts), breaks=100, main="", col="grey80",
#     xlab=expression(Log[10]~"average count"))
#abline(v=log10(0.5), col="blue", lwd=2, lty=2)


#TOP 1000 most variable genes based on IQR
matrix = clean_counts_eset
matrix[,75] = rownames(clean_counts_eset)
rownames(matrix) = 1:11494
iqr_list = apply(matrix[,1:74], 1, IQR)
iqr_list = sort(iqr_list, TRUE)
indices_to_remove = as.numeric(names(iqr_list[1001:length(iqr_list)]))
matrix = matrix[-indices_to_remove,]
top_1000_genes_eset = as.data.frame(matrix)
rownames(top_1000_genes_eset) = top_1000_genes_eset[,75]
top_1000_genes_eset = top_1000_genes_eset[,1:74]

## repeated for the tpm eset
matrix = log_transformed_real_tpm_eset
matrix[,78] = rownames(log_transformed_real_tpm_eset)
rownames(matrix) = 1:20242
iqr_list = apply(matrix[,1:77], 1, IQR)
iqr_list = sort(iqr_list, TRUE)
indices_to_remove = as.numeric(names(iqr_list[1001:length(iqr_list)]))
matrix = matrix[-indices_to_remove,]
top_1000_genes_eset = as.data.frame(matrix)
rownames(top_1000_genes_eset) = top_1000_genes_eset[,78]
top_1000_genes_eset = top_1000_genes_eset[,1:77]


pvcl <- pvclust(data=as.matrix(top_1000_genes_eset), method.hclust="complete", method.dist="correlation",use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)
library(gplots)

heatmap.2(as.matrix(top_1000_genes_eset),Colv = as.dendrogram(pvcl$hclust),dendrogram = "column", scale = "row", trace = "none")
