setwd("~/Documents/E12/QC")
real_tpm_eset <- read.delim("~/Documents/E12/QC/real_tpm_eset.txt")

install.packages("elasticnet")
library(elasticnet)

log2_real_tpm_eset = log2(real_tpm_eset + 1.0)




spc1 <- elasticnet::arrayspc(x=t(as.matrix(log2_real_tpm_eset)), K=1, para=1e3)
myx <- abs(spc1$loadings) > 1e-2


orderedmyx <- order(abs(spc1$loadings), decreasing=TRUE)[1:1000]
top_1000_genes <- rownames(log2_real_tpm_eset)[orderedmyx]

library(pvclust)

top_1000_genes_eset_elasticnet = log2_real_tpm_eset[match(top_1000_genes, rownames(log2_real_tpm_eset)),]

pvcl <- pvclust(data=top_1000_genes_eset_elasticnet, method.hclust="average", method.dist="correlation",use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)

library(gplots)

#write.table(top_1000_genes, "top_1000_genes_elasticnet.txt", sep = "\t")
heatmap.2(as.matrix(top_1000_genes_eset_elasticnet),Colv = as.dendrogram(pvcl$hclust),dendrogram = "column", scale = "row", trace = "none")
