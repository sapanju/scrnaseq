load("~/Downloads/GeneSets-IPA-FinalVersion-EntID.RData")

# Get my gene expression matrix into proper eset format
nbt.data=E12_gsva_eset <- read.delim("~/Documents/E12/Analysis/E12_gsva_eset.txt")
E12_gsva_eset = E12_gsva_eset[,-3]
E12_gsva_eset = E12_gsva_eset[,-1]
E12_gsva_eset = E12_gsva_eset[,-2]
rownames(E12_gsva_eset) = E12_gsva_eset[,1]
E12_gsva_eset = E12_gsva_eset[,-1]
E12_gsva_eset = as.matrix(E12_gsva_eset)

# Get gene sets into proper format to use with GSVA
geneSets = list()
pathways = gSets_IPA_EntID$IPAPathways
pathways = pathways[!duplicated(pathways)]

for (pathway in pathways)
{
  geneSets = append(geneSets, list(pathway=gSets_IPA_EntID[which(gSets_IPA_EntID$IPAPathways == pathway),1]))
  #print(pathway)
  #print(gSets_IPA_EntID[which(gSets_IPA_EntID$IPAPathways == pathway),1])
}

names(geneSets) = pathways

#biocLite("GSVA")
#library(GSVA)
#library(pvclust)

output = gsva(E12_gsva_eset, geneSets, rnaseq=TRUE)$es.obs

# Generate Heatmap
pvcl <- pvclust(data=output, method.hclust="complete", method.dist="correlation",use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)
library(gplots)

heatmap.2(as.matrix(gsvaOutput),Colv = as.dendrogram(pvcl$hclust),dendrogram = "column", scale = "row", trace = "none")

