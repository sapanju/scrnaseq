library(Seurat)
library(dplyr)
library(Matrix)

real_tpm_eset <- read.delim("~/Documents/E12/QC/real_tpm_eset.txt", sep = "\t")
exprs = log2(real_tpm_eset + 1)

colnames(exprs)[42] = "B_LBor2B_C06.a026"
colnames(exprs)[43] = "B_LBor2B_C76.a026"
colnames(exprs)[44] = "B_LBor2B_C96.a026"
colnames(exprs)[38] = "C_EGFRn_TC_RNEASY_1OF2.a026"
colnames(exprs)[39] = "C_EGFRn_TC_RNEASY_2OF2.a026"
colnames(exprs)[40] = "C_EGFRp_TC_RNEASY_1OF2.a026"
colnames(exprs)[41] = "C_EGFRp_TC_RNEASY_2OF2.a026"
colnames(exprs)[49] = "C_TC_200_2OF2.a026" 
colnames(exprs)[50] = "C_TC_RNEASY_1OF2.a026"
colnames(exprs)[51] = "C_TC_RNEASY_2OF2.a026"
colnames(exprs)[45] = "C_SMARTER_NTC.a026"
colnames(exprs)[46] = "C_TC_1000_1OF2.a026"
colnames(exprs)[47] = "C_TC_1000_2OF2.a026" 
colnames(exprs)[48] = "C_TC_200_1OF2.a026"
                 

breast <- new("seurat", raw.data = exprs)
breast <- Setup(breast, do.logNormalize = F, project = "seurat_trying_again")

mito.genes <- grep("^MT-", rownames(breast@data), value = T)
percent.mito <- colSums(expm1(breast@data[mito.genes, ]))/colSums(expm1(breast@data))

ercc.genes <- grep("^ERCC", rownames(breast@data), value = T)
percent.ercc <- colSums(expm1(breast@data[ercc.genes, ]))/colSums(expm1(breast@data))

breast <- AddMetaData(breast, percent.mito, "percent.mito")
breast <- AddMetaData(breast, percent.ercc, "percent.ercc")
VlnPlot(breast, c("nGene", "percent.ercc", "percent.mito"), nCol = 3)


# Remove outlier cells
clean_exprs = exprs[,-82]
clean_exprs = clean_exprs[,-79]
clean_exprs = clean_exprs[,-53]
clean_exprs = clean_exprs[,-52]
clean_exprs = clean_exprs[,-45]
clean_exprs = clean_exprs[,-73]
clean_exprs = clean_exprs[,-54]
clean_exprs = clean_exprs[,-76]
# WB_C75-a026 82
# WB_C68-a026 79
# WB_C04-a026 53
# WB_C03-a026 52
# SMARTER_NTC-a026 45
# WB_C63-a026 73
# WB_C16-a026 54
# WB_C77.a026 76

ave.tpm <- rowMeans(clean_exprs)
keep <- ave.tpm >= 1
sum(keep)
# 10547
clean_exprs = clean_exprs[keep,]

clean_breast <- new("seurat", raw.data = clean_exprs)
clean_breast <- Setup(clean_breast, do.logNormalize = T, min.cells = 1, project = "seurat_trying_again")

mito.genes <- grep("^MT-", rownames(clean_breast@data), value = T)
percent.mito <- colSums(expm1(clean_breast@data[mito.genes, ]))/colSums(expm1(clean_breast@data))

ercc.genes <- grep("^ERCC", rownames(clean_breast@data), value = T)
percent.ercc <- colSums(expm1(clean_breast@data[ercc.genes, ]))/colSums(expm1(clean_breast@data))

clean_breast <- AddMetaData(clean_breast, percent.mito, "percent.mito")
clean_breast <- AddMetaData(clean_breast, percent.ercc, "percent.ercc")
VlnPlot(clean_breast, c("nGene", "percent.ercc", "percent.mito"), nCol = 3)

clean_breast <- MeanVarPlot(clean_breast ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5, do.contour = F)
length(clean_breast@var.genes)
#[1] 3079

clean_breast <- PCA(clean_breast, pc.genes = clean_breast@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)

clean_breast <- ProjectPCA(clean_breast)

PCAPlot(clean_breast, 2, 3)
PCHeatmap(clean_breast, pc.use = 2, cells.use = 100, do.balanced = TRUE)

clean_breast <- JackStraw(clean_breast, num.replicate = 100, do.print = FALSE)
JackStrawPlot(clean_breast, PCs = 1:20)
PCElbowPlot(clean_breast)
PCHeatmap(clean_breast, pc.use = 1:9, cells.use = 100, do.balanced = TRUE, label.columns = F, use.full = FALSE)

clean_breast <- FindClusters(clean_breast, pc.use = 1:5, resolution = 1, print.output = 0, save.SNN = T)

clean_breast <- RunTSNE(clean_breast, dims.use = 1:5)
#note that you can set do.label=T to help label individual clusters
TSNEPlot(clean_breast)

save(clean_breast, file = "~/Documents/E12/seurat/clean_breast.Robj")

VlnPlot(clean_breast, c("KPNA2", "TNFRSF11B", "FDCSP", "HMGB2", "CCL2", "TADA2A", "ANKRD1", "MRPL21", "KBTBD2", "CASK"))
FeaturePlot(clean_breast,  c("KPNA2", "TNFRSF11B", "FDCSP", "HMGB2", "CCL2", "TADA2A", "ANKRD1", "MRPL21", "KBTBD2", "CASK"),cols.use = c("grey","red"))
