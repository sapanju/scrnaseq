

library(edgeR)

load(top_1000_genes_eset)
group = c("G1", "G1", "G1" ,"G2", "G1", "G1", "G1" ,"G2" ,"G2", "G2" ,"G2" ,"G1", "G1" ,"G2", "G2", "G1" ,"G2" ,"G1", "G1" ,"G2", "G1" ,"G2" ,"G2" ,"G1" ,"G1", "G2" ,"G1" ,"G1", "G2", "G1", "G2" ,"G2");

cds <- DGEList( top_1000_genes_eset , group = group )
cds <- calcNormFactors( cds )
# effective library sizes
cds$samples$lib.size * cds$samples$norm.factors

plotMDS.DGEList( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )

cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds , prior.n = 10 )

et = exactTest(cds, prior.count = 10)

resultsByFC <- topTags( et , n = nrow( et$table ) , sort.by = "logFC" )$table
head( resultsByFC)

de.genes <- rownames( resultsByFC )[ resultsByFC$PValue <= 0.05 ]

summary( decideTestsDGE( et , p.value = 0.05 ) ) # the adjusted p-values are used here

hist( resultsTbl.tgw[de.genes.tgw[1:100],"logConc"] , breaks=25 , xlab="Log Concentration" ,
      col="blue" , xlim=c(-18,-6) , ylim=c(0,0.4) , freq=FALSE , main="Tagwise: Top 100" )