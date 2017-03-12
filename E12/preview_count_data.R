e12readcounts <- read.csv("~/Documents/E12/e12readcounts.csv", row.names=1)

genes <- e12readcounts[,1]
cellnames <- e12readcounts[1,]

cellnames <- cellnames[3:87]

counts <- e12readcounts[,3:87]

rownames(counts) <- make.unique(as.character(genes))

counts <- counts[rowSums(counts)>0,]
nGenes <- length(counts[,1])

rowSums(e12readcounts)>0
count(rowSums(e12readcounts)>0)
count.fields(rowSums(e12readcounts)>0)
save(counts,file="counts.RData")
coverage <- colSums(counts)/nGenes

bar.positions <- barplot(coverage,xaxt='n',ylab="Counts per gene")
axis(side=1,labels=cellnames)

axis(side=1,at=c(bar.positions),labels=cellnames, las=2)

counts <- counts[,coverage>1]

cellnames <- cellnames[coverage>1]
coverage <- coverage[coverage>1]
nCells <- length(cellnames)
counts.norm <- t(apply(counts,1,function(x) x/coverage)) # simple normalization method
top.genes <- tail(order(rowSums(counts.norm)),10)
expression <- log2(counts.norm[top.genes,]+0.001) # add a small constant in case there are zeros
library(caroline)
install.packages("caroline")
library(caroline)
violins(as.data.frame(t(expression)),connect=c(),deciles=FALSE,xlab="",ylab="log2 expression")
install.packages("sm")
violins(as.data.frame(t(expression)),connect=c(),deciles=FALSE,xlab="",ylab="log2 expression")
means <- apply(counts.norm,1,mean)
excess.var <- apply(counts,1,var)-means
excess.var[excess.var < 0] <- NA
overdispersion <- excess.var / means^2
hist(log2(overdispersion),main="Variance of read counts is higher than Poisson")

