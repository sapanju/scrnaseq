# Wilcox test

e12_predicted = read.table("e12_predicted.txt", sep = ",")
drug_names = read.table("drug_names.txt", sep = ",", header = FALSE)
top_898_genes_e12_orderedbycluster_cellnames = read.table("top_898_genes_e12_orderedbycluster_cellnames.txt", sep = ",", header = FALSE)

colnames(e12_predicted) = t(top_898_genes_e12_orderedbycluster_cellnames)
rownames(e12_predicted) = t(drug_names)

p = c()
for (drug in 1:length(e12_predicted[,1]))
{
  x = as.numeric(e12_predicted[drug, 1:25])
  y = as.numeric(e12_predicted[drug, 26:length(e12_predicted[drug,])])
  p[drug] = wilcox.test(x,y)$p.value
}