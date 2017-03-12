#knn

source("https://bioconductor.org/biocLite.R")
biocLite("PharmacoGx")
library(PharmacoGx)
library(Hmisc)

load("~/PSets/GDSC1000.RData")
setwd("~/Documents/e09/drugprediction")

load("~/Documents/E12/drugprediction/GDSC1000_auc.RData")
load("~/Documents/E12/drugprediction/GDSC1000_breast_auc.RData")
load("~/Documents/E12/drugprediction/GDSCexpression.RData")
load("~/Documents/E12/drugprediction/GDSCexpression_breast.RData")

annotation = read.table("~/Documents/E12/Analysis/revisedHg19annotationFileSummary.txt", quote="\"")
#top_1000_genes_eset <- read.delim("~/Documents/E12/QC/top_1000_genes_eset.txt")
annotation$V1 = strtrim(annotation$V1, 15) # trim off the decimal places in the ensemble id name
names(annotation)[1] = "EnsembleId"
names(annotation)[2] = "GeneType"
names(annotation)[3] = "GeneName"
gene_names = rownames(top_1000_genes_eset)
feature_names = annotation$EnsembleId[match(gene_names, annotation$GeneName)]
feature_names = paste(feature_names, "_at", sep = "")
annotation$EnsembleId = paste(annotation$EnsembleId, "_at", sep = "")

auc = GDSC1000_breast.auc;
expr = as.matrix(GDSCexpression_breast_allgenes);
expr = expr[,-c(52,51)]; # These rows are all NA
auc = as.matrix(auc[,-c(52,51)]); # remove the corresponding elements in auc response

errors = cbind(rep(0, times = length(auc[,1])), rep(0, times = length(auc[,1])))
rownames(errors) = rownames(auc)
# par(mfrow = c(13,4))
# for (cell in 1:52)
# {
#   plot(expr[cell,], auc[cell,], main = colnames(expr)[cell], xlab = "expr", ylab = "auc");
# }

# correlations = matrix(, nrow(expr), nrow(auc));
# rownames(correlations) = rownames(expr)
# colnames(correlations) = rownames(auc)
# for (i in 1:length(expr[,1]))
# {
#   for (j in 1:length(auc[,1]))
#   {
#     correlations[i,j] = cor(expr[i,], auc[j,], use="pairwise.complete.obs", method = "pearson")
#   }
# }

balanced_drug_list = c();

figs = cut2(1:252, m = 9, onlycuts = TRUE)


for (fig in 1:(length(figs)-1))
{
  png(paste("knn_drugbatch", fig, ".png", sep = ""), width = 1400, height = 800)
  par(mfrow = c(2,5))
  for (drug in figs[fig]:(figs[fig + 1] - 1))
  {
    # SINGLE GENE METHOD FOR DRUG J
    all_inputs = expr
    all_labels = auc[drug,]
    
    all_predicted = c();
    all_valid_labels = c();
    
    # Remove NA's
    toRemove = which(is.na(all_labels))
    if (length(toRemove != 0))
    {
      all_labels = all_labels[-toRemove]
      all_inputs = all_inputs[,-toRemove]
    }
    
    num_sensitive = sum(all_labels <= 0.2)
    num_insensitive = sum(all_labels > 0.2)
    
    ratio = (num_insensitive + 0.0)/num_sensitive
    
    # Only work with drugs that have a balanced class distribution (sensitivity vs insensitivity)
    if (0.2 < ratio && ratio < 1.8)
    {
      balanced_drug_list = c(balanced_drug_list, drug)
      
      cuts = cut2(1:length(all_labels), m = 5, onlycuts = TRUE)
      
      # par(mfrow = c(2,5))
      for (i in 1:(length(cuts)-1))
      {
        start = cuts[i];
        end = cuts[i+1] - 1;
        if (i == (length(cuts)-1))
        {
          end = cuts[i+1];
        }
        
        train_inputs = all_inputs[,-(start:end)]
        train_labels = all_labels[-(start:end)]
        valid_inputs = all_inputs[,start:end]
        valid_labels = all_labels[start:end]
        
        out = knn(t(train_inputs), t(valid_inputs), train_labels, k = 1, prob = FALSE)
        predicted_auc = as.numeric(as.character(out))
        
        mean_sq_err = mean((predicted_auc - valid_labels)^2)
        percent_err = mean(abs(predicted_auc - valid_labels)/valid_labels)
        errors[drug,1] = mean_sq_err
        errors[drug,2] = percent_err
        #predictor = which(corr == max(corr))
        
        # fit <- lm(train_labels ~ train_inputs[predictor,])
        # slope = fit$coefficients[[2]]
        # intercept = fit$coefficients[[1]]
        # 
        # predicted_auc = (valid_inputs[predictor,] * slope) + intercept
        
        #gene_name = annotation$GeneName[which(annotation$EnsembleId == rownames(all_inputs)[predictor])]
        #gene_name = gene_names[which(feature_names == rownames(all_inputs)[predictor])]
        # plot(predicted_auc, valid_labels, main = paste("Fold", i, ":", gene_name))
        # abline(lm(valid_labels ~ predicted_auc))
        
        all_predicted = c(all_predicted, predicted_auc)
        all_valid_labels = c(all_valid_labels, valid_labels)
      }
      
      fit <- lm(all_valid_labels ~ all_predicted)
      slope = fit$coefficients[[2]]
      intercept = fit$coefficients[[1]]
      
      corr = cor(all_predicted, all_valid_labels, use="pairwise.complete.obs", method = "pearson")
      #if (abs(corr) >= 0.6)
      #{
        plot(all_predicted, all_valid_labels, main=rownames(auc)[drug], cex.main = 2)
        corr = cor(all_predicted, all_valid_labels, use="pairwise.complete.obs", method = "pearson")
        abline(intercept, slope)
        intercept = round(intercept, 2)
        slope = round(slope, 2)
        equation = paste("y = ", slope, "x + ", intercept, sep = "")
        mtext(equation, 3, line=-2)
        mtext(round(corr, digits = 2), 3, line=-3)
      #}
    }
   
  }
  dev.off()
}
