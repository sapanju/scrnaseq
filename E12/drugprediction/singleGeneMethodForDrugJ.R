# SINGLE GENE METHOD FOR DRUG J
all_inputs = expr
all_labels = auc[204,]

all_predicted = c();
all_valid_labels = c();

cuts = cut2(1:length(all_labels), m = 5, onlycuts = TRUE)

par(mfrow = c(2,5))
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
  
  corr = c()
  for (singlegene in 1:nrow(train_inputs))
  {
    corr[singlegene] = cor(train_inputs[singlegene,], train_labels, use="pairwise.complete.obs", method = "pearson")
  }
  
  corr = abs(corr)
  predictor = which(corr == max(corr))
  
  fit <- lm(train_labels ~ train_inputs[predictor,])
  slope = fit$coefficients[[2]]
  intercept = fit$coefficients[[1]]
  
  predicted_auc = (valid_inputs[predictor,] * slope) + intercept
  
  gene_name = gene_names[which(feature_names == rownames(all_inputs)[predictor])]
  plot(predicted_auc, valid_labels, main = paste("Fold", i, ":", gene_name))
  
  abline(lm(valid_labels ~ predicted_auc))
  corr = cor(predicted_auc, valid_labels, use="pairwise.complete.obs", method = "pearson")
  intercept = round(intercept, 2)
  slope = round(slope, 2)
  equation = paste("y = ", slope, "x + ", intercept, sep = "")
  #mtext(equation, 3, line=-2)
  mtext(round(corr, digits = 2), 3, line=-3)
  
  all_predicted = c(all_predicted, predicted_auc)
  all_valid_labels = c(all_valid_labels, valid_labels)
}
  