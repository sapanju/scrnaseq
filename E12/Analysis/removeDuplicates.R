# Remove duplicates from a matrix
# NOTE: rownames must be numbers 1...length(matrix)
matrix = real_tpm_eset
rownames(matrix) = 1:20345

duplicates = matrix[which(duplicated(matrix$GeneName)),3]
duplicates = duplicates[!duplicated(duplicates)]
indices_to_remove = c()

for (dup in duplicates)
{
  occurrance_indices = which(matrix$GeneName == dup)
  iqr_list = apply(matrix[occurrance_indices, 5:89], 1, IQR) #1st column is gene name, so start at column 2
  iqr_list = sort(iqr_list, TRUE) # sort descending
  indices_to_remove = append(indices_to_remove, as.numeric(names(iqr_list[2:length(iqr_list)])))
}

# Matrix is the output
matrix = matrix[-indices_to_remove,]