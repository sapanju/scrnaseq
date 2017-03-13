from util import *
import numpy as np
import matplotlib.pyplot as plt


def main():
    # Define parameters
    K = 7
    iters = 200

    # Load data
    auc = np.genfromtxt('GDSC1000_breast_auc.txt', delimiter=',')
    expr = np.genfromtxt('GDSC1000expression_breast.txt', delimiter=',')

    with open('cell_names.txt', 'r') as myfile:
        cell_names = myfile.readlines()
    for line in range(len(cell_names)):
        cell_names[line] = cell_names[line].replace('\n', '').replace('"', '')

    with open('feature_names.txt', 'r') as myfile:
        feature_names = myfile.readlines()
    for line in range(len(feature_names)):
        feature_names[line] = feature_names[line].replace('\n', '').replace('"', '')

    with open('drug_names.txt', 'r') as myfile:
        drug_names = myfile.readlines()
    for line in range(len(drug_names)):
        drug_names[line] = drug_names[line].replace('\n', '').replace('"', '')

    # Remove columns which are just NA (last 2 cell lines)
    expr = np.delete(expr, 51, axis=1)
    expr = np.delete(expr, 50, axis=1)
    auc = np.delete(auc, 51, axis=1)
    auc = np.delete(auc, 50, axis=1)
    cell_names = cell_names[0:50]

    for drug in range(len(drug_names)):
        all_inputs = expr
        all_labels = auc[drug, :]

        # Remove further NA's
        all_inputs = all_inputs[:, ~np.isnan(all_labels)]
        all_labels = all_labels[~np.isnan(all_labels)]

        # Leave one out cross validation
        for iter in range(len(all_inputs[0])):
            fold_identity = np.zeros(len(all_inputs[0]))
            fold_identity[iter] = 1

            train_inputs = all_inputs[:, ~fold_identity]
            train_labels = all_labels[~fold_identity]
            valid_inputs = all_inputs[:, fold_identity]
            valid_labels = all_labels[fold_identity]



if __name__ == '__main__':
    main()
