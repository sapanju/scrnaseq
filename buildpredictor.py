import numpy as np
from sklearn import datasets, linear_model
#import sklearn.cross_validation
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

        # Predicted AUC for each cell line (all for the same drug)
        predicted = np.empty(all_inputs[0].shape)
        mean_sq_err = np.empty(all_inputs[0].shape)

        # Leave one out cross validation
        for i in range(len(all_inputs[0])):
            fold_identity = np.array(np.zeros(len(all_inputs[0])), dtype=bool)
            fold_identity[i] = True

            train_inputs = all_inputs[:, ~fold_identity].T
            train_labels = all_labels[~fold_identity]
            valid_inputs = all_inputs[:, fold_identity].T
            valid_labels = all_labels[fold_identity]

            # Create linear regression object
            regr = linear_model.LinearRegression()

            regr.fit(train_inputs, train_labels)

            #print('Coefficients: \n', regr.coef_)
            # The mean squared error
            print("Drug %d: %s" % (drug, drug_names[drug]))
            print("Mean squared error: %.2f"
                  % np.mean((regr.predict(valid_inputs) - valid_labels) ** 2))
            # Explained variance score: 1 is perfect prediction
            print('Variance score: %.2f' % regr.score(valid_inputs, valid_labels))

            predicted[i] = regr.predict(valid_inputs)

            mean_sq_err[i] = (np.mean(predicted[i] - valid_labels) ** 2)


        # Plot outputs
        plt.scatter(all_labels, predicted, color='black')
        plt.savefig('drug%d.png' % (drug + 1))
        plt.close()
        #plt.plot(valid_inputs, regr.predict(valid_inputs), color='blue', linewidth=3)

        #plt.xticks(())
        #plt.yticks(())

        #plt.show()


if __name__ == '__main__':
    main()
