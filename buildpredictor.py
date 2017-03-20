import numpy as np
from sklearn import datasets, linear_model
from sklearn.neighbors import KNeighborsRegressor
#import sklearn.cross_validation
import matplotlib.pyplot as plt


def load_data():
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

    return auc, expr, cell_names, feature_names, drug_names


def main():

    # Load data
    (auc, expr, cell_names, feature_names, drug_names) = load_data()

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

        # Only work with drugs that have balanced class distribution (AUC < 0.2 is considered sensitive)
        num_sensitive = sum(all_labels <= 0.2)
        num_insensitive = sum(all_labels > 0.2)

        ratio = (num_insensitive + 0.0)/num_sensitive

        if (ratio > 0.2 and ratio < 1.8):
            # Predicted AUC for each cell line (all for the same drug)
            predicted = np.empty(all_inputs[0].shape)
            mean_sq_err = np.empty(all_inputs[0].shape)

            # k = cross_validate_select_k(all_inputs, all_labels, drug_names, drug, False)
            # print "Drug %s: k = %d" % (drug_names[drug], k)

            alpha = cross_validate_select_alpha(all_inputs, all_labels, drug_names, drug, 'lasso', False)
            print "Drug %s: alpha = %3f" % (drug_names[drug], alpha)

            #Leave one out cross validation
            for i in range(len(all_inputs[0])):
                fold_identity = np.array(np.zeros(len(all_inputs[0])), dtype=bool)
                fold_identity[i] = True

                train_inputs = all_inputs[:, ~fold_identity].T
                train_labels = all_labels[~fold_identity]
                valid_inputs = all_inputs[:, fold_identity].T
                valid_labels = all_labels[fold_identity]

                # Create linear regression object
                model = linear_model.Lasso(alpha=alpha)
                # model = linear_model.Ridge(alpha=0.1)
                # model = KNeighborsRegressor(n_neighbors=k, weights='distance')
                model.fit(train_inputs, train_labels)

                #print('Coefficients: \n', model.coef_)

                predicted[i] = model.predict(valid_inputs)

                mean_sq_err[i] = (np.mean(predicted[i] - valid_labels) ** 2)


            # Plot outputs
            fit = np.polyfit(all_labels, predicted, deg=1)
            if (fit[0] > 0.2 and fit[0] < 1.8):
                plt.scatter(all_labels, predicted, color='black')
                plt.plot(all_labels, fit[0] * all_labels + fit[1], color='grey')
                plt.title("%s" % drug_names[drug])
                plt.xlabel("Ground truth AUC")
                plt.ylabel("Predicted AUC")
                plt.savefig('drug%d.png' % (drug + 1))
                plt.close()
                #plt.plot(valid_inputs, model.predict(valid_inputs), color='blue', linewidth=3)

                #plt.xticks(())
                #plt.yticks(())

                #plt.show()


def cross_validate_select_k(all_inputs, all_labels, drug_names, drug, doPlot):
    num_k_to_test = 10
    sum_mean_sq_err = np.empty(num_k_to_test)

    for k in range(0, num_k_to_test):

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

            model = KNeighborsRegressor(n_neighbors=(k+1), weights='distance')
            model.fit(train_inputs, train_labels)

            predicted[i] = model.predict(valid_inputs)

            mean_sq_err[i] = (np.mean(predicted[i] - valid_labels) ** 2)

        sum_mean_sq_err[k] = sum(mean_sq_err)

    if doPlot:
        plt.scatter(range(1, num_k_to_test+1), sum_mean_sq_err)
        plt.title("%s Hyperparameter selection" % drug_names[drug])
        plt.xlabel("K (# of neighbours)")
        plt.ylabel("Sum of Mean Squared Error")
        plt.savefig('drug%d_hyperparam_k.png' % (drug + 1))
        plt.close()

    k_min = np.argmin(sum_mean_sq_err) + 1

    return k_min


def cross_validate_select_alpha(all_inputs, all_labels, drug_names, drug, model_type, doPlot):

    #num_alpha_to_test = 10
    start = 0.01
    end = 1
    test_range = [0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0]
    sum_mean_sq_err = np.empty(len(test_range))

    for alpha in range(len(test_range)):

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

            if model_type == 'lasso':
                model = linear_model.Lasso(alpha=test_range[alpha])
            else:
                model = linear_model.Ridge(alpha=test_range[alpha])

            model.fit(train_inputs, train_labels)

            predicted[i] = model.predict(valid_inputs)

            mean_sq_err[i] = (np.mean(predicted[i] - valid_labels) ** 2)

        sum_mean_sq_err[alpha] = sum(mean_sq_err)

    if (doPlot):
        plt.scatter(test_range, sum_mean_sq_err)
        plt.title("%s Hyperparameter selection" % drug_names[drug])
        plt.xlabel("alpha (regularization strength)")
        plt.ylabel("Sum of Mean Squared Error")
        plt.savefig('drug%d_hyperparam_alpha_%s.png' % (drug + 1, model_type))
        plt.close()

    alpha_min = test_range[np.argmin(sum_mean_sq_err)]

    return alpha_min


if __name__ == '__main__':
    main()
