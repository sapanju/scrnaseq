import numpy as np
from sklearn import datasets, linear_model
from sklearn.neighbors import KNeighborsRegressor
#import sklearn.cross_validation
import matplotlib.pyplot as plt


def load_train_data():
    auc = np.genfromtxt('GDSC1000_breast_auc.txt', delimiter=',')
    expr = np.genfromtxt('GDSC1000expression_breast_logged.txt', delimiter=',')

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


def load_test_data():
    test_data = np.genfromtxt('top_898_genes_e12_orderedbycluster.txt', delimiter=',')

    with open('top_898_genes_e12_orderedbycluster_cellnames.txt', 'r') as myfile:
        test_data_cell_names = myfile.readlines()
    for line in range(len(test_data_cell_names)):
        test_data_cell_names[line] = test_data_cell_names[line].replace('\n', '').replace('"', '')

    return test_data, test_data_cell_names


def main():
    # Load data
    (auc, expr, cell_names, feature_names, drug_names) = load_train_data()
    test_data, test_data_cell_names = load_test_data()

    e12_predicted = np.zeros([len(drug_names), len(test_data[0])])

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
            predicted_lasso = np.empty(all_inputs[0].shape)
            predicted_ridge = np.empty(all_inputs[0].shape)
            predicted_knn = np.empty(all_inputs[0].shape)
            mean_sq_err_lasso = np.empty(all_inputs[0].shape)
            mean_sq_err_ridge = np.empty(all_inputs[0].shape)
            mean_sq_err_knn = np.empty(all_inputs[0].shape)

            e12_predicted_lasso = np.empty([len(drug_names), len(test_data[0])])
            e12_predicted_ridge = np.empty([len(drug_names), len(test_data[0])])
            e12_predicted_knn = np.empty([len(drug_names), len(test_data[0])])

            k = cross_validate_select_k(all_inputs, all_labels, drug_names, drug, False)
            print "Drug %s: k = %d" % (drug_names[drug], k)

            alpha_lasso = cross_validate_select_alpha(all_inputs, all_labels, drug_names, drug, 'lasso', False)
            alpha_ridge = cross_validate_select_alpha(all_inputs, all_labels, drug_names, drug, 'ridge', False)
            # print "Drug %s: alpha = %3f" % (drug_names[drug], alpha)

            #Leave one out cross validation
            for i in range(len(all_inputs[0])):
                fold_identity = np.array(np.zeros(len(all_inputs[0])), dtype=bool)
                fold_identity[i] = True

                train_inputs = all_inputs[:, ~fold_identity].T
                train_labels = all_labels[~fold_identity]
                valid_inputs = all_inputs[:, fold_identity].T
                valid_labels = all_labels[fold_identity]

                # Create linear regression object
                model_lasso = linear_model.Lasso(alpha=alpha_lasso)
                model_ridge = linear_model.Ridge(alpha=alpha_ridge)
                model_knn = KNeighborsRegressor(n_neighbors=k, weights='distance')
                model_lasso.fit(train_inputs, train_labels)
                model_ridge.fit(train_inputs, train_labels)
                model_knn.fit(train_inputs, train_labels)

                #print('Coefficients: \n', model.coef_)

                predicted_lasso[i] = model_lasso.predict(valid_inputs)
                predicted_ridge[i] = model_ridge.predict(valid_inputs)
                predicted_knn[i] = model_knn.predict(valid_inputs)

                mean_sq_err_lasso[i] = (np.mean(predicted_lasso[i] - valid_labels) ** 2)
                mean_sq_err_ridge[i] = (np.mean(predicted_ridge[i] - valid_labels) ** 2)
                mean_sq_err_knn[i] = (np.mean(predicted_knn[i] - valid_labels) ** 2)

            for sample in range(len(test_data[0])):
                e12_predicted_lasso[drug, sample] = model_lasso.predict(test_data[:, sample].reshape(-1, 1).T)
                e12_predicted_ridge[drug, sample] = model_ridge.predict(test_data[:, sample].reshape(-1, 1).T)
                e12_predicted_knn[drug, sample] = model_knn.predict(test_data[:, sample].reshape(-1, 1).T)

            e12_predicted[drug, :] = np.mean(np.array([e12_predicted_lasso[drug, :], e12_predicted_ridge[drug, :],
                                                  e12_predicted_knn[drug, :]]), axis=0)
            #print(e12_predicted[drug, :])
            mean_sq_err = np.mean([np.mean(mean_sq_err_lasso), np.mean(mean_sq_err_ridge), np.mean(mean_sq_err_knn)])

            t = np.arange(1, len(test_data[0])+1)
            # plt.plot(t, e12_predicted[drug, :], 'k--', t, e12_predicted_lasso[drug, :], 'r.', t,
            #                             e12_predicted_ridge[drug, :], 'b.', t, e12_predicted_knn[drug, :], 'g.')
            label_size = 6
            plt.rcParams['xtick.labelsize'] = label_size
            plt.errorbar(t, e12_predicted[drug, :], color='#0058A2', yerr=mean_sq_err, ecolor='#AF2F2A', elinewidth='0.5')
            plt.xticks(t, test_data_cell_names, rotation='vertical', )
            plt.subplots_adjust(bottom=0.25)
            plt.title("%s" % drug_names[drug])
            plt.savefig('drug%d.png' % (drug + 1))
            plt.close()
    np.savetxt('e12_predicted', e12_predicted, delimiter=',')
            # # Plot outputs
            # predicted = np.mean(np.array([predicted_lasso, predicted_ridge, predicted_knn]), axis=0)
            # fit = np.polyfit(all_labels, predicted, deg=1)
            # if fit[0] > 0.2 and fit[0] < 1.8:
            #     plt.scatter(all_labels, predicted, color='black')
            #     plt.plot(all_labels, fit[0] * all_labels + fit[1], color='grey')
            #     plt.title("%s" % drug_names[drug])
            #     plt.xlabel("Ground truth AUC")
            #     plt.ylabel("Predicted AUC")
            #     plt.savefig('drug%d.png' % (drug + 1))
            #     plt.close()


# function to train - output = list of viable drugs, and hyperparameters for ridge, lasso, knn
# function to test - input is output from train, output is auc for each cell under each drug
            # plot auc of each cell - one figure per drug
    # test_data = load_test_data()
    # for i in range(len(test_data[0])):
    #     for drug in range(len(drug_names)):




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
    test_range = [0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 5.0, 10.0]
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
