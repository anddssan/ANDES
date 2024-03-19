#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 10:04:53 2024

@author: rkanjilal
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest

def isolation_forest(input_file):
    dt = pd.read_csv(input_file)
    dt_noinf = dt[~dt.isin([np.nan, np.inf, -np.inf]).any(axis=1)]

    X = dt_noinf.iloc[:, 2:]

    # Normalization
    scaler = StandardScaler()
    X_scaler = scaler.fit_transform(X)

    model_IF = IsolationForest(random_state=0).fit(X_scaler)
    anomaly_scores_if = model_IF.decision_function(X_scaler)

    value = dt_noinf.iloc[:, 0:2].values
    anomaly_scores_if1 = anomaly_scores_if.reshape(anomaly_scores_if.shape[0], 1)

    x_new = np.concatenate((value, anomaly_scores_if1), axis=1)

    datafm = pd.DataFrame(x_new, columns=("Chr", "Pos", "IF_score"))
    return datafm

def oneclass_svm_minibatchtraining(input_file, num_sets, interval_size):
    dt = pd.read_csv(input_file)
    dt_noinf = dt[~dt.isin([np.nan, np.inf, -np.inf]).any(axis=1)]

    X = dt_noinf.iloc[:, 2:]

    # Normalization
    scaler = StandardScaler()
    X_scaler = scaler.fit_transform(X)
    
    sets = [X_scaler[i::num_sets] for i in range(num_sets)]

    clf = OneClassSVM(kernel = 'rbf', gamma='auto')
    
    # Train the One-Class SVM model using mini-batches
    for i, mini_batch in enumerate(sets):
        start_index = i
        mini_batch_subset = mini_batch[start_index::interval_size]  #### sampling in interval_size
        clf.fit(mini_batch_subset)
    
    anomaly_scores_svm = clf.score_samples(X_scaler)

    value = dt_noinf.iloc[:, 0:2].values
    anomaly_scores_svm2 = anomaly_scores_svm.reshape(anomaly_scores_svm.shape[0], 1)

    x_new2 = np.concatenate((value, anomaly_scores_svm2), axis=1)

    datafm2 = pd.DataFrame(x_new2, columns=("Chr", "Pos", "OCSVM_score"))
    return datafm2

def main():
    input_moments = "./M_features.csv"
    input_fda = "./fda_features.csv"
    
    scoreif_mom = isolation_forest(input_moments)
    scoreif_fda = isolation_forest(input_fda)
    
    num_sets = 18
    interval_size = 10
    scoresvm_mom = oneclass_svm_minibatchtraining(input_moments, num_sets, interval_size)
    scoresvm_fda = oneclass_svm_minibatchtraining(input_fda, num_sets, interval_size)
    
    output_moments_if = "./IFscores_M.csv"
    output_fda_if = "./IFscores_F.csv"
    output_moments_svm = "./SVMscores_M.csv"
    output_fda_svm = "./SVMscores_F.csv"
    
    scoreif_mom.to_csv(output_moments_if, index=False)
    scoreif_fda.to_csv(output_fda_if, index=False)
    scoresvm_mom.to_csv(output_moments_svm, index=False)
    scoresvm_fda.to_csv(output_fda_svm, index=False)

main()
