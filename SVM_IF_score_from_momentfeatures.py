#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:19:55 2024

@author: rkanjilal
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest

def perform_one_class_svm(input_file):
    dt = pd.read_csv(input_file)
    dt_noinf = dt[~dt.isin([np.nan, np.inf, -np.inf]).any(1)]

    X = dt_noinf.iloc[:, 2:]

    # Normalization
    scaler = StandardScaler()
    X_scaler = scaler.fit_transform(X)

    clf = OneClassSVM(gamma='auto').fit(X_scaler)
    anomaly_scores_svm = clf.score_samples(X_scaler)

    value = dt_noinf.iloc[:, 0:2].values
    anomaly_scores_svm1 = anomaly_scores_svm.reshape(anomaly_scores_svm.shape[0], 1)

    x_new = np.concatenate((value, anomaly_scores_svm1), axis=1)

    datafm = pd.DataFrame(x_new, columns=("Chr", "Pos", "OCSVM_score"))
    return datafm

def perform_isolation_forest(input_file):
    dt = pd.read_csv(input_file)
    dt_noinf = dt[~dt.isin([np.nan, np.inf, -np.inf]).any(1)]

    X = dt_noinf.iloc[:, 2:]

    # Normalization
    scaler = StandardScaler()
    X_scaler = scaler.fit_transform(X)

    model_IF = IsolationForest(random_state=0).fit(X_scaler)
    anomaly_scores_if = model_IF.decision_function(X_scaler)

    value = dt_noinf.iloc[:, 0:2].values
    anomaly_scores_if1 = anomaly_scores_if.reshape(anomaly_scores_if.shape[0], 1)

    x_new2 = np.concatenate((value, anomaly_scores_if1), axis=1)

    datafm2 = pd.DataFrame(x_new2, columns=("Chr", "Pos", "IF_score"))
    return datafm2

def main():
    input_moments = "moment_features.csv"
    input_fda = "fda_features.csv"
    
    scoresvm_mom = perform_one_class_svm(input_moments)
    scoreif_mom = perform_isolation_forest(input_moments)
    scoreif_fda = perform_isolation_forest(input_fda)
    
    output_moments_svm = "SVM_score_from_moments.csv"
    output_moments_if = "IF_score_from_moments.csv"
    output_fda_if = "IF_score_from_fda.csv"
    
    scoresvm_mom.to_csv(output_moments_svm, index=None)
    scoreif_mom.to_csv(output_moments_if, index=None)
    scoreif_fda.to_csv(output_fda_if, index=None)

# Run the main function
main()

