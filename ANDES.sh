#!/bin/bash

## Command #1
python ./vcf_ss_M_features.py --vcf "$@"

## Command #2
Rscript ./ss_fda_features.R

## Command #3
Rscript ./MD_MF_scores.R

## Command #4
python ./IFSVM_training.py

## Command #5
Rscript ./IFSVM_MF_anomalyscores.R

