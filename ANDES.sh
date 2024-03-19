#!/bin/bash

# run 'chmod +x ANDES.sh' in terminal to give execute permissions to 'ANDES.sh'

# Then run './ANDES.sh' on terminal

# Call the Python script for generating summary statistics and moment features. The user can provide multiple vcf files as 'python ./vcf_ss_M_features.py --vcf CEU21.vcf CEU22.vcf'
python ./vcf_ss_M_features.py --vcf "$1"

# Call the R script for generating FDA features
Rscript ./ss_fda_features.R

# Call the R script for generating MD-M and MD-F anomaly scores
Rscript ./MD_MF_scores.R

# Call the Python script for training IF and SVM from moment and FDA features
python ./IFSVM_training.py

# Call the R script for generating IF-M, IF-F, SVM-M, and SVM-F anomaly scores
Rscript ./IFSVM_MF_anomalyscores.R

