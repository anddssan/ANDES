#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 13:34:07 2024

@author: rkanjilal
"""

#### $ conda install -c conda-forge scikit-allel  or $ pip install scikit-allel[full]  (To install scikit-allel with all optional dependencies via pip)

import argparse
import pandas as pd
import allel
import numpy as np
from scipy.stats import skew, kurtosis
from scipy.spatial.distance import pdist


def unq_count(x):
    dat = np.dtype((np.void, x.dtype.itemsize * x.shape[1]))
    b = np.ascontiguousarray(x).view(dat)
    unq, cnt = np.unique(b, return_counts=True)
    unq = unq.view(x.dtype).reshape(-1, x.shape[1])
    return unq, cnt

def vcf_to_summary_statistics(input_vcf, chr_num):
    callset = allel.read_vcf(input_vcf)

    gt = allel.GenotypeArray(callset['calldata/GT'])
    pos = callset['variants/POS']
    pos = pos.reshape(len(pos), 1)

    dt = np.empty((99, 0))
    for i in range(len(gt)):
        p = gt[i, :, :]
        ind = np.empty((0, 1))
        for j in range(len(p)):
            count = p[j, 0] + p[j, 1]
            count = count.reshape(1, 1)
            ind = np.concatenate([ind, count], axis=0)

        dt = np.concatenate([dt, ind], axis=1)

    arr = np.empty([0, 9])
    l = dt.shape[1]
    ln = (l - 51) + 1
    for i in range(0,ln):
        xx=dt[:,i:i+51]
        
        ### cityblock=> manhattan, pairwise distances between observations in n-dimensional space.
        
        Manhattan_dist1=pdist(xx, 'cityblock')  
        Manhattan_dist1=np.reshape(Manhattan_dist1,(4851,1))
        Manhattan_dist=Manhattan_dist1/51       ###### distance=distance/no. of SNPs
        
        ######### Mean, Var, Skewness, Kurtosis summary statistics computation ############
        
        mean=np.mean(Manhattan_dist)
        mean=mean.reshape(1,1)
        var=np.var(Manhattan_dist)
        var=var.reshape(1,1)
        skew1=skew(Manhattan_dist)
        skew1=np.reshape(skew1,(1,1))
        kurtosis1=kurtosis(Manhattan_dist)
        kurtosis1=np.reshape(kurtosis1,(1,1))
        
        #### finding the position of the middle window
        
        position=pos[i:i+51]
        midwindow_id = 25                  ### midwindow_id = round(len(position) / 2) middle number is 25th because in python the window length 51=> 0 to 50
        midwindow_pos=position[midwindow_id]  
        midwindow_pos=midwindow_pos.reshape(1,1)
        
        ######## frequency of the 1st, 2nd, 3rd and 4th most common string summary statistics computation ##############
        
        un,cn=unq_count(xx)
        if (len(cn)>=4):
            flat=np.sort(cn)[::-1]
            com1=flat[0]/99  # 1st most common   ##flat[-1]/99 ; flat[-2]/99 ; flat[-3]/99 ; flat[-4]/99 
            com1=np.reshape(com1,(1,1))
            com2=flat[1]/99  # 2nd most common
            com2=np.reshape(com2,(1,1))
            com3=flat[2]/99  # 3rd most common
            com3=np.reshape(com3,(1,1))
            com4=flat[3]/99  # 4th most common
            com4=np.reshape(com4,(1,1))
        if(len(cn)==3):
            flat=np.sort(cn)[::-1]
            com1=flat[0]/99  # 1st most common   ##flat[-1]/99 ; flat[-2]/99 ; flat[-3]/99 ; flat[-4]/99 
            com1=np.reshape(com1,(1,1))
            com2=flat[1]/99  # 2nd most common
            com2=np.reshape(com2,(1,1))
            com3=flat[2]/99  # 3rd most common
            com3=np.reshape(com3,(1,1))
            com4=0/99  # 4th most common
            com4=np.reshape(com4,(1,1))
        if(len(cn)==2):
            flat=np.sort(cn)[::-1]
            com1=flat[0]/99  # 1st most common   
            com1=np.reshape(com1,(1,1))
            com2=flat[1]/99  # 2nd most common
            com2=np.reshape(com2,(1,1))
            com3=0/99  # 3rd most common
            com3=np.reshape(com3,(1,1))
            com4=0/99 # 4th most common
            com4=np.reshape(com4,(1,1))
        if(len(cn)==1):
            flat=np.sort(cn)[::-1]
            com1=flat[0]/99  # 1st most common   
            com1=np.reshape(com1,(1,1))
            com2=0/99 # 2nd most common
            com2=np.reshape(com2,(1,1))
            com3=0/99  # 3rd most common
            com3=np.reshape(com3,(1,1))
            com4=0/99 # 4th most common
            com4=np.reshape(com4,(1,1))

        df = np.concatenate((midwindow_pos, mean, var, skew1, kurtosis1, com1, com2, com3, com4), axis=1)
        arr = np.concatenate((arr, df), axis=0)

    chrm = np.ones(len(arr)) * chr_num
    chrm = chrm.reshape(len(chrm), 1)

    chrdata = np.concatenate((chrm, arr), axis=1)
    chr_df = pd.DataFrame(chrdata, columns=('Chr', 'Pos', 'Mean', 'Var', 'Skew', 'Kurtosis',
                                             '1st_commonfreq', '2nd_commonfreq', '3rd_commonfreq', '4th_commonfreq'))

    return chr_df


def summary_statistics_to_moments(input_ss):
    dff = pd.DataFrame()
    for chr_num in range(1, 23): 
        X = input_ss[input_ss["Chr"] == chr_num]
        x = X.values

        array = np.empty([0, 34])
        l = len(x)
        ln = (l - 129) + 1
        for i in range(0, ln):
            pos = x[i:i+129, 1]
            position = pos[64]  # getting a middle position
            position = position.reshape(1, 1)

            xx = x[i:i+129, 2:10]  # Summary Statistics
            mm = np.mean(xx, axis=0)
            m = np.reshape(mm, (1, 8))
            vv = np.var(xx, axis=0)
            v = np.reshape(vv, (1, 8))
            ss = skew(xx, axis=0)
            s = np.reshape(ss, (1, 8))
            kk = kurtosis(xx, axis=0)
            k = np.reshape(kk, (1, 8))
            chrom_num = np.reshape(chr_num, (1, 1))
            df = np.concatenate((chrom_num, position, m, v, s, k), axis=1)
            array = np.concatenate((array, df), axis=0)

        data = pd.DataFrame(array, columns=('Chr', 'Pos', 'Mean of Mean', 'Mean of Var', 'Mean of Skew', 'Mean of Kur', 'Mean of 1stFreq', 'Mean of 2ndFreq', 'Mean of 3rdFreq', 'Mean of 4thFreq', 'Var of Mean', 'Var of Var', 'Var of Skew', 'Var of Kur', 'Var of 1stFreq', 'Var of 2ndFreq', 'Var of 3rdFreq', 'Var of 4thFreq', 'Skew of Mean', 'Skew of var', 'Skew of skew', 'Skew of Kur', 'Skew of 1stFreq', 'Skew of 2ndFreq', 'Skew of 3rdFreq', 'Skew of 4thFreq', 'Kur of Mean', 'Kur of Var', 'Kur of Skew', 'Kur of Kur', 'Kur of 1stFreq', 'Kur of 2ndFreq', 'Kur of 3rdFreq', 'Kur of 4thFreq'))
        dff = pd.concat([dff, data], ignore_index=True)

    return dff

def main():
    parser = argparse.ArgumentParser(description='Process VCF files.')
    parser.add_argument('--output_combined_path', type=str, default='/Users/rkanjilal/Documents/summary_statistics.csv',
                        help='Path to save the combined DataFrame to a CSV file.')

    args = parser.parse_args()

    # Initialize an empty DataFrame to store concatenated results
    combined_df = pd.DataFrame()

    for chr_num in range(1, 23):  # Assuming chromosomes from 1 to 22
        input_vcf_path = f'/Users/rkanjilal/Documents/CEU{chr_num}.vcf'
        
        # Call vcf_to_summary_statistics and get the result DataFrame
        result_df = vcf_to_summary_statistics(input_vcf_path, chr_num)

        # Concatenate the result to the combined DataFrame
        combined_df = pd.concat([combined_df, result_df], ignore_index=True)

    # Save the combined DataFrame to a CSV file
    combined_df.to_csv(args.output_combined_path, index=None)
    
    # Read summary statistics data
    input_ss_path = '/Users/rkanjilal/Documents/summary_statistics.csv'
    output_moments_path = '/Users/rkanjilal/Documents/momemt_features.csv'

    dt = pd.read_csv(input_ss_path)
    
    # Call summary_statistics_to_moments and get the resukts as 32 moment features
    moments_df = summary_statistics_to_moments(dt)

    moments_df.to_csv(output_moments_path, index=None)

if __name__ == "__main__":
    main()



