#!/usr/bin/env python3
#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Calculate the probabilities for R coming from different
                                                countries using Baye's theorem.''')

parser.add_argument('--csv', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to ECDC csv with AMR.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#####FUNCTIONS#####
def apply_bayes(amr_data, outdir):
    '''
    1. Get the year each country reached 1 % R for
    2. Calculate the Bayesian probabilities for spread of R coming from each country
    to all countries.
    '''
    #Select the resistance data
    res_data = amr_data[amr_data['Indicator']=='R - resistant isolates, percentage  ']
    microbe_drug = res_data['Population'].unique() #30 in total
    regions = res_data['RegionName'].unique() #30 in total

    #Get the year each country reached 1 % R, write 0 if not, -1 if no data
    #Need to store when it is not possible to be before also by recording
    #the starting year for all mds. If this is not adjusted for the probabilities
    #will be lowered. Put -1 there.
    resistance_matrix = np.zeros((len(microbe_drug),len(regions)),dtype='int32')
    for i in range(len(microbe_drug)):
        md_data = res_data[res_data['Population']==microbe_drug[i]]
        #Replace '-'
        md_data = md_data.replace({'NumValue':{'-':'nan'}})
        for j in range(len(regions)):
            md_region_data = md_data[md_data['RegionName']==regions[j]]
            #Get year for which the resistance is above 1 %
            r = np.array(md_region_data['NumValue'],dtype='float32')
            excli = np.isnan(r)
            y = np.array(md_region_data['Time'])
            #Remove NaNs
            r = r[~excli]
            y = y[~excli]
            if len(r)<1:#If no data
                resistance_matrix[i,j]=-1
                continue
            if max(r)>1: #Reach 1 %
                above_i = np.where(r>1)[0][0]
                resistance_matrix[i,j]=y[above_i]
            else: #Never reached 1 % - will be represented with zeros
                continue

    #Go through all countries and calculate
    #P(A) = number of times country i is above 1 %/number of possible times
    #P(B) = number of times country j is above 1 %/number of possible times
    above_1_percent = np.zeros(len(regions))
    for i in range(len(regions)):
        region_data = resistance_matrix[:,i]
        above_1_percent[i]=len(np.where(region_data>0)[0])/len(np.where(region_data!=-1)[0])

    #P(B|A)= number of times country j is already above 1 % when country i is above 1 %
    #divided by the number of times country i is above 1 %
    already_above = np.zeros((len(regions),len(regions)))
    for i in range(len(regions)):
        region_i_data = resistance_matrix[:,i]
        #Missing data
        missing_i = np.where(region_i_data!=-1)[0]
        for j in range(len(regions)):
            region_j_data = resistance_matrix[:,j]
            #Missing data
            missing_j = np.where(region_j_data[missing_i]!=-1)[0]
            #Get the year difference in reaching 1 %
            #If >1 % resistance is reached in region j first, this is positive
            ij_diff = region_i_data-region_j_data
            #Get only non-missing data
            ij_diff = ij_diff[missing_i]
            ij_diff = ij_diff[missing_j]
            #Calc prob
            p_b_a = (len(np.where(ij_diff>0)[0])/len(ij_diff))
            p_a = above_1_percent[i]
            p_b = above_1_percent[j]
            #P(A|B)=P(B|A)P(A)/P(B)
            already_above[i,j]=p_b_a*p_a/p_b

    #Visualize Bayesian prob
    for i in range(len(regions)):
        fig,ax = plt.subplots(figsize=(12/2.54, 9/2.54))
        plt.bar(np.arange(len(regions)),already_above[i,:])
        plt.ylim([0,1])
        plt.title(regions[i] + ' ' + str(i))
        plt.xlabel('Country')
        plt.ylabel('Prob.')
        plt.tight_layout()
        plt.savefig(outdir+str(i)+'.png',format='png',dpi=300)
        plt.close()




    return already_above

def viterbi(matrix):
    '''Calculate the most probable path through the matrix using the Viterbi algorithm.
    The i rows represent the regions and the j columns the probabilities that the resistance
    would be spread from each region j to region i (item [i,j])
    '''

    #All regions
    for i in range(matrix.shape[0]-1):
        #Need to be able to start in all possible countries
        prob = np.zeros(matrix.shape[1])
        for j in range(matrix.shape[1]):
            prob[j] = matrix[i,j]*matrix[i+1,j]

        pdb.set_trace()
#####MAIN#####
plt.rcParams.update({'font.size': 6})
args = parser.parse_args()
amr_data = pd.read_csv(args.csv[0])
#'HealthTopic', 'Population', 'Indicator', 'Unit', 'Time', 'RegionCode','RegionName', 'NumValue'
outdir = args.outdir[0]

#Get the transition matrix (prob of spread from all countries to all countries)
try:
    matrix = np.load(outdir+'matrix.npy', allow_pickle=True)
except:
    matrix = apply_bayes(amr_data, outdir)
    np.save(outdir+'matrix.npy', matrix)

#Viterbi
viterbi(matrix)
