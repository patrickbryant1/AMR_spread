#!/usr/bin/env python3
#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

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
    #will be lowered. The start years are therefore saved.
    resistance_matrix = np.zeros((len(microbe_drug),len(regions)),dtype='int32')
    start_times = [] #First recording of data for each microbe_drug combination
    for i in range(len(microbe_drug)):
        md_data = res_data[res_data['Population']==microbe_drug[i]]
        #Replace '-'
        md_data = md_data.replace({'NumValue':{'-':'nan'}})
        #Get min year
        start_times.append(min(md_data['Time']))
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

    #Convert start times to array
    start_times = np.array(start_times)
    #Go through all countries and calculate
    #P(A) = number of times country i is above 1 %/number of possible times
    #P(B) = number of times country j is above 1 %/number of possible times
    above_1_percent = np.zeros(len(regions))
    for i in range(len(regions)):
        region_data = resistance_matrix[:,i]
        above_1_percent[i]=len(np.where(region_data>0)[0])/len(np.where(region_data!=-1)[0])

    #P(B|A)= number of times country j is already above 1 % when country i is above 1 %
    #divided by the number of times country i is above 1 %
    #The number of times both countries are above 1 % at the start of recording are subtracted
    already_above = np.zeros((len(regions),len(regions)))
    for i in range(len(regions)):
        region_i_data = resistance_matrix[:,i]
        #Look at the start times and the possibility of R 1 % being reached somewhere else
        before_possible = region_i_data-start_times
        #Missing data
        missing_i = np.where(before_possible>0)[0]
        print(regions[i],len(missing_i))

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

def graph_traversal(matrix, regions):
    '''Traverse and vidualize the graph using networkx
    '''
    #For directed graphs, explicitly mention create_using=nx.DiGraph,
    #and entry i,j of A corresponds to an edge from i to j.
    G = nx.from_numpy_array(matrix,create_using=nx.DiGraph)
    # use one of the edge properties to control line thickness
    edgewidth = [ d['weight'] for (u,v,d) in G.edges(data=True)]
    # layout
    pos = nx.spring_layout(G, iterations=50)
    labels = {}
    for i in range(30):
        labels[i]=regions[i]
    #pos=nx.circular_layout(G)
    nx.draw_networkx_labels(G,pos=pos,labels=labels)
    nx.draw_networkx_nodes(G, pos, with_labels=True)
    nx.draw_networkx_edges(G, pos, width=edgewidth,)
    plt.show()
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

#Graph
res_data = amr_data[amr_data['Indicator']=='R - resistant isolates, percentage  ']
regions = res_data['RegionName'].unique() #30 in total
graph_traversal(matrix,np.array(regions))
