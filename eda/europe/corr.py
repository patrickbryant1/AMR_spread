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
parser = argparse.ArgumentParser(description = '''Parse and analyze CATH data .''')

parser.add_argument('--csv', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to ECDC csv with AMR.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#####FUNCTIONS#####
def correlate_changes(amr_data, outdir):
    '''Correlate the change in country i with changes in countries i+1,..,n
    x years earlier, where x=0,1
    '''
    #Select the resistance data
    res_data = amr_data[amr_data['Indicator']=='R - resistant isolates, percentage  ']
    years = np.sort(res_data['Time'].unique())
    microbe_drug = res_data['Population'].unique() #30 in total
    regions = res_data['RegionName'].unique() #30 in total

    #md = 'Escherichia coli|Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside)'

    for md in microbe_drug:             #['Acinetobacter spp.|Aminoglycosides']:
        md_data = res_data[res_data['Population']==md]
        #Replace '-'
        md_data = md_data.replace({'NumValue':{'-':'nan'}})
        points = np.array(md_data['NumValue'], dtype='float32')
        points = points[~np.isnan(points)]
        #print(md,points.shape[0])
        fig,ax = plt.subplots(figsize=(12/2.54, 9/2.54))
        for i in range(len(regions)):

            region1 = regions[i]
            md_data_1 = md_data[md_data['RegionName']==region1]
            x = np.array(md_data_1['NumValue'],dtype='float32')
            y=md_data_1['Time']

            #Fix nans
            for xi in range(len(x)):

                if np.isnan(x[xi])==True: #If nan
                    if xi>0:
                        if xi == len(x)-1:
                            x = x[:-1]
                            y = y[:-1]
                        else:
                            x[xi]=(x[xi-1]+x[xi+1])/2
                    else:
                        x[xi]=x[xi+1]

            #Smooth
            x_smoothed = np.zeros(len(x))
            step = 1
            for xi in range(len(x)-step):
                x_smoothed[xi]=np.average(x[xi:xi+step])
            #Set final
            x_smoothed[xi+1:]=x_smoothed[xi]

            #Plot
            plt.plot(y,x,color='mediumseagreen',alpha=0.5,linewidth=1)
            #plt.plot(y,x_smoothed,color='b',alpha=0.5,linewidth=1)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.title(md)
        plt.ylabel('% R')
        plt.tight_layout()
        plt.savefig(outdir+md+'.png',format='png', dpi=300)
        plt.close()
        print(md+'.png')
    pdb.set_trace()


    #'HealthTopic', 'Population', 'Indicator', 'Unit', 'Time', 'RegionCode','RegionName', 'NumValue'

#####MAIN#####
plt.rcParams.update({'font.size': 6})
args = parser.parse_args()
amr_data = pd.read_csv(args.csv[0])
outdir = args.outdir[0]

correlate_changes(amr_data, outdir)
