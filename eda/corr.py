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
        print(md,points.shape[0])
        continue
        for i in range(len(regions)):
            fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
            region1 = regions[i]
            md_data_1 = md_data[md_data['RegionName']==region1]
            max_corr = 0
            min_corr = 0

            for j in range(len(regions)):
                if i == j:
                    continue
                region2 = regions[j] #Comparison region
                md_data_2 = md_data[md_data['RegionName']==region2]
                #Join on year
                merge = pd.merge(md_data_1,md_data_2,on='Time',how='left')
                y = np.array(merge['NumValue_x'],dtype='float32')
                x = np.array(merge['NumValue_y'],dtype='float32')
                #Remove y nans
                x = x[~np.isnan(y)]
                y = y[~np.isnan(y)]
                #Remove x nans
                y = y[~np.isnan(x)]
                x = x[~np.isnan(x)]

                #Correlate
                cs = []
                ps = []
                if len(x)-5<1:
                    continue
                for s in range(len(x)-5):
                    if s == 0:
                        R, p = pearsonr(x,y)
                    else:
                        R, p = pearsonr(x[:-s],y[s:])
                    cs.append(R)
                    ps.append(p)
                #Save max and min corr
                if max(cs)>max_corr:
                    max_cs = cs
                    max_corr = max(cs)
                    max_region = region2
                if min(cs)<min_corr:
                    min_cs = cs
                    min_corr = min(cs)
                    min_region = region2
            #Plot
            plt.plot(np.arange(len(max_cs)),max_cs,label=max_region)
            plt.plot(np.arange(len(min_cs)),min_cs,label=min_region)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.title(region1)
            plt.legend()
            plt.xlabel('Negativ shift (years)')
            plt.ylabel('Pearson R')
            plt.tight_layout()
            plt.savefig(outdir+'combo1/'+region1+'.png',format='png', dpi=300)
    pdb.set_trace()


    #'HealthTopic', 'Population', 'Indicator', 'Unit', 'Time', 'RegionCode','RegionName', 'NumValue'

#####MAIN#####
plt.rcParams.update({'font.size': 6})
args = parser.parse_args()
amr_data = pd.read_csv(args.csv[0])
outdir = args.outdir[0]

correlate_changes(amr_data, outdir)
