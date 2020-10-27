#!/usr/bin/env python3
#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb
#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Parse and analyze CATH data .''')

parser.add_argument('--csv', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to ECDC csv with AMR.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#####FUNCTIONS#####
def correlate_changes(amr_data):
    '''Correlate the change in country i with changes in countries i+1,..,n
    x years earlier, where x=0,1
    '''
    #Select the resistance data
    res_data = amr_data[amr_data['Indicator']=='R - resistant isolates, percentage  ']
    years = np.sort(res_data['Time'].unique())
    microbe_drug = res_data['Population'].unique() #30 in total
    regions = res_data['RegionName'].unique() #30 in total

    md = 'Escherichia coli|Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside)'
    md_data = res_data[res_data['Population']==md]
    #Replace '-'
    md_data = md_data.replace({'NumValue':{'-':'nan'}})
    for i in range(len(regions)):
        region1 = regions[i]
        md_data_1 = md_data[md_data['RegionName']==region1]
        x = np.array(md_data_1['Time'])
        y = np.array(md_data_1['NumValue'],dtype='float32')
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        #Plot
        plt.plot(x,y, label = region1)
    plt.legend()
    plt.show()
    pdb.set_trace()


    #'HealthTopic', 'Population', 'Indicator', 'Unit', 'Time', 'RegionCode','RegionName', 'NumValue'

#####MAIN#####
args = parser.parse_args()
amr_data = pd.read_csv(args.csv[0])
outdir = args.outdir[0]

correlate_changes(amr_data)
