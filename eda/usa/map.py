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
parser = argparse.ArgumentParser(description = '''Analyze AMR strains in the USA.''')

parser.add_argument('--tsv', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to tsv with AMR data.')
parser.add_argument('--outdir', nargs=1, type= str,
                  default=sys.stdin, help = 'Path to output directory. Include /in end')

#####FUNCTIONS#####


#####MAIN#####
plt.rcParams.update({'font.size': 6})
args = parser.parse_args()
amr_data = pd.read_tsv(args.csv[0], sep='\t')
outdir = args.outdir[0]

map_changes(amr_data, outdir)
