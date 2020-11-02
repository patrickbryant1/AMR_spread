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
def map_changes(amr_data, outdir):
    '''Map the AMR changes
    '''
    # (['#Organism Group', 'Strain', 'Serovar', 'Isolate', 'Create date',
    #    'Location', 'Isolation Source', 'Isolation type', 'SNP cluster',
    #    'Min-same', 'Min-diff', 'BioSample', 'Assembly', 'AMR genotypes',
    #    'Collection Date']

    #Map the locations to each state
    locations = amr_data['Location'].unique()


    twolettercodes=['AL', 'AK', 'AS', 'AZ', 'AR', 'CA', 'CO', 'CT',
    'DE', 'DC', 'FM', 'FL', 'GA', 'GU', 'HI', 'ID', 'IL',
    'IN', 'IA', 'KS', 'KY', 'LA', 'ME', 'MH', 'MD', 'MA',
     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 'NV', 'NH', 'NJ',
     'NM', 'NY', 'NC', 'ND', 'MP', 'OH', 'OK', 'OR', 'PW',
     'PA', 'PR', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT',
     'VI', 'VA', 'WA', 'WV', 'WI', 'WY']

    state_names = ['Alabama', 'Alaska', 'American Samoa', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut',
    'Delaware', 'District of Columbia', 'Federated States of Micronesia', 'Florida', 'Georgia', 'Guam', 'Hawaii', 'Idaho',
    'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Marshall Islands', 'Maryland', 'Massachusetts',
    'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey',
    'New Mexico', 'New York', 'North Carolina', 'North Dakota', 'Northern Mariana Islands', 'Ohio', 'Oklahoma', 'Oregon',
    'Palau', 'Pennsylvania', 'Puerto Rico', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah',
    'Vermont', 'Virgin Islands', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming']

    state_keys = {'USA: San Diego, CA':'CA', 'USA:GA':'GA', 'USA:VA':'VA', 'USA:TX':'TX',
       'USA:MD':'MD', 'USA:MI':'MI', 'USA:WY':'WY', 'USA: PA':'PA', 'USA: Menlo Park':'CA',
       'USA:IN':'IN', 'USA:MO':'MO', 'USA: IN':'IN', 'USA: CO':'CO', 'USA: MD':'MD', 'USA:FL':'FL',
       'USA:NJ':'NJ', 'USA:SC';'SC', 'USA:NE':'NE', 'USA:WI':'WI', 'USA: NY';'NY', 'USA: NE':'NE',
       'USA: NJ':'NJ', 'USA:PA':'PA', 'USA: FL':'FL', 'USA:IA':'IA', 'USA:AZ':'AZ', 'USA:TN':'TN',
       'USA:KS':'KS', 'USA:RI':'RI', 'USA:KY':'KY', 'USA: OR':'OR', 'USA: WI':'WI', 'USA: CA':'CA',
       'USA: MN':'MN', 'USA: WV':'WV', 'USA: MO':'MO', 'USA: MI':'MI', 'USA: GA':'GA', 'USA: TN':'TN',
       'USA: ID':'ID', 'USA: TX':'TX', 'USA: OK':'OK', 'USA: UT':'UT', 'USA: VT':'VT', 'USA: IA':'IA',
       'USA: KS':'KS', 'USA: SC':'SC', 'USA: KY':'KY', 'USA:NC':'NC', 'USA:CA':'CA', 'USA: NC':'NC',
       'USA:NY', 'USA:OH', 'USA: Rochester, NY',
       'USA: Imperial County, CA', 'USA:IL', 'USA:SD', 'USA:ME', 'USA:WA',
       'USA: Cambridge, MA, MIT', 'USA:MA', 'USA:ID',
       'USA: University of Florida', 'USA:MN', 'USA:AR', 'USA:CO',
       'USA:NH', 'USA:UT', 'USA:AL', 'USA:LA', 'USA: OH', 'USA: WA',
       'USA:Midwest', 'USA:AK', 'USA:VT', 'USA:OR', 'USA:DE', 'USA:ND',
       'USA:MT', 'USA:WV', 'USA:Boston', 'USA:PR', 'USA:MS',
       'USA: Nebraska', 'USA: Philadelphia', 'USA:OK', 'USA: California',
       'USA:CT', 'USA: Alaska', 'USA: Colorado', 'USA: New York',
       'USA: Mississippi', 'USA: Washington', 'USA:Washington DC',
       'USA: VA', 'USA: AZ', 'USA: MS', 'USA:Illinois', 'USA:Florida',
       'USA: AR', 'USA: WA, Seattle', 'USA: DE',
       'USA: Arkansas Broiler farm', 'USA: San Francisco, CA',
       'USA: Placer County, CA', 'USA: San Diego County, CA', 'USA:DC',
       'USA: Boston', 'USA: New Hampshire,Rindge', 'USA:TN-Nashville',
       'USA:CA-UC Davis', 'USA:TN-Knoxville',
       'USA: Saint Louis, Missouri', 'USA: Alachua, FL', 'USA: Alachua',
       'USA: St. Louis Clyde watershed of Lake Superior',
       'USA: Kaiser Permanente Washington Tacoma lab',
       "USA: Seattle Children's Hospital", 'USA: Minneapolis',
       'USA: Pennsylvania', 'USA: Missouri', 'USA: Virginia',
       'USA: Illinois', 'USA: Wisconsin',
       'USA: Santa Barbara, UCSB animal house',
       'USA: Soldotna landfill, Alaska', 'USA: Soldotna, Alaska',
       'USA: Lower Kenai River, Alaska', 'USA: South',
       'USA:Washington, DC', 'USA:AZ, Flagstaff', 'USA: North Carolina',
       'USA: IL', 'USA: Alachua County, Florida', 'USA: Cambridge, MA',
       'USA: ME', 'USA: Massachusetts, Boston', 'USA: Los Angeles',
       'USA:Tennessee', 'USA: Minnesota', 'USA: Cambridge', 'USA: NH',
       'USA: Florida', 'USA: Alameda County, CA',
       'USA: San Joaquin County, CA', 'USA: Santa Clara County, CA',
       'USA: Berkley, CA', 'USA: Ohio',
       'USA:University of Illinois Medical Center', 'USA: St. Louis, MO',
       'USA:Kansas', 'USA: University of Washington Medical Center',
       'USA: Seattle Harborview Medical Center', 'USA: Maryland',
       'USA: Texas', 'USA: Upper Kenai River, Alaska',
       'USA: Lower Kasilof River, Alaska', 'USA:Virginia',
       'USA: Oklahoma', 'USA: AL', 'USA: Walnut Creek',
       'USA: Gainesville, FL', 'USA:Gainesville, Florida', 'USA:HI',
       'USA: San Mateo County, CA', 'USA: Manteca County, CA',
       'USA: Stockton, CA', 'USA: Orange County, CA',
       'USA: Riverside County, CA', 'USA: Sacramento County, CA',
       'USA: Pittsburgh, Pennsylvania', 'USA: Georgia', 'USA: Northeast',
       'USA:Maine', 'USA: Kansas', 'USA: Massachusetts', 'USA: Michigan',
       'USA: West', 'USA: Sacramento, CA', 'USA: Houston, TX',
       'USA: Ipswich, MA', 'USA: UC Davis Medical Center',
       'USA: UC Davis Medical Center, Davis, Ca', 'USA: AK',
       'USA: Midwest', 'USA: Southeast', 'USA: MA', 'USA: Maywood, IL',
       'USA:University of California at Los Angeles',
       'USA: Oklahoma City, OK', 'USA: Iowa',
       'USA: San Diego VA Medical Center', 'USA: Arkansas', 'USA: Oregon',
       'USA: South Dakota', 'USA: CT', 'USA: Seattle, WA',
       'USA: Pacific Northwest', 'USA:San Francisco',
       'USA:Minneapolis, MN', 'USA:Sacramento CA', 'USA: San Francisco',
       'USA:Michigan', 'USA:Maryland', 'USA: Nevada', 'USA: ND',
       'USA: North Dakota', 'USA: Delaware', 'USA: New York City, NY',
       'USA: Washington,Seattle', 'USA: Indiana', 'USA: Dallas, Texas',
       'USA: Alabama', 'USA: Maine', 'USA:California',
       'USA: Washington DC', 'USA: Burlington', 'USA:Seattle WA',
       'USA:New York City', 'USA: Aurora, CO', 'USA: Salinas, California',
       'USA:Nebraska', 'USA:Ohio', 'USA:Louisiana', 'USA:Cincinnati',
       'USA: Santa Clara, CA',
       'USA: Kaiser Permanente Washington Capitol Hill Urgent Care Clinic',
       'USA: Minneapolis VA Medical Center',
       'USA: Arkansas commercial broiler farm', 'USA: New Hampshire',
       'USA: Connecticut', 'USA: Southeast Michigan', 'USA: Rhode Island',
       'USA: Merced, CA', 'USA: New Mexico', 'USA: Oklahoma City',
       'USA:Ann Arbor MI', 'USA:Jackson MS', 'USA:Dallas TX',
       'USA: NE/KS', 'USA:Missouri,Saint Louis',
       'USA: Monterey, California', 'USA: SOUTH DAKOTA', 'USA:NM',
       'USA: upper Midwest', 'USA: Pennsylvania, Pittsburgh',
       'USA:Minneapolis', 'USA: Seattle', 'USA: District of Columbia',
       'USA: NM', 'USA: Kentucky', 'USA: San Antonio', 'USA: SD',
       'USA:Colorado', 'USA:Texas', 'USA: LA',
       'USA: Missouri, Kansas City', 'USA:Seattle, WA', 'USA:Latah, ID',
       'USA: Houston', 'USA:New York', 'USA: Mississippi river watershed',
       'USA:West', 'USA:Idaho', 'USA:Pennsylvania',
       'USA: Gainesville, Florida', 'USA:Chicago', 'USA:Arizona',
       'USA: Texas Panhandle', 'USA:Walla Walla,WA', 'USA:NV',
       'USA: Davidson County, Tennessee',
       'USA:Ronald Reagan UCLA Medical Center', 'USA: St Louis',
       'USA: Milwaukee', 'USA: Palo Alto, California', 'USA: New Jersey',
       'USA: Santa Clara', 'USA: DC', 'USA:Fort Sam Houston',
       'USA: Alaska, Kasilof River', 'USA: Pittsburgh',
       'USA: Western States', 'USA:Robert Wood Johnson', 'USA:Houston',
       'USA: Baltimore, MD', 'USA: Caguas', 'USA: Nashville, TN',
       'USA: Boston, MA', 'USA: Arizona', 'USA:Atlanta', 'USA:Arkansas',
       'USA: Lafayette', 'USA:Utah', 'USA: Fort Sam, Houston',
       'USA: Curry County, Oregon', 'USA: Alaska, Anchorage mudflats'

    for loc in locations:
        loc = loc.split(':')
        pdb.set_trace()




#####MAIN#####
plt.rcParams.update({'font.size': 6})
args = parser.parse_args()
amr_data = pd.read_csv(args.tsv[0], sep='\t')
outdir = args.outdir[0]

map_changes(amr_data, outdir)
