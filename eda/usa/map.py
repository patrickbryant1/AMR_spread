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
       'USA:NJ':'NJ', 'USA:SC':'SC', 'USA:NE':'NE', 'USA:WI':'WI', 'USA: NY':'NY', 'USA: NE':'NE',
       'USA: NJ':'NJ', 'USA:PA':'PA', 'USA: FL':'FL', 'USA:IA':'IA', 'USA:AZ':'AZ', 'USA:TN':'TN',
       'USA:KS':'KS', 'USA:RI':'RI', 'USA:KY':'KY', 'USA: OR':'OR', 'USA: WI':'WI', 'USA: CA':'CA',
       'USA: MN':'MN', 'USA: WV':'WV', 'USA: MO':'MO', 'USA: MI':'MI', 'USA: GA':'GA', 'USA: TN':'TN',
       'USA: ID':'ID', 'USA: TX':'TX', 'USA: OK':'OK', 'USA: UT':'UT', 'USA: VT':'VT', 'USA: IA':'IA',
       'USA: KS':'KS', 'USA: SC':'SC', 'USA: KY':'KY', 'USA:NC':'NC', 'USA:CA':'CA', 'USA: NC':'NC',
       'USA:NY':'NY', 'USA:OH':'OH', 'USA: Rochester, NY':'NY',
       'USA: Imperial County, CA':'CA', 'USA:IL':'IL', 'USA:SD':'SD', 'USA:ME':'MD', 'USA:WA':'WA',
       'USA: Cambridge, MA, MIT':'MA', 'USA:MA':'MA', 'USA:ID':'ID',
       'USA: University of Florida':'FL', 'USA:MN':'MN', 'USA:AR':'AR', 'USA:CO':'CO',
       'USA:NH':'NH', 'USA:UT':'UT', 'USA:AL':'AL', 'USA:LA':'LA', 'USA: OH':'OH', 'USA: WA':'WA',
       'USA:AK':'AK', 'USA:VT':'VT', 'USA:OR':'OR', 'USA:DE':'DE', 'USA:ND':'ND',
       'USA:MT':'MT', 'USA:WV':'WV', 'USA:Boston':'MA', 'USA:PR':'PR', 'USA:MS':'MS',
       'USA: Nebraska':'NE', 'USA: Philadelphia':'PA', 'USA:OK':'OK', 'USA: California':'CA',
       'USA:CT':'CT', 'USA: Alaska':'AK', 'USA: Colorado':'CO', 'USA: New York':'NY',
       'USA: Mississippi':'MS', 'USA: Washington':'WA', 'USA:Washington DC':'DC',
       'USA: VA':'VA', 'USA: AZ':'AZ', 'USA: MS':'MS', 'USA:Illinois':'IL', 'USA:Florida':'FL',
       'USA: AR':'AR', 'USA: WA, Seattle':'WA', 'USA: DE':'DE',
       'USA: Arkansas Broiler farm':'AK', 'USA: San Francisco, CA':'CA',
       'USA: Placer County, CA':'CA', 'USA: San Diego County, CA':'CA', 'USA:DC':'DC',
       'USA: Boston':'MA', 'USA: New Hampshire,Rindge':'NH', 'USA:TN-Nashville':'TN',
       'USA:CA-UC Davis':'CA', 'USA:TN-Knoxville':'TN',
       'USA: Saint Louis, Missouri':'MO', 'USA: Alachua, FL':'FL', 'USA: Alachua':'FL',
       'USA: St. Louis Clyde watershed of Lake Superior':'MN',
       'USA: Kaiser Permanente Washington Tacoma lab':'WA',
       "USA: Seattle Children's Hospital":'WA', 'USA: Minneapolis':'MN',
       'USA: Pennsylvania':'PA', 'USA: Missouri':'MO', 'USA: Virginia':'VA',
       'USA: Illinois':'IL', 'USA: Wisconsin':'WI',
       'USA: Santa Barbara, UCSB animal house':'CA',
       'USA: Soldotna landfill, Alaska':'AK', 'USA: Soldotna, Alaska':'AK',
       'USA: Lower Kenai River, Alaska':'AK',
       'USA:Washington, DC':'DC', 'USA:AZ, Flagstaff':'AZ', 'USA: North Carolina':'NC',
       'USA: IL':'IL', 'USA: Alachua County, Florida':'FL', 'USA: Cambridge, MA':'MA',
       'USA: ME':'ME', 'USA: Massachusetts, Boston':'MA', 'USA: Los Angeles':'CA',
       'USA:Tennessee':'TN', 'USA: Minnesota':'MN', 'USA: Cambridge':'MA', 'USA: NH':'NH',
       'USA: Florida':'FL', 'USA: Alameda County, CA':'CA',
       'USA: San Joaquin County, CA':'CA', 'USA: Santa Clara County, CA':'CA',
       'USA: Berkley, CA':'CA', 'USA: Ohio':'OH',
       'USA:University of Illinois Medical Center':'IL', 'USA: St. Louis, MO':'MO',
       'USA:Kansas':'KS', 'USA: University of Washington Medical Center':'WA',
       'USA: Seattle Harborview Medical Center':'WA', 'USA: Maryland':'MD',
       'USA: Texas':'TX', 'USA: Upper Kenai River, Alaska':'AK',
       'USA: Lower Kasilof River, Alaska':'AK', 'USA:Virginia':'VA',
       'USA: Oklahoma':'OK', 'USA: AL':'AL', 'USA: Walnut Creek':'CA',
       'USA: Gainesville, FL':'FL', 'USA:Gainesville, Florida':'FL', 'USA:HI':'HI',
       'USA: San Mateo County, CA':'CA', 'USA: Manteca County, CA':'CA',
       'USA: Stockton, CA':'CA', 'USA: Orange County, CA':'CA',
       'USA: Riverside County, CA':'CA', 'USA: Sacramento County, CA':'CA',
       'USA: Pittsburgh, Pennsylvania':'PA', 'USA: Georgia':'GA',
       'USA:Maine':'ME', 'USA: Kansas':'KS', 'USA: Massachusetts':'MA', 'USA: Michigan':'MI',
       'USA: Sacramento, CA':'CA', 'USA: Houston, TX':'TX',
       'USA: Ipswich, MA':'MA', 'USA: UC Davis Medical Center':'CA',
       'USA: UC Davis Medical Center, Davis, Ca':'CA', 'USA: AK':'AK',
       'USA: MA':'MA', 'USA: Maywood, IL':'IL',
       'USA:University of California at Los Angeles':'CA',
       'USA: Oklahoma City, OK':'OK', 'USA: Iowa':'IA',
       'USA: San Diego VA Medical Center':'CA', 'USA: Arkansas':'AR', 'USA: Oregon':'OR',
       'USA: South Dakota':'SD', 'USA: CT':'CT', 'USA: Seattle, WA':'WA',
       'USA:San Francisco':'CA',
       'USA:Minneapolis, MN':'MN', 'USA:Sacramento CA':'CA', 'USA: San Francisco':'CA',
       'USA:Michigan':'MI', 'USA:Maryland':'MD', 'USA: Nevada':'NV', 'USA: ND':'ND',
       'USA: North Dakota':'ND', 'USA: Delaware':'DE', 'USA: New York City, NY':'NY',
       'USA: Washington,Seattle':'WA', 'USA: Indiana':'IN', 'USA: Dallas, Texas':'TX',
       'USA: Alabama':'AL', 'USA: Maine':'ME', 'USA:California':'CA',
       'USA: Washington DC':'DC', 'USA: Burlington':'VT', 'USA:Seattle WA':'WA',
       'USA:New York City':'NY', 'USA: Aurora, CO':'CO', 'USA: Salinas, California':'CA',
       'USA:Nebraska':'NE', 'USA:Ohio':'OH', 'USA:Louisiana':'LA', 'USA:Cincinnati':'OH',
       'USA: Santa Clara, CA':'CA',
       'USA: Kaiser Permanente Washington Capitol Hill Urgent Care Clinic':'DC',
       'USA: Minneapolis VA Medical Center':'MN',
       'USA: Arkansas commercial broiler farm':'AR', 'USA: New Hampshire':'NH',
       'USA: Connecticut':'CT', 'USA: Southeast Michigan':'MI', 'USA: Rhode Island':'RI',
       'USA: Merced, CA':'CA', 'USA: New Mexico':'NM', 'USA: Oklahoma City':'OK',
       'USA:Ann Arbor MI':'MI', 'USA:Jackson MS':'MS', 'USA:Dallas TX':'TX',
       'USA:Missouri,Saint Louis':'MO',
       'USA: Monterey, California':'CA', 'USA: SOUTH DAKOTA':'SD', 'USA:NM':'NM',
       'USA: Pennsylvania, Pittsburgh':'PA',
       'USA:Minneapolis':'MN', 'USA: Seattle':'WA', 'USA: District of Columbia':'DC',
       'USA: NM':'NM', 'USA: Kentucky':'KY', 'USA: San Antonio':'TX', 'USA: SD':'SD',
       'USA:Colorado':'CO', 'USA:Texas':'TX', 'USA: LA':'LA',
       'USA: Missouri, Kansas City':'MO', 'USA:Seattle, WA':'WA', 'USA:Latah, ID':'ID',
       'USA: Houston':'TX', 'USA:New York':'NY', 'USA: Mississippi river watershed':'MS',
        'USA:Idaho':'ID', 'USA:Pennsylvania':'PA',
       'USA: Gainesville, Florida':'FL', 'USA:Chicago':'IL', 'USA:Arizona':'AZ',
       'USA: Texas Panhandle':'TX', 'USA:Walla Walla,WA':'WA', 'USA:NV':'NV',
       'USA: Davidson County, Tennessee':'TN',
       'USA:Ronald Reagan UCLA Medical Center':'CA', 'USA: St Louis':'MO',
       'USA: Milwaukee':'WI', 'USA: Palo Alto, California':'CA', 'USA: New Jersey':'NJ',
       'USA: Santa Clara':'CA', 'USA: DC':'DC', 'USA:Fort Sam Houston':'TX',
       'USA: Alaska, Kasilof River':'AK', 'USA: Pittsburgh':'PA',
       'USA:Robert Wood Johnson':'PA', 'USA:Houston':'TX',
       'USA: Baltimore, MD':'MD', 'USA: Nashville, TN':'TN',
       'USA: Boston, MA':'MA', 'USA: Arizona':'AZ', 'USA:Atlanta':'GA', 'USA:Arkansas':'AR',
       'USA: Lafayette':'GA', 'USA:Utah':'UT', 'USA: Fort Sam, Houston':'TX',
       'USA: Curry County, Oregon':'OR', 'USA: Alaska, Anchorage mudflats':'AK'}

    states = []
    for i in range(len(amr_data)):
        row = amr_data.loc[i]
        loc = row['Location']

        try:
            states.append(state_keys[loc])
        except:
            states.append('NA')
            print(loc)

    #Assign state
    amr_data['State']=states
    pdb.set_trace()




#####MAIN#####
plt.rcParams.update({'font.size': 6})
args = parser.parse_args()
amr_data = pd.read_csv(args.tsv[0], sep='\t')
outdir = args.outdir[0]

map_changes(amr_data, outdir)
