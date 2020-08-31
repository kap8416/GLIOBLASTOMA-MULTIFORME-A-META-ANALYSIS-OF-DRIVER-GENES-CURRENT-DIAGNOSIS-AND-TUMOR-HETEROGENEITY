#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 19:03:49 2020

@author: gabrielemilioherreraoropeza
"""

### This script performs a Mann-Whitney test between the zscores of the 18 driver genes for comparing them GBM age groups by subtypes.

from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats

### --- Making dictionaries with Patients as Keys and Data as Values

onlyfiles = [f for f in listdir('../data/Genes/') if isfile(join('../data/Genes/', f))]

onlyfiles.remove('.DS_Store')

genes = []
Ps = []

for file in onlyfiles:
    file = file.strip('.txt')
    genes.append(file)

genes.sort()

dct = {}

for gene in genes:
    with open('../data/Genes/%s.txt' % gene, 'r') as fileopen:
        fileopen = fileopen.readlines()[3:]
        for line in fileopen:
            line = line.strip('\n')
            line = line.split('\t')
            if not line[0] in dct:
                dct[line[0]] = {}
                dct[line[0]]['genes'] = {}
                dct[line[0]]['classification'] = line[2]
                dct[line[0]]['age'] = line[3]
                dct[line[0]]['days2death'] = line[4]
                dct[line[0]]['mgmt_CH3'] = line[5]
                dct[line[0]]['vital_status'] = line[6]
            if not gene in dct[line[0]]['genes']:
                dct[line[0]]['genes'][gene] = {}
                dct[line[0]]['genes'][gene]['zscore'] = []
            dct[line[0]]['genes'][gene]['zscore'].append(line[1])

gp_ages = [['10', '29'], ['30', '59'], ['60', '89']]
classif = ['N', 'C', 'P', 'M']


for clase in classif:
    dct_data = {}
    YvsA = []
    YvsO = []
    AvsO = []
    for px in dct:
        for gene in genes:
            if not gene in dct_data:
                dct_data[gene] = {}
            if dct[px]['age'] != 'NA':
                    for group in gp_ages:
                        if dct[px]['classification'] == clase:
                            if not str(group) in dct_data[gene]:
                                dct_data[gene][str(group)] = []
                            if int(dct[px]['age']) >= int(group[0]) and int(dct[px]['age']) <= int(group[1]):
                                if len(dct[px]['genes'][gene]['zscore']) == 1:
                                    dct_data[gene][str(group)].append(float(dct[px]['genes'][gene]['zscore'][0]))
                                elif len(dct[px]['genes'][gene]['zscore']) == 2:
                                    dct_data[gene][str(group)].append(float(dct[px]['genes'][gene]['zscore'][0]))
                                    dct_data[gene][str(group)].append(float(dct[px]['genes'][gene]['zscore'][1]))
    for gene in dct_data:
        if dct_data[gene]["['10', '29']"] != [] and dct_data[gene]["['30', '59']"] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]["['10', '29']"], dct_data[gene]["['30', '59']"], alternative = 'two-sided')
            YvsA.append(p)
        elif dct_data[gene]["['10', '29']"] == [] or dct_data[gene]["['30', '59']"] == []:
            YvsA.append('NAN')
        if dct_data[gene]["['10', '29']"] != [] and dct_data[gene]["['60', '89']"] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]["['10', '29']"], dct_data[gene]["['60', '89']"], alternative = 'two-sided')
            YvsO.append(p)
        elif dct_data[gene]["['10', '29']"] == [] or dct_data[gene]["['60', '89']"] == []:
            YvsO.append('NAN')
        if dct_data[gene]["['30', '59']"] != [] and dct_data[gene]["['60', '89']"] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]["['30', '59']"], dct_data[gene]["['60', '89']"], alternative = 'two-sided')
            AvsO.append(p)
        elif dct_data[gene]["['30', '59']"] == [] or dct_data[gene]["['60', '89']"] == []:
            AvsO.append('NAN')
    dct_df = {'gene': genes, 'YvsA': YvsA, 'YvsO': YvsO, 'AvsO': AvsO}
    df = pd.DataFrame.from_dict(data = dct_df)
    df.to_csv('../data/correct/stats_{0}.txt'.format(clase),
              header = True, index = False, sep = '\t')
