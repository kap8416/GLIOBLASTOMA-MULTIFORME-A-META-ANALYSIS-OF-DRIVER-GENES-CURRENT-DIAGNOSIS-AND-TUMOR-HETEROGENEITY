#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 00:23:00 2020

@author: gabrielemilioherreraoropeza
"""
### This script performs a Mann-Whitney test between the zscores of the 18 driver genes for comparing them GBM subtypes by age groups.

from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats

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


for group in gp_ages:
    dct_data = {}
    NvsC = []
    NvsP = []
    NvsM = []
    CvsP = []
    CvsM = []
    PvsM = []
    for px in dct:
        for gene in genes:
            if not gene in dct_data:
                dct_data[gene] = {}
            if dct[px]['age'] != 'NA':
                if int(dct[px]['age']) >= int(group[0]) and int(dct[px]['age']) <= int(group[1]):
                    for clase in classif:
                        if not clase in dct_data[gene]:
                            dct_data[gene][clase] = []
                        if dct[px]['classification'] == clase:
                            if len(dct[px]['genes'][gene]['zscore']) == 1:
                                dct_data[gene][clase].append(float(dct[px]['genes'][gene]['zscore'][0]))
                            elif len(dct[px]['genes'][gene]['zscore']) == 2:
                                dct_data[gene][clase].append(float(dct[px]['genes'][gene]['zscore'][0]))
                                dct_data[gene][clase].append(float(dct[px]['genes'][gene]['zscore'][1]))
    for gene in dct_data:
        if dct_data[gene]['N'] != [] and dct_data[gene]['C'] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]['N'], dct_data[gene]['C'], alternative = 'two-sided')
            NvsC.append(p)
        elif dct_data[gene]['N'] == [] or dct_data[gene]['C'] == []:
            NvsC.append('NAN')
        if dct_data[gene]['N'] != [] and dct_data[gene]['P'] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]['N'], dct_data[gene]['P'], alternative = 'two-sided')
            NvsP.append(p)
        elif dct_data[gene]['N'] == [] or dct_data[gene]['P'] == []:
            NvsP.append('NAN')
        if dct_data[gene]['N'] != [] and dct_data[gene]['M'] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]['N'], dct_data[gene]['M'], alternative = 'two-sided')
            NvsM.append(p)
        elif dct_data[gene]['N'] == [] or dct_data[gene]['M'] == []:
            NvsM.append('NAN')
        if dct_data[gene]['C'] != [] and dct_data[gene]['P'] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]['C'], dct_data[gene]['P'], alternative = 'two-sided')
            CvsP.append(p)
        elif dct_data[gene]['C'] == [] or dct_data[gene]['P'] == []:
            CvsP.append('NAN')
        if dct_data[gene]['C'] != [] and dct_data[gene]['M'] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]['C'], dct_data[gene]['M'], alternative = 'two-sided')
            CvsM.append(p)
        elif dct_data[gene]['C'] == [] or dct_data[gene]['M'] == []:
            CvsM.append('NAN')
        if dct_data[gene]['P'] != [] and dct_data[gene]['M'] != []:
            z, p = scipy.stats.mannwhitneyu(dct_data[gene]['P'], dct_data[gene]['M'], alternative = 'two-sided')
            PvsM.append(p)
        elif dct_data[gene]['P'] == [] or dct_data[gene]['M'] == []:
            PvsM.append('NAN')
    dct_df = {'gene': genes, 'NvsC': NvsC, 'NvsP': NvsP, 'NvsM': NvsM, 'CvsP': CvsP, 'CvsM': CvsM,
              'PvsM': PvsM}
    df = pd.DataFrame.from_dict(data = dct_df)
    df.to_csv('../data/correct/stats_{0}-{1}.txt'.format(group[0], group[1]),
              header = True, index = False, sep = '\t')
