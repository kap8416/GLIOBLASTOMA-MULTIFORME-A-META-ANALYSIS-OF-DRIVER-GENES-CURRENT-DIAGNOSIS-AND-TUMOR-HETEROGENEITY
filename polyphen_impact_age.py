#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 23:49:56 2020

@author: gabrielemilioherreraoropeza
"""

### This script obtains the polyphen impact status of all the genes in GBM by age

import json
import matplotlib.pyplot as plt
import numpy as np
import operator
import pandas as pd

### --- Age Group 10 - 29

'''----------     FIRST PART     ----------'''
gene_dct_cons = {}

with open('../data/10_29.json', 'r') as fileopen:
    mut_dct = json.load(fileopen)

for gene_no in range(len(mut_dct)):
    for element in range(len(mut_dct[gene_no]['consequence'])):
        for transcript in mut_dct[gene_no]['consequence'][element]:
            if mut_dct[gene_no]['consequence'][element]['transcript']['consequence_type'] == 'missense_variant':
                if mut_dct[gene_no]['consequence'][element]['transcript']['is_canonical'] == True:
                    if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] != '':
                        if mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol'] in gene_dct_cons:
                            if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'benign':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'possibly_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'probably_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'unknown':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                        if not mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol'] in gene_dct_cons:
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']] = {}
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'] = []
                            if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'benign':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'possibly_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'probably_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'unknown':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])

'''----------     SECOND PART     ----------'''

genes = []

for gene in gene_dct_cons.keys():
    genes.append(gene)

genes = np.unique(genes)

temporary_dct = {}

for gene in genes:
    if gene in gene_dct_cons:
        temporary_dct[gene] = len(gene_dct_cons[gene]['BE']) + len(gene_dct_cons[gene]['PO']) + len(gene_dct_cons[gene]['PR']) + len(gene_dct_cons[gene]['UN'])

SORT = sorted(temporary_dct.items(), key = operator.itemgetter(1), reverse = True)

'''----------     THIRD PART     ----------'''

genes = []

for group in SORT:
    genes.append(group[0])

dct_gene_len = {}

for gene in genes:
    if gene in gene_dct_cons:
        dct_gene_len[gene] = {}
        for subtype in gene_dct_cons[gene]:
            dct_gene_len[gene][subtype] = len(gene_dct_cons[gene][subtype])

df = pd.DataFrame.from_dict(data = dct_gene_len)
df.to_csv('../data/correct/polyphen_impact_10_29.csv', header = True, index = True)

''' ------------------------------------------------------------------------------------- '''

### --- Age Group 30 - 59

'''----------     FIRST PART     ----------'''
gene_dct_cons = {}

with open('../data/30_59.json', 'r') as fileopen:
    mut_dct = json.load(fileopen)

for gene_no in range(len(mut_dct)):
    for element in range(len(mut_dct[gene_no]['consequence'])):
        for transcript in mut_dct[gene_no]['consequence'][element]:
            if mut_dct[gene_no]['consequence'][element]['transcript']['consequence_type'] == 'missense_variant':
                if mut_dct[gene_no]['consequence'][element]['transcript']['is_canonical'] == True:
                    if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] != '':
                        if mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol'] in gene_dct_cons:
                            if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'benign':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'possibly_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'probably_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'unknown':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                        if not mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol'] in gene_dct_cons:
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']] = {}
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'] = []
                            if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'benign':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'possibly_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'probably_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'unknown':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])

'''----------     SECOND PART     ----------'''

genes = []

for gene in gene_dct_cons.keys():
    genes.append(gene)

genes = np.unique(genes)

temporary_dct = {}

for gene in genes:
    if gene in gene_dct_cons:
        temporary_dct[gene] = len(gene_dct_cons[gene]['BE']) + len(gene_dct_cons[gene]['PO']) + len(gene_dct_cons[gene]['PR']) + len(gene_dct_cons[gene]['UN'])

SORT = sorted(temporary_dct.items(), key = operator.itemgetter(1), reverse = True)

'''----------     THIRD PART     ----------'''

genes = []

for group in SORT:
    genes.append(group[0])

dct_gene_len = {}

for gene in genes:
    if gene in gene_dct_cons:
        dct_gene_len[gene] = {}
        for subtype in gene_dct_cons[gene]:
            dct_gene_len[gene][subtype] = len(gene_dct_cons[gene][subtype])

df = pd.DataFrame.from_dict(data = dct_gene_len)
df.to_csv('../data/correct/polyphen_impact_30_59.csv', header = True, index = True)

''' ------------------------------------------------------------------------------------- '''

### --- Age Group 60 - 89

'''----------     FIRST PART     ----------'''
gene_dct_cons = {}

with open('../data/60_89.json', 'r') as fileopen:
    mut_dct = json.load(fileopen)

for gene_no in range(len(mut_dct)):
    for element in range(len(mut_dct[gene_no]['consequence'])):
        for transcript in mut_dct[gene_no]['consequence'][element]:
            if mut_dct[gene_no]['consequence'][element]['transcript']['consequence_type'] == 'missense_variant':
                if mut_dct[gene_no]['consequence'][element]['transcript']['is_canonical'] == True:
                    if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] != '':
                        if mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol'] in gene_dct_cons:
                            if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'benign':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'possibly_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'probably_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'unknown':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                        if not mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol'] in gene_dct_cons:
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']] = {}
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'] = []
                            gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'] = []
                            if mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'benign':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['BE'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'possibly_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PO'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'probably_damaging':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['PR'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])
                            elif mut_dct[gene_no]['consequence'][element]['transcript']['annotation']['polyphen_impact'] == 'unknown':
                                gene_dct_cons[mut_dct[gene_no]['consequence'][element]['transcript']['gene']['symbol']]['UN'].append(mut_dct[gene_no]['consequence'][element]['transcript']['aa_change'])

'''----------     SECOND PART     ----------'''

genes = []

for gene in gene_dct_cons.keys():
    genes.append(gene)

genes = np.unique(genes)

temporary_dct = {}

for gene in genes:
    if gene in gene_dct_cons:
        temporary_dct[gene] = len(gene_dct_cons[gene]['BE']) + len(gene_dct_cons[gene]['PO']) + len(gene_dct_cons[gene]['PR']) + len(gene_dct_cons[gene]['UN'])

SORT = sorted(temporary_dct.items(), key = operator.itemgetter(1), reverse = True)

'''----------     THIRD PART     ----------'''

genes = []

for group in SORT:
    genes.append(group[0])

dct_gene_len = {}

for gene in genes:
    if gene in gene_dct_cons:
        dct_gene_len[gene] = {}
        for subtype in gene_dct_cons[gene]:
            dct_gene_len[gene][subtype] = len(gene_dct_cons[gene][subtype])

df = pd.DataFrame.from_dict(data = dct_gene_len)
df.to_csv('../data/correct/polyphen_impact_60_89.csv', header = True, index = True)
