#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 23:21:07 2020

@author: gabrielemilioherreraoropeza
"""

### This script gets all the mutations per nucleotide and generates a txt file

import json

nuc_dct = {}

with open('../data/all.json', 'r') as fileopen:
    mut_dct = json.load(fileopen)

for gene_no in range(len(mut_dct)):
    mutation = mut_dct[gene_no]['genomic_dna_change']
    if 'del' in mutation:
        mutation = mutation.split('del')
        mutation = mutation[1]
        for nuc in mutation:
            if not nuc in nuc_dct:
                nuc_dct[nuc] = {}
            if not 'deletions' in nuc_dct[nuc]:
                nuc_dct[nuc]['deletions'] = []
            if not mut_dct[gene_no]['genomic_dna_change'] in nuc_dct[nuc]['deletions']:
                nuc_dct[nuc]['deletions'].append(mut_dct[gene_no]['genomic_dna_change'])
    elif '>' in mutation:
        mutation = mutation.split('>')
        mutation = mutation[0]
        mutation = mutation[-1]
        for nuc in mutation:
            if not nuc in nuc_dct:
                nuc_dct[nuc] = {}
            if not 'sustitutions' in nuc_dct[nuc]:
                nuc_dct[nuc]['sustitutions'] = []
            if not mut_dct[gene_no]['genomic_dna_change'] in nuc_dct[nuc]['sustitutions']:
                nuc_dct[nuc]['sustitutions'].append(mut_dct[gene_no]['genomic_dna_change'])
    elif 'ins' in mutation:
        mutation = mutation.split('ins')
        mutation = mutation[1]
        for nuc in mutation:
            if not nuc in nuc_dct:
                nuc_dct[nuc] = {}
            if not 'insertions' in nuc_dct[nuc]:
                nuc_dct[nuc]['insertions'] = []
            if not mut_dct[gene_no]['genomic_dna_change'] in nuc_dct[nuc]['insertions']:
                nuc_dct[nuc]['insertions'].append(mut_dct[gene_no]['genomic_dna_change'])

nuc_dct_len = {}

for key in nuc_dct:
    nuc_dct_len[key] = {}
    for subtype in nuc_dct[key]:
        nuc_dct_len[key][subtype] = len(nuc_dct[key][subtype])

with open('../data/correct/no_muts_nuc_genes.txt', 'w') as fileopen:
    for nuc in nuc_dct_len:
        for key in nuc_dct_len[nuc]:
            fileopen.write('\t'.join((nuc, key, str(nuc_dct_len[nuc][key]))))
            fileopen.write('\n')
