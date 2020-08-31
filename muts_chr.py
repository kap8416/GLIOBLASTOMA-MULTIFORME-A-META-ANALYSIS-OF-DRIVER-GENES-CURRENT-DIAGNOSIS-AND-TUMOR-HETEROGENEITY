#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 21:12:17 2020

@author: gabrielemilioherreraoropeza
"""

### This script gets all the mutations per chromosome and generates a txt file

import json

chr_dct = {}
files = ['all', '10_29', '30_59', '60_89']

for file in files:
    with open('../data/%s.json' % (file), 'r') as fileopen:
        mut_dct = json.load(fileopen)

    for gene_no in range(len(mut_dct)):
        chromosome = mut_dct[gene_no]['genomic_dna_change']
        chromosome = chromosome.split(':')
        chromosome = chromosome[0]
        if not chromosome in chr_dct:
            chr_dct[chromosome] = {}
        if not mut_dct[gene_no]['mutation_subtype'] in chr_dct[chromosome]:
            chr_dct[chromosome][mut_dct[gene_no]['mutation_subtype']] = []
        if not mut_dct[gene_no]['genomic_dna_change'] in chr_dct[chromosome][mut_dct[gene_no]['mutation_subtype']]:
            chr_dct[chromosome][mut_dct[gene_no]['mutation_subtype']].append(mut_dct[gene_no]['genomic_dna_change'])

    chr_dct_len = {}

    for key in chr_dct:
        chr_dct_len[key] = {}
        for subtype in chr_dct[key]:
            chr_dct_len[key][subtype] = len(chr_dct[key][subtype])

    with open('../data/correct/CHR_no_mutations_%s_NotCONS.txt' % (file), 'w') as fileopen:
        for chromosome in chr_dct_len:
            for key in chr_dct_len[chromosome]:
                fileopen.write('\t'.join((chromosome, key, str(chr_dct_len[chromosome][key]))))
                fileopen.write('\n')
