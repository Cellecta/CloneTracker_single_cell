#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import re

# Load data
data = pd.read_csv('barcode_assignment_umi_5.xls', sep='\t')
data = data.fillna(0)

# Load barcode mappings
# pattern = pd.read_csv('tag.csv', header=None)
# bc14 = pattern[0][:100]
# bc30 = pattern[0][100:]

# bc14_name = pattern[1][:100]
# bc30_name = pattern[1][100:]
    # Load barcode patterns
    #pattern = pd.read_csv('tag.csv', header=None)
bc14 = list(pd.read_csv('Cellecta-CloneTrackerXP-5M-Pool1-BC14-LNGS-300-Library-Design.txt', sep='\t')['BC14 sequence'])
bc30 = list(pd.read_csv('Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt', sep='\t')['BC30 sequence'])
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
bc14 = ["".join(complement.get(base, base) for base in reversed(seq)) for seq in bc14 ]
bc30 = ["".join(complement.get(base, base) for base in reversed(seq)) for seq in bc30 ]
bc14_name = list(pd.read_csv('Cellecta-CloneTrackerXP-5M-Pool1-BC14-LNGS-300-Library-Design.txt', sep='\t')['#BC14 ID'])
bc30_name = list(pd.read_csv('Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt', sep='\t')['#BC30 ID'])
bc14_dict = dict(zip(bc14, bc14_name))
bc30_dict = dict(zip(bc30, bc30_name))

# Functions
def bc_type(bc14_str, bc30_str):
    if bc14_str == 0:
        return "bc14_only"
    if bc30_str == 0:
        return "bc30_only"
    return 'bc14_bc30'

def barcode_name(bc14_str, bc30_str, bc14_dict, bc30_dict, R2_seq):
    if bc14_str != 0 and bc30_str != 0:
        return bc14_dict.get(bc14_str) + "_" + bc30_dict.get(bc30_str)
    if bc14_str == 0:
        return R2_seq[25:39] + "_" + bc30_dict.get(bc30_str)
    return bc14_dict.get(bc14_str) + "_" + R2_seq[-30:]

def final_assigned_barcode_func(barcode_list, umi_count_list):
    if max(umi_count_list) < 5:
        return 'NA'
    if max(umi_count_list) >= 5 and max(umi_count_list) >= (sum(umi_count_list) / 2) and len(umi_count_list) > 1:
        return barcode_list[0]
    if len(umi_count_list) > 1 and (umi_count_list[1] / umi_count_list[0]) >= 0.7 and umi_count_list[0] >= 5:
        return 'multi_barcode'
    if len(umi_count_list) > 1 and (umi_count_list[1] / umi_count_list[0]) < 0.7 and umi_count_list[0] >= 5:
        return barcode_list[0]
    if len(umi_count_list) == 1 and umi_count_list[0] >= 5:
        return barcode_list[0]
    return 'NA'

def final_assigned_type(final_assigned_barcode_str):
    pattern = re.compile(r'^bc14-(\d+)_bc30-(\d+)')
    if final_assigned_barcode_str == "NA":
        return "Undetermined due to low UMI"
    if final_assigned_barcode_str == "multi_barcode":
        return 'multi_barcode'
    if pattern.match(final_assigned_barcode_str):
        return 'One Barcode'
    return 'One Barcode with mutant'

def final_barcode_distribution(bc_type_list):
    counts = Counter(bc_type_list)
    labels, sizes = list(counts.keys()), list(counts.values())
    
    barcode_stat = pd.DataFrame({'Type': labels, 'Count': sizes})
    barcode_stat.to_csv('barcode_stat.xls', sep='\t', index=False)

    plt.figure(figsize=(12, 6))
    plt.pie(sizes, labels=[''] * len(labels), startangle=90)
    plt.title('Barcode Type Distribution')
    plt.legend([f'{label}: {count}' for label, count in zip(labels, sizes)], loc='best')
    plt.savefig('barcode_assignment.png')

def umi_distribution_pie(umi_list):
    bins = np.arange(np.floor(min(umi_list)), np.ceil(max(umi_list)) + 21, 20)
    hist, bin_edges = np.histogram(umi_list, bins=bins)
    percentages = (hist / sum(hist)) * 100

    labels = [f'{int(bin_edges[i])}-{int(bin_edges[i+1])}' for i in range(len(bin_edges) - 1)]
    plt.figure(figsize=(12, 6))
    plt.pie(hist, labels=None, startangle=90, counterclock=False)
    legend_labels = [f'{label}: {hist[i]}' for i, label in enumerate(labels)]
    plt.legend(legend_labels, title='UMI bins', loc='center left', bbox_to_anchor=(0.8, 0.5))
    plt.title('UMI Distribution')
    plt.savefig('umi_distribution.png')

# Processing
data['bc_type'] = data.apply(lambda x: bc_type(x['bc14'], x['bc30']), axis=1)
data_sorted = data.sort_values(by=['cell', 'umi'], ascending=[False, False])
data_sorted['barcode_name'] = data_sorted.apply(lambda x: barcode_name(x['bc14'], x['bc30'], bc14_dict, bc30_dict, x['R2 sequence']), axis=1)

data_umi = data_sorted.groupby('cell')['umi'].apply(list).reset_index(name='umi')
data_barcode = data_sorted.groupby('cell')['barcode_name'].apply(list).reset_index(name='barcode')
data_final = data_umi.merge(data_barcode, on='cell', how='left')

data_final['final_assigned_barcode'] = data_final.apply(lambda x: final_assigned_barcode_func(x['barcode'], x['umi']), axis=1)
data_final['umi_count'] = data_final['umi'].apply(lambda x: x[0])
data_final['barcode_type'] = data_final['final_assigned_barcode'].apply(lambda x: final_assigned_type(x))

# Visualization
final_barcode_distribution(data_final['barcode_type'])
umi_distribution_pie(data_final['umi_count'])

# Save final data
data_final[['cell', 'final_assigned_barcode', 'umi_count', 'barcode_type']].to_csv('barcode_assignment_summary.xls', sep='\t', index=False)
