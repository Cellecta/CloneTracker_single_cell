#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from collections import Counter
import itertools


def find_last_index_greater_than(sorted_desc_list, target_value):
    """
    Find the last index in a descending list where the value is greater than the target value.
    """
    for i in range(len(sorted_desc_list) - 1, -1, -1):
        if sorted_desc_list[i] > target_value:
            return i
    return -1


def bc_search(sequence_list, umi_list, bc14_pattern_list, bc30_pattern_list, umi_cutoff=5):
    """
    Search for barcodes (bc14 and bc30) in the sequence list and return matched barcodes with UMI counts.
    """
    bc14_final_list = []
    bc30_final_list = []

    pos = find_last_index_greater_than(umi_list, umi_cutoff)
    filtered_seq_list = sequence_list[:pos]

    for seq in filtered_seq_list:
        substrings_bc14 = [''.join(seq[i:i+14]) for i in range(len(seq) - 13)]
        bc14_matches = [bc for bc in bc14_pattern_list if bc in substrings_bc14]
        substrings_bc30 = [''.join(seq[i:i+30]) for i in range(len(seq) - 29)]
        bc30_matches = [bc for bc in bc30_pattern_list if bc in substrings_bc30]
        bc14_final_list.append(bc14_matches)
        bc30_final_list.append(bc30_matches)

    for idx, barcodes in enumerate(zip(bc14_final_list, bc30_final_list)):
        if barcodes[0] or barcodes[1]:
            yield filtered_seq_list[idx], barcodes, umi_list[idx]


def main():
    # Load input data
    data = pd.read_csv("cell_umi.xls", sep='\t')
    data['cell'] = data['cell_umi'].str[:16]
    data['umi'] = data['cell_umi'].str[16:]

    # Load barcode whitelist and filter data
    barcodes = pd.read_csv('AD_barcode_cleaned.tsv', sep='\t', header=None)
    data = data[data['cell'].isin(barcodes[0])]

    # Group sequences by cell and count frequencies
    data_seq = data.groupby('cell')['seq'].apply(list).reset_index(name='seq_info')
    data_seq['number'] = data_seq['seq_info'].apply(lambda x: Counter(x))

    # Load barcode patterns
    #pattern = pd.read_csv('tag.csv', header=None)
    bc14 = list(pd.read_csv('Cellecta-CloneTrackerXP-5M-Pool1-BC14-LNGS-300-Library-Design.txt', sep='\t')['BC14 sequence'])
    bc30 = list(pd.read_csv('Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt', sep='\t')['BC30 sequence'])
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bc14 = ["".join(complement.get(base, base) for base in reversed(seq)) for seq in bc14 ]
    bc30 = ["".join(complement.get(base, base) for base in reversed(seq)) for seq in bc30 ]

    # Search for barcodes
    data_seq['bc_search'] = data_seq['number'].apply(
        lambda x: bc_search(
            [key for key, _ in x.most_common()],
            [count for _, count in x.most_common()],
            bc14,
            bc30,
            umi_cutoff=5  # Default UMI cutoff is now 2
        )
    )

    # Write results to the output file
    with open('barcode_assignment_umi_5.xls', 'w') as output:
        output.write('cell\tR2 sequence\tbc14\tbc30\tumi\n')
        for cell, search_results in zip(data_seq['cell'], data_seq['bc_search']):
            for seq, barcodes, umi in search_results:
                bc14_matches = barcodes[0][0] if barcodes[0] else ""
                bc30_matches = barcodes[1][0] if barcodes[1] else ""
                output.write(f"{cell}\t{seq}\t{bc14_matches}\t{bc30_matches}\t{umi}\n")


if __name__ == "__main__":
    main()

