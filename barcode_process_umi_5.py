#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from collections import Counter
import argparse
import os


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

    # Filter sequences based on UMI cutoff
    pos = find_last_index_greater_than(umi_list, umi_cutoff)
    filtered_seq_list = sequence_list[:pos + 1] if pos != -1 else []

    for seq in filtered_seq_list:
        substrings_bc14 = [''.join(seq[i:i + 14]) for i in range(len(seq) - 13)]
        bc14_matches = [bc for bc in bc14_pattern_list if bc in substrings_bc14]
        substrings_bc30 = [''.join(seq[i:i + 30]) for i in range(len(seq) - 29)]
        bc30_matches = [bc for bc in bc30_pattern_list if bc in substrings_bc30]
        bc14_final_list.append(bc14_matches)
        bc30_final_list.append(bc30_matches)

    for idx, barcodes in enumerate(zip(bc14_final_list, bc30_final_list)):
        if barcodes[0] or barcodes[1]:
            yield filtered_seq_list[idx], barcodes, umi_list[idx]


def load_barcodes(filename, reverse_complement=False):
    """
    Load barcode sequences from a file. If reverse_complement is True, return their reverse complements.
    """
    barcodes = pd.read_csv(filename, sep='\t')['BC sequence']
    if reverse_complement:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ["".join(complement.get(base, base) for base in reversed(seq)) for seq in barcodes]
    return list(barcodes)


def main(args):
    # Verify input files
    for file_path in [args.cell_umi, args.bc14_file, args.bc30_file]:
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

    # Load cell UMI data
    data = pd.read_csv(args.cell_umi, sep='\t')
    data['cell'] = data['cell_umi'].str[:16]
    data['umi'] = data['cell_umi'].str[16:]

    # Load barcode whitelist and filter data
    barcodes = pd.read_csv(args.whitelist, sep='\t', header=None)
    data = data[data['cell'].isin(barcodes[0])]

    # Group sequences by cell and count frequencies
    data_seq = data.groupby('cell')['seq'].apply(list).reset_index(name='seq_info')
    data_seq['number'] = data_seq['seq_info'].apply(lambda x: Counter(x))

    # Load barcode patterns
    bc14_patterns = load_barcodes(args.bc14_file, reverse_complement=args.rc)
    bc30_patterns = load_barcodes(args.bc30_file, reverse_complement=args.rc)

    # Search for barcodes
    data_seq['bc_search'] = data_seq['number'].apply(
        lambda x: bc_search(
            [key for key, _ in x.most_common()],
            [count for _, count in x.most_common()],
            bc14_patterns,
            bc30_patterns,
            umi_cutoff=args.umi_cutoff
        )
    )

    # Write results to the output file
    with open(args.output, 'w') as output:
        output.write('cell\tR2_sequence\tbc14\tbc30\tumi\n')
        for cell, search_results in zip(data_seq['cell'], data_seq['bc_search']):
            for seq, barcodes, umi in search_results:
                bc14_matches = barcodes[0][0] if barcodes[0] else ""
                bc30_matches = barcodes[1][0] if barcodes[1] else ""
                output.write(f"{cell}\t{seq}\t{bc14_matches}\t{bc30_matches}\t{umi}\n")

    print(f"Output written to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process cell-UMI data and assign barcodes.")
    parser.add_argument("--cell_umi", required=True, help="Path to the cell_umi file.")
    parser.add_argument("--bc14_file", required=True, help="Path to the BC14 barcode file.")
    parser.add_argument("--bc30_file", required=True, help="Path to the BC30 barcode file.")
    parser.add_argument("--whitelist", required=True, help="Path to the whitelist file.")
    parser.add_argument("--output", required=True, help="Path to the output file.")
    parser.add_argument("--umi_cutoff", type=int, default=5, help="UMI cutoff value (default: 5).")
    parser.add_argument("--rc", action="store_true", help="Apply reverse and complementary transformation to barcodes.")
    args = parser.parse_args()
    main(args)

