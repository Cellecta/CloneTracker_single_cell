#!/usr/bin/env python
# coding: utf-8

import argparse
import os
from collections import Counter

import pandas as pd


def filter_sequences_by_umi(sequence_list, umi_list, umi_cutoff=0):
    """
    Keep all sequences when umi_cutoff <= 0, otherwise keep sequences whose
    per-cell UMI count is >= umi_cutoff. umi_list is expected to be sorted in
    descending order alongside sequence_list.
    """
    if umi_cutoff is None or umi_cutoff <= 0:
        return sequence_list, umi_list

    filtered_pairs = [
        (sequence, umi)
        for sequence, umi in zip(sequence_list, umi_list)
        if umi >= umi_cutoff
    ]

    if not filtered_pairs:
        return [], []

    filtered_seq_list, filtered_umi_list = zip(*filtered_pairs)
    return list(filtered_seq_list), list(filtered_umi_list)


def find_barcode_match(sequence, barcode_set, barcode_length):
    """
    Return the first exact barcode match found in a sliding window search.
    Using set membership is much faster than scanning every whitelist barcode
    against every substring.
    """
    if len(sequence) < barcode_length:
        return ""

    for index in range(len(sequence) - barcode_length + 1):
        candidate = sequence[index:index + barcode_length]
        if candidate in barcode_set:
            return candidate
    return ""


def bc_search(sequence_list, umi_list, bc14_pattern_set, bc30_pattern_set, umi_cutoff=0):
    """
    Search for barcodes (bc14 and bc30) in the sequence list and return matched
    barcodes with per-sequence UMI counts.
    """
    filtered_seq_list, filtered_umi_list = filter_sequences_by_umi(
        sequence_list,
        umi_list,
        umi_cutoff=umi_cutoff,
    )

    for sequence, umi in zip(filtered_seq_list, filtered_umi_list):
        bc14_match = find_barcode_match(sequence, bc14_pattern_set, 14)
        bc30_match = find_barcode_match(sequence, bc30_pattern_set, 30)
        if bc14_match or bc30_match:
            yield sequence, bc14_match, bc30_match, umi


def load_barcodes(filename, reverse_complement=False):
    """
    Load barcode sequences from a file. If reverse_complement is True, return
    their reverse complements.
    """
    barcodes = pd.read_csv(filename, sep='\t')['BC sequence']
    if reverse_complement:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return {
            "".join(complement.get(base, base) for base in reversed(seq))
            for seq in barcodes
        }
    return set(barcodes)


def main(args):
    for file_path in [args.cell_umi, args.bc14_file, args.bc30_file]:
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

    data = pd.read_csv(args.cell_umi, sep='\t')
    data['cell'] = data['cell_umi'].str[:16]
    data['umi'] = data['cell_umi'].str[16:]

    barcodes = pd.read_csv(args.whitelist, sep='\t', header=None)
    data = data[data['cell'].isin(barcodes[0])]

    data_seq = data.groupby('cell')['seq'].apply(list).reset_index(name='seq_info')
    data_seq['number'] = data_seq['seq_info'].apply(lambda x: Counter(x))

    bc14_patterns = load_barcodes(args.bc14_file, reverse_complement=args.rc)
    bc30_patterns = load_barcodes(args.bc30_file, reverse_complement=args.rc)

    data_seq['bc_search'] = data_seq['number'].apply(
        lambda counter: bc_search(
            [key for key, _ in counter.most_common()],
            [count for _, count in counter.most_common()],
            bc14_patterns,
            bc30_patterns,
            umi_cutoff=args.umi_cutoff,
        )
    )

    with open(args.output, 'w') as output:
        output.write('cell\tR2_sequence\tbc14\tbc30\tumi\n')
        for cell, search_results in zip(data_seq['cell'], data_seq['bc_search']):
            for sequence, bc14_match, bc30_match, umi in search_results:
                output.write(f"{cell}\t{sequence}\t{bc14_match}\t{bc30_match}\t{umi}\n")

    print(f"Output written to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process cell-UMI data and assign barcodes.")
    parser.add_argument("--cell_umi", required=True, help="Path to the cell_umi file.")
    parser.add_argument("--bc14_file", required=True, help="Path to the BC14 barcode file.")
    parser.add_argument("--bc30_file", required=True, help="Path to the BC30 barcode file.")
    parser.add_argument("--whitelist", required=True, help="Path to the whitelist file.")
    parser.add_argument("--output", required=True, help="Path to the output file.")
    parser.add_argument(
        "--umi_cutoff",
        type=int,
        default=0,
        help="Minimum per-cell UMI count to keep barcode candidates in the exported table. Use 0 to keep all candidates.",
    )
    parser.add_argument("--rc", action="store_true", help="Apply reverse and complementary transformation to barcodes.")
    args = parser.parse_args()
    main(args)
