#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import argparse
import os
import re


def load_barcode_file(filename, reverse_complement=False):
    """
    Load barcode sequences and IDs from a file. Optionally reverse complement sequences.
    """
    data = pd.read_csv(filename, sep='\t')
    sequences = data['BC sequence'] if 'BC14 sequence' in data else data['BC sequence']
    ids = data['#BC14 ID'] if '#BC14 ID' in data else data['#BC30 ID']

    if reverse_complement:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        sequences = ["".join(complement.get(base, base) for base in reversed(seq)) for seq in sequences]

    return dict(zip(sequences, ids))


def bc_type(bc14_str, bc30_str):
    """Determine barcode type based on presence of bc14 and bc30."""
    if bc14_str == 0:
        return "bc14_only"
    if bc30_str == 0:
        return "bc30_only"
    return 'bc14_bc30'


def barcode_name(bc14_str, bc30_str, bc14_dict, bc30_dict, R2_seq):
    """Generate barcode name based on bc14 and bc30 matches."""
    if bc14_str != 0 and bc30_str != 0:
        return f"{bc14_dict.get(bc14_str, '')}_{bc30_dict.get(bc30_str, '')}"
    if bc14_str == 0:
        return f"{R2_seq[25:39]}_{bc30_dict.get(bc30_str, '')}"
    return f"{bc14_dict.get(bc14_str, '')}_{R2_seq[-30:]}"


def final_assigned_barcode_func(barcode_list, umi_count_list):
    """Assign the final barcode based on UMI counts."""
    if max(umi_count_list) < 5:
        return 'NA'
    if max(umi_count_list) >= (sum(umi_count_list) / 2) and len(umi_count_list) > 1:
        return barcode_list[0]
    if len(umi_count_list) > 1 and (umi_count_list[1] / umi_count_list[0]) >= 0.7 and umi_count_list[0] >= 5:
        return 'multi_barcode'
    if umi_count_list[0] >= 5:
        return barcode_list[0]
    return 'NA'


def final_assigned_type(final_assigned_barcode_str):
    """Determine the type of barcode assignment."""
    pattern = re.compile(r'^bc14-(\d+)_bc30-(\d+)')
    if final_assigned_barcode_str == "NA":
        return "Undetermined due to low UMI"
    if final_assigned_barcode_str == "multi_barcode":
        return 'multi_barcode'
    if pattern.match(final_assigned_barcode_str):
        return 'One Barcode'
    return 'One Barcode with mutant'


def final_barcode_distribution(bc_type_list, output_path):
    """Generate and save barcode distribution pie chart."""
    counts = Counter(bc_type_list)
    labels, sizes = list(counts.keys()), list(counts.values())

    barcode_stat = pd.DataFrame({'Type': labels, 'Count': sizes})
    barcode_stat.to_csv(output_path, sep='\t', index=False)

    plt.figure(figsize=(12, 6))
    plt.pie(sizes, labels=[''] * len(labels), startangle=90)
    plt.title('Barcode Type Distribution')
    plt.legend([f'{label}: {count}' for label, count in zip(labels, sizes)], loc='best')
    plt.savefig(output_path.replace('.tsv', '.png'))


def umi_distribution_pie(umi_list, output_path):
    """Generate and save UMI distribution pie chart."""
    bins = np.arange(np.floor(min(umi_list)), np.ceil(max(umi_list)) + 21, 20)
    hist, bin_edges = np.histogram(umi_list, bins=bins)

    labels = [f'{int(bin_edges[i])}-{int(bin_edges[i + 1])}' for i in range(len(bin_edges) - 1)]
    plt.figure(figsize=(12, 6))
    plt.pie(hist, labels=None, startangle=90, counterclock=False)
    legend_labels = [f'{label}: {hist[i]}' for i, label in enumerate(labels)]
    plt.legend(legend_labels, title='UMI bins', loc='center left', bbox_to_anchor=(0.8, 0.5))
    plt.title('UMI Distribution')
    plt.savefig(output_path)


def main(args):
    # Load input files
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    data = pd.read_csv(args.input, sep='\t').fillna(0)
    bc14_dict = load_barcode_file(args.bc14, reverse_complement=args.rc)
    bc30_dict = load_barcode_file(args.bc30, reverse_complement=args.rc)

    # Process data
    data['bc_type'] = data.apply(lambda x: bc_type(x['bc14'], x['bc30']), axis=1)
    data_sorted = data.sort_values(by=['cell', 'umi'], ascending=[False, False])
    data_sorted['barcode_name'] = data_sorted.apply(
        lambda x: barcode_name(x['bc14'], x['bc30'], bc14_dict, bc30_dict, x['R2_sequence']), axis=1
    )

    # Group and assign final barcodes
    data_umi = data_sorted.groupby('cell')['umi'].apply(list).reset_index(name='umi')
    data_barcode = data_sorted.groupby('cell')['barcode_name'].apply(list).reset_index(name='barcode')
    data_final = data_umi.merge(data_barcode, on='cell', how='left')

    data_final['final_assigned_barcode'] = data_final.apply(
        lambda x: final_assigned_barcode_func(x['barcode'], x['umi']), axis=1
    )
    data_final['umi_count'] = data_final['umi'].apply(len)
    data_final['barcode_type'] = data_final['final_assigned_barcode'].apply(final_assigned_type)

    # Generate visualizations
    final_barcode_distribution(data_final['barcode_type'], args.barcode_stat)
    umi_distribution_pie(data_final['umi_count'], args.umi_pie)

    # Save output
    data_final[['cell', 'final_assigned_barcode', 'umi_count', 'barcode_type']].to_csv(args.output, sep='\t', index=False)
    print(f"Analysis complete. Results saved to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process barcode and UMI data.")
    parser.add_argument("--input", required=True, help="Path to input barcode assignment file.")
    parser.add_argument("--bc14", required=True, help="Path to BC14 barcode file.")
    parser.add_argument("--bc30", required=True, help="Path to BC30 barcode file.")
    parser.add_argument("--rc", action="store_true", help="Apply reverse complementary to barcodes.")
    parser.add_argument("--output", required=True, help="Path to output file.")
    parser.add_argument("--barcode_stat", default="barcode_stat.tsv", help="Path to barcode statistics file.")
    parser.add_argument("--umi_pie", default="umi_distribution.png", help="Path to save UMI pie chart.")
    args = parser.parse_args()
    main(args)
