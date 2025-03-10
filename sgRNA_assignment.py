#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
from collections import Counter
import os

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Process UMI sequences and assign sgRNAs.")
    parser.add_argument("-i", "--input", required=True, help="Input cell UMI file (e.g., cell_umi.xls).")
    parser.add_argument("-o", "--output", required=True, help="Output file for sgRNA assignments (e.g., sgRNA_assignment.csv).")
    parser.add_argument("-t", "--target", required=True, help="Target reference file (e.g., feature_ref.csv).")
    return parser.parse_args()

def find_sgRNA(sequence, sgRNA_dict):
    """
    Find the corresponding sgRNA in the sequence based on the reference dictionary.
    """
    for sgRNA_seq, sgRNA_id in sgRNA_dict.items():
        if sgRNA_seq + 'TTTT' in sequence:
            return sgRNA_id
    return "NA"

def find_most_common_element(seq_list):
    """
    Find the most common sgRNA and its frequency in the list.
    """
    counter = Counter(seq_list)
    most_common_element, frequency = counter.most_common(1)[0]
    return most_common_element, frequency

def find_special_item(counter):
    """
    Identify the dominant sgRNA if it has significantly higher counts than others.
    """
    if len(counter) == 1:
        return next(iter(counter))

    max_item, max_count = counter.most_common(1)[0]
    other_sum = sum(count for item, count in counter.items() if item != max_item)

    return max_item if max_count > other_sum else None

def process_files(input_file, target_file, output_file):
    """ Main processing function """

    # Check if input files exist
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} not found.")
    if not os.path.exists(target_file):
        raise FileNotFoundError(f"Target file {target_file} not found.")

    # Load the input data
    data = pd.read_csv(input_file, sep='\t')
    data['cell'] = data['cell_umi'].str[:16]
    data['umi'] = data['cell_umi'].str[16:]

    # Load sgRNA reference data
    target = pd.read_csv(target_file)
    sgRNA_dict = dict(zip(target['sequence'], target['id']))

    # Assign sgRNA to sequences
    data["sgRNA"] = data['seq'].apply(lambda x: find_sgRNA(x, sgRNA_dict))

    # Filter out unassigned sgRNAs
    data = data[data['sgRNA'] != "NA"]

    # Group by cell barcode and aggregate sgRNA lists
    data_seq = data.groupby('cell')['sgRNA'].apply(list).reset_index(name='seq_info')

    # Extract most common sgRNA and corresponding frequency
    data_seq['seq'] = data_seq['seq_info'].apply(lambda x: find_most_common_element(x)[0])
    data_seq['umi'] = data_seq['seq_info'].apply(lambda x: find_most_common_element(x)[1])
    data_seq['number'] = data_seq['seq_info'].apply(lambda x: len(set(x)))
    data_seq['counter'] = data_seq['seq_info'].apply(lambda x: Counter(x))

    # Assign final sgRNA if one is dominant
    data_seq['sgRNA_assigned'] = data_seq['counter'].apply(lambda x: find_special_item(x))
    data_seq.dropna(inplace=True)

    # Save results to output file
    data_seq.to_csv(output_file, index=False)
    print(f"Processing complete. Results saved to {output_file}")

if __name__ == "__main__":
    args = parse_arguments()
    process_files(args.input, args.target, args.output)
