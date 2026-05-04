#!/usr/bin/env python3

import argparse
import csv
import os
from collections import defaultdict


def reverse_complement_seq(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return str(seq).translate(complement)[::-1]


def filter_sequences_by_umi(sequence_umi_pairs, umi_cutoff):
    """Keep all pairs when umi_cutoff <= 0, otherwise keep pairs above cutoff."""
    if umi_cutoff is None or umi_cutoff <= 0:
        return sequence_umi_pairs
    return [
        (sequence, umi)
        for sequence, umi in sequence_umi_pairs
        if umi >= umi_cutoff
    ]


def find_sgrna_match(sequence, sgrna_set, sgrna_lengths):
    """
    Return the first exact sgRNA match found in a sliding window search.
    Checks for multiple possible lengths if the reference has them.
    """
    for barcode_length in sgrna_lengths:
        if len(sequence) < barcode_length:
            continue

        max_index = len(sequence) - barcode_length + 1
        for index in range(max_index):
            candidate = sequence[index:index + barcode_length]
            if candidate in sgrna_set:
                return candidate
    return ""


def load_sgrnas(filename, reverse_complement=False):
    """
    Load sgRNA sequences from a file. If reverse_complement is True, return
    their reverse complements. Returns a set of sequences and a set of their lengths.
    """
    sgrna_set = set()
    sgrna_lengths = set()
    
    seq_candidates = ["sgRNA sequence", "sgRNA Sequence", "Sequence", "sgRNA", "sequence", "seq"]

    with open(filename, "r") as handle:
        # Check delimiter (tab or comma)
        sample_line = handle.readline()
        delimiter = "\t" if "\t" in sample_line else ","
        handle.seek(0)
        
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        
        seq_col = next((column for column in seq_candidates if column in fieldnames), None)
        if seq_col is None:
             raise ValueError("Cannot find sgRNA sequence column in {0}. Available: {1}".format(filename, fieldnames))

        for row in reader:
            sequence = str(row.get(seq_col, "")).strip()
            if not sequence:
                continue
            if reverse_complement:
                sequence = reverse_complement_seq(sequence)
            sgrna_set.add(sequence)
            sgrna_lengths.add(len(sequence))

    # sort lengths descending so we match longer ones first if they overlap
    return sgrna_set, sorted(list(sgrna_lengths), reverse=True)


def load_whitelist(whitelist_path):
    whitelist = set()
    with open(whitelist_path, "r") as handle:
        for line in handle:
            barcode = line.strip()
            if barcode:
                whitelist.add(barcode)
    return whitelist


def load_cell_sequence_counts(cell_umi_path, whitelist_set):
    """
    Aggregate one count per surviving cell+UMI best sequence.
    Returns {cell: {sequence: umi_count}}.
    """
    cell_sequence_counts = defaultdict(lambda: defaultdict(int))

    with open(cell_umi_path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cell_umi", "seq"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing columns in {0}: {1}".format(cell_umi_path, sorted(missing)))

        for row in reader:
            cell_umi = str(row["cell_umi"]).strip()
            sequence = str(row["seq"]).strip()

            if len(cell_umi) < 17 or not sequence:
                continue

            cell = cell_umi[:16]
            if cell not in whitelist_set:
                continue

            cell_sequence_counts[cell][sequence] += 1

    return cell_sequence_counts


def iter_sgrna_hits(cell_sequence_counts, sgrna_patterns, sgrna_lengths, umi_cutoff):
    for cell in sorted(cell_sequence_counts):
        sequence_counts = cell_sequence_counts[cell]
        ranked_pairs = sorted(
            sequence_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )
        for sequence, umi in filter_sequences_by_umi(ranked_pairs, umi_cutoff):
            sgrna_match = find_sgrna_match(sequence, sgrna_patterns, sgrna_lengths)
            if sgrna_match:
                yield cell, sequence, sgrna_match, umi


def main(args):
    for file_path in [args.cell_umi, args.sgrna_file, args.whitelist]:
        if not os.path.isfile(file_path):
            raise FileNotFoundError("File not found: {0}".format(file_path))

    whitelist_set = load_whitelist(args.whitelist)
    cell_sequence_counts = load_cell_sequence_counts(args.cell_umi, whitelist_set)

    sgrna_patterns, sgrna_lengths = load_sgrnas(args.sgrna_file, reverse_complement=args.rc)

    with open(args.output, "w") as output:
        output.write("cell\tR2_sequence\tsgrna\tumi\n")
        for cell, sequence, sgrna_match, umi in iter_sgrna_hits(
            cell_sequence_counts,
            sgrna_patterns,
            sgrna_lengths,
            args.umi_cutoff,
        ):
            output.write(
                "{0}\t{1}\t{2}\t{3}\n".format(
                    cell,
                    sequence,
                    sgrna_match,
                    umi,
                )
            )

    print("Output written to {0}".format(args.output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process cell-UMI data and assign sgRNAs.")
    parser.add_argument("--cell_umi", required=True, help="Path to the cell_umi file.")
    parser.add_argument("--sgrna_file", required=True, help="Path to the sgRNA reference file.")
    parser.add_argument("--whitelist", required=True, help="Path to the whitelist file.")
    parser.add_argument("--output", required=True, help="Path to the output file.")
    parser.add_argument(
        "--umi_cutoff",
        type=int,
        default=0,
        help="Minimum per-cell UMI count to keep sgRNA candidates in the exported table. Use 0 to keep all candidates.",
    )
    parser.add_argument("--rc", action="store_true", help="Apply reverse and complementary transformation to sgRNAs.")
    args = parser.parse_args()
    main(args)
