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


def find_barcode_match(sequence, barcode_set, barcode_length):
    """
    Return the first exact barcode match found in a sliding window search.
    Using set membership is much faster than scanning every whitelist barcode
    against every substring.
    """
    if len(sequence) < barcode_length:
        return ""

    max_index = len(sequence) - barcode_length + 1
    for index in range(max_index):
        candidate = sequence[index:index + barcode_length]
        if candidate in barcode_set:
            return candidate
    return ""


def load_barcodes(filename, reverse_complement=False):
    """
    Load barcode sequences from a file. If reverse_complement is True, return
    their reverse complements.
    """
    barcode_set = set()
    with open(filename, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "BC sequence" not in reader.fieldnames:
            raise ValueError("Missing 'BC sequence' column in {0}".format(filename))

        for row in reader:
            sequence = str(row.get("BC sequence", "")).strip()
            if not sequence:
                continue
            if reverse_complement:
                sequence = reverse_complement_seq(sequence)
            barcode_set.add(sequence)

    return barcode_set


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


def iter_barcode_hits(cell_sequence_counts, bc14_patterns, bc30_patterns, umi_cutoff):
    for cell in sorted(cell_sequence_counts):
        sequence_counts = cell_sequence_counts[cell]
        ranked_pairs = sorted(
            sequence_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )
        for sequence, umi in filter_sequences_by_umi(ranked_pairs, umi_cutoff):
            bc14_match = find_barcode_match(sequence, bc14_patterns, 14)
            bc30_match = find_barcode_match(sequence, bc30_patterns, 30)
            if bc14_match or bc30_match:
                yield cell, sequence, bc14_match, bc30_match, umi


def main(args):
    for file_path in [args.cell_umi, args.bc14_file, args.bc30_file, args.whitelist]:
        if not os.path.isfile(file_path):
            raise FileNotFoundError("File not found: {0}".format(file_path))

    whitelist_set = load_whitelist(args.whitelist)
    cell_sequence_counts = load_cell_sequence_counts(args.cell_umi, whitelist_set)

    bc14_patterns = load_barcodes(args.bc14_file, reverse_complement=args.rc)
    bc30_patterns = load_barcodes(args.bc30_file, reverse_complement=args.rc)

    with open(args.output, "w") as output:
        output.write("cell\tR2_sequence\tbc14\tbc30\tumi\n")
        for cell, sequence, bc14_match, bc30_match, umi in iter_barcode_hits(
            cell_sequence_counts,
            bc14_patterns,
            bc30_patterns,
            args.umi_cutoff,
        ):
            output.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                    cell,
                    sequence,
                    bc14_match,
                    bc30_match,
                    umi,
                )
            )

    print("Output written to {0}".format(args.output))


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
