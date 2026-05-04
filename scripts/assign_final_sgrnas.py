#!/usr/bin/env python3
# coding: utf-8

import argparse
import csv
import os
import re
from collections import Counter, defaultdict

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def reverse_complement_seq(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return str(seq).translate(complement)[::-1]


def load_sgrna_file(filename, reverse_complement=False):
    seq_candidates = ["sgRNA sequence", "sgRNA Sequence", "Sequence", "sgRNA", "sequence", "seq"]
    id_candidates = ["sgRNA ID", "sgRNA_ID", "ID", "id"]

    with open(filename, "r") as handle:
        sample_line = handle.readline()
        delimiter = "\t" if "\t" in sample_line else ","
        handle.seek(0)
        
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = reader.fieldnames or []

        seq_col = next((column for column in seq_candidates if column in fieldnames), None)
        id_col = next((column for column in id_candidates if column in fieldnames), None)

        if seq_col is None:
            raise ValueError("Cannot find sgRNA sequence column in {0}".format(filename))

        if id_col is None:
            raise ValueError("Cannot find sgRNA ID column in {0}".format(filename))

        result = {}
        for row in reader:
            sequence = str(row.get(seq_col, "")).strip()
            sgrna_id = str(row.get(id_col, "")).strip()
            if not sequence:
                continue
            if reverse_complement:
                sequence = reverse_complement_seq(sequence)
            result[sequence] = sgrna_id
    return result


def is_missing(value):
    return str(value).strip() in ["0", "", "NA", "None"]


def final_assigned_sgrna_func(sgrna_list, umi_count_list, min_total_umi=3, min_top_umi=3):
    if not umi_count_list:
        return "NA"

    total_umi = sum(umi_count_list)
    top_umi = umi_count_list[0]

    if total_umi < min_total_umi or top_umi < min_top_umi:
        return "NA"

    valid_sgrnas = [sgrna for sgrna in sgrna_list if sgrna != "NA"]
    if not valid_sgrnas:
        return "NA"

    if len(umi_count_list) == 1:
        return sgrna_list[0]

    if top_umi >= (total_umi / 2.0):
        return sgrna_list[0]

    if umi_count_list[1] / float(top_umi) >= 0.7:
        return "multi_sgrna"

    return sgrna_list[0]


def final_assigned_type(final_assigned_sgrna_str):
    value = str(final_assigned_sgrna_str).strip()

    if value == "NA":
        return "Undetermined due to low UMI"

    if value == "multi_sgrna":
        return "multi_sgrna"

    return "One sgRNA"


def final_sgrna_distribution(sgrna_type_list, output_path):
    counts = Counter(sgrna_type_list)
    labels = list(counts.keys())
    sizes = list(counts.values())

    with open(output_path, "w") as handle:
        handle.write("Type\tCount\n")
        for label, size in zip(labels, sizes):
            handle.write("{0}\t{1}\n".format(label, size))

    if plt is None or not sizes:
        return

    plt.figure(figsize=(10, 6))
    plt.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=90)
    plt.title("sgRNA Type Distribution")
    plt.tight_layout()
    plt.savefig(output_path.replace(".tsv", ".png"), dpi=300)
    plt.close()


def umi_distribution_pie(umi_list, output_path):
    umi_list = [value for value in umi_list if value is not None]

    if plt is None or not umi_list:
        return

    min_umi = min(umi_list) if umi_list else 0
    max_umi = max(umi_list) if umi_list else 0
    start = (min_umi // 20) * 20
    end = ((max_umi // 20) + 1) * 20
    bucket_counts = Counter()
    for value in umi_list:
        bucket_start = ((value - start) // 20) * 20 + start
        bucket_end = bucket_start + 20
        bucket_counts[(bucket_start, bucket_end)] += 1

    nonzero = [
        (count, bucket_start, bucket_end)
        for (bucket_start, bucket_end), count in sorted(bucket_counts.items())
        if count > 0
    ]
    if not nonzero:
        return

    labels = ["{0}-{1}".format(int(start), int(end)) for _, start, end in nonzero]
    sizes = [count for count, _, _ in nonzero]

    plt.figure(figsize=(10, 6))
    plt.pie(sizes, labels=labels, startangle=90)
    plt.title("UMI Distribution")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def format_list(values):
    if isinstance(values, list):
        return ";".join(str(value) for value in values)
    return "NA"


def load_assignments(input_path, sgrna_dict):
    grouped = defaultdict(list)

    with open(input_path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_cols = ["cell", "umi", "sgrna", "R2_sequence"]
        missing_cols = [column for column in required_cols if column not in (reader.fieldnames or [])]
        if missing_cols:
            raise ValueError("Missing columns: {0}".format(missing_cols))

        for row in reader:
            cell = str(row["cell"]).strip()
            if not cell:
                continue

            umi = int(row["umi"])
            sgrna_seq = str(row["sgrna"]).strip()
            
            sgrna_id = sgrna_dict.get(sgrna_seq, sgrna_seq)
            if is_missing(sgrna_seq):
                sgrna_id = "NA"
                
            grouped[cell].append((umi, sgrna_id))

    return grouped


def summarize_cells(grouped, min_total_umi, min_top_umi):
    rows = []
    for cell in sorted(grouped):
        ranked = sorted(grouped[cell], key=lambda item: (-item[0], item[1]))
        umi_values = [umi for umi, _barcode in ranked]
        sgrna_values = [sgrna for _umi, sgrna in ranked]
        final_sgrna = final_assigned_sgrna_func(
            sgrna_values,
            umi_values,
            min_total_umi=min_total_umi,
            min_top_umi=min_top_umi,
        )
        rows.append(
            {
                "cell": cell,
                "sgrna": sgrna_values,
                "umi": umi_values,
                "final_assigned_sgrna": final_sgrna,
                "umi_count": umi_values[0] if umi_values else 0,
                "sgrna_type": final_assigned_type(final_sgrna),
            }
        )
    return rows


def write_debug_csv(rows, output_path):
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "cell",
            "sgrna",
            "umi",
            "final_assigned_sgrna",
            "umi_count",
            "sgrna_type",
        ])
        for row in rows:
            writer.writerow([
                row["cell"],
                format_list(row["sgrna"]),
                format_list(row["umi"]),
                row["final_assigned_sgrna"],
                row["umi_count"],
                row["sgrna_type"],
            ])


def write_summary(rows, output_path):
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cell", "final_assigned_sgrna", "umi_count", "sgrna_type"])
        for row in rows:
            writer.writerow([
                row["cell"],
                row["final_assigned_sgrna"],
                row["umi_count"],
                row["sgrna_type"],
            ])


def write_cell_sgrna_table(rows, output_path):
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "cell_barcode",
            "sgrnas",
            "sgrna_umis",
        ])
        for row in rows:
            writer.writerow([
                "{0}-1".format(row["cell"]),
                format_list(row["sgrna"]),
                format_list(row["umi"]),
            ])


def main(args):
    if not os.path.exists(args.input):
        raise FileNotFoundError(args.input)

    sgrna_dict = load_sgrna_file(args.sgrna_file, reverse_complement=args.rc)

    grouped = load_assignments(args.input, sgrna_dict)
    rows = summarize_cells(
        grouped,
        min_total_umi=args.assignment_min_total_umi,
        min_top_umi=args.assignment_min_top_umi,
    )

    if args.debug_csv:
        write_debug_csv(rows, args.debug_csv)

    final_sgrna_distribution(
        [row["sgrna_type"] for row in rows],
        args.sgrna_stat,
    )
    umi_distribution_pie(
        [row["umi_count"] for row in rows],
        args.umi_pie,
    )

    write_summary(rows, args.output)
    write_cell_sgrna_table(rows, args.cell_sgrna_table)

    print("Analysis complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process sgRNA and UMI data.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--sgrna_file", required=True)
    parser.add_argument("--rc", action="store_true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--sgrna_stat", default="sgrna_stat.tsv")
    parser.add_argument("--umi_pie", default="umi_distribution.png")
    parser.add_argument("--debug_csv", default=None)
    parser.add_argument("--cell_sgrna_table", default="cell_sgrna_table.tsv")
    parser.add_argument(
        "--assignment_min_total_umi",
        type=int,
        default=3,
        help="Minimum total sgRNA-supporting UMIs required for a final assignment.",
    )
    parser.add_argument(
        "--assignment_min_top_umi",
        type=int,
        default=3,
        help="Minimum top sgRNA UMI count required for a final assignment.",
    )
    args = parser.parse_args()
    main(args)
