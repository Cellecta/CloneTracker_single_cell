#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import re
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def reverse_complement_seq(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return str(seq).translate(complement)[::-1]


def load_barcode_file(filename, reverse_complement=False):
    data = pd.read_csv(filename, sep="\t").fillna("")

    seq_candidates = ["BC14 sequence", "BC30 sequence", "BC sequence"]
    id_candidates = ["#BC14 ID", "#BC30 ID", "BC14 ID", "BC30 ID", "BC ID"]

    seq_col = next((c for c in seq_candidates if c in data.columns), None)
    id_col = next((c for c in id_candidates if c in data.columns), None)

    if seq_col is None:
        raise ValueError(f"Cannot find barcode sequence column in {filename}")

    if id_col is None:
        raise ValueError(f"Cannot find barcode ID column in {filename}")

    sequences = data[seq_col].astype(str).str.strip()
    ids = data[id_col].astype(str).str.strip()

    if reverse_complement:
        sequences = sequences.apply(reverse_complement_seq)

    return dict(zip(sequences, ids))


def bc_type(bc14_str, bc30_str):
    bc14 = str(bc14_str)
    bc30 = str(bc30_str)

    if bc14 in ["0", "", "NA"] and bc30 not in ["0", "", "NA"]:
        return "bc30_only"

    if bc30 in ["0", "", "NA"] and bc14 not in ["0", "", "NA"]:
        return "bc14_only"

    if bc14 in ["0", "", "NA"] and bc30 in ["0", "", "NA"]:
        return "no_barcode"

    return "bc14_bc30"


def is_missing(x):
    return str(x) in ["0", "", "NA"]


def barcode_name(bc14_str, bc30_str, bc14_dict, bc30_dict, r2_seq):
    r2_seq = str(r2_seq)

    bc14_missing = is_missing(bc14_str)
    bc30_missing = is_missing(bc30_str)

    left = r2_seq[25:39] if len(r2_seq) >= 39 else r2_seq
    right = r2_seq[-30:] if len(r2_seq) >= 30 else r2_seq

    if not bc14_missing and not bc30_missing:
        return f"{bc14_dict.get(bc14_str, bc14_str)}_{bc30_dict.get(bc30_str, bc30_str)}"

    if bc14_missing and not bc30_missing:
        return f"{left}_{bc30_dict.get(bc30_str, bc30_str)}"

    if bc30_missing and not bc14_missing:
        return f"{bc14_dict.get(bc14_str, bc14_str)}_{right}"

    return "NA"


def final_assigned_barcode_func(barcode_list, umi_count_list, min_total_umi=3, min_top_umi=3):
    if not umi_count_list:
        return "NA"

    total_umi = sum(umi_count_list)
    top_umi = umi_count_list[0]

    if total_umi < min_total_umi or top_umi < min_top_umi:
        return "NA"

    valid_barcodes = [barcode for barcode in barcode_list if barcode != "NA"]
    if len(valid_barcodes) == 0:
        return "NA"

    if len(umi_count_list) == 1:
        return barcode_list[0]

    if top_umi >= (total_umi / 2.0):
        return barcode_list[0]

    if umi_count_list[1] / float(top_umi) >= 0.7:
        return "multi_barcode"

    return barcode_list[0]


def final_assigned_type(final_assigned_barcode_str):
    value = str(final_assigned_barcode_str).strip()
    pattern = re.compile(r"bc14-\d+_bc30-\d+", re.IGNORECASE)

    if value == "NA":
        return "Undetermined due to low UMI"

    if value == "multi_barcode":
        return "multi_barcode"

    if pattern.fullmatch(value):
        return "One Barcode"

    return "One Barcode with mutant"


def final_barcode_distribution(bc_type_list, output_path):
    counts = Counter(bc_type_list)
    labels = list(counts.keys())
    sizes = list(counts.values())

    barcode_stat = pd.DataFrame({"Type": labels, "Count": sizes})
    barcode_stat.to_csv(output_path, sep="\t", index=False)

    plt.figure(figsize=(10, 6))
    plt.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=90)
    plt.title("Barcode Type Distribution")
    plt.tight_layout()
    plt.savefig(output_path.replace(".tsv", ".png"), dpi=300)
    plt.close()


def umi_distribution_pie(umi_list, output_path):
    umi_list = [x for x in umi_list if pd.notna(x)]

    if not umi_list:
        return

    bins = np.arange(min(umi_list), max(umi_list) + 20, 20)
    hist, bin_edges = np.histogram(umi_list, bins=bins)

    labels = [
        f"{int(bin_edges[i])}-{int(bin_edges[i + 1])}"
        for i in range(len(bin_edges) - 1)
    ]

    plt.figure(figsize=(10, 6))
    plt.pie(hist, labels=labels, startangle=90)
    plt.title("UMI Distribution")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main(args):
    if not os.path.exists(args.input):
        raise FileNotFoundError(args.input)

    data = pd.read_csv(args.input, sep="\t").fillna("")

    required_cols = ["cell", "umi", "bc14", "bc30", "R2_sequence"]
    missing_cols = [column for column in required_cols if column not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing columns: {missing_cols}")

    bc14_dict = load_barcode_file(args.bc14, reverse_complement=args.rc)
    bc30_dict = load_barcode_file(args.bc30, reverse_complement=args.rc)

    data["bc_type"] = data.apply(lambda row: bc_type(row["bc14"], row["bc30"]), axis=1)

    data_sorted = data.sort_values(by=["cell", "umi"], ascending=[True, False]).copy()
    data_sorted["barcode_name"] = data_sorted.apply(
        lambda row: barcode_name(
            row["bc14"],
            row["bc30"],
            bc14_dict,
            bc30_dict,
            row["R2_sequence"],
        ),
        axis=1,
    )

    data_final = (
        data_sorted.groupby("cell")
        .agg({
            "umi": list,
            "barcode_name": list,
        })
        .reset_index()
    )
    data_final.rename(columns={"barcode_name": "barcode"}, inplace=True)

    data_final["final_assigned_barcode"] = data_final.apply(
        lambda row: final_assigned_barcode_func(
            row["barcode"],
            row["umi"],
            min_total_umi=args.assignment_min_total_umi,
            min_top_umi=args.assignment_min_top_umi,
        ),
        axis=1,
    )

    data_final["umi_count"] = data_final["umi"].apply(lambda values: values[0] if len(values) > 0 else 0)
    data_final["barcode_type"] = data_final["final_assigned_barcode"].apply(final_assigned_type)

    if args.debug_csv:
        data_final.to_csv(args.debug_csv, index=False)

    final_barcode_distribution(data_final["barcode_type"], args.barcode_stat)
    umi_distribution_pie(data_final["umi_count"].tolist(), args.umi_pie)

    data_final[[
        "cell",
        "final_assigned_barcode",
        "umi_count",
        "barcode_type",
    ]].to_csv(args.output, sep="\t", index=False)

    def format_list(value):
        if isinstance(value, list):
            return ";".join(map(str, value))
        return "NA"

    cell_table = pd.DataFrame({
        "cell_barcode": data_final["cell"].astype(str) + "-1",
        "clonetracker_barcodes": data_final["barcode"].apply(format_list),
        "clonetracker_barcode_umis": data_final["umi"].apply(format_list),
    })
    cell_table.to_csv(args.cell_barcode_table, sep="\t", index=False)

    print("Analysis complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process barcode and UMI data.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--bc14", required=True)
    parser.add_argument("--bc30", required=True)
    parser.add_argument("--rc", action="store_true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--barcode_stat", default="barcode_stat.tsv")
    parser.add_argument("--umi_pie", default="umi_distribution.png")
    parser.add_argument("--debug_csv", default=None)
    parser.add_argument("--cell_barcode_table", default="cell_clonetracker_barcode_table.tsv")
    parser.add_argument(
        "--assignment_min_total_umi",
        type=int,
        default=3,
        help="Minimum total barcode-supporting UMIs required for a final assignment.",
    )
    parser.add_argument(
        "--assignment_min_top_umi",
        type=int,
        default=3,
        help="Minimum top barcode UMI count required for a final assignment.",
    )
    args = parser.parse_args()
    main(args)
