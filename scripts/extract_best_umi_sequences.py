#!/usr/bin/env python3

import argparse
import gzip
from collections import defaultdict


def extract_umi_from_header(header_str):
    """Extract the concatenated cell barcode and UMI from a FASTQ header."""
    parts = header_str.split("_")
    if len(parts) < 3:
        raise ValueError("Unexpected header format: {0}".format(header_str))
    return parts[1] + parts[2].split(" ", 1)[0]


def iter_fastq_sequences(input_file):
    """Yield FASTQ header and sequence using a lightweight 4-line parser."""
    with gzip.open(input_file, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError("Incomplete FASTQ record encountered in {0}".format(input_file))

            yield header.rstrip("\n\r"), sequence.rstrip("\n\r")


def process_fastq(input_file, output_file):
    umi_sequences = defaultdict(lambda: defaultdict(int))

    for header, sequence in iter_fastq_sequences(input_file):
        umi = extract_umi_from_header(header)
        umi_sequences[umi][sequence] += 1

    with open(output_file, "w") as output:
        output.write("cell_umi\tseq\treads_count\n")
        for umi in sorted(umi_sequences):
            sequence_counts = umi_sequences[umi]
            best_sequence, best_count = max(
                sequence_counts.items(),
                key=lambda item: (item[1], item[0]),
            )
            output.write("{0}\t{1}\t{2}\n".format(umi, best_sequence, best_count))

    print("Processing complete. Results saved to {0}".format(output_file))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process an NGS FASTQ file to extract UMI sequences.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file (gzipped).")
    parser.add_argument("-o", "--output", required=True, help="Output file for UMI sequences.")
    args = parser.parse_args()
    process_fastq(args.input, args.output)
