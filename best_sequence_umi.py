import argparse
from Bio import SeqIO
from collections import defaultdict
import gzip

# Function to extract UMI from the record description
def extract_umi_from_header(header_str):
    parts = header_str.split('_')
    return parts[1] + parts[2].split(' ')[0]

# Main function
def main(args):
    # Dictionary to store sequences and their frequencies for each UMI
    umi_sequences = defaultdict(lambda: defaultdict(int))
    
    # Parse the FASTQ file and count sequences for each UMI
    with gzip.open(args.input, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Extract UMI from the header
            umi = extract_umi_from_header(record.description)
            
            # Count occurrences of each sequence for the given UMI
            umi_sequences[umi][record.seq] += 1

    # Identify the most frequent sequence for each UMI
    best_sequences = {
        umi: max(sequences, key=sequences.get)
        for umi, sequences in umi_sequences.items()
    }

    # Write the results to a tab-delimited file
    with open(args.output, "w") as output:
        # Write the header
        output.write("cell_umi\tseq\treads_count\n")
        
        # Write UMI, best sequence, and read counts
        for umi, best_sequence in best_sequences.items():
            counts = umi_sequences[umi][best_sequence]
            output.write(f"{umi}\t{best_sequence}\t{counts}\n")

    print(f"Output written to {args.output}")

# Argument parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a FASTQ file to extract UMIs and determine the most frequent sequences.")
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the input FASTQ file (gzip compressed)."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path to the output TSV file."
    )
    args = parser.parse_args()
    main(args)
