import argparse
from Bio import SeqIO
from collections import defaultdict
import gzip

def extract_umi_from_header(header_str):
    """ Extracts UMI from FASTQ header assuming a specific format. """
    try:
        parts = header_str.split('_')
        return parts[1] + parts[2].split(' ')[0]
    except IndexError:
        raise ValueError(f"Unexpected header format: {header_str}")

def process_fastq(input_file, output_file):
    # Dictionary to store sequences and their frequencies for each UMI
    umi_sequences = defaultdict(lambda: defaultdict(int))

    # Read and process the FASTQ file
    with gzip.open(input_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            umi = extract_umi_from_header(record.description)
            umi_sequences[umi][str(record.seq)] += 1

    # Determine the most frequent sequence for each UMI
    best_sequences = {
        umi: max(sequences, key=sequences.get) for umi, sequences in umi_sequences.items()
    }

    # Write results to output file
    with open(output_file, 'w') as output:
        output.write('cell_umi\tseq\treads_count\n')
        for umi, best_sequence in best_sequences.items():
            counts = umi_sequences[umi][best_sequence]
            output.write(f"{umi}\t{best_sequence}\t{counts}\n")

    print(f"Processing complete. Results saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process an NGS FASTQ file to extract UMI sequences.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file (gzipped).")
    parser.add_argument("-o", "--output", required=True, help="Output file for UMI sequences.")

    args = parser.parse_args()

    # Run the processing function with provided arguments
    process_fastq(args.input, args.output)
