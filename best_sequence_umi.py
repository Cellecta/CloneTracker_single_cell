from Bio import SeqIO
from collections import defaultdict
import gzip

fastq_file_path = "AD_FBP1_S7_R2_001_extracted.fastq.gz"

# Dictionary to store sequences and their frequencies for each UMI
umi_sequences = defaultdict(lambda: defaultdict(int))
def extract_umi_from_header(header_str):
    return header_str.split('_')[1]+header_str.split('_')[2].split(' ')[0]

# Iterate over each record in the FASTQ file
with gzip.open(fastq_file_path, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        # Extract UMI from the header or wherever it's located in your specific data
        umi = extract_umi_from_header(record.description)
        
        # Count occurrences of each sequence for the given UMI
        umi_sequences[umi][record.seq] += 1

# Find the best sequence (highest frequency) for each UMI
best_sequences = {umi: max(sequences, key=sequences.get) for umi, sequences in umi_sequences.items()}
output = open("cell_umi.xls", 'w')
output.write('cell_umi'+'\t'+'seq'+'\t'+'reads_count'+'\n')
for umi, best_sequence in best_sequences.items():
    counts = umi_sequences[umi][best_sequence]
    output.write(str(umi)+'\t'+str(best_sequence)+'\t'+str(counts)+'\n')
output.close()