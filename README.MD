---
title: "Cellecta Clonetracker Barcode Single Cell Data Processing"
author: "Dongfang Hu"
date: "2024-12-15"
---

# Introduction

This document provides a step-by-step workflow for processing Cellecta CloneTracker Barcode single-cell data and assigning CloneTracker barcodes to individual cells using shell and Python scripts. The goal is to ensure clarity and reproducibility in sequencing data processing and barcode analysis.

Cellecta has designed 10X-compatible single-cell NGS primers to amplify the CloneTracker barcode fragment, facilitating the assignment of CloneTracker barcodes to 10X cell barcodes. After sequencing, each sample will generate two NGS libraries:

Gene Expression Library – containing gene expression profiling data.
CloneTracker Barcode Library – containing the CloneTracker barcode data.
This workflow will guide you through processing these libraries to accurately link CloneTracker barcodes with corresponding 10X cell barcodes.

# Input Data

The input data for this workflow includes:

- **FASTQ files from Gene Expression Profiling NGS data:**
  - `/path/to/10X_Gene_Expression_profiling_fastqs/Sample1_GEX/Sample1_GEX_S3_R1_001.fastq.gz`
  - `/path/to/10X_Gene_Expression_profiling_fastqs/Sample1_GEX/Sample1_GEX_S3_R2_001.fastq.gz`
- **FASTQ files from the Cellecta Clonetracker barcode NGS library:**
  - `/path/to/10X_cellecta_cloneracker_barcode_NGS_fastqs/Sample1_CTB/Sample1_CTB_S7_R1_001.fastq.gz`
  - `/path/to/10X_cellecta_cloneracker_barcode_NGS_fastqs/Sample1_CTB/Sample1_CTB_S7_R2_001.fastq.gz`
- **Reference barcode sequences provided by Cellecta:**
  - `Cellecta-CloneTrackerXP-5M-Pool1-BC14-LNGS-300-Library-Design.txt`
  - `Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt`

# Step 1: Load Reference Barcode Sequences

Analyze Gene Expression Profiling NGS data and extract cell barcodes using `cellranger count`:

```bash
/mnt/project/Pipeline/Software/cellranger-8.0.0/cellranger count \
    --id=Sample1 \
    --fastqs=/path/to/10X_Gene_Expression_profiling_fastqs/Sample1_GEX/ \
    --transcriptome=/mnt/project/Pipeline/Reference/10XGenomics/refdata-gex-GRCh38-2024-A \
    --create-bam true \
    --include-introns false

zcat /Path/to/cellranger_count_output/Sample1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > /Path/to/cellranger_count_output/Sample1/Sample1_barcode.xls

sed 's/-.*//' /Path/to/cellranger_count_output/Sample1/Sample1_barcode.xls > /Path/to/cellranger_count_output/Sample1/Sample1_barcode_cleaned.tsv
```

# Step 2: Process Cellecta CloneTracker Barcode NGS Library Using UMI-Tools

Use `umi_tools` to extract and match 10X cell barcodes:

### Shell Script

```bash
#!/bin/bash

FILEIN1=/path/to/10X_cellecta_cloneracker_barcode_NGS_fastqs/Sample1_CTB/Sample1_CTB_S7_R1_001.fastq.gz
FILEOUT1=/path/to/10X_cellecta_cloneracker_barcode_NGS_fastqs/`basename ${FILEIN1} .fastq.gz`_extracted.fastq.gz
FILEIN2=/path/to/10X_cellecta_cloneracker_barcode_NGS_fastqs/Sample1_CTB/Sample1_CTB_S7_R2_001.fastq.gz
FILEOUT2=/path/to/10X_cellecta_cloneracker_barcode_NGS_fastqs/`basename ${FILEIN2} .fastq.gz`_extracted.fastq.gz

WHITELIST=/Path/to/cellranger_count_output/Sample1/Sample1_barcode_cleaned.tsv

umi_tools extract \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
    --stdin $FILEIN1 \
    --stdout $FILEOUT1 \
    --read2-in $FILEIN2 \
    --read2-out=$FILEOUT2 \
    --whitelist=$WHITELIST
```

# Step 3: Identify the Best Sequence for Each UMI

Run the following Python script to select the best-quality sequence for each UMI:

```bash
python best_sequence_umi.py
```

# Step 4: Assign Clonetracker Barcodes to 10X Cell Barcodes

Use this script to map Clonetracker barcodes to cell barcodes:

```bash
python barcode_process_umi_5.py
```

# Step 5: Final Barcode Assignment Using UMI Distribution

Perform the final barcode assignment based on UMI distribution and criteria:

```bash
python umi_distribution_analysis_umi_5.py
```

# Conclusion

This workflow provides a systematic approach to process Cellecta Clonetracker single-cell data, assign barcodes, and analyze barcode distributions. Scripts and parameters can be adjusted based on dataset-specific requirements. By following these steps, you can ensure accurate and reproducible results for your analysis.
