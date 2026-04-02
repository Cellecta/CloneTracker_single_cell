<p align="center">
  <img src="./logo_head.png" alt="Cellecta" width="100%">
</p>

# Cellecta Single-Cell Pipeline

This repository provides a customer-facing workspace for running Cellecta
single-cell analysis workflows, with a current focus on CloneTracker barcode
assignment, QC, and end-to-end processing from FASTQ inputs.

CloneTracker is Cellecta's lentiviral barcode platform for labeling cells with
stable, heritable, sequenceable barcodes so researchers can follow clonal
behavior over time. For single-cell applications, CloneTracker XP libraries are
designed so barcodes can be detected in RNA-seq workflows, making it possible
to connect clonal identity with transcriptional state in single-cell studies.

For 10x Chromium 3' workflows, CloneTracker XP barcodes are positioned farther
upstream than the standard short 3' gene-expression read typically reaches. In
practice, this means barcode recovery generally requires an additional targeted
PCR amplification step using a barcode-specific primer before downstream
library indexing and sequencing. This repository provides the downstream
analysis workflow for processing those recovered single-cell CloneTracker
barcode reads.

Product information:

- CloneTracker barcode libraries and kits:
  https://cellecta.com/collections/clonetracker-barcode-libraries-and-kits
- Example CloneTracker XP single-cell compatible library:
  https://cellecta.com/products/clonetracker-xp-50m-barcode-3-library-with-rfp-puro
- Cellecta manual section: 10x Chromium 3' Detection of Barcodes/sgRNA in
  pScribe Libraries:
  https://manuals.cellecta.com/crispr-pooled-lentiviral-sgrna-libraries/v3a/en/topic/10x-chromium-3-detection-of-barcodes-sgrna-in-pscribe-libraries
- CloneTracker single-cell analysis workflow repository:
  https://github.com/Cellecta/CloneTracker_single_cell

## What This Repository Supports

- `cellranger count` orchestration for gene expression libraries
- CloneTracker barcode assignment from barcode FASTQs
- QC and reporting that combine CellRanger and CloneTracker outputs
- WSL- and Docker-friendly execution patterns

## Repository Layout

```text
configs/                          Example and future runtime configuration
docs/                             Usage notes and architecture details
examples/                         Example sample sheets
scripts/                          Convenience launchers and helper scripts
src/cellecta_sc_pipeline/         Installable Python package
tests/                            Lightweight test/config scaffolding
```

## Installable Commands

After installing with `pip install .`, the package exposes:

- `cellecta-full-pipeline`
- `cellecta-clonetracker-batch`
- `cellecta-scrnaseq-qc`

Convenience wrappers are also available under `scripts/`:

- `scripts/run_full_pipeline.py`
- `scripts/run_clonetracker_batch.py`
- `scripts/run_scrnaseq_qc.py`

## Workflow Summary

The full pipeline expects a sample sheet with:

- `sample`
- `gex_fastq_dir`
- `clonetracker_fastq_dir`

It then performs:

1. `cellranger count` for each gene expression library
2. CloneTracker barcode assignment using the filtered CellRanger barcodes
3. QC and reporting that join GEX outputs and CloneTracker assignments

The packaged workflow automatically uses these helper scripts from `scripts/`
unless you override them explicitly:

- `extract_best_umi_sequences.py`
- `process_barcode_umis.py`
- `assign_final_barcodes.py`

## Why This Assignment Method

This pipeline uses a conservative CloneTracker assignment strategy on purpose.
For targeted single-cell CloneTracker or sgRNA libraries, an extra PCR
enrichment step can recover multiple weak barcode-related signals, partial
matches, or low-level background molecules in the same cell. Methods designed
for more permissive feature detection can therefore overcall multiple barcodes
per cell and produce assignments that are difficult to interpret biologically.

The current workflow is designed to prioritize reliable per-cell barcode calls.
It keeps one best sequence per cell+UMI, searches for supported barcode
components, ranks barcode candidates by UMI support, and assigns a final
barcode only when there is enough evidence for a confident call. This approach
typically gives fewer but more interpretable assignments than permissive
multi-feature calling and is better suited for CloneTracker and Cellecta-style
targeted libraries.

## Barcode Orientation

Current barcode-orientation note:

- The current pipeline is designed to reverse-complement the input BC14 and BC30 reference barcode sequences first during barcode matching and final assignment.
- In practice, this means the barcode design files passed with `--bc14-file` and `--bc30-file` are currently interpreted through the pipeline's `--rc` path.

For command details, see [docs/full_pipeline_usage.md](docs/full_pipeline_usage.md).

## Quick Start

```bash
python -m pip install --upgrade pip
python -m pip install .
```

Dry-run the full pipeline with the bundled example sample sheet:

```bash
cellecta-full-pipeline \
  --samples-csv examples/pipeline_samples.csv \
  --pipeline-root /path/to/run_output \
  --transcriptome /path/to/refdata-gex-GRCh38-2024-A \
  --bc14-file /path/to/Cellecta-CloneTrackerXP-50M-BC14-LNGS-300-Library-Design.txt \
  --bc30-file /path/to/Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt \
  --barcode-search-umi-cutoff 0 \
  --assignment-min-total-umi 3 \
  --assignment-min-top-umi 3 \
  --dry-run
```

The checked-in example sample sheets currently target the `Control` dataset
available in this workspace.

## Docker

The included `Dockerfile` installs the package and helper scripts, but it does
not bundle 10x Genomics `cellranger`. That makes the image suitable for:

- CloneTracker barcode processing
- QC/report generation
- Full-pipeline dry runs
- Full-pipeline runs where `cellranger` is mounted or provided separately

Typical mounted locations:

- `/work/input`
- `/work/output`
- `/work/refs`

Example build:

```bash
docker build -t cellecta-sc-pipeline .
```

Example run:

```bash
docker run --rm \
  -v /path/to/input:/work/input \
  -v /path/to/output:/work/output \
  -v /path/to/refs:/work/refs \
  cellecta-sc-pipeline \
  --samples-csv /work/input/pipeline_samples.csv \
  --pipeline-root /work/output/run1 \
  --transcriptome /work/refs/refdata-gex-GRCh38-2024-A \
  --bc14-file /work/refs/Cellecta-CloneTrackerXP-50M-BC14-LNGS-300-Library-Design.txt \
  --bc30-file /work/refs/Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt
```
