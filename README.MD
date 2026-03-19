# Cellecta Single-Cell Pipeline Prototype

This repository is a customer-facing example workspace for Cellecta single-cell
analysis projects. The long-term target is a unified package that supports:

- CloneTracker barcode assignment and QC
- CRISPR single-cell workflows
- Shared reporting, configuration, and helper utilities
- Docker-based execution for reproducible customer runs

## Current layout

```text
configs/                          Example and future runtime configuration
docs/                             Notes about architecture and pipeline design
examples/                         Small example inputs and usage patterns
scripts/                          Convenience launchers and helper scripts
src/cellecta_sc_pipeline/         Reusable Python package
tests/                            Future automated tests
```

## Python package layout

```text
src/cellecta_sc_pipeline/
  pipelines/
    clonetracker/
    crispr/
  shared/
```

## Installable CLI commands

After installing the repo with `pip install .`, the package exposes:

- `cellecta-full-pipeline`
- `cellecta-clonetracker-batch`
- `cellecta-scrnaseq-qc`

The root-level launchers remain available for convenience during prototyping:

- `run_pipeline.py`
- `clonetracker_batch.py`
- `scRNAseq_QC.py`

## Full pipeline concept

The end-to-end pipeline expects one sample sheet with:

- `sample`
- `gex_fastq_dir`
- `clonetracker_fastq_dir`

It then runs:

1. `cellranger count` for gene expression FASTQs
2. CloneTracker barcode assignment for barcode FASTQs
3. Final QC/reporting that joins the CellRanger and CloneTracker outputs

If the repo root contains `best_sequence_umi.py`, `barcode_process_umi_5.py`,
and `barcode_final_assignment_50M.py`, those helper scripts are used by default.

See [docs/full_pipeline_usage.md](docs/full_pipeline_usage.md) for the current
command shape.

## Local installation

```bash
python -m pip install --upgrade pip
python -m pip install .
```

## Docker plan

A starter `Dockerfile` is included so this repo can evolve into a containerized
customer delivery. The current image installs the Python package and helper
scripts, but it does not bundle 10x Genomics `cellranger`.

That means the Docker image is ready today for:

- CloneTracker barcode processing
- QC/report generation
- Full-pipeline dry runs
- Full-pipeline runs where `cellranger` is provided separately

Typical mounted paths should look like:

- `/work/input`
- `/work/output`
- `/work/refs`
- `/work/config`

## WSL-focused usage

The primary target runtime is now WSL/Linux. In practice that means:

- run the package from an Ubuntu shell in WSL
- use Linux-style mount paths like `/mnt/d/...`
- prefer Docker commands launched from WSL rather than PowerShell
- pass Linux paths to `cellranger`, FASTQ folders, references, and output roots

Example build:

```bash
docker build -t cellecta-sc-pipeline .
```

Example full-pipeline run:

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

Example CloneTracker-only run from the same image:

```bash
docker run --rm \
  --entrypoint cellecta-clonetracker-batch \
  -v /path/to/input:/work/input \
  -v /path/to/output:/work/output \
  cellecta-sc-pipeline \
  --samples_csv /work/input/generated_clonetracker_samples.csv \
  --out_root /work/output/clonetracker \
  --bc14_file /work/input/Cellecta-CloneTrackerXP-50M-BC14-LNGS-300-Library-Design.txt \
  --bc30_file /work/input/Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt
```

Example WSL launch:

```bash
python -m pip install .

cellecta-full-pipeline \
  --samples-csv /mnt/d/Project/Clonetracker/103951/examples/pipeline_samples.csv \
  --pipeline-root /mnt/d/Project/Clonetracker/103951/test_run_wsl \
  --cellranger-bin /mnt/d/Project/Software/cellranger-10.0.0/cellranger \
  --transcriptome /mnt/d/Project/Clonetracker/103951/data/refdata-gex-GRCh38-2024-A \
  --bc14-file /mnt/d/Project/Clonetracker/103951/data/Cellecta-CloneTrackerXP-50M-BC14-LNGS-300-Library-Design.txt \
  --bc30-file /mnt/d/Project/Clonetracker/103951/data/Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt \
  --barcode-search-umi-cutoff 0 \
  --assignment-min-total-umi 3 \
  --assignment-min-top-umi 3
```

For a WSL-first Docker launch, use [scripts/run_docker_wsl.sh](/d:/Project/Clonetracker/103951/scripts/run_docker_wsl.sh).

The checked-in example sample sheet currently references only the `Control`
dataset, because that is the data present under `data/` in this workspace.
