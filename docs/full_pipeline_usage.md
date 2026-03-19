# Full Pipeline Usage

The prototype pipeline can be launched from either the installed CLI or the
root-level convenience wrapper. Because `best_sequence_umi.py`,
`barcode_process_umi_5.py`, and `barcode_final_assignment_50M.py` are present
in the repo root, the pipeline uses them automatically unless you override the
paths.

The checked-in example CSVs currently reference the `Control` sample only,
because that is the dataset available in the repository `data/` folder.

## Installed CLI

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

## Wrapper script

```bash
python run_pipeline.py \
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

## Docker

The starter Docker image uses `cellecta-full-pipeline` as its entrypoint.

```bash
docker build -t cellecta-sc-pipeline .

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

Note: the current image does not bundle `cellranger`, so full end-to-end runs
still require you to provide `cellranger` separately or use `--skip-cellranger`
when reusing existing outputs.

## WSL

For WSL, prefer Linux-style paths from your Ubuntu shell:

```bash
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

For a Docker wrapper from WSL, use `scripts/run_docker_wsl.sh`.

## UMI thresholds

There are two separate UMI thresholds in the pipeline:

- `--barcode-search-umi-cutoff`
- `--assignment-min-total-umi`
- `--assignment-min-top-umi`

`--barcode-search-umi-cutoff` controls which barcode candidates are retained in
the exported per-cell barcode table before final assignment. Use `0` to keep
all recovered barcode candidates and UMIs for each cell.

`--assignment-min-total-umi` controls the minimum total barcode-supporting UMIs
required to make a final CloneTracker barcode assignment.

`--assignment-min-top-umi` controls the minimum UMI count of the top barcode
required to make a final CloneTracker barcode assignment.

Recommended current setting for your goal:

- `--barcode-search-umi-cutoff 0`
- `--assignment-min-total-umi 3`
- `--assignment-min-top-umi 3`

## Optional overrides

Use these only if you want to point at non-default helper scripts:

- `--best-sequence-umi-py`
- `--barcode-process-py`
- `--final-assignment-py`

## Expected stages

1. `cellranger count` for each sample using `gex_fastq_dir`
2. CloneTracker barcode assignment using `clonetracker_fastq_dir`
3. Final QC/report generation that joins GEX and CloneTracker outputs

## Output layout

```text
<pipeline-root>/
  analysis/
  cellranger/
  clonetracker/
  configs/
```

## Sample sheet columns

- `sample`
- `gex_fastq_dir`
- `clonetracker_fastq_dir`
