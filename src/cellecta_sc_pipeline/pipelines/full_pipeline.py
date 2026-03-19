#!/usr/bin/env python3
import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Iterable


SOURCE_REPO_ROOT = Path(__file__).resolve().parents[3]
REQUIRED_SAMPLE_COLUMNS = {"sample", "gex_fastq_dir", "clonetracker_fastq_dir"}


def resolve_local_helper(filename: str) -> Path:
    """Resolve helper scripts from common source and installed-package locations."""
    candidates = [
        Path.cwd() / filename,
        SOURCE_REPO_ROOT / filename,
        Path(sys.executable).resolve().parent.parent / filename,
        Path('/app') / filename,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[1]


BEST_SEQUENCE_DEFAULT = resolve_local_helper("best_sequence_umi.py")
BARCODE_PROCESS_DEFAULT = resolve_local_helper("barcode_process_umi_5.py")
FINAL_ASSIGNMENT_DEFAULT = resolve_local_helper("barcode_final_assignment_50M.py")


def local_tool_arg(parser: argparse.ArgumentParser, flag: str, default_path: Path, help_text: str) -> None:
    """Register a tool-path argument, requiring it only when no local default exists."""
    kwargs = {"help": help_text}
    if default_path.exists():
        kwargs["default"] = str(default_path)
        kwargs["required"] = False
    else:
        kwargs["required"] = True
    parser.add_argument(flag, **kwargs)


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI for the end-to-end single-cell processing workflow."""
    parser = argparse.ArgumentParser(
        description="Run the Cellecta single-cell pipeline from FASTQs to final QC analysis."
    )
    parser.add_argument("--samples-csv", required=True, help="Sample sheet with per-sample FASTQ directories")
    parser.add_argument("--pipeline-root", required=True, help="Root output directory for the full pipeline run")
    parser.add_argument("--cellranger-bin", default="cellranger", help="Path to the cellranger executable")
    parser.add_argument("--transcriptome", required=True, help="CellRanger transcriptome reference")
    parser.add_argument("--create-bam", default="true", choices=["true", "false"], help="Pass through to cellranger count --create-bam")
    parser.add_argument("--include-introns", default="false", choices=["true", "false"], help="Pass through to cellranger count --include-introns")
    parser.add_argument("--bc14-file", required=True, help="BC14 reference file")
    parser.add_argument("--bc30-file", required=True, help="BC30 reference file")
    local_tool_arg(parser, "--best-sequence-umi-py", BEST_SEQUENCE_DEFAULT, "Path to best_sequence_umi.py")
    local_tool_arg(parser, "--barcode-process-py", BARCODE_PROCESS_DEFAULT, "Path to barcode_process_umi_5.py")
    local_tool_arg(parser, "--final-assignment-py", FINAL_ASSIGNMENT_DEFAULT, "Path to barcode_final_assignment_50M.py")
    parser.add_argument("--bc-pattern", default="CCCCCCCCCCCCCCCCNNNNNNNNNNNN", help="umi_tools extract barcode pattern")
    parser.add_argument(
        "--barcode-search-umi-cutoff",
        "--umi-cutoff",
        dest="barcode_search_umi_cutoff",
        type=int,
        default=0,
        help="Minimum per-cell UMI count to keep barcode candidates before final assignment. Use 0 to keep all.",
    )
    parser.add_argument("--assignment-min-total-umi", type=int, default=3, help="Minimum total barcode-supporting UMIs required for a final assignment")
    parser.add_argument("--assignment-min-top-umi", type=int, default=3, help="Minimum top barcode UMI count required for a final assignment")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=8000)
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--max-mt-pct", type=float, default=15)
    parser.add_argument("--skip-cellranger", action="store_true", help="Skip cellranger count and reuse existing outputs")
    parser.add_argument("--skip-clonetracker", action="store_true", help="Skip CloneTracker barcode assignment and reuse existing outputs")
    parser.add_argument("--skip-qc", action="store_true", help="Skip final QC/report generation")
    parser.add_argument("--force", action="store_true", help="Force rerun of steps where supported")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them")
    return parser


def read_samples_csv(path: Path) -> list[dict]:
    """Load and validate the required sample metadata used by downstream stages."""
    if not path.exists():
        raise FileNotFoundError(f"samples.csv not found: {path}")

    with path.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError("samples.csv does not contain any sample rows")

    for column in REQUIRED_SAMPLE_COLUMNS:
        if column not in rows[0].keys():
            raise ValueError(f"samples.csv missing required column: {column}")

    for index, row in enumerate(rows, start=1):
        for column in REQUIRED_SAMPLE_COLUMNS:
            if not row.get(column, "").strip():
                raise ValueError(f"samples.csv row {index} has empty value for '{column}'")
    return rows


def require_existing_path(path: Path, description: str) -> None:
    """Raise a clear error when a required file or directory is missing."""
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {path}")


def validate_pipeline_inputs(rows: list[dict], args: argparse.Namespace) -> None:
    """Check user-provided references and per-sample FASTQ directories up front."""
    require_existing_path(Path(args.transcriptome), "Transcriptome reference")
    require_existing_path(Path(args.bc14_file), "BC14 reference file")
    require_existing_path(Path(args.bc30_file), "BC30 reference file")

    if not args.skip_cellranger:
        cellranger_bin = Path(args.cellranger_bin)
        if args.cellranger_bin != "cellranger" and not cellranger_bin.exists():
            raise FileNotFoundError(f"cellranger executable not found: {cellranger_bin}")

    for row in rows:
        sample = row["sample"].strip()
        require_existing_path(Path(row["gex_fastq_dir"].strip()), f"[{sample}] gex_fastq_dir")
        require_existing_path(
            Path(row["clonetracker_fastq_dir"].strip()),
            f"[{sample}] clonetracker_fastq_dir",
        )


def expected_cellranger_barcodes(pipeline_root: Path, sample: str) -> Path:
    """Return the filtered barcode whitelist expected from CellRanger count."""
    return (
        pipeline_root
        / "cellranger"
        / f"{sample}_GEX"
        / "outs"
        / "filtered_feature_bc_matrix"
        / "barcodes.tsv.gz"
    )


def validate_reused_cellranger_outputs(rows: list[dict], pipeline_root: Path) -> None:
    """Ensure reuse mode has the CellRanger outputs required downstream."""
    for row in rows:
        sample = row["sample"].strip()
        require_existing_path(
            expected_cellranger_barcodes(pipeline_root, sample),
            f"[{sample}] existing CellRanger barcodes for --skip-cellranger",
        )


def validate_reused_clonetracker_outputs(rows: list[dict], pipeline_root: Path) -> Path:
    """Ensure reuse mode has the CloneTracker outputs needed by QC."""
    clonetracker_root = pipeline_root / "clonetracker"
    for row in rows:
        sample = row["sample"].strip()
        require_existing_path(
            clonetracker_root / sample / f"{sample}_barcode_assignment_summary.tsv",
            f"[{sample}] existing CloneTracker summary for --skip-clonetracker",
        )
        require_existing_path(
            clonetracker_root / sample / f"{sample}_cell_clonetracker_barcode_table.tsv",
            f"[{sample}] existing CloneTracker barcode table for --skip-clonetracker",
        )
    return clonetracker_root


def command_to_text(cmd: Iterable[str]) -> str:
    """Render a subprocess command list into a readable log string."""
    return " ".join(str(part) for part in cmd)


def run_command(cmd: list[str], *, cwd: Path | None = None, dry_run: bool = False) -> None:
    """Log and optionally execute a subprocess command."""
    print("\n[CMD] " + command_to_text(cmd), flush=True)
    if dry_run:
        return
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def write_clonetracker_samples_csv(rows: list[dict], output_path: Path, cellranger_root: Path) -> None:
    """Generate the sample sheet expected by the batch CloneTracker runner."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fastq_dir", "cellranger_barcodes_gz"])
        writer.writeheader()
        for row in rows:
            sample = row["sample"].strip()
            writer.writerow({
                "sample": sample,
                "fastq_dir": row["clonetracker_fastq_dir"].strip(),
                # CloneTracker consumes the filtered CellRanger barcodes for barcode assignment.
                "cellranger_barcodes_gz": str(
                    expected_cellranger_barcodes(cellranger_root.parent, sample)
                ),
            })


def run_cellranger_for_sample(row: dict, args: argparse.Namespace, cellranger_root: Path) -> None:
    """Run `cellranger count` for one sample and verify the expected barcode output exists."""
    sample = row["sample"].strip()
    sample_id = f"{sample}_GEX"
    sample_out = cellranger_root / sample_id
    sample_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(args.cellranger_bin), "count",
        f"--id={sample_id}",
        f"--fastqs={row['gex_fastq_dir'].strip()}",
        f"--transcriptome={args.transcriptome}",
        f"--create-bam={args.create_bam}",
        f"--include-introns={args.include_introns}",
    ]
    run_command(cmd, cwd=cellranger_root, dry_run=args.dry_run)

    if not args.dry_run:
        expected = sample_out / "outs" / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"
        if not expected.exists():
            raise FileNotFoundError(f"Expected CellRanger output not found: {expected}")


def run_clonetracker(rows: list[dict], args: argparse.Namespace, pipeline_root: Path) -> Path:
    """Run the batch CloneTracker assignment stage for all samples."""
    generated_csv = pipeline_root / "configs" / "generated_clonetracker_samples.csv"
    cellranger_root = pipeline_root / "cellranger"
    write_clonetracker_samples_csv(rows, generated_csv, cellranger_root)

    clonetracker_out = pipeline_root / "clonetracker"
    cmd = [
        sys.executable,
        "-m",
        "cellecta_sc_pipeline.pipelines.clonetracker.batch",
        "--samples_csv", str(generated_csv),
        "--out_root", str(clonetracker_out),
        "--bc14_file", str(args.bc14_file),
        "--bc30_file", str(args.bc30_file),
        "--best_sequence_umi_py", str(args.best_sequence_umi_py),
        "--barcode_process_py", str(args.barcode_process_py),
        "--final_assignment_py", str(args.final_assignment_py),
        "--bc_pattern", str(args.bc_pattern),
        "--barcode_search_umi_cutoff", str(args.barcode_search_umi_cutoff),
        "--assignment_min_total_umi", str(args.assignment_min_total_umi),
        "--assignment_min_top_umi", str(args.assignment_min_top_umi),
    ]
    if args.force:
        cmd.append("--force")
    run_command(cmd, cwd=Path.cwd(), dry_run=args.dry_run)
    return clonetracker_out


def run_qc_for_sample(sample: str, args: argparse.Namespace, pipeline_root: Path, clonetracker_root: Path) -> None:
    """Launch the downstream QC/reporting step for a single processed sample."""
    gex_outs = pipeline_root / "cellranger" / f"{sample}_GEX" / "outs"
    qc_out = pipeline_root / "analysis" / sample
    summary_path = clonetracker_root / sample / f"{sample}_barcode_assignment_summary.tsv"
    cell_barcode_table_path = clonetracker_root / sample / f"{sample}_cell_clonetracker_barcode_table.tsv"

    if not args.dry_run:
        require_existing_path(gex_outs, f"[{sample}] CellRanger outs directory for QC")
        require_existing_path(summary_path, f"[{sample}] CloneTracker summary for QC")
        require_existing_path(
            cell_barcode_table_path,
            f"[{sample}] CloneTracker barcode table for QC",
        )

    cmd = [
        sys.executable,
        "-m",
        "cellecta_sc_pipeline.shared.scrnaseq_qc",
        "--input", str(gex_outs),
        "--output", str(qc_out),
        "--sample-name", sample,
        "--min-genes", str(args.min_genes),
        "--max-genes", str(args.max_genes),
        "--min-counts", str(args.min_counts),
        "--max-mt-pct", str(args.max_mt_pct),
        "--clonetracker-summary", str(summary_path),
        "--clonetracker-umi", str(cell_barcode_table_path),
    ]
    run_command(cmd, cwd=Path.cwd(), dry_run=args.dry_run)


def main() -> None:
    """Coordinate CellRanger, CloneTracker, and QC stages for all requested samples."""
    args = build_parser().parse_args()
    pipeline_root = Path(args.pipeline_root)
    pipeline_root.mkdir(parents=True, exist_ok=True)

    rows = read_samples_csv(Path(args.samples_csv))
    validate_pipeline_inputs(rows, args)

    cellranger_root = pipeline_root / "cellranger"
    if not args.skip_cellranger:
        # First stage: generate filtered feature-barcode matrices for each GEX library.
        for row in rows:
            run_cellranger_for_sample(row, args, cellranger_root)
    elif not args.dry_run:
        validate_reused_cellranger_outputs(rows, pipeline_root)

    if args.skip_clonetracker:
        clonetracker_root = (
            pipeline_root / "clonetracker"
            if args.dry_run
            else validate_reused_clonetracker_outputs(rows, pipeline_root)
        )
    else:
        # Second stage: assign CloneTracker barcodes using the CellRanger barcode whitelist.
        clonetracker_root = run_clonetracker(rows, args, pipeline_root)

    if not args.skip_qc:
        # Final stage: combine expression outputs and barcode assignments into QC summaries.
        for row in rows:
            run_qc_for_sample(row["sample"].strip(), args, pipeline_root, clonetracker_root)

    print("\nPipeline finished.")
    print(f"Pipeline root: {pipeline_root}")
    print(f"Source repo root: {SOURCE_REPO_ROOT}")


if __name__ == "__main__":
    main()
