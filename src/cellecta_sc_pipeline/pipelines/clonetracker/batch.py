#!/usr/bin/env python3
import argparse
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


SOURCE_REPO_ROOT = Path(__file__).resolve().parents[4]


def resolve_local_helper(filename: str) -> Path:
    """Resolve helper scripts from common source and installed-package locations."""
    candidates = [
        Path.cwd() / filename,
        Path.cwd() / "scripts" / filename,
        SOURCE_REPO_ROOT / filename,
        SOURCE_REPO_ROOT / "scripts" / filename,
        Path(sys.executable).resolve().parent.parent / filename,
        Path(sys.executable).resolve().parent.parent / "scripts" / filename,
        Path('/app') / filename,
        Path('/app') / "scripts" / filename,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[1]


BEST_SEQUENCE_DEFAULT = resolve_local_helper("extract_best_umi_sequences.py")
BARCODE_PROCESS_DEFAULT = resolve_local_helper("process_barcode_umis.py")
FINAL_ASSIGNMENT_DEFAULT = resolve_local_helper("assign_final_barcodes.py")


def local_tool_arg(parser: argparse.ArgumentParser, flag: str, default_path: Path, help_text: str) -> None:
    kwargs = {"help": help_text}
    if default_path.exists():
        kwargs["default"] = str(default_path)
        kwargs["required"] = False
    else:
        kwargs["required"] = True
    parser.add_argument(flag, **kwargs)


def run(cmd: List[str], *, cwd: Optional[Path] = None) -> None:
    print("\n[CMD] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def read_samples_csv(path: Path) -> List[dict]:
    with path.open("r", newline="") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError("samples.csv does not contain any sample rows")

    required = {"sample", "fastq_dir", "cellranger_barcodes_gz"}
    for index, row in enumerate(rows, start=1):
        missing = required - set(key for key, value in row.items() if key is not None)
        if missing:
            raise ValueError(
                f"samples.csv row {index} missing columns: {sorted(missing)}"
            )
        for key in required:
            if not row.get(key, "").strip():
                raise ValueError(
                    f"samples.csv row {index} has empty value for '{key}'"
                )
    return rows


def make_whitelist(barcodes_gz: Path, whitelist_out: Path) -> None:
    whitelist_out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(barcodes_gz, "rt") as fin, whitelist_out.open("w") as fout:
        for line in fin:
            barcode = line.strip()
            if not barcode:
                continue
            fout.write(barcode.split("-", 1)[0] + "\n")


def merge_gz_members(inputs: List[Path], out_gz: Path) -> None:
    out_gz.parent.mkdir(parents=True, exist_ok=True)
    with out_gz.open("wb") as writer:
        for file_path in inputs:
            with file_path.open("rb") as reader:
                shutil.copyfileobj(reader, writer)


def find_fastqs(fastq_dir: Path, read: str) -> List[Path]:
    return sorted(fastq_dir.glob(f"*{read}*.fastq.gz"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch CloneTracker barcode assignment for multiple samples."
    )
    parser.add_argument("--samples_csv", required=True, help="CSV with sample,fastq_dir,cellranger_barcodes_gz")
    parser.add_argument("--out_root", required=True, help="Output root folder (per-sample subfolders created)")
    parser.add_argument("--bc14_file", required=True, help="BC14 reference file")
    parser.add_argument("--bc30_file", required=True, help="BC30 reference file")
    local_tool_arg(parser, "--best_sequence_umi_py", BEST_SEQUENCE_DEFAULT, "Path to extract_best_umi_sequences.py")
    local_tool_arg(parser, "--barcode_process_py", BARCODE_PROCESS_DEFAULT, "Path to process_barcode_umis.py")
    local_tool_arg(parser, "--final_assignment_py", FINAL_ASSIGNMENT_DEFAULT, "Path to assign_final_barcodes.py")
    parser.add_argument("--bc_pattern", default="CCCCCCCCCCCCCCCCNNNNNNNNNNNN", help="umi_tools --bc-pattern")
    parser.add_argument(
        "--barcode_search_umi_cutoff",
        "--umi_cutoff",
        dest="barcode_search_umi_cutoff",
        type=int,
        default=0,
        help="Minimum per-cell UMI count to keep barcode candidates before final assignment. Use 0 to keep all.",
    )
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
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs")
    return parser


def process_sample(sample: str, fastq_dir: Path, barcodes_gz: Path, out_root: Path, args: argparse.Namespace) -> None:
    sample_out = out_root / sample
    sample_out.mkdir(parents=True, exist_ok=True)

    whitelist = sample_out / f"{sample}_barcode_whitelist.tsv"
    merged_r1 = sample_out / f"{sample}_merged_R1.fastq.gz"
    merged_r2 = sample_out / f"{sample}_merged_R2.fastq.gz"
    extracted_r1 = sample_out / f"{sample}_extracted_R1.fastq.gz"
    extracted_r2 = sample_out / f"{sample}_extracted_R2.fastq.gz"
    cell_umi_tsv = sample_out / f"{sample}_cell_umi.tsv"
    assign_umi_tsv = sample_out / f"{sample}_barcode_assignment_umi.tsv"
    summary_tsv = sample_out / f"{sample}_barcode_assignment_summary.tsv"
    cell_barcode_table_tsv = sample_out / f"{sample}_cell_clonetracker_barcode_table.tsv"

    print(f"\n========== Processing sample: {sample} ==========")

    if not fastq_dir.exists():
        raise FileNotFoundError(f"[{sample}] fastq_dir not found: {fastq_dir}")
    if not barcodes_gz.exists():
        raise FileNotFoundError(f"[{sample}] cellranger_barcodes_gz not found: {barcodes_gz}")

    r1_files = find_fastqs(fastq_dir, "R1")
    r2_files = find_fastqs(fastq_dir, "R2")
    if not r1_files:
        raise FileNotFoundError(f"[{sample}] No R1 fastqs found in {fastq_dir}")
    if not r2_files:
        raise FileNotFoundError(f"[{sample}] No R2 fastqs found in {fastq_dir}")

    if whitelist.exists() and not args.force:
        print(f"[SKIP] whitelist exists: {whitelist}")
    else:
        print(f"[MAKE] whitelist: {whitelist}")
        make_whitelist(barcodes_gz, whitelist)

    if merged_r1.exists() and merged_r2.exists() and not args.force:
        print(f"[SKIP] merged fastqs exist: {merged_r1}, {merged_r2}")
    else:
        print(f"[MERGE] R1 files: {len(r1_files)} -> {merged_r1}")
        merge_gz_members(r1_files, merged_r1)
        print(f"[MERGE] R2 files: {len(r2_files)} -> {merged_r2}")
        merge_gz_members(r2_files, merged_r2)

    if extracted_r1.exists() and extracted_r2.exists() and not args.force:
        print(f"[SKIP] extracted fastqs exist: {extracted_r1}, {extracted_r2}")
    else:
        run([
            "umi_tools", "extract",
            f"--bc-pattern={args.bc_pattern}",
            "--stdin", str(merged_r1),
            "--stdout", str(extracted_r1),
            "--read2-in", str(merged_r2),
            "--read2-out", str(extracted_r2),
            "--whitelist", str(whitelist),
        ], cwd=sample_out)

    if cell_umi_tsv.exists() and not args.force:
        print(f"[SKIP] cell_umi exists: {cell_umi_tsv}")
    else:
        run([
            sys.executable,
            str(args.best_sequence_umi_py),
            "-i", str(extracted_r2),
            "-o", str(cell_umi_tsv),
        ], cwd=sample_out)

    if assign_umi_tsv.exists() and not args.force:
        print(f"[SKIP] barcode assignment exists: {assign_umi_tsv}")
    else:
        run([
            sys.executable,
            str(args.barcode_process_py),
            "--cell_umi", str(cell_umi_tsv),
            "--bc14_file", str(args.bc14_file),
            "--bc30_file", str(args.bc30_file),
            "--whitelist", str(whitelist),
            "--output", str(assign_umi_tsv),
            "--umi_cutoff", str(args.barcode_search_umi_cutoff),
            "--rc",
        ], cwd=sample_out)

    if summary_tsv.exists() and cell_barcode_table_tsv.exists() and not args.force:
        print(f"[SKIP] summary exists: {summary_tsv}")
    else:
        run([
            sys.executable,
            str(args.final_assignment_py),
            "--input", str(assign_umi_tsv),
            "--bc14", str(args.bc14_file),
            "--bc30", str(args.bc30_file),
            "--output", str(summary_tsv),
            "--cell_barcode_table", str(cell_barcode_table_tsv),
            "--assignment_min_total_umi", str(args.assignment_min_total_umi),
            "--assignment_min_top_umi", str(args.assignment_min_top_umi),
            "--rc",
        ], cwd=sample_out)

    print(f"[DONE] {sample} -> {summary_tsv}")


def main() -> None:
    args = build_parser().parse_args()
    rows = read_samples_csv(Path(args.samples_csv))
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    for row in rows:
        process_sample(
            sample=row["sample"].strip(),
            fastq_dir=Path(row["fastq_dir"].strip()),
            barcodes_gz=Path(row["cellranger_barcodes_gz"].strip()),
            out_root=out_root,
            args=args,
        )

    print("\nAll samples completed.")


if __name__ == "__main__":
    main()
