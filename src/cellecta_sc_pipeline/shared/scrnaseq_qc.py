#!/usr/bin/env python3
import argparse
import base64
import shutil
from datetime import datetime
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
except ModuleNotFoundError as exc:
    plt = None
    np = None
    pd = None
    sc = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None


CELLECTA_GREEN = "#2E8B57"
CELLECTA_GRID = "#E5E5E5"

if plt is not None:
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=[CELLECTA_GREEN])


METRIC_DESCRIPTIONS = {
    "sample": "Sample identifier",
    "cells_total": "Total number of cells before QC filtering",
    "cells_final": "Number of cells remaining after QC filtering",
    "median_counts": "Median UMI counts per cell",
    "median_genes": "Median number of detected genes per cell",
    "median_pct_mt": "Median percentage of mitochondrial transcripts",
    "cells_one_barcode": "Cells with confident single CloneTracker barcode assignment",
    "cells_multiple_barcodes": "Cells containing multiple CloneTracker barcodes",
    "cells_undetermined_low_umi": (
        "Cells without confident barcode assignment due to low UMI support"
    ),
    "cells_one_barcode_with_mutant": (
        "Cells with one barcode assignment that contains partial or mutated sequence"
    ),
    "fraction_cells_with_one_barcode_assignment": (
        "Fraction of filtered cells with confident single barcode assignment"
    ),
    "n_unique_clones": (
        "Number of unique CloneTracker barcodes among single-barcode cells"
    ),
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Cellecta CloneTracker single-cell QC pipeline"
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-name")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=8000)
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--max-mt-pct", type=float, default=15)
    parser.add_argument("--clonetracker-summary")
    parser.add_argument("--clonetracker-umi")
    return parser


def find_filtered_matrix_dir(input_dir: Path) -> Path:
    input_dir = Path(input_dir)
    if input_dir.name == "filtered_feature_bc_matrix":
        return input_dir

    matches = list(input_dir.rglob("filtered_feature_bc_matrix"))
    if not matches:
        raise FileNotFoundError("filtered_feature_bc_matrix not found")
    if len(matches) > 1:
        raise ValueError("Multiple filtered_feature_bc_matrix directories found")
    return matches[0]


def require_existing_file(path: str | Path, description: str) -> Path:
    """Normalize and validate a required file path."""
    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"{description} not found: {resolved}")
    return resolved


def add_qc_gene_flags(adata) -> None:
    genes = adata.var_names
    adata.var["mt"] = genes.str.startswith(("MT-", "mt-"))
    adata.var["ribo"] = genes.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = genes.str.match(r"^HB[ABDEGMQZ]")


def plot_qc_violin(adata, outdir: Path, sample: str) -> None:
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        multi_panel=True,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(outdir / f"{sample}.qc_violin.png", dpi=150)
    plt.close()


def plot_clonetracker_types(adata, outdir: Path, sample: str) -> None:
    counts = adata.obs["clonetracker_barcode_type"].value_counts()
    plt.figure(figsize=(10, 6))
    counts.plot(kind="bar", color=CELLECTA_GREEN, edgecolor="black")
    plt.ylabel("Cells")
    plt.xticks(rotation=30, ha="right")
    plt.grid(axis="y", color=CELLECTA_GRID)
    plt.tight_layout()
    plt.savefig(outdir / f"{sample}.clonetracker_barcode_type_barplot.png", dpi=150)
    plt.close()


def plot_clone_size_distribution(adata, outdir: Path, sample: str):
    total_cells = adata.n_obs
    mask = adata.obs["clonetracker_barcode_type"].isin(
        ["One Barcode", "One Barcode with mutant"]
    )
    assigned = adata.obs.loc[mask, "clonetracker_final_barcode"].dropna()
    clone_sizes = assigned.value_counts()
    top50 = clone_sizes.head(50)

    table = top50.reset_index()
    table.columns = ["Clone_barcode", "Cell_count"]
    table["Fraction_of_total_cells"] = table["Cell_count"] / total_cells
    table.to_csv(outdir / f"{sample}.clone_size_top50.tsv", sep="\t", index=False)

    plt.figure(figsize=(12, 5))
    top50.plot(kind="bar", color=CELLECTA_GREEN)
    plt.ylabel("Cells")
    plt.grid(axis="y", color=CELLECTA_GRID)
    plt.tight_layout()
    plt.savefig(outdir / f"{sample}.clone_size_top50.png", dpi=150)
    plt.close()

    return table


def add_clonetracker(adata, summary_file: Path):
    adata.obs["cell_nosuffix"] = adata.obs_names.str.replace(r"-\d+$", "", regex=True)

    summary = pd.read_csv(summary_file, sep="\t").rename(
        columns={
            "final_assigned_barcode": "clonetracker_final_barcode",
            "umi_count": "clonetracker_umi_count",
            "barcode_type": "clonetracker_barcode_type",
        }
    )
    summary = summary.set_index("cell")

    adata.obs["clonetracker_final_barcode"] = adata.obs["cell_nosuffix"].map(
        summary["clonetracker_final_barcode"]
    )
    adata.obs["clonetracker_umi_count"] = adata.obs["cell_nosuffix"].map(
        summary["clonetracker_umi_count"]
    )
    adata.obs["clonetracker_barcode_type"] = adata.obs["cell_nosuffix"].map(
        summary["clonetracker_barcode_type"]
    )

    missing = adata.obs["clonetracker_barcode_type"].isna()
    adata.obs.loc[missing, "clonetracker_barcode_type"] = (
        "Undetermined due to low UMI"
    )
    return adata


def compute_clone_stats(adata) -> dict:
    obs = adata.obs
    total_cells = len(obs)
    type_counts = obs["clonetracker_barcode_type"].value_counts()

    one = type_counts.get("One Barcode", 0)
    multi = type_counts.get("multi_barcode", 0)
    undetermined = type_counts.get("Undetermined due to low UMI", 0)
    one_mutant = type_counts.get("One Barcode with mutant", 0)

    assigned = obs.loc[
        obs["clonetracker_barcode_type"] == "One Barcode",
        "clonetracker_final_barcode",
    ]

    return {
        "cells_one_barcode": int(one),
        "cells_multiple_barcodes": int(multi),
        "cells_undetermined_low_umi": int(undetermined),
        "cells_one_barcode_with_mutant": int(one_mutant),
        "fraction_cells_with_one_barcode_assignment": (
            (one + one_mutant) / total_cells if total_cells else 0.0
        ),
        "n_unique_clones": assigned.nunique(),
    }


def embed_png(path: Path) -> str:
    if not Path(path).exists():
        return ""
    with open(path, "rb") as handle:
        data = base64.b64encode(handle.read()).decode()
    return f'<img src="data:image/png;base64,{data}" width="700">'


def write_html(sample: str, outdir: Path, summary: dict, clone_table, args) -> None:
    violin = embed_png(outdir / f"{sample}.qc_violin.png")
    typebar = embed_png(outdir / f"{sample}.clonetracker_barcode_type_barplot.png")
    clones = embed_png(outdir / f"{sample}.clone_size_top50.png")

    rows = []
    for metric, value in summary.items():
        if isinstance(value, float):
            value = f"{value:.4f}"
        rows.append(
            {
                "Metric": metric,
                "Value": value,
                "Description": METRIC_DESCRIPTIONS.get(metric, ""),
            }
        )
    summary_table = pd.DataFrame(rows).to_html(index=False)

    qc_filter_text = f"""
<div class="info-panel">

<b>Cell quality filtering</b><br><br>

Cells were filtered using the following criteria:

<ul>
<li>Minimum detected genes per cell: {args.min_genes}</li>
<li>Maximum detected genes per cell: {args.max_genes}</li>
<li>Minimum UMI counts per cell: {args.min_counts}</li>
<li>Maximum mitochondrial transcript percentage: {args.max_mt_pct}%</li>
</ul>

Cells not meeting these thresholds were removed prior to downstream analysis.

</div>
"""

    barcode_definition_text = """
<div class="info-panel">

<b>CloneTracker barcode assignment</b><br><br>

CloneTracker barcodes were assigned to cells based on UMI support from detected
BC14-BC30 barcode pairs.

Barcode assignment categories:

<ul>

<li><b>One Barcode</b><br>
A single barcode dominates the UMI counts within a cell
(top barcode >= 50% of total barcode UMIs and top UMI >= 3).</li>

<li><b>multi_barcode</b><br>
Two barcodes have similar support within a cell
(second barcode UMI >= 70% of the top barcode UMI),
suggesting possible barcode collision or multiplet.</li>

<li><b>One Barcode with mutant</b><br>
A barcode containing partial or mutated sequence detected during barcode recovery.</li>

<li><b>Undetermined due to low UMI</b><br>
Cells lacking sufficient UMI evidence for confident barcode assignment
(total barcode UMIs < 3 or top barcode UMI < 3).</li>

</ul>

</div>
"""

    clone_table_html = ""
    if clone_table is not None:
        clone_table_html = clone_table.to_html(index=False, float_format="%.4f")

    html = f"""
<html>

<head>

<title>Cellecta CloneTracker Barcode Assignment and QC Report</title>

<style>

body {{
font-family: Arial, Helvetica, sans-serif;
margin: 40px;
color: #333333;
}}

.header-banner {{
background-color: {CELLECTA_GREEN};
color: white;
padding: 20px;
border-radius: 6px;
}}

.info-panel {{
margin-top: 20px;
padding: 15px;
background-color: #F6FBF7;
border-left: 5px solid {CELLECTA_GREEN};
}}

table {{
border-collapse: collapse;
width: 85%;
}}

th {{
background-color: {CELLECTA_GREEN};
color: white;
padding: 8px;
}}

td {{
border: 1px solid #E5E5E5;
padding: 6px;
}}

tr:nth-child(even) {{
background-color: #F6FBF7;
}}

h2 {{
color: {CELLECTA_GREEN};
margin-top: 30px;
}}

</style>

</head>

<body>

<div class="header-banner">

<h1>Cellecta CloneTracker Barcode Assignment and QC Report</h1>

Single Cell Dataset Analysis

</div>

<div class="info-panel">

<b>Sample:</b> {sample}<br>
<b>Generated:</b> {datetime.now()}

</div>

<h2>Cell Quality Filtering</h2>

{qc_filter_text}

<h2>Summary</h2>

{summary_table}

<h2>CloneTracker Barcode Assignment</h2>

{barcode_definition_text}

<h2>QC Metrics</h2>

{violin}

<h2>CloneTracker Barcode Types</h2>

{typebar}

<h2>Clone Size Distribution</h2>

{clones}

<h2>Top 50 Clone Table</h2>

{clone_table_html}

</body>
</html>
"""

    (outdir / f"{sample}.QC_report.html").write_text(html)


def write_readme(outdir: Path, sample: str) -> None:
    text = f"""
Cellecta CloneTracker Barcode Assignment QC

Main report
-----------
{sample}.QC_report.html

Summary
-------
qc_summary.tsv

QC metrics
----------
{sample}.cell_qc_metrics.tsv

Filtered dataset
----------------
{sample}.filtered.h5ad

CloneTracker outputs
--------------------
{sample}.clonetracker_obs.tsv: cell metadata with CloneTracker assignments
{sample}.cell_clonetracker_barcode_table.tsv: CloneTracker barcode UMI table

Clone size analysis
-------------------
{sample}.clone_size_top50.tsv
{sample}.clone_size_top50.png

QC figures
----------
{sample}.qc_violin.png
{sample}.clonetracker_barcode_type_barplot.png
"""
    (outdir / "README.txt").write_text(text)


def main() -> None:
    args = build_parser().parse_args()

    if IMPORT_ERROR is not None:
        raise ModuleNotFoundError(
            "Required dependencies for scRNAseq_QC.py are not installed. "
            "Install matplotlib, numpy, pandas, and scanpy to run QC."
        ) from IMPORT_ERROR

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    matrix_dir = find_filtered_matrix_dir(input_dir)
    sample = args.sample_name if args.sample_name else matrix_dir.parent.name

    print("Processing:", sample)

    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", make_unique=True)
    cells_total_raw = adata.n_obs

    add_qc_gene_flags(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True)

    adata.obs.to_csv(output_dir / f"{sample}.cell_qc_metrics.tsv", sep="\t")
    plot_qc_violin(adata, output_dir, sample)

    keep = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes)
        & (adata.obs["n_genes_by_counts"] <= args.max_genes)
        & (adata.obs["total_counts"] >= args.min_counts)
        & (adata.obs["pct_counts_mt"] <= args.max_mt_pct)
    )
    adata = adata[keep].copy()

    clone_table = None
    if bool(args.clonetracker_summary) != bool(args.clonetracker_umi):
        raise ValueError(
            "Provide both --clonetracker-summary and --clonetracker-umi, or neither."
        )

    if args.clonetracker_summary and args.clonetracker_umi:
        clonetracker_summary = require_existing_file(
            args.clonetracker_summary,
            "CloneTracker summary",
        )
        clonetracker_umi = require_existing_file(
            args.clonetracker_umi,
            "CloneTracker barcode UMI table",
        )
        shutil.copy2(
            clonetracker_umi,
            output_dir / f"{sample}.cell_clonetracker_barcode_table.tsv",
        )
        adata = add_clonetracker(adata, clonetracker_summary)
        adata.obs[
            [column for column in adata.obs.columns if column.startswith("clonetracker_")]
        ].to_csv(output_dir / f"{sample}.clonetracker_obs.tsv", sep="\t")
        plot_clonetracker_types(adata, output_dir, sample)
        clone_table = plot_clone_size_distribution(adata, output_dir, sample)

    adata.write(output_dir / f"{sample}.filtered.h5ad")

    summary = {
        "sample": sample,
        "cells_total": cells_total_raw,
        "cells_final": adata.n_obs,
        "median_counts": float(np.median(adata.obs["total_counts"])),
        "median_genes": float(np.median(adata.obs["n_genes_by_counts"])),
        "median_pct_mt": float(np.median(adata.obs["pct_counts_mt"])),
    }
    if "clonetracker_barcode_type" in adata.obs:
        summary.update(compute_clone_stats(adata))

    pd.DataFrame([summary]).to_csv(output_dir / "qc_summary.tsv", sep="\t", index=False)

    if clone_table is not None:
        write_html(sample, output_dir, summary, clone_table, args)

    write_readme(output_dir, sample)
    print("QC finished")


if __name__ == "__main__":
    main()
