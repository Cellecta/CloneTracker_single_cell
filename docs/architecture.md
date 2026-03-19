# Prototype Architecture

## Goals

- Keep CloneTracker and CRISPR workflows separated where assay logic differs.
- Share QC, reporting, and file-parsing utilities where possible.
- Preserve simple command-line entry points while organizing code as a package.

## Planned modules

- `pipelines.clonetracker`
  CloneTracker-specific barcode extraction, assignment, and summarization.
- `pipelines.crispr`
  Future CRISPR-specific guide assignment and analysis logic.
- `shared`
  Shared QC, plotting, reporting, and common helper functions.

## Prototype rule of thumb

Assay-specific code should live under the corresponding pipeline folder.
Anything that would reasonably be reused by both CRISPR and CloneTracker should
move into `shared`.
