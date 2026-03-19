#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  cat <<'EOF'
Usage:
  scripts/run_docker_wsl.sh <image> <input_dir> <output_dir> <refs_dir> <samples_csv> [extra pipeline args...]

Example:
  scripts/run_docker_wsl.sh \
    cellecta-sc-pipeline \
    /mnt/d/Project/Clonetracker/103951/examples \
    /mnt/d/Project/Clonetracker/103951/test_run_wsl \
    /mnt/d/Project/Clonetracker/103951/data \
    pipeline_samples.csv \
    --skip-cellranger
EOF
  exit 1
fi

IMAGE="$1"
INPUT_DIR="$(realpath "$2")"
OUTPUT_DIR="$(realpath "$3")"
REFS_DIR="$(realpath "$4")"
SAMPLES_CSV="$5"
shift 5

mkdir -p "$OUTPUT_DIR"

docker run --rm \
  -v "${INPUT_DIR}:/work/input" \
  -v "${OUTPUT_DIR}:/work/output" \
  -v "${REFS_DIR}:/work/refs" \
  "$IMAGE" \
  --samples-csv "/work/input/${SAMPLES_CSV}" \
  --pipeline-root /work/output/run1 \
  --transcriptome /work/refs/refdata-gex-GRCh38-2024-A \
  --bc14-file /work/refs/Cellecta-CloneTrackerXP-50M-BC14-LNGS-300-Library-Design.txt \
  --bc30-file /work/refs/Cellecta-CloneTrackerXP-50M-BC30-LNGS-300-Library-Design.txt \
  "$@"
