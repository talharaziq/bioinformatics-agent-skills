#!/usr/bin/env bash
# Run manifest validation followed by fastp -> FastQC -> Salmon quantification.
# Usage: scripts/run_salmon_pipeline.sh manifest.tsv salmon_index output_dir [threads]
set -euo pipefail
if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 manifest.tsv salmon_index output_dir [threads]" >&2
  exit 2
fi
THREADS="${4:-4}"
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then echo "threads must be a positive integer" >&2; exit 2; fi
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
python -m bioskills.cli validate-manifest "$1" --check-fastq >/dev/null
python -m bioskills.cli run-salmon --manifest "$1" --index "$2" --outdir "$3" --threads "$THREADS"
