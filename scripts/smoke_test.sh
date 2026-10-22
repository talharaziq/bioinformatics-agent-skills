#!/usr/bin/env bash
# Run automated tests and core CLI smoke tests using included synthetic data.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cd "$ROOT"
python -m unittest discover -s tests -v
python scripts/validate_skills.py
bash -n scripts/run_salmon_pipeline.sh
python -m bioskills.cli sequence-qc examples/sequences/toy.fasta --json "$TMP/sequence_qc.json"
python -m bioskills.cli expression-qc --matrix examples/expression/counts.tsv --metadata examples/expression/metadata.tsv --outdir "$TMP/expression-qc"
python -m bioskills.cli ora --genes examples/pathways/query_genes.txt --gmt examples/pathways/toy.gmt --universe examples/pathways/universe.txt --output "$TMP/ora.tsv" --min-size 2
python -m bioskills.cli validate-manifest examples/expression/manifest.tsv --check-fastq --json "$TMP/manifest.json"
python - "$TMP" <<'PY'
import json, sys
from pathlib import Path
root = Path(sys.argv[1])
assert json.loads((root / "sequence_qc.json").read_text())["records"] == 3
assert json.loads((root / "expression-qc" / "expression_qc.json").read_text())["samples"] == 4
assert len((root / "ora.tsv").read_text().splitlines()) >= 2
assert json.loads((root / "manifest.json").read_text())["n_samples"] == 4
print("All core smoke tests passed.")
PY
