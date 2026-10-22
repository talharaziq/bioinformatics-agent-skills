# Bioinformatics Agent Skills

[![Tests](https://github.com/talharaziq/Advanced-Bioinformatics/actions/workflows/test.yml/badge.svg)](https://github.com/talharaziq/Advanced-Bioinformatics/actions/workflows/test.yml)

**A practical, executable collection of agent skills for common bioinformatics tasks.** Each skill includes the open `SKILL.md` format, decision rules, quality gates, and helper scripts. The core command-line tools use Python's standard library, so basic sequence QC, matrix QC, manifest validation, and local pathway over-representation analysis do not require a large Python stack.

This project is built to work both with coding agents that understand filesystem-based skills and directly from a terminal. Results are deterministic where the underlying method is deterministic, tabular, and designed to be inspected—not treated as unquestionable biological conclusions.

## Included skills

| Skill | What it actually runs | Useful outputs |
|---|---|---|
| `sequence-analysis` | Streaming FASTA/FASTQ parser; supports gzip; checks record structure; reports length distribution, GC%, ambiguity, and Phred+33 quality summaries | JSON summary; optional per-record TSV |
| `geo-dataset-exploration` | Fetches a GEO Series through GEOparse; exports series/sample metadata and optionally a single-platform VALUE matrix | `series_metadata.json`, `samples.tsv`, optional `series_matrix.tsv` |
| `bulk-rnaseq` | Validates a sample sheet; orchestrates fastp, FastQC, and Salmon; provides DESeq2 scripts for raw gene-count or Salmon/tximport data | Trimming/QC logs, quantification, normalized counts, differential-expression table, report |
| `single-cell-rnaseq` | Optional Scanpy pipeline for 10x matrix directories, 10x H5, or AnnData `.h5ad` files | Cell QC table, processed `.h5ad`, summary JSON, optional UMAP plot |
| `pathway-enrichment` | Local GMT gene-set over-representation analysis using a hypergeometric upper tail and Benjamini–Hochberg FDR | Ranked enrichment TSV with overlaps, p-values and FDR |

## Quick start (core tools)

Requires Python 3.10 or newer. Install from this directory:

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

Smoke-test the install with included toy data:

```bash
bio-skills sequence-qc examples/sequences/toy.fasta --json results/sequence_qc.json
bio-skills expression-qc --matrix examples/expression/counts.tsv --metadata examples/expression/metadata.tsv --outdir results/expression-qc
bio-skills ora --genes examples/pathways/query_genes.txt --gmt examples/pathways/toy.gmt --universe examples/pathways/universe.txt --output results/ora.tsv --min-size 2
bio-skills validate-manifest examples/expression/manifest.tsv --check-fastq
```

Equivalent direct invocation without an editable install: `PYTHONPATH=src python -m bioskills.cli sequence-qc examples/sequences/toy.fasta`.

## Install optional workflows

### GEO retrieval

```bash
python -m pip install -e '.[geo]'
bio-skills geo-fetch GSEXXXXX --outdir results/GSEXXXXX
# Optional: also export the platform's VALUE expression matrix (may be large)
bio-skills geo-fetch GSEXXXXX --outdir results/GSEXXXXX --include-matrix
```

GEO is a public archive of heterogeneous studies. The tool exports metadata; it does not guarantee that every series has a comparable expression matrix, clean condition labels, compatible platforms, or enough replicates. Inspect the returned `samples.tsv` and original paper before modelling. The optional matrix export refuses to merge multi-platform series into one table.

### Bulk RNA-seq (fastp + FastQC + Salmon)

Install these external executables using their official installation guides: `fastp`, `FastQC`, and `Salmon`. Build a Salmon index from the exact transcriptome reference you intend to quantify. Validate your manifest first:

```bash
bio-skills validate-manifest examples/expression/manifest.tsv --check-fastq
```

Then run on your own real FASTQ data and Salmon index:

```bash
scripts/run_salmon_pipeline.sh /path/to/manifest.tsv /path/to/salmon_index results/rnaseq 8
```

The sample sheet is tab-separated and must contain `sample_id`, `condition`, `fastq_r1`, and `fastq_r2`. Set `fastq_r2` to an empty field for single-end reads. FASTQ paths are resolved relative to the manifest location. Extra metadata columns (such as `batch`, donor, or paired-subject identifiers) are preserved in `run_manifest.tsv` for downstream R models. The pipeline never downloads a reference or guesses a genome build. Existing per-sample Salmon result folders are not overwritten.

For gene-level differential-expression testing, provide the exact transcript-to-gene table used for the transcriptome reference (`transcript_id`, `gene_id`):

```bash
Rscript scripts/tximport_deseq2.R results/rnaseq/run_manifest.tsv /path/to/tx2gene.tsv results/rnaseq condition control treatment
```

Required R packages: Bioconductor `tximport` and `DESeq2`. For an already-created raw integer gene count matrix, use:

```bash
Rscript scripts/deseq2_counts.R counts.tsv metadata.tsv results/deseq2 condition control treatment
# Optional batch-adjusted model, only if batch is not confounded with condition:
Rscript scripts/deseq2_counts.R counts.tsv metadata.tsv results/deseq2 condition control treatment batch
```

DESeq2 expects raw, non-negative integer counts—not TPM, FPKM, CPM, log-transformed data, or an arbitrary normalized matrix. The scripts check sample alignment, duplicate IDs, integer counts, and design-matrix rank. They cannot rescue inadequate replication, confounded study design, or mismatched references. Inspect QC before interpreting adjusted p-values.

### Single-cell RNA-seq (optional Scanpy)

```bash
python -m pip install -e '.[single-cell]'
bio-skills single-cell /path/to/filtered_feature_bc_matrix --outdir results/single-cell
# Also accepts a 10x .h5 or existing AnnData .h5ad file.
```

Thresholds are exposed as CLI arguments (`--min-genes`, `--min-cells`, `--max-pct-mito`). The default mitochondrial prefix is human/mouse-style `MT-`; adapt this for organisms with a different gene-symbol convention. This is a generic preprocessing workflow, not a substitute for doublet detection, sample-aware integration, batch-effect assessment, or cell-type annotation.

### Local pathway enrichment

Provide a newline-delimited gene list, a GMT gene-set file, and preferably a biologically appropriate tested-gene universe (for RNA-seq, often all genes eligible for the differential test after filtering):

```bash
bio-skills ora --genes hits.txt --gmt gene_sets.gmt --universe tested_genes.txt --output ora.tsv
```

ORA uses exact identifier matching and a hypergeometric upper-tail test. It does not map aliases or Ensembl versions automatically. Use the same identifier namespace in all files. Terms are tested only if their size within the universe is between `--min-size` and `--max-size`; FDR is Benjamini–Hochberg across all qualifying gene sets in the GMT, including zero-overlap sets (zero-overlap rows are omitted from the result table). The result depends strongly on the universe and gene-set collection.

## Skill installation and usage by agents

Each folder in `skills/` contains a valid `SKILL.md` with YAML frontmatter and task-specific procedures. Point your compatible agent at the **whole repository** so the referenced command-line utilities and the sibling `src/` package are present. For Claude Code, one option is to copy a chosen skill directory into `.claude/skills/`; if you do that, also install this project's Python package (`python -m pip install -e /path/to/this/project`) because the helper wrappers import its CLI. Other agents can read the same `SKILL.md` and use the documented commands.

Validate skill metadata with:

```bash
python scripts/validate_skills.py
```

## Outputs, limitations, and responsible use

- The included sample datasets are **synthetic fixtures**, not evidence for any biological hypothesis.
- The utilities do not send data to external services, except `geo-fetch`, which downloads the accession explicitly requested by the user.
- Keep raw data immutable; write outputs to a separate directory; preserve command logs, references, software versions, and sample metadata.
- Never infer that a missing file, a failed download, or an empty output means no biological effect.
- Validate genome build, annotation release, gene identifier version, read layout, library preparation, experimental design, and biological replicates before the analysis.
- Treat exploratory QC, differential-expression tests, ORA, and pathway labels as distinct steps. A significant term or gene is not proof of mechanism or causality.

## Development and tests

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python scripts/validate_skills.py
bash -n scripts/run_salmon_pipeline.sh
# End-to-end smoke test of all dependency-light commands
bash scripts/smoke_test.sh
```

Core test suite has no large bioinformatics runtime dependencies. GEO, Scanpy, external read-processing/quantification tools and R/Bioconductor are optional and have their own prerequisite checks. For scientific assumptions and citations, see [`docs/methods-and-references.md`](docs/methods-and-references.md) and the relevant skill folder.

## Roadmap

Near-term contributions should prioritize actual reproducible fixtures and validated outputs (for example, a small public GEO series with a documented license, a tiny Salmon index built from a pinned reference, or a multi-donor single-cell dataset), rather than adding more prompt-only skill descriptions. See [`docs/validation-plan.md`](docs/validation-plan.md).
