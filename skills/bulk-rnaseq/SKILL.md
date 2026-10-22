---
name: bulk-rnaseq
description: Plan, validate, and run reproducible bulk RNA-seq workflows from FASTQ files through trimming, quality control, Salmon transcript quantification, and DESeq2 differential expression. Use when asked to analyze bulk RNA-seq, validate a sample sheet, quantify reads, compare conditions, or audit DE results.
---

# Bulk RNA-seq workflow

## Start with the design, not the aligner

1. Ask/inspect what organism, reference/annotation release, strandedness, read layout, library preparation, groups, biological replicates, pairing, batch, and contrast are intended. Never choose a reference build by guessing.
2. Validate the manifest. It must be a tab-separated file containing `sample_id`, `condition`, `fastq_r1`, `fastq_r2`; an empty `fastq_r2` means single-end. Paths resolve relative to the manifest.
3. Confirm raw FASTQs exist and first records are well formed. Preserve the original files. Review FastQC/MultiQC, read depth, adapter content, duplication, quality, and contamination signals before setting trim options.
4. Run fastp + FastQC + Salmon against a pre-built Salmon index for the intended transcriptome. Keep logs, tool versions, index provenance, transcriptome FASTA checksum, and annotation release.
5. Use the generated Salmon manifest and the exact transcript-to-gene map used for the index with `tximport_deseq2.R`. For an existing gene count matrix use `deseq2_counts.R` instead.
6. Inspect PCA/ sample-level QC, library sizes, sample correlations, outlier candidates, and design rank. A batch term is appropriate only if estimable and scientifically justified.
7. Report effect size (log2 fold change), raw p-value, adjusted p-value, base mean, filtering method, contrast, tested universe, and limitations. Use FDR (usually BH) rather than unadjusted p-values for genome-wide claims.

## Validate sample sheet

```bash
bio-skills validate-manifest /path/to/manifest.tsv --check-fastq
```

Manifest paths are resolved relative to the manifest file. Validation catches missing files, duplicate IDs, unsafe folder names, malformed first FASTQ records when requested, and single-condition designs. A one-replicate group produces a warning, not an assertion of valid inferential design.

## Run quantification

Required command-line tools: `fastp`, `FastQC`, `Salmon`. Build the Salmon index yourself and record the reference provenance.

```bash
scripts/run_salmon_pipeline.sh /path/to/manifest.tsv /path/to/salmon_index results/rnaseq 8
```

The pipeline saves per-sample logs, FastQC output, fastp HTML/JSON, and Salmon quantification under the output directory. Existing Salmon result folders trigger a stop to prevent accidental overwrite.

## Differential expression

Salmon transcript quantifications:

```bash
Rscript scripts/tximport_deseq2.R results/rnaseq/run_manifest.tsv tx2gene.tsv results/rnaseq condition control treatment
```

Raw gene-count matrix:

```bash
Rscript scripts/deseq2_counts.R counts.tsv metadata.tsv results/deseq2 condition control treatment
```

Optional batch: append the metadata column name after `treatment` only when the design is full-rank. Required packages are Bioconductor `tximport` and/or `DESeq2` as appropriate.

## Statistical and reproducibility guardrails

- DESeq2 uses a negative-binomial model for raw counts; do not supply TPM/FPKM/CPM/log values or rounded normalized values.
- Do not drop inconvenient samples solely because they change significance. Use documented QC reasons and preserve an audit trail.
- Do not treat technical replicates as independent biological replicates.
- Paired/repeated-measures studies may need blocking or subject-level terms; the generic two-group design is not sufficient for every experiment.
- Condition and batch must not be perfectly confounded. The scripts reject a non-full-rank design but cannot validate all scientific assumptions.
- The example manifest and reads are synthetic. Do not treat their results as biological evidence.
