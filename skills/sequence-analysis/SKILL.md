---
name: sequence-analysis
description: Quality-check, summarize, and validate nucleotide DNA/RNA FASTA or FASTQ files, including gzip-compressed reads. Use when asked to inspect sequencing reads, read length, GC content, ambiguous bases, FASTQ quality, sequence formats, or produce a concise sequence QC report.
---

# Sequence analysis and QC

## When to use this skill

Use for first-pass validation and descriptive QC of `.fa`, `.fasta`, `.fna`, `.faa`, `.fq`, `.fastq`, and compressed `.gz` equivalents. Do not use this summary as a replacement for adapter/quality trimming, contamination screening, read mapping, or variant calling.

## Procedure

1. Confirm the file exists, format, compression, source, and expected sequence type. Do not infer the organism or reference from filenames alone.
2. Run the streaming parser. It validates headers and FASTQ sequence/quality lengths while calculating summaries without loading all reads into memory.
3. Read the JSON output and inspect the optional per-record TSV for outliers. Verify the inferred format; use `--format` explicitly for unusual extensions.
4. For FASTQ, interpret the quality summary as Phred+33 only. Check instrument and preprocessing provenance if older encodings may be involved.
5. Report anomalies and limitations. Do not silently edit raw files. Keep raw input immutable and store reports separately.

## Commands

From the project root after `python -m pip install -e .`:

```bash
bio-skills sequence-qc reads.fastq.gz --json results/reads_qc.json --records-tsv results/read_metrics.tsv
bio-skills sequence-qc contigs.fasta --format fasta --json results/assembly_qc.json
```

Without installation, use `PYTHONPATH=src python -m bioskills.cli sequence-qc reads.fastq.gz`.

A helper wrapper is available at `skills/sequence-analysis/scripts/run_qc.py`.

## What is reported

Record count; total bases; min/max/mean/median length; exact length distribution; GC percentage; non-ACGT/ambiguous base count and percentage; base counts; for FASTQ, mean ASCII-33 base-quality score. Optional TSV reports each record's ID, length, GC%, ambiguous bases, and mean Phred+33 quality.

## Quality gates and pitfalls

- Stop on malformed FASTA/FASTQ records rather than producing a partial silent report.
- Check read length, GC%, ambiguity, quality, and length distribution against the library protocol and expected organism; there is no universal pass/fail threshold.
- A mean base quality can hide a poor 3′ tail. Use FastQC/MultiQC for per-cycle diagnostics on production datasets.
- GC outliers may reflect biology, contamination, or library bias; do not label them contaminants without evidence.
- The parser treats symbols outside `ACGTU` as ambiguous for the ambiguity metric; this is descriptive and may need adapting for protein sequences.
- The program assumes Phred+33 for FASTQ score conversion and records that assumption in the JSON.

## Expected outputs

A machine-readable JSON report and optionally a per-record TSV. Keep the exact command and package version with the results. Use synthetic fixtures from `examples/sequences/` to verify installation.
