---
name: geo-dataset-exploration
description: Retrieve and audit public NCBI Gene Expression Omnibus (GEO) GSE series metadata, sample characteristics, platform information, and optionally a single-platform expression matrix. Use when a user provides a GSE accession, asks to explore a public GEO study, or needs a sample manifest before reanalysis.
---

# GEO dataset exploration

## Scope

This skill retrieves metadata from a specifically requested GEO Series using GEOparse. It does not guarantee a complete, correctly annotated expression matrix or perform differential expression by default.

## Procedure

1. Resolve the exact GSE accession and inspect the study title, summary, organism, platform(s), sample count, publication, and supplementary-file notes.
2. Fetch series metadata into a new output directory. Never overwrite user files or reuse an accession inferred only from a vague description without explaining the choice.
3. Inspect `samples.tsv`. The `characteristics_ch1` fields are often inconsistent, use free text, and may encode treatment, time, tissue, sex, genotype, or batch. Normalize them into explicit metadata columns only after checking source records and the publication.
4. Compare sample identifiers, sample count, platform count, organism, and experimental design with the source paper and GEO web record.
5. Export the optional VALUE matrix only for a series with exactly one platform. The utility refuses to combine multi-platform data into a single matrix.
6. For differential expression, determine whether the data are raw counts, normalized microarray intensities, or another value type. Select a method appropriate to the assay; do not feed arbitrary GEO `VALUE` data into DESeq2.
7. Document exclusions, repeated measures, biological replicates, paired designs, batch effects, and any sample-label assumptions before interpreting results.

## Run

Install the optional client and retrieve metadata:

```bash
python -m pip install -e '.[geo]'
bio-skills geo-fetch GSEXXXXX --outdir results/GSEXXXXX
```

For a single-platform series where GEOparse can produce a VALUE matrix:

```bash
bio-skills geo-fetch GSEXXXXX --outdir results/GSEXXXXX --include-matrix
```

Use the wrapper `skills/geo-dataset-exploration/scripts/fetch_geo.py` for direct execution from the repository.

## Expected outputs

- `series_metadata.json`: accession, title, summary, study design, PMID field if present, sample/platform counts.
- `samples.tsv`: sample identifier, title, source, platform, characteristics JSON, and raw metadata JSON.
- `series_matrix.tsv`: optional matrix created only when one platform is present; values may still require background correction, normalization, probeset-to-gene mapping, or study-specific processing.
- `download/`: GEOparse's local downloaded/cached material.

## Quality gates and pitfalls

- GEO series may contain multiple platforms, superseries/subseries, cell lines, technical replicates, paired samples, or mixed study arms.
- The archive's label is not proof of biological sample identity; verify against publication and supplementary files.
- Avoid comparing normalized microarray values to RNA-seq counts as if they were the same measurement scale.
- Confirm licensing, publication restrictions, and human-data governance before re-distribution.
- The script requires network access and an optional third-party package; network/server errors should be reported rather than guessed around.
