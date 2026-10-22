# Methods and references

This file records upstream documentation for the methods and formats used by the project. Confirm versions when you run production analyses; this repository deliberately does not pin a specific reference genome or pathway database for every organism/study.

## Agent Skills format

- [Agent Skills overview and SKILL.md structure](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/overview)
- [Open format specification maintained in Microsoft's skills repository](https://github.com/microsoft/skills/blob/main/.github/docs/skills.md)

## Sequence files and QC

- [NCBI FASTA format help](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
- [FASTQ format overview (Wellcome Sanger Institute)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_help.html)
- [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

The sequence summary parser uses a streaming approach and computes GC as `(G+C)/total bases`; ambiguity is defined as every symbol outside A/C/G/T/U. Quality conversion assumes Phred+33 for FASTQ; older encodings must be investigated from provenance.

## GEO

- [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
- [GEOquery Bioconductor documentation](https://bioconductor.org/packages/GEOquery/)
- [GEOparse documentation](https://geoparse.readthedocs.io/)

GEO is heterogeneous: series design, platform, data processing, and measurement scales vary. `series_matrix.tsv` is named to describe a matrix export, not to promise raw counts.

## Bulk RNA-seq

- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Love, Huber & Anders (2014), Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://doi.org/10.1186/s13059-014-0550-8)
- [tximport Bioconductor package](https://bioconductor.org/packages/tximport/)
- [Salmon documentation](https://salmon.readthedocs.io/)
- [fastp documentation](https://github.com/OpenGene/fastp)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

DESeq2 models raw counts with a negative-binomial generalized linear model and uses dispersion estimation/shrinkage and multiple-testing correction. The provided script uses `DESeq2::results` and BH-adjusted `padj`; it does not automatically do a full contrast sensitivity analysis or perform outlier replacement beyond DESeq2 defaults. Salmon output requires transcript-to-gene mapping matched to the reference.

## Single-cell RNA-seq

- [Scanpy documentation](https://scanpy.readthedocs.io/)
- [Scanpy PBMC 3k tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html)
- [AnnData documentation](https://anndata.readthedocs.io/)

Single-cell thresholds are dataset-specific. Donor-aware inference, doublet detection, ambient RNA correction, integration, and annotation require deliberate additional steps.

## Pathway over-representation

ORA computes a one-sided hypergeometric test of overlap against an explicit gene universe. Multiple-testing correction is Benjamini–Hochberg. Reference for multiple-testing adjustment: [Benjamini & Hochberg (1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). The method is sensitive to gene universe, term sizes, identifier mapping, database choice, and redundancy among related gene sets.
