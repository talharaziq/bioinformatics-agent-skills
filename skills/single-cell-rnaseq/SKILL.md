---
name: single-cell-rnaseq
description: Perform cautious single-cell RNA-seq preprocessing with Scanpy on a 10x matrix directory, 10x H5 file, or AnnData H5AD object. Use when asked for single-cell QC, filtering, normalization, highly variable genes, PCA, neighborhood graphs, Leiden clusters, UMAP, or a saved processed AnnData object.
---

# Single-cell RNA-seq preprocessing

## Procedure

1. Inspect organism, chemistry, cell input, sample/donor identifiers, sequencing depth, experimental conditions, and whether input counts are raw. Never merge donors or conditions into one anonymous sample without recording sample IDs.
2. Load the 10x matrix directory, 10x H5, or `.h5ad` file. Confirm dimensions, feature identifiers, count provenance, and that genes/cells are unique.
3. Calculate QC metrics and review distributions for genes per cell, total counts, mitochondrial fraction, and possible empty droplets/doublets. The default `MT-` mitochondrial symbol prefix may be wrong for non-human/mouse organisms or non-standard annotations.
4. Choose filtering thresholds based on the actual data and sample. Defaults are starting points only, not universal rules. Record thresholds and cell counts before/after each filter.
5. Preserve a raw-count layer, then normalize total counts and apply `log1p`; identify highly variable genes, scale, compute PCA/neighbors, and calculate UMAP. Leiden clustering is attempted when its optional dependencies are available.
6. Save the processed `.h5ad`, QC table, summary JSON, plot if available, environment versions, and parameter record.
7. Before biological interpretation, assess sample/donor effects, doublets, ambient RNA, cell-cycle effects, and integration strategy. Differential expression should account for biological replication (often sample-level pseudobulk); cells from one donor are not independent replicates.

## Run

```bash
python -m pip install -e '.[single-cell]'
bio-skills single-cell /path/to/filtered_feature_bc_matrix --outdir results/scRNA
# Set data-dependent thresholds explicitly when warranted:
bio-skills single-cell input.h5 --outdir results/scRNA --min-genes 300 --min-cells 5 --max-pct-mito 15
```

The helper is `skills/single-cell-rnaseq/scripts/run_scanpy.py`.

## Outputs

`cell_qc.tsv` records pre-filter QC metrics plus a retained/not-retained flag; `processed.h5ad` contains a raw-count layer and normalized/log-transformed expression; `single_cell_summary.json`; and `umap.png` when plotting dependencies succeed. UMAP is a visualization, not evidence that clusters are discrete cell types.

## Limitations

This skill does not perform doublet detection, SoupX/ambient RNA correction, integration across studies, automated cell type annotation, trajectory inference, or donor-aware pseudobulk differential testing. Do not present a cluster label as ground truth. Report filtering losses and assess parameter sensitivity.
