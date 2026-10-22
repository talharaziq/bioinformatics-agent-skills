#!/usr/bin/env Rscript
# Differential expression from a raw gene-by-sample integer count matrix.
# Usage: Rscript scripts/deseq2_counts.R counts.tsv metadata.tsv outdir condition_col reference test [batch_col]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6 || length(args) > 7) {
  stop("Usage: Rscript scripts/deseq2_counts.R counts.tsv metadata.tsv outdir condition_col reference test [batch_col]", call. = FALSE)
}
counts_path <- args[[1]]; metadata_path <- args[[2]]; outdir <- args[[3]]
condition_col <- args[[4]]; reference <- args[[5]]; test <- args[[6]]
batch_col <- if (length(args) == 7) args[[7]] else NULL
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 from Bioconductor before running this script.", call. = FALSE)
if (!identical(make.names(condition_col), condition_col)) stop("Use an R-safe metadata column name for condition_col (letters, digits, dot, underscore).", call. = FALSE)
if (!is.null(batch_col) && !identical(make.names(batch_col), batch_col)) stop("Use an R-safe metadata column name for batch_col.", call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Read identifiers as character to preserve leading zeros and versioned IDs.
counts_df <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
if (ncol(counts_df) < 3) stop("Counts table needs gene_id plus at least two samples.", call. = FALSE)
gene_ids <- trimws(counts_df[[1]])
if (anyNA(gene_ids) || any(!nzchar(gene_ids))) stop("Empty gene_id values found.", call. = FALSE)
if (anyDuplicated(gene_ids)) stop("Duplicate gene IDs found. Resolve duplicates before DE.", call. = FALSE)
count_mat <- as.matrix(counts_df[, -1, drop = FALSE]); rownames(count_mat) <- gene_ids
if (anyDuplicated(colnames(count_mat))) stop("Duplicate sample columns in count matrix.", call. = FALSE)
suppressWarnings(storage.mode(count_mat) <- "numeric")
if (anyNA(count_mat) || any(!is.finite(count_mat)) || any(count_mat < 0)) stop("Counts must be finite and non-negative.", call. = FALSE)
if (any(abs(count_mat - round(count_mat)) > 1e-8)) stop("DESeq2 requires raw integer counts; fractional/normalized values were detected.", call. = FALSE)
count_mat <- round(count_mat)

meta <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
if (!"sample_id" %in% names(meta)) stop("Metadata must contain sample_id.", call. = FALSE)
meta$sample_id <- trimws(meta$sample_id)
if (anyNA(meta$sample_id) || any(!nzchar(meta$sample_id))) stop("Metadata contains an empty sample_id.", call. = FALSE)
if (anyDuplicated(meta$sample_id)) stop("Duplicate sample_id in metadata.", call. = FALSE)
if (!condition_col %in% names(meta)) stop(paste("Metadata lacks condition column:", condition_col), call. = FALSE)
if (anyNA(meta[[condition_col]]) || any(!nzchar(trimws(meta[[condition_col]])))) stop("Missing condition labels in metadata.", call. = FALSE)
if (!setequal(colnames(count_mat), meta$sample_id)) stop("Count matrix sample columns and metadata sample_id values do not match exactly.", call. = FALSE)
meta <- meta[match(colnames(count_mat), meta$sample_id), , drop = FALSE]
meta[[condition_col]] <- trimws(meta[[condition_col]])
col_data <- meta[, condition_col, drop = FALSE]
rownames(col_data) <- meta$sample_id
col_data[[condition_col]] <- factor(col_data[[condition_col]])
if (!all(c(reference, test) %in% levels(col_data[[condition_col]]))) stop("Reference/test must match observed condition labels exactly.", call. = FALSE)
col_data[[condition_col]] <- relevel(col_data[[condition_col]], ref = reference)
if (sum(col_data[[condition_col]] == reference) < 2 || sum(col_data[[condition_col]] == test) < 2) warning("One or both groups have fewer than two samples; estimates/inference are likely unreliable.")
if (is.null(batch_col)) design_formula <- as.formula(paste("~", condition_col)) else {
  if (!batch_col %in% names(meta)) stop(paste("Metadata lacks batch column:", batch_col), call. = FALSE)
  if (anyNA(meta[[batch_col]]) || any(!nzchar(trimws(meta[[batch_col]])))) stop("Missing batch labels in metadata.", call. = FALSE)
  meta[[batch_col]] <- trimws(meta[[batch_col]])
  col_data[[batch_col]] <- factor(meta[[batch_col]])
  design_formula <- as.formula(paste("~", batch_col, "+", condition_col))
}
model_matrix <- model.matrix(design_formula, col_data)
if (qr(model_matrix)$rank < ncol(model_matrix)) stop("Design matrix is not full rank; inspect confounded condition/batch labels or insufficient replication.", call. = FALSE)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat, colData = col_data, design = design_formula)
dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
if (nrow(dds) == 0) stop("No genes remain after the low-count filter.", call. = FALSE)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
res <- DESeq2::results(dds, contrast = c(condition_col, test, reference), alpha = 0.05)
res_df <- as.data.frame(res); res_df$gene_id <- rownames(res_df); res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
write.table(res_df, file = file.path(outdir, "differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
norm <- as.data.frame(DESeq2::counts(dds, normalized = TRUE)); norm <- cbind(gene_id = rownames(norm), norm)
write.table(norm, file = file.path(outdir, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tryCatch({
  vsd <- DESeq2::vst(dds, blind = FALSE)
  pca <- DESeq2::plotPCA(vsd, intgroup = condition_col, returnData = TRUE)
  write.table(cbind(sample_id = rownames(pca), pca), file = file.path(outdir, "pca_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}, error = function(e) message("PCA export skipped: ", conditionMessage(e)))
summary_lines <- c("# DESeq2 differential expression", "", paste0("- Reference: `", reference, "`"), paste0("- Test: `", test, "`"), paste0("- Design: `", deparse(design_formula), "`"), paste0("- Samples: ", ncol(dds)), paste0("- Genes tested: ", nrow(dds)), paste0("- FDR < 0.05: ", sum(res$padj < 0.05, na.rm = TRUE)), "", "Adjusted p-values use the Benjamini-Hochberg procedure. Inspect PCA, sample-level QC, design assumptions, and effect sizes before interpreting hits.", "")
writeLines(summary_lines, file.path(outdir, "report.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Finished. Results: ", file.path(outdir, "differential_expression.tsv"))
