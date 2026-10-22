#!/usr/bin/env Rscript
# Transcript quantifications from Salmon -> tximport -> gene-level DESeq2.
# Usage: Rscript scripts/tximport_deseq2.R manifest.tsv tx2gene.tsv outdir condition_col reference test [batch_col]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6 || length(args) > 7) stop("Usage: Rscript scripts/tximport_deseq2.R manifest.tsv tx2gene.tsv outdir condition_col reference test [batch_col]", call. = FALSE)
manifest_path <- args[[1]]; tx2gene_path <- args[[2]]; outdir <- args[[3]]
condition_col <- args[[4]]; reference <- args[[5]]; test <- args[[6]]
batch_col <- if (length(args) == 7) args[[7]] else NULL
for (pkg in c("tximport", "DESeq2")) if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("Missing Bioconductor package:", pkg), call. = FALSE)
if (!identical(make.names(condition_col), condition_col)) stop("Use an R-safe metadata column name for condition_col.", call. = FALSE)
if (!is.null(batch_col) && !identical(make.names(batch_col), batch_col)) stop("Use an R-safe metadata column name for batch_col.", call. = FALSE)
manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
tx2gene <- read.delim(tx2gene_path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
if (!all(c("sample_id", condition_col) %in% names(manifest))) stop("Manifest must contain sample_id and the requested condition column.", call. = FALSE)
if (!all(c("transcript_id", "gene_id") %in% names(tx2gene))) stop("tx2gene TSV must have transcript_id and gene_id columns.", call. = FALSE)
manifest$sample_id <- trimws(manifest$sample_id)
if (anyNA(manifest$sample_id) || any(!nzchar(manifest$sample_id))) stop("Empty sample_id in manifest.", call. = FALSE)
if (anyDuplicated(manifest$sample_id)) stop("Duplicate sample_id in manifest.", call. = FALSE)
if (anyDuplicated(tx2gene$transcript_id)) stop("tx2gene transcript_id values must be unique.", call. = FALSE)
if (any(!nzchar(tx2gene$transcript_id)) || any(!nzchar(tx2gene$gene_id))) stop("Empty transcript_id/gene_id in tx2gene.", call. = FALSE)
quant_files <- file.path(outdir, "salmon", manifest$sample_id, "quant.sf")
if (!all(file.exists(quant_files))) stop(paste("Missing Salmon quant.sf for:", paste(manifest$sample_id[!file.exists(quant_files)], collapse = ", ")), call. = FALSE)
names(quant_files) <- manifest$sample_id
if (anyNA(manifest[[condition_col]]) || any(!nzchar(trimws(manifest[[condition_col]])))) stop("Missing condition labels in manifest.", call. = FALSE)
manifest[[condition_col]] <- trimws(manifest[[condition_col]])
col_data <- manifest[, condition_col, drop = FALSE]
rownames(col_data) <- manifest$sample_id
col_data[[condition_col]] <- factor(col_data[[condition_col]])
if (!all(c(reference, test) %in% levels(col_data[[condition_col]]))) stop("Reference/test labels not present in the condition column.", call. = FALSE)
col_data[[condition_col]] <- relevel(col_data[[condition_col]], ref = reference)
if (sum(col_data[[condition_col]] == reference) < 2 || sum(col_data[[condition_col]] == test) < 2) warning("One or both groups have fewer than two samples; inference may be unreliable.")
if (is.null(batch_col)) design <- as.formula(paste("~", condition_col)) else {
  if (!batch_col %in% names(manifest)) stop(paste("Requested batch column missing:", batch_col), call. = FALSE)
  if (anyNA(manifest[[batch_col]]) || any(!nzchar(trimws(manifest[[batch_col]])))) stop("Missing batch labels in manifest.", call. = FALSE)
  manifest[[batch_col]] <- trimws(manifest[[batch_col]])
  col_data[[batch_col]] <- factor(manifest[[batch_col]])
  design <- as.formula(paste("~", batch_col, "+", condition_col))
}
mm <- model.matrix(design, col_data)
if (qr(mm)$rank < ncol(mm)) stop("Design is not full rank; likely condition/batch confounding.", call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
txi <- tximport::tximport(quant_files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")], countsFromAbundance = "lengthScaledTPM")
dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = as.data.frame(col_data), design = design)
dds <- dds[rowSums(DESeq2::counts(dds)) >= 10, ]
if (!nrow(dds)) stop("No genes remain after low-count filtering.", call. = FALSE)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
res <- DESeq2::results(dds, contrast = c(condition_col, test, reference), alpha = 0.05)
res_df <- as.data.frame(res); res_df$gene_id <- rownames(res_df); res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]
write.table(res_df, file.path(outdir, "differential_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
norm <- as.data.frame(DESeq2::counts(dds, normalized = TRUE)); norm <- cbind(gene_id = rownames(norm), norm)
write.table(norm, file.path(outdir, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tryCatch({
  vsd <- DESeq2::vst(dds, blind = FALSE)
  pca <- DESeq2::plotPCA(vsd, intgroup = condition_col, returnData = TRUE)
  write.table(cbind(sample_id = rownames(pca), pca), file.path(outdir, "pca_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}, error = function(e) message("PCA export skipped: ", conditionMessage(e)))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
writeLines(c("# Salmon + tximport + DESeq2", "", paste0("Design: `", deparse(design), "`"), paste0("Samples: ", ncol(dds)), paste0("Genes tested: ", nrow(dds)), paste0("FDR < 0.05: ", sum(res$padj < 0.05, na.rm = TRUE)), "", "Confirm that the provided transcript-to-gene mapping matches the exact transcriptome used to build the Salmon index."), file.path(outdir, "report.md"))
message("Finished. Results: ", file.path(outdir, "differential_expression.tsv"))
