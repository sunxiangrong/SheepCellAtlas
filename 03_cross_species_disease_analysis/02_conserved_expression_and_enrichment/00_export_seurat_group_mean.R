#!/usr/bin/env Rscript
# Adapted from scpagwas-xinong/scripts/export_sheep_seurat_pseudobulk.R.
# Export the mean of the selected RNA layer per group. This is an expression
# summary for correlation, NOT a replicate-aware raw-count pseudobulk for DE.
# Require an explicit layer: do not silently normalize or substitute counts.
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop('Usage: Rscript 00_export_seurat_group_mean.R <input.rds> <group_column> <layer> <output.tsv>')
}
input <- args[[1]]
group_column <- args[[2]]
layer <- args[[3]]
output <- args[[4]]
sidecars <- c(output, paste0(output, '.groups.tsv'), paste0(output, '.session.txt'))
if (any(file.exists(sidecars))) stop('At least one output or sidecar already exists')
object <- readRDS(input)
if (!inherits(object, 'Seurat')) stop('Input must be a Seurat object')
if (!'RNA' %in% Assays(object)) stop('Missing RNA assay')
if (!group_column %in% colnames(object[[]])) stop('Missing grouping metadata column')
if (!layer %in% Layers(object[['RNA']])) stop('Requested exact RNA layer is absent; join split layers explicitly first')
matrix <- GetAssayData(object, assay = 'RNA', layer = layer)
metadata <- object[[]]
if (!all(colnames(matrix) %in% rownames(metadata))) stop('Expression cells missing from metadata')
groups <- as.character(metadata[colnames(matrix), group_column])
valid <- !is.na(groups) & nzchar(trimws(groups))
if (!any(valid)) stop('No cells have valid group labels')
matrix <- matrix[, valid, drop = FALSE]
groups <- groups[valid]
labels <- sort(unique(groups))
if ('gene' %in% labels) stop('Group name gene conflicts with output identifier column')
# vapply preserves the gene-by-group dimensions even for a single group.
average <- vapply(labels, function(label) {
  Matrix::rowMeans(matrix[, groups == label, drop = FALSE])
}, numeric(nrow(matrix)))
colnames(average) <- labels
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(gene = rownames(matrix), average, check.names = FALSE), output,
            sep = '\t', quote = FALSE, row.names = FALSE)
write.table(data.frame(group = labels, n_cells = as.integer(table(factor(groups, levels = labels))),
                       assay = 'RNA', layer = layer, aggregation = 'arithmetic_mean_selected_layer'),
            paste0(output, '.groups.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(output, '.session.txt'))
