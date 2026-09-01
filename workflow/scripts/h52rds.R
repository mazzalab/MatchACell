#!/usr/bin/env Rscript

# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Utility — H5AD -> RDS converter
# ════════════════════════════════════════════════════════════════════════
#
# Converts an AnnData .h5ad file into a Seurat object serialized as .rds --
# the reverse of workflow/scripts/rds2h5.R.
#
# Reads via zellkonverter::readH5AD(reader = "R") -- the same pure-R h5ad
# parser addmodulescore_tool.R already relies on -- rather than a custom
# rhdf5 parser, so this reads any well-formed .h5ad (scanpy-written or
# rds2h5.R-written), not just this pipeline's own output. Requires the
# matchacell-rule-r conda env (Seurat, SingleCellExperiment, zellkonverter),
# same as addmodulescore_tool.R -- not the plain rds2h5.R deps.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
})

# =====================
# CLI
# =====================
option_list <- list(
  make_option(c("-i", "--input"),  type = "character", help = "Input .h5ad file"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory (a 'rds' subfolder is created inside it)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Both --input and --output are required. Use --help for usage.")
}

# =====================
# Output path
# =====================
rds_dir <- file.path(opt$output, "rds")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
basename_noext <- sub("\\.h5ad$", "", basename(opt$input), ignore.case = TRUE)
rds_path <- file.path(rds_dir, paste0(basename_noext, ".rds"))

# =====================
# Read input
# =====================
message("Reading ", opt$input, " ...")
sce <- zellkonverter::readH5AD(opt$input, reader = "R")

assay_names <- SummarizedExperiment::assayNames(sce)
data_layer   <- if ("X" %in% assay_names) "X" else assay_names[1]
counts_layer <- if ("counts" %in% assay_names) "counts" else data_layer
message(sprintf(
  "Assays found: %s -> counts = '%s', data = '%s'",
  paste(assay_names, collapse = ", "), counts_layer, data_layer
))

# =====================
# SCE -> Seurat (colData -> meta.data, reducedDims -> reductions, both
# handled automatically by as.Seurat)
# =====================
seurat_obj <- Seurat::as.Seurat(sce, counts = counts_layer, data = data_layer)

# =====================
# Write .rds
# =====================
message("Writing ", rds_path, " ...")
saveRDS(seurat_obj, rds_path)
message("Done -> ", rds_path)
