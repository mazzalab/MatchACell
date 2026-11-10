#!/usr/bin/env Rscript

# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 1 (optional) — doublet detection with scDblFinder
# ════════════════════════════════════════════════════════════════════════
#
# Reads the raw counts of a staged .h5ad (.X, or a layer such as "counts" for
# normalized input), runs scDblFinder (Germain et al., F1000Research 2021) and
# writes one row per cell: cell, doublet_score, doublet_class. Step 1 then
# removes (or flags) the doublets as one more QC criterion.
#
# Doublets only form within one capture (10x lane/channel): with several
# captures in one file, give the .obs column naming them (--samples-col), so
# scDblFinder simulates doublets within each capture only.
#
# The .h5ad is read with rhdf5 following AnnData's on-disk format (CSR/CSC or
# dense matrices, string or categorical .obs columns), so this needs neither
# Python nor zellkonverter.

suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Staged .h5ad"),
  make_option(c("-o", "--output"), type = "character", help = "Output CSV"),
  make_option("--layer", type = "character", default = "",
              help = "Layer with raw counts; empty uses .X [default: .X]"),
  make_option("--samples-col", type = "character", default = "",
              help = ".obs column with the capture of each cell (optional)"),
  make_option("--dbr", type = "double", default = NA,
              help = "Expected doublet rate; default: scDblFinder's (1%% per 1,000 cells)"),
  make_option("--threads", type = "integer", default = 1),
  make_option("--seed", type = "integer", default = 0)
)
opt <- parse_args(OptionParser(option_list = option_list))

read_matrix <- function(file, path) {
  # genes x cells dgCMatrix from an AnnData matrix stored at `path`
  attrs <- h5readAttributes(file, path)
  enc <- attrs[["encoding-type"]]
  if (!is.null(enc) && enc %in% c("csr_matrix", "csc_matrix")) {
    shape <- as.integer(attrs[["shape"]])            # n_obs, n_vars
    rd <- function(x) h5read(file, paste0(path, "/", x), bit64conversion = "double")
    data <- as.numeric(rd("data"))
    indices <- as.integer(rd("indices"))
    indptr <- as.numeric(rd("indptr"))
    if (enc == "csr_matrix") {
      # rows are cells: read as the CSC of its transpose, genes x cells
      return(sparseMatrix(i = indices, p = indptr, x = data,
                          dims = c(shape[2], shape[1]), index1 = FALSE))
    }
    m <- sparseMatrix(i = indices, p = indptr, x = data,
                      dims = c(shape[1], shape[2]), index1 = FALSE)
    return(t(m))
  }
  # dense: HDF5 is row-major, so rhdf5 already returns genes x cells
  as(h5read(file, path), "CsparseMatrix")
}

read_index <- function(file, group) {
  idx <- h5readAttributes(file, group)[["_index"]]
  as.character(h5read(file, paste0(group, "/", idx)))
}

read_obs_column <- function(file, col) {
  path <- paste0("obs/", col)
  info <- h5ls(file, recursive = TRUE)
  is_group <- any(info$group == paste0("/", path))
  if (is_group) {                                     # categorical
    codes <- as.integer(h5read(file, paste0(path, "/codes")))
    cats <- as.character(h5read(file, paste0(path, "/categories")))
    out <- cats[codes + 1L]
    out[codes < 0] <- NA
    return(out)
  }
  as.character(h5read(file, path))
}

file <- opt$input
mpath <- if (nzchar(opt$layer)) paste0("layers/", opt$layer) else "X"
counts <- read_matrix(file, mpath)
colnames(counts) <- read_index(file, "obs")
rownames(counts) <- make.unique(read_index(file, "var"))
message(sprintf("[scDblFinder] %d cells x %d genes from %s", ncol(counts), nrow(counts), mpath))

nz <- counts@x[seq_len(min(length(counts@x), 1e6))]
if (any(nz != round(nz))) {
  stop("The matrix at '", mpath, "' is not raw counts (non-integer values). ",
       "Point --layer at the raw counts.")
}

sce <- SingleCellExperiment(list(counts = counts))
args <- list(sce = sce)
if (nzchar(opt$`samples-col`)) {
  args$samples <- read_obs_column(file, opt$`samples-col`)
  if (anyNA(args$samples)) stop("--samples-col ", opt$`samples-col`, " has missing values")
  message(sprintf("[scDblFinder] %d captures in '%s'", length(unique(args$samples)),
                  opt$`samples-col`))
  args$BPPARAM <- if (opt$threads > 1) MulticoreParam(opt$threads, RNGseed = opt$seed)
                  else SerialParam(RNGseed = opt$seed)
}
if (!is.na(opt$dbr)) args$dbr <- opt$dbr

set.seed(opt$seed)
sce <- do.call(scDblFinder, args)

out <- data.frame(cell = colnames(sce),
                  doublet_score = sce$scDblFinder.score,
                  doublet_class = as.character(sce$scDblFinder.class))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
write.csv(out, opt$output, row.names = FALSE)
n <- sum(out$doublet_class == "doublet")
message(sprintf("[scDblFinder] %d doublets of %d cells (%.1f%%) -> %s",
                n, nrow(out), 100 * n / nrow(out), opt$output))
