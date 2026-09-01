#!/usr/bin/env Rscript

# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Utility — RDS -> H5AD converter
# ════════════════════════════════════════════════════════════════════════
#
# Converts a serialized R object (.rds) holding a Seurat or
# SingleCellExperiment object into an AnnData .h5ad file, the format the
# rest of the pipeline reads (addmodulescore_tool.R's zellkonverter
# consumer, score_genes.py, celltypist_tool.py, ...).
#
# Implementation note: SeuratDisk (the "usual" Seurat -> h5ad bridge) was
# tried first and dropped -- it chokes on Seurat v5 Assay5 objects and on
# spatial/FOV data ("guess_dtype: unknown type") and is unmaintained. This
# writes the .h5ad directly with rhdf5, following AnnData's on-disk
# dataframe/categorical/dense-array encoding (see
# https://anndata.readthedocs.io/en/latest/fileformat-prose.html), which
# keeps the dependency footprint to what the pipeline's R env already has
# (Seurat, rhdf5, SingleCellExperiment, optparse) and was round-trip
# verified against both python's anndata.read_h5ad and
# zellkonverter::readH5AD(reader = "R") -- the actual reader
# addmodulescore_tool.R uses.
#
# Only dense arrays are written (no CSR/CSC sparse encoding): fine for
# targeted panels (Xenium, ...) with a few hundred genes: not meant for
# full-transcriptome assays with tens of thousands of genes x many cells.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SingleCellExperiment)
  library(rhdf5)
})

# =====================
# CLI
# =====================
option_list <- list(
  make_option(c("-i", "--input"),  type = "character", help = "Input .rds file (Seurat or SingleCellExperiment object)"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory (a 'h5' subfolder is created inside it)"),
  make_option(c("-a", "--assay"),  type = "character", default = NA, help = "Assay to export [default: object's DefaultAssay]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Both --input and --output are required. Use --help for usage.")
}

# =====================
# rhdf5 helpers -- AnnData on-disk encoding
# =====================
# rhdf5::h5writeDataset writes an R matrix of dim (a, b) as an HDF5 dataset
# of shape (b, a) -- a true transpose, reconciling R's column-major layout
# with HDF5/numpy's row-major one. So to land on-disk in AnnData's expected
# (n_obs, n_var) shape, every dense matrix below must be *fed in* as
# (n_var, n_obs) -- i.e. Seurat's native gene x cell layout, unchanged.

write_dense <- function(mat, h5loc, name) {
  h5writeDataset(mat, h5loc, name)
  did <- H5Dopen(h5loc, name)
  h5writeAttribute("array", did, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", did, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  H5Dclose(did)
}

write_index <- function(names_vec, h5loc) {
  h5writeDataset(names_vec, h5loc, "_index")
  did <- H5Dopen(h5loc, "_index")
  h5writeAttribute("string-array", did, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", did, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  H5Dclose(did)
}

write_string_col <- function(vec, h5loc, name) {
  h5writeDataset(as.character(vec), h5loc, name)
  did <- H5Dopen(h5loc, name)
  h5writeAttribute("string-array", did, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", did, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  H5Dclose(did)
}

write_array_col <- function(vec, h5loc, name) {
  h5writeDataset(vec, h5loc, name)
  did <- H5Dopen(h5loc, name)
  h5writeAttribute("array", did, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", did, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  H5Dclose(did)
}

write_categorical_col <- function(fac, h5loc, name) {
  sub <- H5Gcreate(h5loc, name)
  h5writeAttribute("categorical", sub, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", sub, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute(0L, sub, "ordered", asScalar = TRUE)
  cats <- levels(fac)
  codes <- as.integer(fac) - 1L
  codes[is.na(codes)] <- -1L
  h5writeDataset(cats, sub, "categories")
  cid <- H5Dopen(sub, "categories")
  h5writeAttribute("string-array", cid, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", cid, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  H5Dclose(cid)
  h5writeDataset(codes, sub, "codes")
  coid <- H5Dopen(sub, "codes")
  h5writeAttribute("array", coid, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", coid, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  H5Dclose(coid)
  H5Gclose(sub)
}

# cols: named list of R vectors (factor / character / numeric / integer / logical)
write_dataframe <- function(h5loc, name, index_vec, cols) {
  g <- H5Gcreate(h5loc, name)
  h5writeAttribute("dataframe", g, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.2.0", g, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("_index", g, "_index", variableLengthString = TRUE, asScalar = TRUE)
  col_names <- names(cols)
  if (length(col_names) == 0) {
    # rhdf5 chokes on a zero-length *character* attribute; anndata itself
    # writes an empty column-order as a plain empty (float64) array when
    # there are no extra columns, so match that instead of forcing string.
    h5writeAttribute(double(0), g, "column-order")
  } else {
    h5writeAttribute(col_names, g, "column-order", variableLengthString = TRUE)
  }
  write_index(index_vec, g)
  for (cn in col_names) {
    v <- cols[[cn]]
    if (is.factor(v)) {
      write_categorical_col(v, g, cn)
    } else if (is.character(v)) {
      write_string_col(v, g, cn)
    } else if (is.logical(v)) {
      write_array_col(as.integer(v), g, cn)
    } else {
      write_array_col(v, g, cn)
    }
  }
  H5Gclose(g)
}

write_dict_group <- function(h5loc, name) {
  g <- H5Gcreate(h5loc, name)
  h5writeAttribute("dict", g, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
  h5writeAttribute("0.1.0", g, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)
  g
}

# =====================
# Read input
# =====================
message("Reading ", opt$input, " ...")
obj <- readRDS(opt$input)

if (is(obj, "SingleCellExperiment")) {
  message("Input is a SingleCellExperiment -> converting to Seurat ...")
  has_logcounts <- "logcounts" %in% SummarizedExperiment::assayNames(obj)
  obj <- Seurat::as.Seurat(obj, counts = "counts", data = if (has_logcounts) "logcounts" else "counts")
} else if (!is(obj, "Seurat")) {
  stop("Unsupported object class '", paste(class(obj), collapse = ", "), "' -- expected Seurat or SingleCellExperiment.")
}

assay_to_use <- if (!is.na(opt$assay)) opt$assay else Seurat::DefaultAssay(obj)
message("Using assay: ", assay_to_use)

assay_obj <- obj[[assay_to_use]]
avail_layers <- SeuratObject::Layers(assay_obj)

counts_mat <- as.matrix(SeuratObject::LayerData(assay_obj, layer = "counts"))  # genes x cells
data_layer <- if ("data" %in% avail_layers) "data" else "counts"
data_mat   <- as.matrix(SeuratObject::LayerData(assay_obj, layer = data_layer))  # genes x cells

var_names <- rownames(assay_obj)
obs_names <- colnames(obj)
message(sprintf("Exporting %d genes x %d cells (X = '%s' layer)", length(var_names), length(obs_names), data_layer))

# =====================
# obs columns from meta.data
# =====================
obs_cols <- as.list(obj@meta.data)

# =====================
# Reductions -> obsm (dims x cells on write, per the transpose note above)
# =====================
reduction_names <- names(obj@reductions)

# =====================
# Output path
# =====================
h5_dir <- file.path(opt$output, "h5")
dir.create(h5_dir, recursive = TRUE, showWarnings = FALSE)
basename_noext <- sub("\\.rds$", "", basename(opt$input), ignore.case = TRUE)
h5ad_path <- file.path(h5_dir, paste0(basename_noext, ".h5ad"))
unlink(h5ad_path)

# =====================
# Write .h5ad
# =====================
message("Writing ", h5ad_path, " ...")
h5createFile(h5ad_path)
fid <- H5Fopen(h5ad_path)
h5writeAttribute("anndata", fid, "encoding-type", variableLengthString = TRUE, asScalar = TRUE)
h5writeAttribute("0.1.0", fid, "encoding-version", variableLengthString = TRUE, asScalar = TRUE)

write_dense(data_mat, fid, "X")

layers_g <- write_dict_group(fid, "layers")
write_dense(counts_mat, layers_g, "counts")
H5Gclose(layers_g)

write_dataframe(fid, "obs", obs_names, obs_cols)
write_dataframe(fid, "var", var_names, list())

obsm_g <- write_dict_group(fid, "obsm")
for (red in reduction_names) {
  emb <- Seurat::Embeddings(obj, reduction = red)  # cells x dims
  write_dense(t(emb), obsm_g, paste0("X_", red))   # dims x cells on write -> (n_obs, n_dims) on disk
}
H5Gclose(obsm_g)

for (nm in c("varm", "obsp", "varp")) {
  g <- write_dict_group(fid, nm)
  H5Gclose(g)
}

uns_g <- write_dict_group(fid, "uns")
H5Gclose(uns_g)

H5Fclose(fid)
message("Done -> ", h5ad_path)
