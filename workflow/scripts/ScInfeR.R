#!/usr/bin/env Rscript

# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: ScInfeR (reference mapping)
# ════════════════════════════════════════════════════════════════════════
#
# Wraps swainasish/ScInfeR (github.com/swainasish/ScInfeR) into the same
# CLI/input shape as the other annotators (addmodulescore_tool.R, scanvi_tool.py):
# Step 1's clustered_multi_resolution.rds + MatchA_Verdict.txt in, an
# annotated .rds out for rds2h5.R to convert.
#
#### WARNING !!
## If `remotes::install_github("swainasish/ScInfeR")` fails on git credentials,
## run `gitcreds::gitcreds_delete()` first, then retry the install.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ScInfeR)
  library(openxlsx)
  library(ggplot2)
})

# =====================
# CLI
# =====================
option_list <- list(
  make_option(c("-i", "--input"),    type = "character", help = "Input Seurat object (clustered_multi_resolution.rds from Step 1)"),
  make_option(c("-r", "--reference"), type = "character", help = "Annotated single-cell reference Seurat object (.rds)"),
  make_option(c("--ref_label_col"),   type = "character", default = "cell_type", help = "Column in the reference metadata holding cell type labels [default: %default]"),
  make_option(c("--tissue"),          type = "character", default = "", help = "Tissue type passed to ScInfeR::fetch_markerset() to add curated signature markers on top of the reference-derived ones; empty skips it"),
  make_option(c("--num_marker_per_ct"), type = "integer", default = 10, help = "Markers per cell type drawn from the reference [default: %default]"),
  make_option(c("--own_weightage"),   type = "double",  default = 0.5, help = "Weight given to the user's own markers vs. the reference-derived ones [default: %default]"),
  make_option(c("--n_neighbor"),      type = "integer", default = 10, help = "Neighbors used by ScInfeR's prediction graph [default: %default]"),
  make_option(c("-v", "--verdict"),  type = "character", help = "Verdict file"),
  make_option(c("-o", "--output"),   type = "character", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

# =====================
# Extract best res
# =====================
extract_best_res <- function(verdict_file) {
  txt <- paste(readLines(verdict_file, warn = FALSE), collapse = "\n")
  m <- regmatches(txt, regexpr("Recommended resolution\\s*:\\s*\\S+", txt))
  if (length(m) == 0) {
    stop(sprintf("'Recommended resolution' not found in %s", verdict_file))
  }
  trimws(sub(".*:", "", m))
}

leiden_col <- extract_best_res(opt$verdict)

# =====================
# Read input (.rds, already a Seurat object -- see workflow/scripts/h52rds.R)
# =====================
message("Reading ", opt$input, " ...")
seurat_obj <- readRDS(opt$input)
seurat_obj[[leiden_col]] <- as.character(seurat_obj[[leiden_col, drop = TRUE]])
Idents(seurat_obj) <- leiden_col

# =====================
# Reference single-cell dataset
# =====================
message("Reading reference ", opt$reference, " ...")
reference_obj <- readRDS(opt$reference)

if (!opt$ref_label_col %in% colnames(reference_obj[[]])) {
  stop(sprintf("Column '%s' not found in reference metadata.", opt$ref_label_col))
}

umap_candidates <- c("X_umap", "umap", "UMAP")
ref_umap_reduction <- Find(function(r) r %in% names(reference_obj@reductions), umap_candidates)
if (is.null(ref_umap_reduction)) {
  stop(
    "No UMAP reduction found in reference among ", paste(umap_candidates, collapse = ", "),
    " -- ScInfeR::get_marker_from_ref_matrix() needs one."
  )
}

ref_assay <- Seurat::DefaultAssay(reference_obj)
ref_counts <- SeuratObject::LayerData(reference_obj[[ref_assay]], layer = "counts")

reference_marker <- ScInfeR::get_marker_from_ref_matrix(
  exp_mat = ref_counts,
  annotations = reference_obj[[opt$ref_label_col, drop = TRUE]],
  umap_cor = Seurat::Embeddings(reference_obj, reduction = ref_umap_reduction),
  num_marker_per_ct = opt$num_marker_per_ct
)

# =====================
# Optional curated signature markers, merged with the reference-derived ones
# =====================
if (nzchar(opt$tissue)) {
  message("Fetching curated marker set for tissue: ", opt$tissue)
  scinfer_marker_signature <- ScInfeR::fetch_markerset(tissue_type = opt$tissue)
  marker_combined <- rbind(scinfer_marker_signature, reference_marker)
} else {
  marker_combined <- reference_marker
}

# =====================
# Output folders
# =====================
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =====================
# Predict cell types on the query object
# =====================
# ScInfeR::predict_celltype_scRNA_seurat()'s `reduction` arg defaults to
# "umap" -- after h52rds.R's SCE -> Seurat conversion the actual name is
# "X_umap" (see addmodulescore_tool.R's note on the same quirk), so detect
# it the same way before predicting, not just for the plot below.
query_umap_reduction <- Find(function(r) r %in% names(seurat_obj@reductions), umap_candidates)
if (is.null(query_umap_reduction)) {
  query_umap_reduction <- names(seurat_obj@reductions)[1]
  warning(
    "No UMAP reduction found among ", paste(umap_candidates, collapse = ", "),
    " -- falling back to '", query_umap_reduction, "' (may not be a 2D UMAP layout)."
  )
}

query_assay <- Seurat::DefaultAssay(seurat_obj)
message("Running ScInfeR::predict_celltype_scRNA_seurat() ...")
scinfer_pred <- ScInfeR::predict_celltype_scRNA_seurat(
  s_object = seurat_obj,
  group_annt = seurat_obj[[leiden_col, drop = TRUE]],
  ct_marker_df = marker_combined,
  subtype_present = FALSE,
  subtype_info = FALSE,
  assay_name = query_assay,
  slot_name = "counts",
  reduction = query_umap_reduction,
  own_weightage = opt$own_weightage,
  n_neighbor = opt$n_neighbor
)
seurat_obj$cell_type_pred <- scinfer_pred$celltype

# =====================
# UMAP
# =====================
p_annot <- DimPlot(seurat_obj, group.by = "cell_type_pred", reduction = query_umap_reduction) +
  ggtitle("ScInfeR annotation (cell_type_pred)")
ggsave(
  file.path(output_dir, paste0(leiden_col, "_cell_type_pred_umap.png")),
  p_annot, width = 7, height = 5, dpi = 150
)

# =====================
# Cluster -> cell type mapping table
# =====================
message("Creating cluster-to-celltype mapping table...")

count_mat <- table(
  cluster = seurat_obj[[leiden_col, drop = TRUE]],
  cell_type_pred = seurat_obj$cell_type_pred
)
count_cols <- colnames(count_mat)
df_unified <- data.frame(cluster = rownames(count_mat), unclass(count_mat), check.names = FALSE)

percentages <- as.matrix(count_mat)
percentages <- percentages / rowSums(percentages) * 100

df_unified$dominant_cell_type <- count_cols[apply(percentages, 1, which.max)]
df_unified$dominant_pct <- apply(percentages, 1, max)

openxlsx::write.xlsx(
  df_unified,
  file.path(output_dir, "cluster_to_celltype_mapping.xlsx")
)

# =====================
# Save the annotated Seurat object -- workflow/rules/scinfer.smk runs
# rds2h5.R on this right after, same as addmodulescore.smk does.
# =====================
rds_dir <- file.path(output_dir, "rds")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(rds_dir, "scinfer_annotated.rds")
saveRDS(seurat_obj, rds_path)
message("Saved annotated Seurat object -> ", rds_path)
