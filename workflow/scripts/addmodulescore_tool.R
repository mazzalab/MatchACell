#!/usr/bin/env Rscript

# ════════════════════════════════════════════════════════════════════════
#  MatchACell · Step 2 — Cell Type Annotation: AddModuleScore (Seurat)
# ════════════════════════════════════════════════════════════════════════
#
# Same control-gene-binning algorithm as workflow/scripts/score_genes.py
# (Scanpy's sc.tl.score_genes), run instead through Seurat::AddModuleScore.
# Kept as a separate R script (rather than a Python re-implementation) so
# the two independently confirm each other rather than sharing one
# implementation's blind spots.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

# =====================
# CLI
# =====================
option_list <- list(
  make_option(c("-i", "--input"),   type = "character", help = "Input Seurat object (clustered_multi_resolution.rds from Step 1)"),
  make_option(c("-a", "--annot"),   type = "character", help = "Excel file containing gene signatures"),
  make_option(c("-t", "--thr"),     type = "double",     help = "Threshold on module score"),
  make_option(c("-v", "--verdict"), type = "character", help = "Verdict file"),
  make_option(c("-o", "--output"),  type = "character", help = "Output directory")
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
# Read input (.rds, already a Seurat object -- see workflow/scripts/h52rds.R,
# which did the same zellkonverter/as.Seurat conversion addmodulescore_tool.R
# used to do here itself)
# =====================
message("Reading ", opt$input, " ...")
seurat_obj <- readRDS(opt$input)
seurat_obj[[leiden_col]] <- as.character(seurat_obj[[leiden_col, drop = TRUE]])
Idents(seurat_obj) <- leiden_col

# =====================
# Gene signatures
# =====================
gene_universe <- rownames(seurat_obj)
sheets <- openxlsx::getSheetNames(opt$annot)

signature_list <- list()
for (sheet in sheets) {
  df <- openxlsx::read.xlsx(opt$annot, sheet = sheet, colNames = FALSE)
  genes <- trimws(as.character(df[[1]]))
  genes <- genes[genes %in% gene_universe]
  if (length(genes) > 0) {
    signature_list[[sheet]] <- genes
  }
}

if (length(signature_list) == 0) {
  stop("No signature in ", opt$annot, " matched a gene in the input data.")
}

# AddModuleScore bins all genes into `nbin` expression bins and samples
# `ctrl` control genes per bin *without replacement* -- on a small/targeted
# gene panel (HVG-filtered subsets, Xenium panels, ...) the default ctrl=100
# can exceed a bin's actual size and crash with "cannot take a sample larger
# than the population". Cap ctrl to the average bin size instead.
n_bin <- 24
ctrl_size <- max(1, min(100, floor(length(gene_universe) / n_bin)))
message(sprintf(
  "Gene universe: %d genes, %d bins -> ctrl=%d control genes/bin",
  length(gene_universe), n_bin, ctrl_size
))

# =====================
# Output folders
# =====================
output_dir <- opt$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
violin_dir <- file.path(output_dir, "violin")
dir.create(violin_dir, recursive = TRUE, showWarnings = FALSE)

# =====================
# AddModuleScore per signature
# =====================
score_names <- character(0)

for (sig_name in names(signature_list)) {
  message("Scoring signature: ", sig_name)

  seurat_obj <- Seurat::AddModuleScore(
    seurat_obj,
    features = list(signature_list[[sig_name]]),
    name = sig_name,
    nbin = n_bin,
    ctrl = ctrl_size
  )

  # AddModuleScore always suffixes the requested `name` with "1" for a
  # single feature list -- rename back to the plain signature name.
  generated_col <- paste0(sig_name, "1")
  colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == generated_col] <- sig_name
  score_names <- c(score_names, sig_name)

  global_median <- median(seurat_obj[[sig_name, drop = TRUE]])

  p <- VlnPlot(seurat_obj, features = sig_name, group.by = leiden_col, pt.size = 0) +
    geom_hline(
      yintercept = global_median, color = "red", linetype = "dashed"
    ) +
    ggtitle(sig_name) +
    ylab("Module score") +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))

  ggsave(
    file.path(violin_dir, paste0(sig_name, "_violin.png")),
    p, width = 8, height = 5, dpi = 300
  )
}

# =====================
# Heatmap (mean score per cluster)
# =====================
df_scores <- seurat_obj@meta.data[, c(leiden_col, score_names), drop = FALSE]

df_mean <- df_scores %>%
  group_by(.data[[leiden_col]]) %>%
  summarise(across(all_of(score_names), mean), .groups = "drop")

heat_mat <- as.matrix(df_mean[, score_names, drop = FALSE])
rownames(heat_mat) <- df_mean[[leiden_col]]

png(
  file.path(output_dir, "Heatmap_Signature_Validation.png"),
  width = 1600, height = 1000, res = 150
)
pheatmap::pheatmap(
  heat_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Validation of Cell Type Identity via Module Scores"
)
dev.off()

# =====================
# UMAP
# =====================
umap_reduction <- if ("UMAP" %in% names(seurat_obj@reductions)) "UMAP" else names(seurat_obj@reductions)[1]

p_umap <- DimPlot(seurat_obj, group.by = leiden_col, reduction = umap_reduction)
ggsave(
  file.path(output_dir, paste0(leiden_col, "_leiden_umap.png")),
  p_umap, width = 6, height = 5, dpi = 150
)

for (sig_name in score_names) {
  p <- FeaturePlot(seurat_obj, features = sig_name, reduction = umap_reduction) +
    scale_color_viridis_c()
  ggsave(
    file.path(output_dir, paste0(sig_name, "_umap.png")),
    p, width = 6, height = 5, dpi = 150
  )
}

# =====================
# Mapping automatico e Creazione Tabella di Report
# =====================
message("Creating annotation summary table with safety thresholds...")

final_annotations <- character(nrow(df_mean))
highest_scores <- numeric(nrow(df_mean))

for (i in seq_len(nrow(df_mean))) {
  row_scores <- unlist(df_mean[i, score_names])
  valid <- row_scores[row_scores >= opt$thr]

  if (length(valid) > 0) {
    valid <- sort(valid, decreasing = TRUE)
    final_annotations[i] <- paste(names(valid), collapse = ", ")
    highest_scores[i] <- valid[1]
  } else {
    final_annotations[i] <- "Unknown"
    highest_scores[i] <- max(row_scores)
  }
}

df_report <- as.data.frame(df_mean)
df_report$final_annotation <- final_annotations
df_report$highest_score <- highest_scores
colnames(df_report)[colnames(df_report) == leiden_col] <- "cluster"

openxlsx::write.xlsx(
  df_report,
  file.path(output_dir, "cluster_annotation_summary.xlsx")
)

message(
  "Summary table saved to ",
  file.path(output_dir, "cluster_annotation_summary.xlsx")
)

# =====================
# Applica formalmente la mappatura all'oggetto per le analisi a valle
# =====================
mapping_vec <- setNames(final_annotations, as.character(df_mean[[leiden_col]]))
seurat_obj$cell_type_pred <- unname(
  mapping_vec[as.character(seurat_obj[[leiden_col, drop = TRUE]])]
)

# =====================
# Save the annotated Seurat object -- workflow/rules/addmodulescore.smk runs
# rds2h5.R on this right after, as its own pure-R rhdf5 writer (unlike
# zellkonverter::writeH5AD, which has no R-only backend and falls back to a
# basilisk-managed python env) to land the final .h5ad in this method's
# output folder.
# =====================
rds_dir <- file.path(output_dir, "rds")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
rds_path <- file.path(rds_dir, "addmodulescore_annotated.rds")
saveRDS(seurat_obj, rds_path)
message("Saved annotated Seurat object -> ", rds_path)
