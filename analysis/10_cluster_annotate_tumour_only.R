# ============================================================
# 10_cluster_annotate_tumour_only.R
# ============================================================
# Thesis Methods section: 4.1.6 (clustering and marker-based annotation)
#
# Reads : CD8_TumourOnly_scVI.rds, written by
#         09_reconstruct_tumour_only_seurat.R. That object already carries the "scVI_snn" neighbour graph
#         (FindNeighbors, reduction = "scVI", dims = 1:30, k.param = 40).
# Writes: CD8_TumourOnly_scVI.rds (re-saved with the resolution-0.27 clusters)
#         tumour_only_scvi_outputs/markers_allclusters_wilcox.csv
#         tumour_only_scvi_outputs/markers_top10_percluster.csv
#         tumour_only_scvi_outputs/DotPlot_canonical_markers.png
#
# PROVENANCE: this script was not saved as a .R file during the original
# analysis. It is reconstructed from the R command history (.Rhistory) of the
# session, so the commands below are the ones that were actually run, with
# repeated and abandoned attempts removed. It has not been re-run here.
#
# Known problems in the original session, left uncorrected (see README):
#   - The clusters were first computed at resolution 0.25 (7 clusters) and
#     FindAllMarkers was run on that solution before re-clustering at 0.27.
#     Only the final resolution-0.27 marker tables are reported in the thesis;
#     the earlier ones were overwritten in place.
#   - A metadata column named "cluster_res025" was assigned the resolution-0.27
#     identities. The name is misleading; the contents are the 0.27 solution.
#   - One DotPlot call referenced "cluster_res027" before that column existed,
#     and one print() call referenced an object name that was never defined.
#     Both errored in the original session and were re-run afterwards.
#
# Input data are not included in this repository. Set the DATA_ROOT
# environment variable, or edit the fallback below. See README.md.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "data")
WORK_DIR  <- file.path(DATA_ROOT, "scVI/13 Merged 2-3-5-10-11-12")
setwd(WORK_DIR)

CD8_TumourOnly_scVI <- readRDS("CD8_TumourOnly_scVI.rds")

# -----------------------------
# 1) Clustering at resolution 0.27
#    The "scVI_snn" graph already exists on the object.
# -----------------------------
CD8_TumourOnly_scVI <- FindClusters(
  CD8_TumourOnly_scVI,
  graph.name = "scVI_snn",
  resolution = 0.27,
  algorithm  = 1   # 1 = Louvain
)

Idents(CD8_TumourOnly_scVI) <- "seurat_clusters"

# -----------------------------
# 2) Sanity checks
# -----------------------------
obj <- CD8_TumourOnly_scVI
DefaultAssay(obj) <- "RNA"

cat("Cells:", ncol(obj), " Genes:", nrow(obj), "\n")
cat("Reductions:", paste(names(obj@reductions), collapse = ", "), "\n")

if ("scVI_snn_res.0.27" %in% colnames(obj@meta.data)) {
  Idents(obj) <- obj$scVI_snn_res.0.27
  cat("Using Idents = scVI_snn_res.0.27\n")
} else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  Idents(obj) <- obj$seurat_clusters
  cat("Using Idents = seurat_clusters\n")
} else {
  stop("No clustering column found. Expected scVI_snn_res.0.27 or seurat_clusters.")
}

Idents(obj) <- factor(Idents(obj), levels = sort(unique(as.character(Idents(obj)))))
obj$cluster_res027 <- Idents(obj)

cat("Clusters:", length(levels(Idents(obj))), "\n")
print(table(Idents(obj)))

# -----------------------------
# 3) Cluster markers
# -----------------------------
marker_args <- list(
  object          = obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

cat("Running FindAllMarkers...\n")
markers <- do.call(FindAllMarkers, marker_args)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

out_dir <- "tumour_only_scvi_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(markers, file.path(out_dir, "markers_allclusters_wilcox.csv"), row.names = FALSE)
write.csv(top10,   file.path(out_dir, "markers_top10_percluster.csv"),   row.names = FALSE)
cat("Markers saved to:", out_dir, "\n")

# -----------------------------
# 4) Canonical marker DotPlot (annotation aid)
#    NOTE: this list is for visual inspection only. It is NOT the gene panel
#    used for Prolif-Tex / NonProlif-Tex assignment; see Methods 4.1.7 and
#    the panels/ folder.
# -----------------------------
canonical <- c(
  # Exhaustion / Tex
  "PDCD1", "LAG3", "HAVCR2", "TOX", "ENTPD1", "CTLA4", "TIGIT",
  # Proliferation
  "MKI67", "TOP2A", "HMGB2", "TYMS", "STMN1",
  # Effector cytotoxicity
  "NKG7", "GZMB", "PRF1", "IFNG", "GNLY",
  # Memory-like
  "IL7R", "CCR7", "LTB", "TCF7",
  # Activation
  "CD69", "HLA-DRA"
)

canonical_present <- canonical[canonical %in% rownames(obj)]
cat("Canonical genes present:", length(canonical_present), "/", length(canonical), "\n")

p_dot <- DotPlot(obj, features = canonical_present, group.by = "cluster_res027") +
  RotatedAxis() +
  ggtitle("Canonical markers across clusters (res=0.27)")

ggsave(file.path(out_dir, "DotPlot_canonical_markers.png"),
       p_dot, width = 14, height = 6, dpi = 200)

# -----------------------------
# 5) Cluster plot on the Python UMAP
# -----------------------------
DimPlot(
  CD8_TumourOnly_scVI,
  reduction = "umap",
  group.by  = "seurat_clusters",
  label     = TRUE
)

stopifnot(identical(
  rownames(Embeddings(CD8_TumourOnly_scVI, "umap")),
  colnames(CD8_TumourOnly_scVI)
))

# -----------------------------
# 6) Save
# -----------------------------
saveRDS(CD8_TumourOnly_scVI, file = "CD8_TumourOnly_scVI.rds")
