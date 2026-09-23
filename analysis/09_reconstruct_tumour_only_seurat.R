# ============================================================
# 09_reconstruct_tumour_only_seurat.R
# ============================================================
# Thesis Methods section: 4.1.5 (import of the scVI export into Seurat)
#
# Reads : the tumour-only Matrix Market export written by the final block of
#         07_merge_scvi_export.ipynb (counts.mtx, features.tsv, barcodes.tsv,
#         metadata.csv, scvi_embedding.csv, umap_coords.csv)
# Writes: CD8_TumourOnly_scVI.rds
#
# Notes:
#  - The scVI latent space and the UMAP coordinates are imported from Python.
#    Neither is recomputed here.
#  - This script builds the "scVI_snn" neighbour graph (k.param = 40, 30 scVI
#    dimensions) that 10_cluster_annotate_tumour_only.R depends on. It must be
#    run before that script.
#  - The FindClusters call below uses resolution 0.50. That solution is NOT
#    the one reported in the thesis. It was superseded by resolution 0.27,
#    applied in script 10. The line is retained because it was part of the
#    original script and because the FindNeighbors step above it is required.
#  - Input data are not included in this repository. Set the DATA_ROOT
#    environment variable, or edit the fallback below. See README.md.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

DATA_ROOT  <- Sys.getenv("DATA_ROOT", unset = "data")
EXPORT_DIR <- file.path(DATA_ROOT,
                        "scVI/13 Merged 2-3-5-10-11-12",
                        "Reconstruction of CD8_TumourOnly",
                        "Seurat Object Reconstruction files")
setwd(EXPORT_DIR)


# ============================================================
# Load tumour-only atlas export (counts + metadata + scVI + UMAP)
# and reproduce the *same* scVI-UMAP layout from Python in Seurat.
# ============================================================

# -----------------------------
# 1) Load counts + create object
# -----------------------------
m <- ReadMtx(
  mtx = "counts.mtx",
  features = "features.tsv",
  cells = "barcodes.tsv",
  feature.column = 2   # use gene_name (2nd col)
)

CD8_TumourOnly_scVI <- CreateSeuratObject(
  counts = m,
  assay  = "RNA",
  project = "CD8_TumourOnly",
  min.cells = 0,
  min.features = 0
)

# -----------------------------
# 2) Attach metadata (aligned to cell names)
# -----------------------------
meta <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)

# align + keep same order as Seurat object
meta <- meta[colnames(CD8_TumourOnly_scVI), , drop = FALSE]
CD8_TumourOnly_scVI <- AddMetaData(CD8_TumourOnly_scVI, meta)

# -----------------------------
# 3) QC helper + normalization
# -----------------------------
DefaultAssay(CD8_TumourOnly_scVI) <- "RNA"

# optional: percent.mt if you have pct_counts_mt in metadata
if ("pct_counts_mt" %in% colnames(CD8_TumourOnly_scVI@meta.data)) {
  CD8_TumourOnly_scVI$percent.mt <- as.numeric(CD8_TumourOnly_scVI$pct_counts_mt)
}

CD8_TumourOnly_scVI <- NormalizeData(
  CD8_TumourOnly_scVI,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

# (Optional) sanity checks
print(dim(CD8_TumourOnly_scVI))
if ("Tissue" %in% colnames(CD8_TumourOnly_scVI@meta.data)) {
  print(table(CD8_TumourOnly_scVI$Tissue))
}

# -----------------------------
# 4) scVI latent -> Seurat reduction
# -----------------------------
scvi <- read.csv("scvi_embedding.csv", row.names = 1, check.names = FALSE)
scvi <- scvi[colnames(CD8_TumourOnly_scVI), , drop = FALSE]

scvi <- as.matrix(sapply(scvi, as.numeric))
rownames(scvi) <- colnames(CD8_TumourOnly_scVI)
colnames(scvi) <- paste0("scVI_", seq_len(ncol(scvi)))

CD8_TumourOnly_scVI[["scVI"]] <- CreateDimReducObject(
  embeddings = scvi,
  key = "scVI_",
  assay = DefaultAssay(CD8_TumourOnly_scVI)
)

# -----------------------------
# 5) Import the *exact same* UMAP coordinates from Python
#    (do NOT recompute UMAP in Seurat)
# -----------------------------
umap <- read.csv("umap_coords.csv", row.names = 1, check.names = FALSE)
umap <- umap[colnames(CD8_TumourOnly_scVI), , drop = FALSE]

umap <- as.matrix(sapply(umap, as.numeric))
rownames(umap) <- colnames(CD8_TumourOnly_scVI)
colnames(umap) <- c("UMAP_1", "UMAP_2")

CD8_TumourOnly_scVI[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(CD8_TumourOnly_scVI)
)

# -----------------------------
# 6) Neighbors + clustering in scVI space
# -----------------------------
nd <- ncol(Embeddings(CD8_TumourOnly_scVI, "scVI"))

set.seed(1)
CD8_TumourOnly_scVI <- FindNeighbors(
  CD8_TumourOnly_scVI,
  reduction = "scVI",
  dims = 1:nd,
  k.param = 40,
  graph.name = "scVI_snn"
)

CD8_TumourOnly_scVI <- FindClusters(
  CD8_TumourOnly_scVI,
  graph.name = "scVI_snn",
  resolution = 0.50,
  algorithm = 1  # 1=Louvain; 4=Leiden if available
)

Idents(CD8_TumourOnly_scVI) <- "seurat_clusters"

# -----------------------------
# 7) Plot: same UMAP layout as Python
# -----------------------------
DimPlot(
  CD8_TumourOnly_scVI,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

# (Optional) quickly verify identical cell ordering in the reduction
stopifnot(identical(
  rownames(Embeddings(CD8_TumourOnly_scVI, "umap")),
  colnames(CD8_TumourOnly_scVI)
))

# -----------------------------
# 8) Save Seurat object
# -----------------------------
saveRDS(CD8_TumourOnly_scVI, file = "CD8_TumourOnly_scVI.rds")