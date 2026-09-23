# ============================================================
# 08_reconstruct_full_atlas_seurat.R
# ============================================================
# Thesis Methods section: 4.1.5 (import of the scVI export into Seurat)
#
# Reads : the full-atlas Matrix Market export written by 07_merge_scvi_export.ipynb
#         (counts.mtx, features.tsv, barcodes.tsv, metadata.csv,
#          scvi_embedding.csv, umap_coords.csv)
# Writes: CD8_Full_scVI.rds
#
# Notes:
#  - The scVI latent space and the UMAP coordinates are imported from Python.
#    Neither is recomputed here.
#  - The clustering performed below (resolution 0.25) is NOT reported in the
#    thesis. Methods 4.1.6 describes clustering of the tumour-only atlas only.
#    It is retained because it was part of the original script.
#  - Input data are not included in this repository. Set the DATA_ROOT
#    environment variable, or edit the fallback below. See README.md.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# ============================================================
# 0) Set working directory to the FULL atlas export folder
#    (contains counts.mtx, features.tsv, barcodes.tsv, metadata.csv,
#     scvi_embedding.csv, umap_coords.csv)
# ============================================================
DATA_ROOT  <- Sys.getenv("DATA_ROOT", unset = "data")
EXPORT_DIR <- file.path(DATA_ROOT,
                        "scVI/13 Merged 2-3-5-10-11-12",
                        "Reconstruction of CD8_Full")
setwd(EXPORT_DIR)

# -----------------------------
# 1) Load counts + create object
# -----------------------------
m <- ReadMtx(
  mtx = "counts.mtx",
  features = "features.tsv",
  cells = "barcodes.tsv",
  feature.column = 2   # use gene_name (2nd col)
)

CD8_Full_scVI <- CreateSeuratObject(
  counts = m,
  assay  = "RNA",
  project = "CD8_Full_scVI",
  min.cells = 0,
  min.features = 0
)

# -----------------------------
# 2) Attach metadata (aligned to cell names)
# -----------------------------
meta <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)

# align + keep same order as Seurat object
meta <- meta[colnames(CD8_Full_scVI), , drop = FALSE]
CD8_Full_scVI <- AddMetaData(CD8_Full_scVI, meta)

# -----------------------------
# 3) QC helper + normalization
# -----------------------------
DefaultAssay(CD8_Full_scVI) <- "RNA"

# optional: percent.mt if you have pct_counts_mt in metadata
if ("pct_counts_mt" %in% colnames(CD8_Full_scVI@meta.data)) {
  CD8_Full_scVI$percent.mt <- as.numeric(CD8_Full_scVI$pct_counts_mt)
}

CD8_Full_scVI <- NormalizeData(
  CD8_Full_scVI,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

# sanity checks
print(dim(CD8_Full_scVI))
if ("Tissue" %in% colnames(CD8_Full_scVI@meta.data)) {
  print(table(CD8_Full_scVI$Tissue))
}
if ("compartment" %in% colnames(CD8_Full_scVI@meta.data)) {
  print(table(CD8_Full_scVI$compartment))
}

# -----------------------------
# 4) scVI latent -> Seurat reduction
# -----------------------------
scvi <- read.csv("scvi_embedding.csv", row.names = 1, check.names = FALSE)
scvi <- scvi[colnames(CD8_Full_scVI), , drop = FALSE]

scvi <- as.matrix(sapply(scvi, as.numeric))
rownames(scvi) <- colnames(CD8_Full_scVI)
colnames(scvi) <- paste0("scVI_", seq_len(ncol(scvi)))

CD8_Full_scVI[["scVI"]] <- CreateDimReducObject(
  embeddings = scvi,
  key = "scVI_",
  assay = DefaultAssay(CD8_Full_scVI)
)

# -----------------------------
# 5) Import the *exact same* UMAP coordinates from Python
#    (do NOT recompute UMAP in Seurat)
# -----------------------------
umap <- read.csv("umap_coords.csv", row.names = 1, check.names = FALSE)
umap <- umap[colnames(CD8_Full_scVI), , drop = FALSE]

umap <- as.matrix(sapply(umap, as.numeric))
rownames(umap) <- colnames(CD8_Full_scVI)
colnames(umap) <- c("UMAP_1", "UMAP_2")

CD8_Full_scVI[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(CD8_Full_scVI)
)

# -----------------------------
# 6) Neighbors + clustering in scVI space
# -----------------------------
nd <- ncol(Embeddings(CD8_Full_scVI, "scVI"))

set.seed(1)
CD8_Full_scVI <- FindNeighbors(
  CD8_Full_scVI,
  reduction = "scVI",
  dims = 1:nd,
  k.param = 40,
  graph.name = "scVI_snn"
)

# Pick your resolution (start with 0.25 to mirror what you used)
CD8_Full_scVI <- FindClusters(
  CD8_Full_scVI,
  graph.name = "scVI_snn",
  resolution = 0.25,
  algorithm = 1  # 1=Louvain; 4=Leiden if available
)

Idents(CD8_Full_scVI) <- "seurat_clusters"

# -----------------------------
# 7) Plots: same UMAP layout as Python
# -----------------------------
# Clusters
DimPlot(
  CD8_Full_scVI,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

# Cross-compartment view (use whichever field exists: Tissue or compartment)
if ("Tissue" %in% colnames(CD8_Full_scVI@meta.data)) {
  print(DimPlot(CD8_Full_scVI, reduction = "umap", group.by = "Tissue"))
}
if ("compartment" %in% colnames(CD8_Full_scVI@meta.data)) {
  print(DimPlot(CD8_Full_scVI, reduction = "umap", group.by = "compartment"))
}

# Optional: within-dataset or patient checks (only if fields exist)
if ("Dataset" %in% colnames(CD8_Full_scVI@meta.data)) {
  print(DimPlot(CD8_Full_scVI, reduction = "umap", group.by = "Dataset"))
}
if ("PatientID" %in% colnames(CD8_Full_scVI@meta.data)) {
  # can be huge; keep for spot checks on subsets
  message("PatientID exists (not plotting by default; too many levels).")
}

# Verify identical cell ordering in the reduction
stopifnot(identical(
  rownames(Embeddings(CD8_Full_scVI, "umap")),
  colnames(CD8_Full_scVI)
))

# -----------------------------
# 8) Save Seurat object
# -----------------------------
saveRDS(CD8_Full_scVI, file = "CD8_Full_scVI.rds")
message("✅ Saved: CD8_Full_scVI.rds")