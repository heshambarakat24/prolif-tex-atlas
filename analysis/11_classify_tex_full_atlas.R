# ============================================================
# 11_classify_tex_full_atlas.R
# ============================================================
# Thesis Methods sections: 4.1.7 (score-based assignment) and
#                          4.1.8 (cross-compartment prevalence)
#
# Requires: CD8_Full_scVI in the R session, produced by
#           08_reconstruct_full_atlas_seurat.R. This script does not load it
#           from disk; run 08 first in the same session.
# Adds    : ExhaustionScore, ProlifScore, Exh_z, Prolif_z and TexClass to the
#           object's meta.data. Downstream scripts 12-14 depend on TexClass.
#
# Gene panels: the exhaustion and proliferation gene lists used here are the
# ones reported in Methods 4.1.7 and deposited in panels/. They are NOT the
# canonical marker list used for cluster DotPlots in script 10.
#
# Notes:
#  - Percentile thresholds (Q75 exhaustion, Q75/Q25 proliferation) are computed
#    on TUMOUR cells only and then applied to all compartments, so that the
#    cross-compartment comparison uses a common cut-off.
#  - The third category is named "Unlabelled" here; the thesis calls it "Other".
#  - Chunk 5 passes meta.data through dplyr group_by/mutate and assigns the
#    result back to the Seurat object. dplyr returns a tibble and drops row
#    names, which Seurat uses to track cell identity. Verify that
#    rownames(CD8_Full_scVI@meta.data) still holds cell barcodes before reuse.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
})


# ============================================================
# CHUNK 1 — Safety checks (ensure required metadata fields exist)
# Purpose: Fail early if the object lacks required columns for this workflow.
# ============================================================
stopifnot(exists("CD8_Full_scVI"))
stopifnot(inherits(CD8_Full_scVI, "Seurat"))
stopifnot(all(c("Dataset", "Tissue") %in% colnames(CD8_Full_scVI@meta.data)))

message("✅ Seurat object found: CD8_Full_scVI")
message("✅ Required metadata columns present: Dataset, Tissue")


# ============================================================
# CHUNK 2 — Define curated gene signatures
# Purpose: Provide gene sets exactly as stated in the Methods section.
# ============================================================
exhaustion_genes <- c("PDCD1","LAYN","HAVCR2","CTLA4","TOX","TIGIT","LAG3","CXCL13","ENTPD1")
prolif_genes     <- c("MKI67","TOP2A","PCNA","HMGB2","CDC20","UBE2C","BIRC5","AURKB","CCNB1")


# ============================================================
# CHUNK 3 — Optional sanity check of gene presence in the object
# Purpose: Warn if many genes are missing (common if var_names are Ensembl IDs).
# ============================================================
gene_universe <- rownames(CD8_Full_scVI[["RNA"]])
exh_missing   <- setdiff(exhaustion_genes, gene_universe)
pro_missing   <- setdiff(prolif_genes, gene_universe)

if (length(exh_missing) > 0) message("⚠️ Missing exhaustion genes: ", paste(exh_missing, collapse = ", "))
if (length(pro_missing) > 0) message("⚠️ Missing proliferation genes: ", paste(pro_missing, collapse = ", "))

if (length(exh_missing) >= 5 || length(pro_missing) >= 5) {
  message("⚠️ Many genes are missing. If your features are Ensembl IDs, you must map to gene symbols first.")
}


# ============================================================
# CHUNK 4 — Compute module scores (only if not already present)
# Purpose: Add exhaustion/proliferation module scores to meta.data.
# Notes: AddModuleScore appends a numeric suffix (typically '1').
# ============================================================
DefaultAssay(CD8_Full_scVI) <- "RNA"

if (!any(grepl("^ExhaustionScore", colnames(CD8_Full_scVI@meta.data)))) {
  CD8_Full_scVI <- AddModuleScore(
    CD8_Full_scVI,
    features = list(exhaustion_genes),
    name = "ExhaustionScore"
  )
  message("✅ Added ExhaustionScore module score.")
} else {
  message("ℹ️ ExhaustionScore already present. Skipping AddModuleScore.")
}

if (!any(grepl("^ProlifScore", colnames(CD8_Full_scVI@meta.data)))) {
  CD8_Full_scVI <- AddModuleScore(
    CD8_Full_scVI,
    features = list(prolif_genes),
    name = "ProlifScore"
  )
  message("✅ Added ProlifScore module score.")
} else {
  message("ℹ️ ProlifScore already present. Skipping AddModuleScore.")
}

# Detect the actual module score column names (first match is fine)
exh_col <- grep("^ExhaustionScore", colnames(CD8_Full_scVI@meta.data), value = TRUE)[1]
pro_col <- grep("^ProlifScore",     colnames(CD8_Full_scVI@meta.data), value = TRUE)[1]

message("Using exhaustion score column: ", exh_col)
message("Using proliferation score column: ", pro_col)


# ============================================================
# CHUNK 5 — Z-normalize scores within Dataset
# Purpose: Reduce dataset/platform differences in score magnitude as per Methods.
# Output columns: Exh_z, Prolif_z
# ============================================================
md <- CD8_Full_scVI@meta.data

md <- md %>%
  group_by(Dataset) %>%
  mutate(
    Exh_z    = as.numeric(scale(.data[[exh_col]])),
    Prolif_z = as.numeric(scale(.data[[pro_col]]))
  ) %>%
  ungroup()

CD8_Full_scVI@meta.data <- md
rm(md)

message("✅ Added z-normalized columns: Exh_z and Prolif_z")


# ============================================================
# CHUNK 6 — Compute tumour-based percentile thresholds
# Purpose: Define Q75/Q25 cutoffs using ONLY tumour cells (merged cohort).
# IMPORTANT: Adjust tumour label spelling if yours differs.
# ============================================================
md <- CD8_Full_scVI@meta.data

tum_mask <- md$Tissue %in% c("Tumour", "Tumor")
stopifnot(any(tum_mask))  # if this fails, your Tissue naming differs

exh_q75 <- quantile(md$Exh_z[tum_mask], probs = 0.75, na.rm = TRUE)
pro_q75 <- quantile(md$Prolif_z[tum_mask], probs = 0.75, na.rm = TRUE)
pro_q25 <- quantile(md$Prolif_z[tum_mask], probs = 0.25, na.rm = TRUE)

message("Tumour thresholds:")
message("  Exh_z Q75    = ", round(exh_q75, 3))
message("  Prolif_z Q75 = ", round(pro_q75, 3))
message("  Prolif_z Q25 = ", round(pro_q25, 3))


# ============================================================
# CHUNK 7 — Classify cells into Prolif-Tex / Non-Prolif-Tex / Unlabelled
# Purpose: Apply percentile-based rules at single-cell level (independent of clusters).
# Output column: TexClass
# ============================================================
md$TexClass <- "Unlabelled"
md$TexClass[ md$Exh_z >= exh_q75 & md$Prolif_z >= pro_q75 ] <- "Prolif-Tex"
md$TexClass[ md$Exh_z >= exh_q75 & md$Prolif_z <= pro_q25 ] <- "Non-Prolif-Tex"

md$TexClass <- factor(md$TexClass, levels = c("Prolif-Tex","Non-Prolif-Tex","Unlabelled"))

CD8_Full_scVI@meta.data <- md
rm(md)

message("✅ Added classification column: TexClass")
print(table(CD8_Full_scVI$TexClass, useNA = "ifany"))


# ============================================================
# CHUNK 8 — (Recommended) Stacked barplot: composition within labelled Tex only
# Purpose: For each Tissue compartment, show Prolif vs Non-Prolif among labelled Tex cells.
# Interpretation: "Of Tex cells we confidently classify, how do they split by Tissue?"
# ============================================================
df_comp_labelled <- CD8_Full_scVI@meta.data %>%
  filter(TexClass %in% c("Prolif-Tex","Non-Prolif-Tex")) %>%
  count(Tissue, TexClass, name = "n") %>%
  group_by(Tissue) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

p1 <- ggplot(df_comp_labelled, aes(x = Tissue, y = freq, fill = TexClass)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Distribution of Prolif-Tex vs Non-Prolif-Tex across Tissue compartments",
    subtitle = "Composition within labelled Tex cells (tumour-derived thresholds; dataset-wise z-normalised scores)",
    x = "Tissue compartment",
    y = "Fraction within labelled Tex",
    fill = NULL
  ) +
  theme_classic(base_size = 12)

print(p1)


# ============================================================
# CHUNK 9 — Optional stacked barplot: include Unlabelled (full composition)
# Purpose: Show full breakdown per Tissue including Unlabelled cells.
# Interpretation: "What fraction of all CD8 in each Tissue are classified vs unclassified?"
# ============================================================
df_comp_all <- CD8_Full_scVI@meta.data %>%
  count(Tissue, TexClass, name = "n") %>%
  group_by(Tissue) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

p2 <- ggplot(df_comp_all, aes(x = Tissue, y = freq, fill = TexClass)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Prolif-Tex / Non-Prolif-Tex classification across Tissue compartments",
    subtitle = "Includes Unlabelled cells (tumour-derived thresholds; dataset-wise z-normalised scores)",
    x = "Tissue compartment",
    y = "Fraction of all CD8 cells",
    fill = NULL
  ) +
  theme_classic(base_size = 12)

print(p2)


# ============================================================
# CHUNK 10 — Optional "prevalence" plot: fraction of all CD8 that are Prolif/Non-Prolif Tex
# Purpose: Estimate prevalence among all CD8 (not just composition of labelled Tex).
# Interpretation: "In each Tissue, what % of all CD8 are Prolif-Tex vs Non-Prolif-Tex?"
# ============================================================
df_prev <- CD8_Full_scVI@meta.data %>%
  mutate(TexClass2 = ifelse(TexClass %in% c("Prolif-Tex","Non-Prolif-Tex"), as.character(TexClass), "Other")) %>%
  count(Tissue, TexClass2, name = "n") %>%
  group_by(Tissue) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  filter(TexClass2 != "Other")

p3 <- ggplot(df_prev, aes(x = Tissue, y = freq, fill = TexClass2)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Prevalence of Prolif-Tex vs Non-Prolif-Tex across Tissue compartments",
    subtitle = "Fraction among all CD8 cells (tumour-derived thresholds; dataset-wise z-normalised scores)",
    x = "Tissue compartment",
    y = "Prevalence among all CD8",
    fill = NULL
  ) +
  theme_classic(base_size = 12)

print(p3)


# ============================================================
# CHUNK 11 — Save results back into the Seurat object
# Purpose: Persist TexClass + score columns for downstream analyses.
# ============================================================
saveRDS(CD8_Full_scVI, file = "CD8_Full_scVI_with_TexClass.rds")
message("✅ Saved: CD8_Full_scVI_with_TexClass.rds")