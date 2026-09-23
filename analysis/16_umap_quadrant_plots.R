# ============================================================
# 16_umap_quadrant_plots.R
# ============================================================
# Thesis Methods section: 4.1.7 (score-based assignment)
# Produces: Figure 6B (exhaustion-proliferation quadrant plot)
#           Figure 6C (Prolif-Tex / NonProlif-Tex on the tumour-only UMAP)
#
# Requires: CD8_TumourOnly_scVI in the R session, produced by
#           09_reconstruct_tumour_only_seurat.R. Not loaded from disk here.
#
# Gene panels: the nine-gene exhaustion and nine-gene proliferation panels of
# Methods 4.1.7, as deposited in panels/. These are the same panels used by
# scripts 11, 15 and 17.
#
# Note: AddModuleScore selects control genes at random and no seed is set in
# this script, so a re-run shifts a small number of borderline cells between
# configurations. Directions and thresholds are unaffected. See README.
#
# Input data are not included in this repository. See README.md.
# ============================================================

# ============================================================
# TWO MAIN THESIS FIGURES:
# 1) UMAP: Prolif-Tex vs Non-Prolif-Tex on tumour-only scVI atlas
# 2) Exhaustion vs Proliferation quadrant plot
#
# Purpose:
# - Score tumour CD8+ T cells for exhaustion and proliferation
# - Z-normalise scores within Dataset
# - Define:
#     Prolif-Tex     = Exhaustion >= Q75 AND Proliferation >= Q75
#     Non-prolif-Tex = Exhaustion >= Q75 AND Proliferation <= Q25
# - Generate:
#     Figure 1: UMAP overlay of Prolif-Tex and Non-prolif-Tex
#     Figure 2: Exhaustion vs Proliferation quadrant plot
#
# Assumptions:
# - Object: CD8_TumourOnly_scVI
# - Existing UMAP reduction is already present in the Seurat object
# - RNA assay is present
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 0) Input object
# ------------------------------------------------------------

obj <- CD8_TumourOnly_scVI
DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------
# 1) Check object structure
# ------------------------------------------------------------

cat("\nAvailable assays:\n")
print(Assays(obj))

cat("\nAvailable reductions:\n")
print(Reductions(obj))

if (!"umap" %in% Reductions(obj)) {
  stop("UMAP reduction not found in object. Please compute/import UMAP first.")
}

if (!"Dataset" %in% colnames(obj@meta.data)) {
  stop("Dataset column not found in metadata. It is required for within-dataset z-normalisation.")
}

# ------------------------------------------------------------
# 2) Curated signatures
# ------------------------------------------------------------

exhaustion_genes <- c(
  "PDCD1", "LAYN", "HAVCR2", "CTLA4", "TOX",
  "TIGIT", "LAG3", "CXCL13", "ENTPD1"
)

proliferation_genes <- c(
  "MKI67", "TOP2A", "PCNA", "HMGB2", "CDC20",
  "UBE2C", "BIRC5", "AURKB", "CCNB1"
)

exhaustion_genes_present <- exhaustion_genes[exhaustion_genes %in% rownames(obj)]
proliferation_genes_present <- proliferation_genes[proliferation_genes %in% rownames(obj)]

cat("\nExhaustion genes present:\n")
print(exhaustion_genes_present)

cat("\nProliferation genes present:\n")
print(proliferation_genes_present)

if (length(exhaustion_genes_present) < 3) {
  stop("Too few exhaustion genes found in object.")
}

if (length(proliferation_genes_present) < 3) {
  stop("Too few proliferation genes found in object.")
}

# ------------------------------------------------------------
# 3) Remove previous module-score columns if the block was run before
# ------------------------------------------------------------

old_score_cols <- grep(
  pattern = "^(ExhaustionScore|ProliferationScore)[0-9]+$",
  x = colnames(obj@meta.data),
  value = TRUE
)

if (length(old_score_cols) > 0) {
  cat("\nRemoving old score columns:\n")
  print(old_score_cols)
  obj@meta.data[, old_score_cols] <- NULL
}

# ------------------------------------------------------------
# 4) Add module scores
# ------------------------------------------------------------

set.seed(123)

obj <- AddModuleScore(
  object = obj,
  features = list(exhaustion_genes_present),
  assay = "RNA",
  name = "ExhaustionScore"
)

obj <- AddModuleScore(
  object = obj,
  features = list(proliferation_genes_present),
  assay = "RNA",
  name = "ProliferationScore"
)

# AddModuleScore creates:
# ExhaustionScore1
# ProliferationScore1

# ------------------------------------------------------------
# 5) Extract metadata and z-normalise scores within Dataset
# ------------------------------------------------------------

meta <- obj@meta.data %>%
  mutate(
    cell_id = rownames(obj@meta.data),
    ExhaustionScore_raw = ExhaustionScore1,
    ProliferationScore_raw = ProliferationScore1,
    Dataset = as.character(Dataset)
  ) %>%
  group_by(Dataset) %>%
  mutate(
    ExhaustionScore_z = as.numeric(scale(ExhaustionScore_raw)),
    ProliferationScore_z = as.numeric(scale(ProliferationScore_raw))
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 6) Define pooled thresholds across the tumour-only atlas
# ------------------------------------------------------------

exh_q75  <- quantile(meta$ExhaustionScore_z, probs = 0.75, na.rm = TRUE)
prol_q75 <- quantile(meta$ProliferationScore_z, probs = 0.75, na.rm = TRUE)
prol_q25 <- quantile(meta$ProliferationScore_z, probs = 0.25, na.rm = TRUE)

cat("\nThresholds used:\n")
cat("Exhaustion Q75 =", round(exh_q75, 4), "\n")
cat("Proliferation Q75 =", round(prol_q75, 4), "\n")
cat("Proliferation Q25 =", round(prol_q25, 4), "\n")

# ------------------------------------------------------------
# 7) Classify cells
# ------------------------------------------------------------

meta <- meta %>%
  mutate(
    TexClass = case_when(
      ExhaustionScore_z >= exh_q75 & ProliferationScore_z >= prol_q75 ~ "Prolif-Tex",
      ExhaustionScore_z >= exh_q75 & ProliferationScore_z <= prol_q25 ~ "Non-prolif-Tex",
      TRUE ~ "Other"
    ),
    TexClass = factor(
      TexClass,
      levels = c("Other", "Non-prolif-Tex", "Prolif-Tex")
    )
  )

cat("\nCell counts by class:\n")
print(table(meta$TexClass, useNA = "ifany"))

cat("\nCell counts by class and Dataset:\n")
print(table(meta$Dataset, meta$TexClass))

# ------------------------------------------------------------
# 8) Add scores and classification back to working object only
#    No new RDS object is saved.
# ------------------------------------------------------------

obj$ExhaustionScore_z <- meta$ExhaustionScore_z[match(colnames(obj), meta$cell_id)]
obj$ProliferationScore_z <- meta$ProliferationScore_z[match(colnames(obj), meta$cell_id)]
obj$TexClass <- meta$TexClass[match(colnames(obj), meta$cell_id)]

# ============================================================
# FIGURE 1: UMAP Prolif-Tex vs Non-prolif-Tex
# ============================================================

# ------------------------------------------------------------
# 9) Extract UMAP coordinates
# ------------------------------------------------------------

umap_df <- Embeddings(obj, reduction = "umap") %>%
  as.data.frame()

colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell_id <- rownames(umap_df)

umap_df <- umap_df %>%
  left_join(
    meta %>% select(cell_id, TexClass),
    by = "cell_id"
  )

# ------------------------------------------------------------
# 10) Split into background and highlighted cells
# ------------------------------------------------------------

bg_df <- umap_df %>%
  filter(TexClass == "Other")

fg_df <- umap_df %>%
  filter(TexClass %in% c("Prolif-Tex", "Non-prolif-Tex")) %>%
  mutate(
    TexClass = factor(
      TexClass,
      levels = c("Prolif-Tex", "Non-prolif-Tex")
    )
  )

# ------------------------------------------------------------
# 11) Plot UMAP
# ------------------------------------------------------------

p_tex_umap <- ggplot() +
  geom_point(
    data = bg_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey85",
    size = 0.2,
    alpha = 0.5
  ) +
  geom_point(
    data = fg_df,
    aes(x = UMAP_1, y = UMAP_2, color = TexClass),
    size = 0.28,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "Prolif-Tex" = "#D95F02",
      "Non-prolif-Tex" = "#1F78B4"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  ) +
  labs(
    title = "UMAP: Prolif-Tex vs Non-prolif-Tex",
    subtitle = "Tumour-only scVI atlas; other CD8 T cells shown in grey",
    x = "UMAP1",
    y = "UMAP2",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12)
  )

print(p_tex_umap)

# ------------------------------------------------------------
# 12) Save UMAP figure
# ------------------------------------------------------------

ggsave(
  filename = "FIGURE_UMAP_ProlifTex_vs_NonProlifTex_tumour_only.png",
  plot = p_tex_umap,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "FIGURE_UMAP_ProlifTex_vs_NonProlifTex_tumour_only.pdf",
  plot = p_tex_umap,
  width = 8,
  height = 6
)

# ============================================================
# FIGURE 2: Exhaustion vs Proliferation quadrant plot
# ============================================================

# ------------------------------------------------------------
# 13) Prepare scatter plot dataframe
# ------------------------------------------------------------

scatter_df <- meta %>%
  select(
    cell_id,
    Dataset,
    ExhaustionScore_z,
    ProliferationScore_z,
    TexClass
  ) %>%
  mutate(
    plot_order = case_when(
      TexClass == "Other" ~ 1,
      TexClass == "Non-prolif-Tex" ~ 2,
      TexClass == "Prolif-Tex" ~ 3
    ),
    TexClass = factor(
      TexClass,
      levels = c("Other", "Non-prolif-Tex", "Prolif-Tex")
    )
  ) %>%
  arrange(plot_order)

# ------------------------------------------------------------
# 14) Plot exhaustion vs proliferation
# ------------------------------------------------------------

p_exh_prol <- ggplot(
  scatter_df,
  aes(x = ExhaustionScore_z, y = ProliferationScore_z)
) +
  geom_point(
    aes(color = TexClass),
    size = 0.18,
    alpha = 0.45
  ) +
  geom_vline(
    xintercept = exh_q75,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_hline(
    yintercept = prol_q75,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_hline(
    yintercept = prol_q25,
    linetype = "dotted",
    linewidth = 0.45,
    color = "grey35"
  ) +
  scale_color_manual(
    values = c(
      "Other" = "grey80",
      "Non-prolif-Tex" = "#1F78B4",
      "Prolif-Tex" = "#D95F02"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  ) +
  labs(
    title = "Exhaustion vs Proliferation",
    subtitle = "Quadrants define Prolif-Tex (Q75/Q75) and Non-prolif-Tex (Q75/Q25)",
    x = "Exhaustion score (z, within dataset)",
    y = "Proliferation score (z, within dataset)",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12)
  )

print(p_exh_prol)

# ------------------------------------------------------------
# 15) Save exhaustion vs proliferation figure
# ------------------------------------------------------------

ggsave(
  filename = "FIGURE_Exhaustion_vs_Proliferation_quadrants.png",
  plot = p_exh_prol,
  width = 8,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = "FIGURE_Exhaustion_vs_Proliferation_quadrants.pdf",
  plot = p_exh_prol,
  width = 8,
  height = 5.5
)

cat("\nDone. Generated and saved:\n")
cat("1) FIGURE_UMAP_ProlifTex_vs_NonProlifTex_tumour_only.png/pdf\n")
cat("2) FIGURE_Exhaustion_vs_Proliferation_quadrants.png/pdf\n")