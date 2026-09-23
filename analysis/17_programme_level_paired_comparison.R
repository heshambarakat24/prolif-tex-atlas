## ============================================================
## 17_programme_level_paired_comparison.R
## ============================================================
## Thesis Methods section: 4.1.9   |   Produces: Figure 8
##
## Requires: CD8_TumourOnly_scVI in the R session, produced by
##           09_reconstruct_tumour_only_seurat.R. Not loaded from disk here.
## Writes  : Programme_Module_Scoring_ProlifTex_vs_NonProlifTex_v2_panels417/
##           Programme_level_statistics.csv, Programme_gene_coverage.csv,
##           Patient_level_programme_scores_*.csv, FunctionalConfig_*.csv,
##           Figure_Definition_axes.*, Figure_Programme_level_differences.*,
##           Figure_Combined_definition_and_programme_level.*
##
## What it does:
##  - Scores fifteen gene programmes with AddModuleScore on the tumour-only
##    object, z-normalised within Dataset. Programme gene lists are those
##    reported in Methods 4.1.9 and deposited in panels/.
##  - Assigns Prolif-Tex and NonProlif-Tex using the nine-gene exhaustion and
##    nine-gene proliferation panels of Methods 4.1.7. The twelve-gene
##    Exhaustion programme is compared but does not define the configurations.
##  - Builds a patient-level paired table, keeping patients with at least
##    10 cells in both configurations.
##  - Reports median paired differences with bootstrap 95% confidence
##    intervals (2000 iterations), paired Wilcoxon signed-rank tests with
##    exact = FALSE, Benjamini-Hochberg adjustment, and paired rank-biserial
##    correlation as the effect size.
##
## Input data are not included in this repository. See README.md.
## ============================================================

## ------------------------------------------------------------
## Programme-level comparison: Prolif-Tex vs NonProlif-Tex
## ------------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

set.seed(123)

## ------------------------------------------------------------
## 0. Output folder
## ------------------------------------------------------------

out_dir <- "Programme_Module_Scoring_ProlifTex_vs_NonProlifTex"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------
## 1. Object setup
## ------------------------------------------------------------

obj <- CD8_TumourOnly_scVI
DefaultAssay(obj) <- "RNA"

meta_cols <- colnames(obj@meta.data)

patient_col <- if ("PatientID_std" %in% meta_cols) {
  "PatientID_std"
} else if ("PatientID" %in% meta_cols) {
  "PatientID"
} else {
  stop("No PatientID or PatientID_std column found.")
}

if (!"Dataset" %in% meta_cols) {
  obj$Dataset <- "Dataset_1"
}

dataset_col <- "Dataset"

## ------------------------------------------------------------
## 2. Define gene programmes
## ------------------------------------------------------------

programmes_raw <- list(
  
  ## Definition axes
  Cycling = c(
    "MKI67", "TOP2A", "PCNA", "HMGB2", "CDC20",
    "UBE2C", "BIRC5", "AURKB", "CCNB1"
  ),
  
  Exhaustion = c(
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX",
    "LAYN", "ENTPD1", "CXCL13", "EOMES", "NR4A1", "NR4A2"
  ),
  
  ## Interpretable non-defining programmes
  Glycolysis = c(
    "HK2", "GAPDH", "ALDOA", "ENO1", "LDHA",
    "SLC2A3", "PFKP", "PGK1", "PKM", "TPI1"
  ),
  
  OXPHOS = c(
    "ATP5F1A", "ATP5F1B", "COX5A", "COX5B", "COX7A2",
    "UQCRC1", "UQCRQ", "NDUFA8", "NDUFA9", "NDUFB8", "SDHB"
  ),
  
  p53_checkpoint_stress = c(
    "TP53", "CDKN1A", "GADD45A", "GADD45B", "DDIT3",
    "BBC3", "PMAIP1", "BAX", "RRM2B", "CHEK1",
    "WEE1", "CLSPN", "RAD17"
  ),
  
  mTOR_MYC = c(
    "MTOR", "AKT1", "EIF4EBP1", "RPS6KB1", "TSC2",
    "RPTOR", "RHEB", "RPS6", "EIF4E", "MYC", "NPM1"
  ),
  
  Cytotoxicity = c(
    "GZMB", "GNLY", "PRF1", "NKG7", "FGFBP2",
    "GZMH", "GZMA", "KLRD1", "KLRK1", "SPON2"
  ),
  
  TypeI_IFN_ISG = c(
    "IFIT1", "IFIT3", "ISG15", "MX1", "OAS1",
    "IFI6", "IFI27", "IFIT2", "IRF7", "RSAD2"
  ),
  
  TRM = c(
    "ITGAE", "CXCR6", "CD69", "ZNF683", "RUNX3",
    "CXCR3", "CXCR4", "ITGA1", "ITGAL"
  ),
  
  Apoptosis_death_response = c(
    "BCL2", "BCL2L1", "BAX", "BCL2A1", "CASP3",
    "FAS", "FASLG", "BID", "BAD", "BAK1", "BBC3", "PMAIP1"
  ),
  
  Checkpoint_receptors = c(
    "PDCD1", "HAVCR2", "LAG3", "TIGIT",
    "ENTPD1", "LAYN", "CTLA4", "LGALS9"
  ),
  
  TCF7_memory = c(
    "TCF7", "CCR7", "LEF1", "IL7R", "BCL2",
    "SELL", "CXCR3", "KLF2"
  ),
  
  TCR_signalling = c(
    "CD3D", "CD3E", "CD3G", "TRAC", "LCK",
    "ZAP70", "LAT", "PTPRC", "VAV1", "ITK"
  ),
  
  NFkB = c(
    "NFKB1", "RELA", "TNFAIP3", "BIRC3",
    "IKBKB", "NFKBIA", "REL", "TRAF1", "TRAF6"
  ),
  
  Costimulation = c(
    "CD28", "ICOS", "CD27", "TNFRSF9",
    "CD40LG", "TNFRSF4", "TNFRSF18"
  )
)

programme_labels <- c(
  Cycling = "Cycling",
  Exhaustion = "Exhaustion",
  Glycolysis = "Glycolysis",
  OXPHOS = "OXPHOS",
  p53_checkpoint_stress = "p53 / checkpoint stress",
  mTOR_MYC = "mTOR / MYC",
  Cytotoxicity = "Cytotoxicity",
  TypeI_IFN_ISG = "Type-I IFN / ISG",
  TRM = "TRM",
  Apoptosis_death_response = "Apoptosis / death response",
  Checkpoint_receptors = "Checkpoint receptors",
  TCF7_memory = "TCF7-memory",
  TCR_signalling = "TCR signalling",
  NFkB = "NF-κB",
  Costimulation = "Costimulation"
)

## ------------------------------------------------------------
## 3. Keep only genes present in object
## ------------------------------------------------------------

features_obj <- rownames(obj)
feature_map <- setNames(features_obj, toupper(features_obj))

programmes_present <- lapply(programmes_raw, function(g) {
  g_upper <- toupper(g)
  present <- intersect(g_upper, names(feature_map))
  unname(feature_map[present])
})

coverage_tbl <- data.frame(
  Programme = names(programmes_raw),
  n_input_genes = lengths(programmes_raw),
  n_present_genes = lengths(programmes_present),
  present_genes = sapply(programmes_present, paste, collapse = ";")
)

write.csv(
  coverage_tbl,
  file = file.path(out_dir, "Programme_gene_coverage.csv"),
  row.names = FALSE
)

programmes_present <- programmes_present[lengths(programmes_present) >= 3]

## ------------------------------------------------------------
## 4. Add module scores
## ------------------------------------------------------------

for (prog in names(programmes_present)) {
  
  message("Scoring programme: ", prog)
  
  tmp_prefix <- paste0("MODULE_", prog)
  
  obj <- AddModuleScore(
    object = obj,
    features = list(programmes_present[[prog]]),
    assay = "RNA",
    name = tmp_prefix
  )
  
  generated_col <- paste0(tmp_prefix, "1")
  final_col <- paste0(prog, "_score")
  
  obj[[final_col]] <- obj@meta.data[[generated_col]]
  obj@meta.data[[generated_col]] <- NULL
}

## ------------------------------------------------------------
## 5. Z-normalise each programme within dataset
## ------------------------------------------------------------

z_by_group <- function(x, group) {
  ave(x, group, FUN = function(v) {
    s <- sd(v, na.rm = TRUE)
    m <- mean(v, na.rm = TRUE)
    
    if (is.na(s) || s == 0) {
      rep(NA_real_, length(v))
    } else {
      (v - m) / s
    }
  })
}

for (prog in names(programmes_present)) {
  
  score_col <- paste0(prog, "_score")
  z_col <- paste0(prog, "_z")
  
  obj[[z_col]] <- z_by_group(
    obj@meta.data[[score_col]],
    obj@meta.data[[dataset_col]]
  )
}

## ------------------------------------------------------------
## 6. Define Prolif-Tex and NonProlif-Tex
## ------------------------------------------------------------

## Classification follows Methods 4.1.7, which uses a nine-gene exhaustion
## panel. The Exhaustion PROGRAMME scored above has twelve genes and is used
## for the programme-level comparison only, not for assigning configurations.
## The Cycling programme above is now identical to the 4.1.7 proliferation
## panel, so it serves directly as the proliferation axis.

exhaustion_genes_417 <- c(
  "PDCD1", "LAYN", "HAVCR2", "CTLA4", "TOX",
  "TIGIT", "LAG3", "CXCL13", "ENTPD1"
)
exhaustion_genes_417 <- exhaustion_genes_417[exhaustion_genes_417 %in% rownames(obj)]

obj <- AddModuleScore(
  object   = obj,
  features = list(exhaustion_genes_417),
  assay    = "RNA",
  name     = "Exh417"
)
obj$Exh417_z <- z_by_group(obj@meta.data$Exh4171, obj@meta.data[[dataset_col]])

exh_z <- obj@meta.data$Exh417_z
prolif_z <- obj@meta.data$Cycling_z

exh_q75 <- quantile(exh_z, 0.75, na.rm = TRUE)
prolif_q75 <- quantile(prolif_z, 0.75, na.rm = TRUE)
prolif_q25 <- quantile(prolif_z, 0.25, na.rm = TRUE)

obj$FunctionalConfig_programme_final <- dplyr::case_when(
  exh_z >= exh_q75 & prolif_z >= prolif_q75 ~ "Prolif-Tex",
  exh_z >= exh_q75 & prolif_z <= prolif_q25 ~ "NonProlif-Tex",
  TRUE ~ "Other"
)

threshold_tbl <- data.frame(
  threshold = c("Exhaustion Q75", "Proliferation Q75", "Proliferation Q25"),
  value = c(exh_q75, prolif_q75, prolif_q25)
)

write.csv(
  threshold_tbl,
  file = file.path(out_dir, "FunctionalConfig_thresholds.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(obj$FunctionalConfig_programme_final)),
  file = file.path(out_dir, "FunctionalConfig_cell_counts.csv"),
  row.names = FALSE
)

## ------------------------------------------------------------
## 7. Build patient-level paired table
## ------------------------------------------------------------

min_cells_per_patient_group <- 10

z_cols <- paste0(names(programmes_present), "_z")

meta <- obj@meta.data
meta$CellID <- rownames(meta)

keep_cols <- c(
  "CellID",
  patient_col,
  dataset_col,
  "FunctionalConfig_programme_final",
  z_cols
)

meta2 <- meta[, keep_cols, drop = FALSE]

colnames(meta2)[colnames(meta2) == patient_col] <- "PatientID"
colnames(meta2)[colnames(meta2) == dataset_col] <- "Dataset"
colnames(meta2)[colnames(meta2) == "FunctionalConfig_programme_final"] <- "FunctionalConfig"

meta_long <- meta2 %>%
  filter(FunctionalConfig %in% c("Prolif-Tex", "NonProlif-Tex")) %>%
  mutate(PatientKey = paste(Dataset, PatientID, sep = "__")) %>%
  pivot_longer(
    cols = all_of(z_cols),
    names_to = "Programme_col",
    values_to = "score_z"
  ) %>%
  mutate(
    Programme = sub("_z$", "", Programme_col),
    FunctionalConfig_short = recode(
      FunctionalConfig,
      "Prolif-Tex" = "ProlifTex",
      "NonProlif-Tex" = "NonProlifTex"
    )
  )

patient_programme_scores <- meta_long %>%
  group_by(PatientKey, Dataset, PatientID, FunctionalConfig_short, Programme) %>%
  summarise(
    n_cells = n(),
    median_score = median(score_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_patient_group)

paired_scores <- patient_programme_scores %>%
  pivot_wider(
    names_from = FunctionalConfig_short,
    values_from = c(median_score, n_cells)
  ) %>%
  filter(
    !is.na(median_score_ProlifTex),
    !is.na(median_score_NonProlifTex)
  ) %>%
  mutate(
    paired_diff = median_score_ProlifTex - median_score_NonProlifTex
  )

write.csv(
  patient_programme_scores,
  file = file.path(out_dir, "Patient_level_programme_scores_long.csv"),
  row.names = FALSE
)

write.csv(
  paired_scores,
  file = file.path(out_dir, "Patient_level_programme_scores_paired.csv"),
  row.names = FALSE
)

## ------------------------------------------------------------
## 8. Statistics
## ------------------------------------------------------------

boot_median_ci <- function(x, B = 2000) {
  x <- x[is.finite(x)]
  
  if (length(x) < 3) {
    return(c(NA_real_, NA_real_))
  }
  
  boot_meds <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  as.numeric(quantile(boot_meds, probs = c(0.025, 0.975), na.rm = TRUE))
}

wilcox_paired_safe <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) return(NA_real_)
  if (all(x == y)) return(1)
  
  tryCatch(
    wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

rank_biserial_paired <- function(d) {
  d <- d[is.finite(d)]
  d <- d[d != 0]
  
  n <- length(d)
  if (n < 1) return(NA_real_)
  
  r <- rank(abs(d))
  W_pos <- sum(r[d > 0])
  W_neg <- sum(r[d < 0])
  W_total <- n * (n + 1) / 2
  
  (W_pos - W_neg) / W_total
}

stats_tbl <- paired_scores %>%
  group_by(Programme) %>%
  summarise(
    n_patients = n(),
    median_diff = median(paired_diff, na.rm = TRUE),
    ci_low = boot_median_ci(paired_diff)[1],
    ci_high = boot_median_ci(paired_diff)[2],
    p_value = wilcox_paired_safe(
      median_score_ProlifTex,
      median_score_NonProlifTex
    ),
    rank_biserial = rank_biserial_paired(paired_diff),
    .groups = "drop"
  ) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    direction = ifelse(
      median_diff >= 0,
      "Prolif-Tex higher",
      "NonProlif-Tex higher"
    ),
    sig_label = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01 ~ "**",
      q_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Programme_label = programme_labels[Programme]
  )

write.csv(
  stats_tbl,
  file = file.path(out_dir, "Programme_level_statistics.csv"),
  row.names = FALSE
)

## ------------------------------------------------------------
## 9. Clean plotting setup
## ------------------------------------------------------------

cols <- c(
  "Prolif-Tex higher" = "#B2182B",
  "NonProlif-Tex higher" = "#2166AC"
)

theme_pub <- theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(face = "bold", size = 25, hjust = 0),
    plot.subtitle = element_text(size = 15, hjust = 0),
    plot.caption = element_text(size = 12, hjust = 0),
    plot.title.position = "plot",
    axis.title = element_text(face = "bold", size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 15),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    plot.margin = margin(t = 25, r = 55, b = 25, l = 25)
  )

get_xlim <- function(df) {
  rng <- range(c(df$ci_low, df$ci_high, df$median_diff, 0), na.rm = TRUE)
  pad <- diff(rng) * 0.08
  if (!is.finite(pad) || pad == 0) pad <- 0.05
  c(rng[1] - pad, rng[2] + pad)
}

## ------------------------------------------------------------
## 10. Plot A — definition axes
## ------------------------------------------------------------

definition_order_top <- c("Cycling", "Exhaustion")

stats_definition <- stats_tbl %>%
  filter(Programme %in% definition_order_top) %>%
  mutate(
    y_num = match(Programme, rev(definition_order_top))
  ) %>%
  arrange(y_num)

p_definition <- ggplot(stats_definition, aes(x = median_diff, y = y_num)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0,
    linewidth = 0.9,
    color = "grey40"
  ) +
  geom_point(
    aes(fill = direction),
    shape = 21,
    size = 5.5,
    color = "black",
    stroke = 0.8
  ) +
  geom_text(
    aes(label = sig_label),
    nudge_y = 0.25,
    size = 7,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(
    breaks = stats_definition$y_num,
    labels = stats_definition$Programme_label,
    expand = expansion(mult = c(0.20, 0.35))
  ) +
  coord_cartesian(
    xlim = get_xlim(stats_definition),
    clip = "off"
  ) +
  labs(
    title = "Definition axes",
    subtitle = "Cycling is expected by definition; exhaustion verifies that both groups remain within the exhausted compartment",
    x = "Median paired difference: Prolif-Tex − NonProlif-Tex",
    y = NULL
  ) +
  theme_pub

## ------------------------------------------------------------
## 11. Plot B — non-defining programme differences
## ------------------------------------------------------------

nondef_order_top <- c(
  "Glycolysis",
  "OXPHOS",
  "p53_checkpoint_stress",
  "mTOR_MYC",
  "Cytotoxicity",
  "TypeI_IFN_ISG",
  "TRM",
  "Apoptosis_death_response",
  "Checkpoint_receptors",
  "TCF7_memory",
  "TCR_signalling",
  "NFkB",
  "Costimulation"
)

stats_nondef <- stats_tbl %>%
  filter(Programme %in% nondef_order_top) %>%
  mutate(
    y_num = match(Programme, rev(nondef_order_top))
  ) %>%
  arrange(y_num)

p_programmes <- ggplot(stats_nondef, aes(x = median_diff, y = y_num)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0,
    linewidth = 0.9,
    color = "grey40"
  ) +
  geom_point(
    aes(fill = direction),
    shape = 21,
    size = 5.5,
    color = "black",
    stroke = 0.8
  ) +
  geom_text(
    aes(label = sig_label),
    nudge_y = 0.25,
    size = 7,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(
    breaks = stats_nondef$y_num,
    labels = stats_nondef$Programme_label,
    expand = expansion(mult = c(0.05, 0.13))
  ) +
  coord_cartesian(
    xlim = get_xlim(stats_nondef),
    clip = "off"
  ) +
  labs(
    title = "Programme-level differences between Prolif-Tex and NonProlif-Tex",
    subtitle = "Patient-level paired medians; lines show bootstrap 95% CI; stars show BH-adjusted paired Wilcoxon significance",
    x = "Median paired difference: Prolif-Tex − NonProlif-Tex",
    y = NULL,
    caption = "*** q < 0.001; ** q < 0.01; * q < 0.05. Full statistics exported separately."
  ) +
  theme_pub

## ------------------------------------------------------------
## 12. Save plots
## ------------------------------------------------------------

ggsave(
  filename = file.path(out_dir, "Figure_Definition_axes.pdf"),
  plot = p_definition,
  width = 14,
  height = 5,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Figure_Definition_axes.png"),
  plot = p_definition,
  width = 14,
  height = 5,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Figure_Programme_level_differences.pdf"),
  plot = p_programmes,
  width = 16,
  height = 9.5,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Figure_Programme_level_differences.png"),
  plot = p_programmes,
  width = 16,
  height = 9.5,
  units = "in",
  dpi = 600
)

## Optional combined figure if patchwork is installed
if (requireNamespace("patchwork", quietly = TRUE)) {
  
  library(patchwork)
  
  p_combined <- p_definition / p_programmes +
    plot_layout(heights = c(0.32, 0.68))
  
  ggsave(
    filename = file.path(out_dir, "Figure_Combined_definition_and_programmes.pdf"),
    plot = p_combined,
    width = 16,
    height = 14,
    units = "in"
  )
  
  ggsave(
    filename = file.path(out_dir, "Figure_Combined_definition_and_programmes.png"),
    plot = p_combined,
    width = 16,
    height = 14,
    units = "in",
    dpi = 600
  )
}

## ------------------------------------------------------------
## 13. Return updated object to environment
## ------------------------------------------------------------

CD8_TumourOnly_scVI <- obj

## Show final plots
p_definition
p_programmes