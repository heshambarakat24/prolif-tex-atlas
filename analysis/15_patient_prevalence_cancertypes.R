# ============================================================
# 15_patient_prevalence_cancertypes.R
# ============================================================
# Thesis Methods section: 4.1.8 (patient-level prevalence across cancer types)
# Produces: Figures 7B and 7C
#
# Requires: CD8_TumourOnly_scVI in the R session, produced by
#           09_reconstruct_tumour_only_seurat.R. Not loaded from disk here.
# Writes  : Patient_level_Tex_prevalence_*.csv, CancerType_*.csv,
#           CancerType_patient_level_prevalence_main.png and related figures.
#
# Original file header marked this VERSION v1.0, FROZEN FOR THESIS FIGURES,
# dated 2026-04-08, with an instruction not to modify it in place. That header
# is retained below.
#
# Notes:
#  - Gene panels are the nine exhaustion and nine proliferation genes reported
#    in Methods 4.1.7 and deposited in panels/.
#  - Module scores and percentile thresholds are recomputed here on the
#    TUMOUR-ONLY object. Script 11 computes them on the FULL atlas. Both follow
#    Methods 4.1.8, which recalculates per analysis object, but because
#    z-normalisation runs within Dataset over different cell sets, the numeric
#    thresholds printed by the two scripts differ. This is expected.
#  - Significance brackets are annotated manually rather than with
#    stat_pvalue_manual(); see the original note below.
#  - Input data are not included in this repository. See README.md.
# ============================================================

# ============================================================
# VERSION: v1.0 (FROZEN FOR THESIS FIGURES)
# DATE: 2026-04-08
#
# DESCRIPTION:
# Patient-level prevalence of Prolif-Tex and Non-Prolif-Tex
# across cancer types from the tumour-only scVI Seurat object.
#
# Main logic:
# 1) Compute exhaustion and proliferation module scores
# 2) Z-normalise scores within dataset
# 3) Classify cells using pooled percentile thresholds
# 4) Summarise prevalence at patient level within cancer type
# 5) Keep cancer types with >= 5 patients for the main figure
# 6) Plot patient-level distributions across cancer types
# 7) Run Kruskal-Wallis + pairwise Wilcoxon tests within each class
# 8) Export tables and save clean + annotated figures
#
# NOTES:
# - Uses tumour-only scVI Seurat object
# - Patient is the unit of analysis
# - Main figure is descriptive + patient-level
# - Full pairwise statistics are exported to CSV
# - Manual bracket annotation is used instead of stat_pvalue_manual()
#   because ggpubr produced unstable errors in this workflow
#
# IMPORTANT:
# This version is frozen for thesis figures.
# Do not modify directly. Create a new version instead.
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(scales)
library(stringr)
library(forcats)

# ------------------------------------------------------------
# 0) INPUT OBJECT
# ------------------------------------------------------------
obj <- CD8_TumourOnly_scVI

# ------------------------------------------------------------
# 1) BASIC CHECKS
# ------------------------------------------------------------
cat("\nMetadata columns:\n")
print(colnames(obj@meta.data))

cat("\nTissue levels:\n")
print(table(obj$Tissue, useNA = "ifany"))

cat("\nCancer types:\n")
print(sort(unique(obj$CancerType)))

cat("\nNumber of unique patients:\n")
print(length(unique(obj$PatientID)))

# ------------------------------------------------------------
# 2) CURATED SIGNATURES
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

# ------------------------------------------------------------
# 3) ADD MODULE SCORES
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 4) EXTRACT METADATA
# ------------------------------------------------------------
meta <- obj@meta.data %>%
  mutate(
    ExhaustionScore_raw = ExhaustionScore1,
    ProliferationScore_raw = ProliferationScore1,
    Dataset = as.character(Dataset),
    PatientID = as.character(PatientID),
    CancerType = as.character(CancerType),
    Tissue = as.character(Tissue)
  )

# ------------------------------------------------------------
# 5) OPTIONAL: confirm tumour-only object
# ------------------------------------------------------------
cat("\nTissue distribution in object:\n")
print(table(meta$Tissue, useNA = "ifany"))

# ------------------------------------------------------------
# 6) Z-NORMALISE WITHIN DATASET
# ------------------------------------------------------------
meta <- meta %>%
  group_by(Dataset) %>%
  mutate(
    ExhaustionScore_z = as.numeric(scale(ExhaustionScore_raw)),
    ProliferationScore_z = as.numeric(scale(ProliferationScore_raw))
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 7) CLASSIFY CELLS USING POOLED THRESHOLDS
# ------------------------------------------------------------
exh_q75  <- quantile(meta$ExhaustionScore_z, probs = 0.75, na.rm = TRUE)
prol_q75 <- quantile(meta$ProliferationScore_z, probs = 0.75, na.rm = TRUE)
prol_q25 <- quantile(meta$ProliferationScore_z, probs = 0.25, na.rm = TRUE)

cat("\nThresholds used:\n")
cat("Exhaustion Q75 =", round(exh_q75, 4), "\n")
cat("Proliferation Q75 =", round(prol_q75, 4), "\n")
cat("Proliferation Q25 =", round(prol_q25, 4), "\n")

meta <- meta %>%
  mutate(
    TexClass = case_when(
      ExhaustionScore_z >= exh_q75 & ProliferationScore_z >= prol_q75 ~ "Prolif-Tex",
      ExhaustionScore_z >= exh_q75 & ProliferationScore_z <= prol_q25 ~ "Non-Prolif-Tex",
      TRUE ~ "Other"
    )
  )

cat("\nCell counts by class:\n")
print(table(meta$TexClass, useNA = "ifany"))

# ------------------------------------------------------------
# 8) PATIENT-LEVEL PREVALENCE WITHIN CANCER TYPE
# ------------------------------------------------------------
patient_summary <- meta %>%
  group_by(CancerType, PatientID) %>%
  summarise(
    total_cells = n(),
    prolif_tex_cells = sum(TexClass == "Prolif-Tex"),
    nonprolif_tex_cells = sum(TexClass == "Non-Prolif-Tex"),
    prolif_tex_prop = prolif_tex_cells / total_cells,
    nonprolif_tex_prop = nonprolif_tex_cells / total_cells,
    .groups = "drop"
  )

# ------------------------------------------------------------
# 9) PATIENT COUNTS PER CANCER TYPE
# ------------------------------------------------------------
n_patients_df <- patient_summary %>%
  group_by(CancerType) %>%
  summarise(
    n_patients = n_distinct(PatientID),
    median_cells_per_patient = median(total_cells),
    min_cells_per_patient = min(total_cells),
    max_cells_per_patient = max(total_cells),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients), CancerType)

cat("\nPatient counts per cancer type:\n")
print(n_patients_df)

# ------------------------------------------------------------
# 10) FILTER MAIN FIGURE TO CANCER TYPES WITH >= 5 PATIENTS
# ------------------------------------------------------------
eligible_cancers <- n_patients_df %>%
  filter(n_patients >= 5)

excluded_cancers <- n_patients_df %>%
  filter(n_patients < 5)

cat("\nCancer types kept in main figure (n >= 5):\n")
print(eligible_cancers %>% select(CancerType, n_patients))

cat("\nCancer types excluded from main figure (n < 5):\n")
print(excluded_cancers %>% select(CancerType, n_patients))

patient_summary_main <- patient_summary %>%
  semi_join(eligible_cancers %>% select(CancerType), by = "CancerType")

patient_summary_supp <- patient_summary %>%
  semi_join(excluded_cancers %>% select(CancerType), by = "CancerType")

# ------------------------------------------------------------
# 11) LONG FORMAT FOR PLOTTING / TESTING
# ------------------------------------------------------------
plot_df <- patient_summary_main %>%
  select(CancerType, PatientID, total_cells, prolif_tex_prop, nonprolif_tex_prop) %>%
  pivot_longer(
    cols = c(prolif_tex_prop, nonprolif_tex_prop),
    names_to = "Class",
    values_to = "Prevalence"
  ) %>%
  mutate(
    Class = recode(
      Class,
      prolif_tex_prop = "Prolif-Tex",
      nonprolif_tex_prop = "Non-Prolif-Tex"
    )
  )

# ------------------------------------------------------------
# 12) ORDER CANCER TYPES BY MEDIAN PREVALENCE WITHIN CLASS
#     Use Prolif-Tex ordering for both figures for consistency
# ------------------------------------------------------------
order_df <- plot_df %>%
  filter(Class == "Prolif-Tex") %>%
  group_by(CancerType) %>%
  summarise(med_prev = median(Prevalence), .groups = "drop") %>%
  left_join(eligible_cancers %>% select(CancerType, n_patients), by = "CancerType") %>%
  mutate(label = paste0(CancerType, " (n=", n_patients, ")")) %>%
  arrange(med_prev)

label_map <- eligible_cancers %>%
  select(CancerType, n_patients) %>%
  mutate(label = paste0(CancerType, " (n=", n_patients, ")"))

plot_df <- plot_df %>%
  left_join(label_map, by = "CancerType") %>%
  mutate(
    label = factor(label, levels = order_df$label)
  )

# ------------------------------------------------------------
# 13) SUMMARY TABLES FOR DISPLAY / EXPORT
# ------------------------------------------------------------
summary_stats <- plot_df %>%
  group_by(Class, CancerType, label) %>%
  summarise(
    n_patients = n(),
    median_prevalence = median(Prevalence),
    mean_prevalence = mean(Prevalence),
    IQR_prevalence = IQR(Prevalence),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 14) OMNIBUS KRUSKAL-WALLIS TESTS WITHIN EACH CLASS
# ------------------------------------------------------------
kw_results <- plot_df %>%
  group_by(Class) %>%
  group_modify(~{
    kt <- kruskal.test(Prevalence ~ label, data = .x)
    tibble(
      statistic = unname(kt$statistic),
      df = unname(kt$parameter),
      p_value = kt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  )

cat("\nKruskal-Wallis results:\n")
print(kw_results)

# ------------------------------------------------------------
# 15) PAIRWISE WILCOXON TESTS WITHIN EACH CLASS
# ------------------------------------------------------------
pairwise_results <- map_dfr(unique(plot_df$Class), function(this_class) {
  
  tmp <- plot_df %>% filter(Class == this_class)
  
  pw <- pairwise.wilcox.test(
    x = tmp$Prevalence,
    g = tmp$label,
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  pw_mat <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE)
  colnames(pw_mat) <- c("group1", "group2", "p_adj")
  
  pw_mat %>%
    filter(!is.na(p_adj)) %>%
    mutate(Class = this_class) %>%
    select(Class, group1, group2, p_adj)
})

group_medians <- plot_df %>%
  group_by(Class, label) %>%
  summarise(median_prev = median(Prevalence), .groups = "drop")

pairwise_results <- pairwise_results %>%
  left_join(group_medians, by = c("Class", "group1" = "label")) %>%
  rename(group1_median = median_prev) %>%
  left_join(group_medians, by = c("Class", "group2" = "label")) %>%
  rename(group2_median = median_prev) %>%
  mutate(
    abs_median_diff = abs(group1_median - group2_median),
    significant = p_adj < 0.05
  ) %>%
  arrange(Class, p_adj)

cat("\nTop pairwise results:\n")
print(head(pairwise_results, 20))

# ------------------------------------------------------------
# 16) SELECT SIGNIFICANT COMPARISONS FOR OPTIONAL ANNOTATED PLOT
#     Too many brackets destroy readability, so keep only:
#     - BH-adjusted p < 0.05
#     - top 8 per class by adjusted p-value
# ------------------------------------------------------------
pairwise_sig_for_plot <- pairwise_results %>%
  filter(significant) %>%
  group_by(Class) %>%
  arrange(p_adj, desc(abs_median_diff), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

max_y_df <- plot_df %>%
  group_by(Class) %>%
  summarise(max_y = max(Prevalence, na.rm = TRUE), .groups = "drop")

pairwise_sig_for_plot <- pairwise_sig_for_plot %>%
  left_join(max_y_df, by = "Class") %>%
  group_by(Class) %>%
  mutate(
    y.position = max_y + seq(0.03, by = 0.03, length.out = n()),
    p_label = case_when(
      p_adj < 1e-4 ~ format(p_adj, scientific = TRUE, digits = 2),
      TRUE ~ sprintf("%.4f", p_adj)
    )
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 17) MAIN CLEAN FIGURE
# ------------------------------------------------------------
kw_text <- kw_results %>%
  mutate(
    label_txt = paste0(
      Class, ": Kruskal-Wallis BH-adjusted p = ",
      ifelse(
        p_adj_BH < 1e-4,
        format(p_adj_BH, scientific = TRUE, digits = 2),
        sprintf("%.4f", p_adj_BH)
      )
    )
  ) %>%
  pull(label_txt)

subtitle_text <- paste(
  "Patient is the unit of analysis. Only cancer types with n >= 5 patients are shown.",
  "Boxplots summarize patient-level prevalence; dots indicate individual patients.",
  paste(kw_text, collapse = " | "),
  sep = "\n"
)

p_main <- ggplot(plot_df, aes(x = label, y = Prevalence)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    fill = "grey88",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.16,
    height = 0,
    size = 1.6,
    alpha = 0.8,
    colour = "grey20"
  ) +
  facet_wrap(~ Class, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Patient-level prevalence of Prolif-Tex and Non-Prolif-Tex across cancer types",
    subtitle = subtitle_text,
    x = "Cancer type",
    y = "Fraction of tumour CD8 cells"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_main)

# ------------------------------------------------------------
# 18) OPTIONAL ANNOTATED FIGURE WITH TOP SIGNIFICANT PAIRWISE TESTS
#     REMOVED: ggpubr::stat_pvalue_manual() was unstable here.
#     Use chunk 20 manual annotation instead.
# ------------------------------------------------------------

# ------------------------------------------------------------
# 19) OPTIONAL SEPARATE FIGURES PER CLASS
#     Sometimes cleaner for papers / thesis figures
# ------------------------------------------------------------
plot_single_class <- function(df, this_class, subtitle_suffix = NULL) {
  
  tmp <- df %>% filter(Class == this_class)
  
  kw_p <- kw_results %>%
    filter(Class == this_class) %>%
    pull(p_adj_BH)
  
  kw_p_txt <- ifelse(
    kw_p < 1e-4,
    format(kw_p, scientific = TRUE, digits = 2),
    sprintf("%.4f", kw_p)
  )
  
  subtitle_text <- paste(
    "Patient-level prevalence across cancer types (n >= 5 patients per cancer type).",
    paste0("Kruskal-Wallis BH-adjusted p = ", kw_p_txt),
    subtitle_suffix
  )
  
  ggplot(tmp, aes(x = label, y = Prevalence)) +
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA,
      fill = "grey88",
      colour = "black",
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.16,
      height = 0,
      size = 1.7,
      alpha = 0.8,
      colour = "grey20"
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(
      title = paste0(this_class, " prevalence by cancer type"),
      subtitle = subtitle_text,
      x = "Cancer type",
      y = "Fraction of tumour CD8 cells"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

p_prolif <- plot_single_class(plot_df, "Prolif-Tex")
p_nonprolif <- plot_single_class(plot_df, "Non-Prolif-Tex")

print(p_prolif)
print(p_nonprolif)

# ------------------------------------------------------------
# 20) OPTIONAL ANNOTATED SINGLE-CLASS FIGURES (MANUAL BRACKETS)
#     Compact version: same result, robust to empty annotation sets
# ------------------------------------------------------------
add_manual_annotations <- function(base_plot, class_name, plot_df, pairwise_sig_for_plot) {
  
  ann <- pairwise_sig_for_plot %>%
    filter(Class == class_name) %>%
    transmute(
      group1 = as.character(group1),
      group2 = as.character(group2),
      p_label = as.character(p_label)
    )
  
  if (nrow(ann) == 0) {
    message("No retained pairwise annotations for ", class_name, ".")
    return(base_plot)
  }
  
  x_levels <- levels(plot_df %>% filter(Class == class_name) %>% pull(label))
  x_map <- setNames(seq_along(x_levels), x_levels)
  
  ymax <- max(plot_df$Prevalence[plot_df$Class == class_name], na.rm = TRUE)
  tick <- ymax * 0.012
  text_offset <- ymax * 0.035
  
  ann <- ann %>%
    mutate(
      x1 = unname(x_map[group1]),
      x2 = unname(x_map[group2]),
      xmid = (x1 + x2) / 2
    ) %>%
    filter(!is.na(x1) & !is.na(x2)) %>%
    arrange(x2 - x1, x1) %>%
    mutate(
      y = ymax + seq(
        from = ymax * 0.07,
        by   = ymax * 0.07,
        length.out = n()
      )
    )
  
  base_plot +
    geom_segment(
      data = ann,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3
    ) +
    geom_segment(
      data = ann,
      aes(x = x1, xend = x1, y = y - tick, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3
    ) +
    geom_segment(
      data = ann,
      aes(x = x2, xend = x2, y = y - tick, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3
    ) +
    geom_text(
      data = ann,
      aes(x = xmid, y = y + text_offset, label = p_label),
      inherit.aes = FALSE,
      size = 3
    ) +
    coord_cartesian(
      ylim = c(0, max(ann$y + text_offset) * 1.03),
      clip = "off"
    )
}

p_prolif_annot <- add_manual_annotations(
  base_plot = p_prolif,
  class_name = "Prolif-Tex",
  plot_df = plot_df,
  pairwise_sig_for_plot = pairwise_sig_for_plot
)

p_nonprolif_annot <- add_manual_annotations(
  base_plot = p_nonprolif,
  class_name = "Non-Prolif-Tex",
  plot_df = plot_df,
  pairwise_sig_for_plot = pairwise_sig_for_plot
)

print(p_prolif_annot)
print(p_nonprolif_annot)

# ------------------------------------------------------------
# 21) EXPORT TABLES
# ------------------------------------------------------------
write.csv(patient_summary, "Patient_level_Tex_prevalence_all_cancers.csv", row.names = FALSE)
write.csv(patient_summary_main, "Patient_level_Tex_prevalence_main_n_ge_5.csv", row.names = FALSE)
write.csv(patient_summary_supp, "Patient_level_Tex_prevalence_excluded_n_lt_5.csv", row.names = FALSE)

write.csv(n_patients_df, "CancerType_patient_counts.csv", row.names = FALSE)
write.csv(eligible_cancers, "CancerType_included_main_figure_n_ge_5.csv", row.names = FALSE)
write.csv(excluded_cancers, "CancerType_excluded_main_figure_n_lt_5.csv", row.names = FALSE)

write.csv(summary_stats, "CancerType_patient_level_prevalence_summary_stats.csv", row.names = FALSE)
write.csv(kw_results, "CancerType_Kruskal_Wallis_results.csv", row.names = FALSE)
write.csv(pairwise_results, "CancerType_pairwise_Wilcoxon_BH_results.csv", row.names = FALSE)
write.csv(pairwise_sig_for_plot, "CancerType_top_significant_pairwise_for_plot.csv", row.names = FALSE)

# ------------------------------------------------------------
# 22) SAVE FIGURES
# ------------------------------------------------------------
ggsave("CancerType_patient_level_prevalence_main.png", p_main, width = 14, height = 10, dpi = 300)

ggsave("ProlifTex_patient_level_prevalence_by_cancer_type.png", p_prolif, width = 13, height = 7, dpi = 300)
ggsave("NonProlifTex_patient_level_prevalence_by_cancer_type.png", p_nonprolif, width = 13, height = 7, dpi = 300)

ggsave("ProlifTex_patient_level_prevalence_by_cancer_type_annotated.png", p_prolif_annot, width = 14, height = 8, dpi = 300)
ggsave("NonProlifTex_patient_level_prevalence_by_cancer_type_annotated.png", p_nonprolif_annot, width = 14, height = 8, dpi = 300)