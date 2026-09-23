# ============================================================
# 14_prevalence_plot_final.R
# ============================================================
# Thesis Methods section: 4.1.8   |   Produces: Figure 7A
#
# Requires: patient_tex in the R session, produced by script 12.
#
# This file contains two alternative layouts for the same comparison:
#   FIX A - all pairwise Wilcoxon comparisons across tissue compartments,
#           BH-adjusted, faceted by TexClass.  THIS IS FIGURE 7A.
#   FIX B - a single panel showing only Tumour-versus-other contrasts.
#           Not used in the thesis; retained as it was part of the file.
#
# Figure 7A's legend describes pairwise comparisons across compartments within
# each class, which corresponds to Fix A.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggpubr)
})

# ============================================================
# CHUNK 1 — Safety checks + consistent Tissue ordering
# Purpose: Ensure required objects/columns exist and Tissue is ordered.
# ============================================================
stopifnot(exists("patient_tex"))
stopifnot(all(c("PatientID","Tissue","TexClass","prevalence") %in% colnames(patient_tex)))

tissue_levels <- c("Blood","Normal","Peritumour","Tumour")
patient_tex$Tissue <- factor(patient_tex$Tissue, levels = tissue_levels)

# Ensure TexClass is clean
patient_tex$TexClass <- factor(patient_tex$TexClass, levels = c("Prolif-Tex","Non-Prolif-Tex"))

# Utility: compute y-positions within each TexClass so brackets don't overlap inside facets
add_y_positions_by_group <- function(pw_tbl, df, step = 0.05, pad = 0.02) {
  pw_tbl %>%
    group_by(TexClass) %>%
    arrange(group1, group2) %>%
    mutate(
      y.position = max(df$prevalence[df$TexClass == first(TexClass)], na.rm = TRUE) +
        pad + step * row_number()
    ) %>%
    ungroup()
}

# ============================================================
# FIX A — Facet by TexClass (cleanest)
# Purpose: Show all pairwise tissue p-values separately for each Tex class.
# ============================================================

# --- Compute BH-adjusted pairwise Wilcoxon tests within each TexClass
pw_byclass_A <- compare_means(
  prevalence ~ Tissue,
  data = patient_tex,
  group.by = "TexClass",
  method = "wilcox.test",
  p.adjust.method = "BH"
)

# --- Add y-positions (within each facet)
pw_byclass_A <- add_y_positions_by_group(pw_byclass_A, patient_tex, step = 0.05, pad = 0.02)

# --- Plot (faceted)
p_fixA <- ggplot(patient_tex, aes(x = Tissue, y = prevalence)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.6) +
  stat_pvalue_manual(
    pw_byclass_A,
    label = "p.format",
    tip.length = 0.01,
    hide.ns = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.30))
  ) +
  facet_wrap(~TexClass, nrow = 1) +
  labs(
    title = "Patient-level prevalence of Prolif-Tex and Non-Prolif-Tex",
    subtitle = "All pairwise Wilcoxon comparisons across Tissue within each class (BH-adjusted)",
    x = "Tissue compartment",
    y = "Fraction of CD8 cells"
  ) +
  theme_classic(base_size = 12)

print(p_fixA)

# Optional: print the p-value table used for Fix A
message("=== Fix A: Pairwise p-values by TexClass (BH-adjusted) ===")
print(pw_byclass_A %>% select(TexClass, group1, group2, p, p.adj, p.format, p.signif))


# ============================================================
# FIX B — One panel, but only key comparisons (Tumour vs others)
# Purpose: Avoid p-value clutter while keeping the most important contrasts.
# Note: Full pairwise table can be reported in Supplementary (Fix A table).
# ============================================================

key_comparisons <- list(
  c("Blood", "Tumour"),
  c("Normal", "Tumour"),
  c("Peritumour", "Tumour")
)

# For Fix B we compute p-values for ONLY those comparisons, within each TexClass
pw_byclass_B <- compare_means(
  prevalence ~ Tissue,
  data = patient_tex %>% filter(Tissue %in% tissue_levels),  # keep ordered
  group.by = "TexClass",
  method = "wilcox.test",
  comparisons = key_comparisons,
  p.adjust.method = "BH"
)

# Add y-positions (stagger within class; then we’ll dodge by class visually)
pw_byclass_B <- pw_byclass_B %>%
  group_by(TexClass) %>%
  arrange(group1, group2) %>%
  mutate(
    y.position = max(patient_tex$prevalence[patient_tex$TexClass == first(TexClass)], na.rm = TRUE) +
      0.02 + 0.06 * row_number()
  ) %>%
  ungroup()

p_fixB <- ggplot(patient_tex, aes(x = Tissue, y = prevalence, fill = TexClass)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(color = TexClass),
    size = 1.8,
    alpha = 0.6,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    )
  ) +
  stat_pvalue_manual(
    pw_byclass_B,
    label = "p.format",
    tip.length = 0.01,
    hide.ns = FALSE
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.30))
  ) +
  labs(
    title = "Patient-level prevalence of Prolif-Tex and Non-Prolif-Tex",
    subtitle = "Key contrasts shown (Tumour vs others); Wilcoxon within each class (BH-adjusted)",
    x = "Tissue compartment",
    y = "Fraction of CD8 cells",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12)

print(p_fixB)

# Optional: print the p-value table used for Fix B
message("=== Fix B: Key p-values by TexClass (BH-adjusted) ===")
print(pw_byclass_B %>% select(TexClass, group1, group2, p, p.adj, p.format, p.signif))
