# ============================================================
# 13_patient_pvalues_compartments.R
# ============================================================
# Thesis Methods section: 4.1.8 (pairwise tests across tissue compartments)
#
# Requires: patient_tex in the R session, produced by script 12.
# Writes  : pairwise_pvalues_Tissue_combinedTex.csv
#           pairwise_pvalues_Tissue_byTexClass.csv
#
# Notes:
#  - Pairwise Wilcoxon rank-sum tests, Benjamini-Hochberg adjusted.
#  - This is an earlier version of the compartment prevalence figure. The
#    version used for Figure 7A is script 14.
#  - Four trailing lines of the original file referenced objects from the
#    spatial analysis (model_all, cor_all) and were removed; they belonged to
#    Methods 4.1.12 and errored here. See README.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggpubr)   # for compare_means() + stat_pvalue_manual()
})

# ============================================================
# CHUNK 1 — Safety checks
# Purpose: Ensure required objects and columns exist.
# ============================================================
stopifnot(exists("patient_tex"), exists("patient_tex_total"))

stopifnot(all(c("PatientID","Tissue","TexClass","prevalence") %in% colnames(patient_tex)))
stopifnot(all(c("PatientID","Tissue","prevalence") %in% colnames(patient_tex_total)))

# enforce Tissue order (edit if you prefer a different biological order)
tissue_levels <- c("Blood","Normal","Peritumour","Tumour")
patient_tex$Tissue       <- factor(patient_tex$Tissue, levels = tissue_levels)
patient_tex_total$Tissue <- factor(patient_tex_total$Tissue, levels = tissue_levels)

# ============================================================
# CHUNK 2 — Define all pairwise Tissue comparisons
# Purpose: Generate ALL combinations automatically (no manual list).
# ============================================================
tissues_present <- levels(droplevels(patient_tex_total$Tissue))
pairwise_list <- combn(tissues_present, 2, simplify = FALSE)

# helper: create y-positions for p-value brackets
make_y_positions <- function(df, step = 0.05, pad = 0.02) {
  ymax <- max(df$prevalence, na.rm = TRUE)
  base <- ymax + pad
  base + step * seq_along(pairwise_list)
}

# ============================================================
# CHUNK 3 — Pairwise tests for COMBINED Tex prevalence
# Purpose: Compute Wilcoxon tests for all Tissue combinations.
# ============================================================
pw_total <- compare_means(
  prevalence ~ Tissue,
  data = patient_tex_total,
  method = "wilcox.test",
  p.adjust.method = "BH"
)

# Add bracket positions
pw_total <- pw_total %>%
  mutate(
    y.position = make_y_positions(patient_tex_total, step = 0.05, pad = 0.02),
    label = paste0("p=", p.format)  # nice formatted p for plotting
  )

message("=== Pairwise tissue p-values (COMBINED Tex) ===")
print(pw_total %>% select(group1, group2, p, p.adj, p.format, p.signif))

# ============================================================
# CHUNK 4 — Plot A: COMBINED Tex prevalence + all pairwise p-values
# Purpose: Show patient-level Tex prevalence and annotate all p-values.
# ============================================================
pA <- ggplot(patient_tex_total, aes(x = Tissue, y = prevalence)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  stat_pvalue_manual(
    pw_total,
    label = "p.format",     # or use "p.signif"
    tip.length = 0.01,
    hide.ns = FALSE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0.02, 0.25))) +
  labs(
    title = "Patient-level prevalence of exhausted CD8 T cells",
    subtitle = "Combined Prolif-Tex and Non-Prolif-Tex (pairwise Wilcoxon; BH-adjusted table printed to console)",
    x = "Tissue compartment",
    y = "Fraction of CD8 cells"
  ) +
  theme_classic(base_size = 12)

print(pA)

# ============================================================
# CHUNK 5 — Pairwise tests for EACH TexClass separately
# Purpose: Compute pairwise Tissue p-values within Prolif-Tex and Non-Prolif-Tex.
# ============================================================
pw_byclass <- compare_means(
  prevalence ~ Tissue,
  data = patient_tex,
  group.by = "TexClass",
  method = "wilcox.test",
  p.adjust.method = "BH"
)

# Add bracket positions PER class (so brackets don't overlap between classes)
# We stagger them: Prolif-Tex higher, Non-Prolif-Tex slightly lower.
y_pos_base <- make_y_positions(patient_tex, step = 0.05, pad = 0.02)

pw_byclass <- pw_byclass %>%
  group_by(TexClass) %>%
  mutate(
    y.position = y_pos_base +
      ifelse(first(TexClass) == "Prolif-Tex", 0.10, 0.02),
    label = paste0("p=", p.format)
  ) %>%
  ungroup()

message("=== Pairwise tissue p-values (BY TexClass) ===")
print(pw_byclass %>% select(TexClass, group1, group2, p, p.adj, p.format, p.signif))

# ============================================================
# CHUNK 6 — Plot B: Prolif-Tex vs Non-Prolif-Tex + all pairwise p-values
# Purpose: Visualize each TexClass and annotate all tissue pairwise tests.
# ============================================================

pB <- ggplot(patient_tex, aes(x = Tissue, y = prevalence, fill = TexClass)) +
  
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
    pw_byclass,
    label = "p.format",
    tip.length = 0.01,
    hide.ns = FALSE
  ) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.28))
  ) +
  
  labs(
    title = "Patient-level prevalence of Prolif-Tex and Non-Prolif-Tex",
    subtitle = "Pairwise Wilcoxon across Tissue within each class (BH-adjusted)",
    x = "Tissue compartment",
    y = "Fraction of CD8 cells",
    fill = NULL,
    color = NULL
  ) +
  
  theme_classic(base_size = 12)

print(pB)

# ============================================================
# CHUNK 7 — Save p-value tables (optional)
# Purpose: Persist statistics for thesis reporting.
# ============================================================
write.csv(pw_total,   file = "pairwise_pvalues_Tissue_combinedTex.csv", row.names = FALSE)
write.csv(pw_byclass, file = "pairwise_pvalues_Tissue_byTexClass.csv",  row.names = FALSE)

message("✅ Saved p-value tables:")
message("  - pairwise_pvalues_Tissue_combinedTex.csv")
message("  - pairwise_pvalues_Tissue_byTexClass.csv")
