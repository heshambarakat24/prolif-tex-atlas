# ==============================================================================
# 22_spatial_pressure_subtypes.R
# ==============================================================================
# Thesis Methods section: 4.1.12 (spatial analysis of microenvironmental
#                                 pressure subtypes)
#
# Reads : the eight CK-prefixed 10x Visium samples of GSE301720, each as a
#         directory of Space Ranger output under
#         <DATA_ROOT>/GSE301720_Seurat_Ready/. Not included in this repository;
#         see README.md for the accession.
# Writes: Spatial_pressure_subtypes_2/ - per-sample Spearman correlations,
#         per-sample linear-model coefficients, cross-sample Wilcoxon tests
#         and the associated figures.
#
# What it does:
#  - Builds one Seurat object per sample with the Spatial assay, retaining
#    spots with more than 10 counts and more than 10 detected features, and
#    normalises each sample independently with SCTransform.
#  - Defines T-cell-enriched regions as the upper 30% of the per-sample T-cell
#    score distribution, and drops samples with fewer than 30 such spots.
#  - Scores proliferation and three pressure subtypes (IFN/STAT1, NF-kB,
#    cytokine/JAK-STAT3) with AddModuleScore. Gene lists are those reported in
#    Methods 4.1.12 and deposited in panels/.
#  - Within T-cell-enriched spots of each sample, computes Spearman
#    correlations between the proliferation score and each pressure score, and
#    fits Proliferation ~ IFN + NFkB + Cytokine.
#  - Tests the cross-sample distributions of correlations and coefficients
#    against zero with one-sample Wilcoxon signed-rank tests.
#
# Scope: eight samples, one cancer type, p-values are nominal and uncorrected.
# All analyses are at Visium spot level; T-cell-enriched spots are spatial
# neighbourhoods, not purified T-cell measurements.
# ==============================================================================

# ==============================================================================
# SPATIAL SUBTYPE PRESSURE ANALYSIS
# ------------------------------------------------------------------------------
# Purpose:
#   Split "pressure" into biologically grounded subtypes and test which one is
#   most consistently associated with reduced proliferative signal in
#   T-cell-enriched regions.
#
# Main questions:
#   1) Which subtype (IFN, NFkB, Cytokine-feedback) shows the strongest inverse
#      association with Prolif?
#   2) In a multivariable model, which subtype has the strongest negative effect?
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(future)
})

# ------------------------------------------------------------------------------
# 0) SETTINGS
# ------------------------------------------------------------------------------

plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

DATA_ROOT <- Sys.getenv("DATA_ROOT", unset = "data")
base_dir  <- file.path(DATA_ROOT, "GSE301720_Seurat_Ready")
out_dir  <- file.path(base_dir, "Spatial_pressure_subtypes_2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
samples <- samples[grepl("^CK", samples)]

min_count_threshold <- 10
min_feature_threshold <- 10
tcell_quantile <- 0.70
min_spots_for_analysis <- 30

# ------------------------------------------------------------------------------
# 1) SIGNATURES
# ------------------------------------------------------------------------------

# T-cell-enriched neighborhood proxy
tcell_sig <- c("CD3D", "CD3E", "CD3G", "TRAC", "CD2", "CD8A", "CD8B")

# Proliferative competence
prolif_sig <- c(
  "MKI67", "PCNA",
  "MCM2", "MCM5",
  "TOP2A",
  "BIRC5", "UBE2C", "AURKB", "CDC20", "CCNB1"
)

# Pressure subtypes
ifn_sig <- c("CXCL9", "CXCL10", "STAT1", "IRF1", "IFITM1", "ISG15")

nfkb_sig <- c("NFKBIA", "TNFAIP3", "BIRC3", "RELA", "NFKB1")

cytokine_sig <- c("SOCS3", "SOCS1", "IL6ST", "JAK1", "STAT3")

# ------------------------------------------------------------------------------
# 2) HELPERS
# ------------------------------------------------------------------------------

score_module <- function(obj, genes, name) {
  genes_use <- intersect(genes, rownames(obj))
  if (length(genes_use) >= 2) {
    obj <- AddModuleScore(obj, features = list(genes_use), name = name)
  } else if (length(genes_use) == 1) {
    # single gene fallback
    obj[[paste0(name, "1")]] <- FetchData(obj, vars = genes_use)[, 1]
    warning("Module ", name, " has 1 detected gene; using raw expression of ", genes_use)
  } else {
    obj[[paste0(name, "1")]] <- 0
    warning("Module ", name, " had 0 detected genes; filled with 0.")
  }
  obj
}

safe_cor_test <- function(x, y, min_n = 10) {
  ok <- complete.cases(x, y)
  if (sum(ok) < min_n) {
    return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

safe_wilcox_zero <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", signif(p, 3))
}

fit_multivariable <- function(df) {
  d <- df %>%
    select(Prolif1, IFN1, NFkB1, Cytokine1) %>%
    filter(complete.cases(.))
  
  if (nrow(d) < 20) {
    return(tibble(
      beta_ifn = NA_real_, p_ifn = NA_real_,
      beta_nfkb = NA_real_, p_nfkb = NA_real_,
      beta_cytokine = NA_real_, p_cytokine = NA_real_,
      model_r2 = NA_real_, n_model = nrow(d)
    ))
  }
  
  fit <- lm(Prolif1 ~ IFN1 + NFkB1 + Cytokine1, data = d)
  sm <- summary(fit)
  coefs <- sm$coefficients
  
  get_coef <- function(term, col = 1) {
    if (!(term %in% rownames(coefs))) return(NA_real_)
    coefs[term, col]
  }
  
  tibble(
    beta_ifn = get_coef("IFN1", 1),
    p_ifn = get_coef("IFN1", 4),
    beta_nfkb = get_coef("NFkB1", 1),
    p_nfkb = get_coef("NFkB1", 4),
    beta_cytokine = get_coef("Cytokine1", 1),
    p_cytokine = get_coef("Cytokine1", 4),
    model_r2 = sm$r.squared,
    n_model = nrow(d)
  )
}

# ------------------------------------------------------------------------------
# 3) PROCESS SAMPLES
# ------------------------------------------------------------------------------

processed_objects <- list()
spot_tables <- list()
cor_results <- list()
model_results <- list()

for (s in samples) {
  cat("\n====================================================\n")
  cat("Processing:", s, "\n")
  cat("====================================================\n")
  
  data_dir <- file.path(base_dir, s)
  spatial_dir <- file.path(data_dir, "spatial")
  
  if (!dir.exists(spatial_dir)) {
    message("Skipping ", s, " (no spatial folder)")
    next
  }
  
  counts <- Read10X(data.dir = data_dir)
  obj <- CreateSeuratObject(counts = counts, assay = "Spatial", project = s)
  
  img <- Read10X_Image(spatial_dir)
  DefaultAssay(img) <- "Spatial"
  obj[["slice1"]] <- img
  
  obj <- subset(
    obj,
    subset = nCount_Spatial > min_count_threshold & nFeature_Spatial > min_feature_threshold
  )
  
  if (ncol(obj) == 0) {
    message("Skipping ", s, " (no spots remained after QC)")
    next
  }
  
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
  
  # Score modules
  obj <- score_module(obj, tcell_sig, "TcellNeighborhood")
  obj <- score_module(obj, prolif_sig, "Prolif")
  obj <- score_module(obj, ifn_sig, "IFN")
  obj <- score_module(obj, nfkb_sig, "NFkB")
  obj <- score_module(obj, cytokine_sig, "Cytokine")
  
  thresh <- as.numeric(quantile(obj$TcellNeighborhood1, tcell_quantile, na.rm = TRUE))
  obj$Tcell_region <- ifelse(obj$TcellNeighborhood1 >= thresh, "High", "Low")
  
  df <- FetchData(
    obj,
    vars = c(
      "TcellNeighborhood1",
      "Prolif1",
      "IFN1",
      "NFkB1",
      "Cytokine1",
      "Tcell_region"
    )
  ) %>%
    rownames_to_column("spot_id") %>%
    mutate(sample = s)
  
  df_tcell <- df %>% filter(Tcell_region == "High")
  
  if (nrow(df_tcell) < min_spots_for_analysis) {
    message("Skipping ", s, " (too few T-cell-enriched spots)")
    next
  }
  
  # Per-sample correlations
  res_ifn <- safe_cor_test(df_tcell$Prolif1, df_tcell$IFN1, min_n = 20)
  res_nfkb <- safe_cor_test(df_tcell$Prolif1, df_tcell$NFkB1, min_n = 20)
  res_cyt <- safe_cor_test(df_tcell$Prolif1, df_tcell$Cytokine1, min_n = 20)
  
  cor_tbl <- tibble(
    sample = s,
    n_spots_total = nrow(df),
    n_spots_tcellhigh = nrow(df_tcell),
    
    rho_prolif_ifn = res_ifn$rho,
    p_prolif_ifn = res_ifn$p,
    
    rho_prolif_nfkb = res_nfkb$rho,
    p_prolif_nfkb = res_nfkb$p,
    
    rho_prolif_cytokine = res_cyt$rho,
    p_prolif_cytokine = res_cyt$p
  )
  
  model_tbl <- fit_multivariable(df_tcell) %>%
    mutate(sample = s)
  
  processed_objects[[s]] <- obj
  spot_tables[[s]] <- df_tcell
  cor_results[[s]] <- cor_tbl
  model_results[[s]] <- model_tbl
  
  write_csv(df_tcell, file.path(out_dir, paste0(s, "_tcellhigh_spot_table.csv")))
  write_csv(cor_tbl, file.path(out_dir, paste0(s, "_subtype_correlations.csv")))
  write_csv(model_tbl, file.path(out_dir, paste0(s, "_multivariable_model.csv")))
}

# ------------------------------------------------------------------------------
# 4) COMBINE RESULTS
# ------------------------------------------------------------------------------

cor_all <- bind_rows(cor_results)
model_all <- bind_rows(model_results)

write_csv(cor_all, file.path(out_dir, "ALL_samples_subtype_correlations.csv"))
write_csv(model_all, file.path(out_dir, "ALL_samples_multivariable_models.csv"))

# ------------------------------------------------------------------------------
# 5) CROSS-SAMPLE TESTS
# ------------------------------------------------------------------------------

p_ifn <- safe_wilcox_zero(cor_all$rho_prolif_ifn)
p_nfkb <- safe_wilcox_zero(cor_all$rho_prolif_nfkb)
p_cyt <- safe_wilcox_zero(cor_all$rho_prolif_cytokine)

p_beta_ifn <- safe_wilcox_zero(model_all$beta_ifn)
p_beta_nfkb <- safe_wilcox_zero(model_all$beta_nfkb)
p_beta_cyt <- safe_wilcox_zero(model_all$beta_cytokine)

cross_sample_stats <- tibble(
  metric = c(
    "Prolif vs IFN",
    "Prolif vs NFkB",
    "Prolif vs Cytokine",
    "Beta IFN (multivariable)",
    "Beta NFkB (multivariable)",
    "Beta Cytokine (multivariable)"
  ),
  median_effect = c(
    median(cor_all$rho_prolif_ifn, na.rm = TRUE),
    median(cor_all$rho_prolif_nfkb, na.rm = TRUE),
    median(cor_all$rho_prolif_cytokine, na.rm = TRUE),
    median(model_all$beta_ifn, na.rm = TRUE),
    median(model_all$beta_nfkb, na.rm = TRUE),
    median(model_all$beta_cytokine, na.rm = TRUE)
  ),
  wilcox_p_value = c(
    p_ifn, p_nfkb, p_cyt,
    p_beta_ifn, p_beta_nfkb, p_beta_cyt
  )
)

write_csv(cross_sample_stats, file.path(out_dir, "CROSS_SAMPLE_STATS.csv"))

# ------------------------------------------------------------------------------
# 6) MAIN PLOT: CORRELATION BETWEEN PRESSURE SUBTYPES AND PROLIFERATION
# ------------------------------------------------------------------------------

cor_metric_labels <- c(
  rho_prolif_cytokine = "Cytokine/JAK–STAT3",
  rho_prolif_ifn      = "IFN/STAT1",
  rho_prolif_nfkb     = "NF-κB"
)

cor_metric_levels <- unname(cor_metric_labels)

plot_cor <- cor_all %>%
  select(sample, all_of(names(cor_metric_labels))) %>%
  pivot_longer(
    cols = -sample,
    names_to = "metric_raw",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric_raw, !!!cor_metric_labels),
    metric = factor(metric, levels = cor_metric_levels)
  )

p_values_cor <- c(
  "Cytokine/JAK–STAT3" = p_cyt,
  "IFN/STAT1"          = p_ifn,
  "NF-κB"              = p_nfkb
)

label_df <- tibble(
  metric = factor(names(p_values_cor), levels = cor_metric_levels),
  label = vapply(p_values_cor, fmt_p, character(1)),
  y = max(plot_cor$value, na.rm = TRUE) + 0.05
)

p_main <- ggplot(plot_cor, aes(x = metric, y = value)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.9) +
  geom_text(
    data = label_df,
    aes(x = metric, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Microenvironmental pressure subtypes show axis-specific inverse associations\nwith proliferative signal",
    x = NULL,
    y = "Spearman rho: proliferation score vs pressure subtype score"
  ) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "PRESSURE_SUBTYPE_COMPARISON.png"),
  p_main,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(out_dir, "PRESSURE_SUBTYPE_COMPARISON.pdf"),
  p_main,
  width = 9,
  height = 6
)


# ------------------------------------------------------------------------------
# 7) SECONDARY PLOT: MULTIVARIABLE MODEL BETAS
# ------------------------------------------------------------------------------

beta_metric_labels <- c(
  beta_cytokine = "Cytokine/JAK–STAT3",
  beta_ifn      = "IFN/STAT1",
  beta_nfkb     = "NF-κB"
)

beta_metric_levels <- unname(beta_metric_labels)

plot_beta <- model_all %>%
  select(sample, all_of(names(beta_metric_labels))) %>%
  pivot_longer(
    cols = -sample,
    names_to = "metric_raw",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric_raw, !!!beta_metric_labels),
    metric = factor(metric, levels = beta_metric_levels)
  )

p_values_beta <- c(
  "Cytokine/JAK–STAT3" = p_beta_cyt,
  "IFN/STAT1"          = p_beta_ifn,
  "NF-κB"              = p_beta_nfkb
)

label_beta <- tibble(
  metric = factor(names(p_values_beta), levels = beta_metric_levels),
  label = vapply(p_values_beta, fmt_p, character(1)),
  y = max(plot_beta$value, na.rm = TRUE) + 0.05
)

p_beta <- ggplot(plot_beta, aes(x = metric, y = value)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.2, alpha = 0.9) +
  geom_text(
    data = label_beta,
    aes(x = metric, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Multivariable contributions of pressure subtypes to proliferative signal",
    x = NULL,
    y = "Linear model beta"
  ) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "PRESSURE_SUBTYPE_MULTIVARIABLE_BETAS.png"),
  p_beta,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(out_dir, "PRESSURE_SUBTYPE_MULTIVARIABLE_BETAS.pdf"),
  p_beta,
  width = 9,
  height = 6
)

# ------------------------------------------------------------------------------
# 9) INTERPRETATION NOTES
# ------------------------------------------------------------------------------

interpretation_tbl <- tibble(
  conclusion = c(
    "This analysis remains contextual, not cell-intrinsic.",
    "Spots are mixed-cell; T-cell-enriched does not mean pure T-cell.",
    "Pressure was split into predefined biological subtypes: IFN, NFkB, and cytokine-feedback.",
    "The main goal is to identify which subtype most consistently opposes proliferative signal.",
    "The main interpretation should be based on whichever pressure subtype shows the most consistent negative association with proliferative signal across samples."
  )
)

write_csv(interpretation_tbl, file.path(out_dir, "INTERPRETATION_NOTES.csv"))

# ------------------------------------------------------------------------------
# 10) SUMMARY
# ------------------------------------------------------------------------------

cat("\n====================================================\n")
cat("Pressure subtype analysis completed.\n")
cat("====================================================\n")
cat("Outputs written to:\n", normalizePath(out_dir), "\n")
cat("\nCross-sample stats:\n")
print(cross_sample_stats)
cat("\nInterpretation rule:\n")
cat(cat("  Strongest support = the same pressure subtype is consistently negative in both correlation and multivariable analyses.\n"))
cat("  Weak support = mixed directions or no subtype stands out.\n")