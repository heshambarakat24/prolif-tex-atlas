# ===========================================================
# 20_residual_gene_characterisation.R
# ===========================================================
# Thesis Methods section: 4.1.11 (retained-gene characterisation)
#
# Reads : Residual_pipeline_outputs/Residual_pipeline_obj.rds
#         Residual_pipeline_outputs/ResidualA_<label>_up_genes.csv
#         Residual_pipeline_outputs/ResidualB_<label>_up_genes.csv
# Writes: Residual_pipeline_outputs/Script3_residual_characterization_<label>/
#
# Subsets the object on CMP_B, reusing the labels assigned in script 18. It
# does not re-derive configurations.
#
# Characterises top residual genes within the exhausted-cell space: descriptive
# module scores (G1/S, Checkpoint, Mitotic, Inflammatory, Metabolic,
# Exhaustion - gene lists as in Methods 4.1.11, deposited in panels/), a linear
# model of each gene on configuration with mitotic score and patient baseline,
# and correlation structure among retained genes.
#
# gene_set_label near the top of the script selects the gene-set size. The
# thesis reports the top 100 residual upregulated genes per contrast, so set
# gene_set_label <- "top100" to reproduce the reported analysis.
#
# Minimum 20 cells per patient-configuration group; alpha 0.05.
#
# Paths are relative to the working directory. See README.md.
# ===========================================================

# ===========================================================
# SCRIPT 3 — Residual survivor characterization in exhausted cells
#
# Goal:
#   Characterize top-ranked residual genes within the exhausted-state
#   space only (Prolif-Tex vs NonProlif-Tex), without rerunning DE.
#
# Purpose:
#   This script does not claim discovery of a new transcriptional
#   program. Instead, it tests whether residual genes show coherent
#   module-like behavior or heterogeneous residual structure.
#
# Main questions:
#   1) How do top residual genes behave within exhausted cells?
#   2) Are they still associated with proliferative state after
#      accounting for mitotic score and patient baseline?
#   3) Do stringent residual genes co-vary as a coherent module?
#
# Inputs:
#   Residual_pipeline_outputs/
#     Residual_pipeline_obj.rds
#     ResidualA_topXX_up_genes.csv
#     ResidualB_topXX_up_genes.csv
#
# Required metadata in saved object:
#   CMP_B
#   PatientID
#
# Main outputs:
#   1) Patient-state summaries
#   2) Residual gene screening summaries
#   3) Heatmap-long tables for plotting
#   4) Correlation matrix for stringent residual survivors
#
# Notes:
#   - Heuristic labels are descriptive only
#   - Correlation / heatmap outputs are used to assess whether
#     residual genes form a coherent module
# ===========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(broom)
  library(tibble)
})

# -----------------------------------------------------------
# 0) SETTINGS
# -----------------------------------------------------------
pipeline_dir <- "Residual_pipeline_outputs"

# CHANGE THIS ONLY:
gene_set_label <- "top100"   # use "top20" or "top100"

outdir <- file.path(
  pipeline_dir,
  paste0("Script3_residual_characterization_", gene_set_label)
)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

min_cells_per_group <- 20
alpha_sig <- 0.05

# -----------------------------------------------------------
# 1) Load inputs
# -----------------------------------------------------------
obj_path <- file.path(pipeline_dir, "Residual_pipeline_obj.rds")
topA_path <- file.path(pipeline_dir, paste0("ResidualA_", gene_set_label, "_up_genes.csv"))
topB_path <- file.path(pipeline_dir, paste0("ResidualB_", gene_set_label, "_up_genes.csv"))

if (!file.exists(obj_path)) stop("Missing object: ", obj_path)
if (!file.exists(topA_path)) stop("Missing top gene file: ", topA_path)
if (!file.exists(topB_path)) stop("Missing top gene file: ", topB_path)

obj <- readRDS(obj_path)
topA <- read_csv(topA_path, show_col_types = FALSE)
topB <- read_csv(topB_path, show_col_types = FALSE)

DefaultAssay(obj) <- "RNA"

cat("\nRunning Script 3 with gene_set_label = ", gene_set_label, "\n", sep = "")

# -----------------------------------------------------------
# 2) Subset to exhausted-state interpretation space
# -----------------------------------------------------------
tex_obj <- subset(
  obj,
  subset = CMP_B %in% c("Prolif-Tex", "NonProlif-Tex")
)

tex_obj$CMP_B <- factor(tex_obj$CMP_B, levels = c("NonProlif-Tex", "Prolif-Tex"))

cat("\nCells in tex_obj:\n")
print(table(tex_obj$CMP_B, useNA = "ifany"))

required_meta <- c("CMP_B", "PatientID")
missing_meta <- setdiff(required_meta, colnames(tex_obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata in saved object: ", paste(missing_meta, collapse = ", "))
}

# -----------------------------------------------------------
# 3) Define descriptive biological axes
# -----------------------------------------------------------
g1s_genes <- c(
  "CCND1","CCND2","CCNE1",
  "E2F1","E2F2","E2F3",
  "CDK2",
  "CDC6","MCM2","MCM3","MCM5","MCM6"
)

checkpoint_genes <- c(
  "ATR","CHEK1","WEE1","CLSPN",
  "CDC25C","RAD17","HUS1","RFC4",
  "GADD45A","GADD45B","RRM2B"
)

mitotic_genes <- c(
  "CDK1","CCNB1","CCNB2",
  "FOXM1","PLK1",
  "AURKA","AURKB",
  "BUB1","BUB1B","CDC20"
)

inflammatory_genes <- c(
  "NFKB1","RELA","TNFAIP3","NFKBIA",
  "STAT3","IL6ST","SOCS3",
  "IRF1","IFITM1",
  "CXCL10","CCL5","BIRC3"
)

metabolic_genes <- c(
  "MYC","NPM1",
  "SLC2A1","HK2","LDHA",
  "NDUFA1","COX5A","ATP5F1B",
  "HADHA","ACADM","CPT1A"
)

exhaustion_genes <- c(
  "PDCD1","HAVCR2","LAG3","TIGIT",
  "TOX","ENTPD1","CXCL13",
  "NR4A1","NR4A2","NR4A3"
)

module_list_raw <- list(
  G1S = g1s_genes,
  Checkpoint = checkpoint_genes,
  Mitotic = mitotic_genes,
  Inflammatory = inflammatory_genes,
  Metabolic = metabolic_genes,
  Exhaustion = exhaustion_genes
)

filter_genes_present <- function(seurat_obj, gene_vector) {
  unique(gene_vector[gene_vector %in% rownames(seurat_obj)])
}

module_list <- lapply(module_list_raw, function(x) filter_genes_present(tex_obj, x))

module_presence_report <- bind_rows(lapply(names(module_list_raw), function(nm) {
  tibble(
    module = nm,
    input_gene = module_list_raw[[nm]],
    present = module_list_raw[[nm]] %in% module_list[[nm]]
  )
}))

write_csv(module_presence_report, file.path(outdir, "module_gene_presence_report.csv"))

# -----------------------------------------------------------
# 4) Recompute module scores
# -----------------------------------------------------------
tex_obj <- ScaleData(tex_obj, verbose = FALSE)

for (nm in names(module_list)) {
  if (length(module_list[[nm]]) >= 2) {
    tex_obj <- AddModuleScore(
      object = tex_obj,
      features = list(module_list[[nm]]),
      name = paste0(nm, "_Score")
    )
  } else {
    warning("Module ", nm, " has fewer than 2 detected genes and was not scored.")
  }
}

required_score_meta <- c(
  "G1S_Score1", "Checkpoint_Score1", "Mitotic_Score1",
  "Inflammatory_Score1", "Metabolic_Score1", "Exhaustion_Score1"
)

missing_scores <- setdiff(required_score_meta, colnames(tex_obj@meta.data))
if (length(missing_scores) > 0) {
  stop("Missing required module score columns: ", paste(missing_scores, collapse = ", "))
}

# -----------------------------------------------------------
# 5) Build exhausted-only metadata table
# -----------------------------------------------------------
meta_exhausted <- tex_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(
    State = factor(ifelse(CMP_B == "Prolif-Tex", "Prolif", "NonProlif"),
                   levels = c("NonProlif", "Prolif")),
    G1S_score = G1S_Score1,
    Checkpoint_score = Checkpoint_Score1,
    Mitotic_score = Mitotic_Score1,
    Inflammatory_score = Inflammatory_Score1,
    Metabolic_score = Metabolic_Score1,
    Exhaustion_score = Exhaustion_Score1
  ) %>%
  select(
    cell_id, PatientID, CMP_B, State,
    G1S_score, Checkpoint_score, Mitotic_score,
    Inflammatory_score, Metabolic_score, Exhaustion_score
  )

module_patient_state_summary <- meta_exhausted %>%
  group_by(PatientID, State) %>%
  summarise(
    G1S_score = mean(G1S_score, na.rm = TRUE),
    Checkpoint_score = mean(Checkpoint_score, na.rm = TRUE),
    Mitotic_score = mean(Mitotic_score, na.rm = TRUE),
    Inflammatory_score = mean(Inflammatory_score, na.rm = TRUE),
    Metabolic_score = mean(Metabolic_score, na.rm = TRUE),
    Exhaustion_score = mean(Exhaustion_score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

write_csv(module_patient_state_summary, file.path(outdir, "module_patient_state_summary.csv"))

# -----------------------------------------------------------
# 6) Helper functions
# -----------------------------------------------------------
safe_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

safe_wilcox_paired <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
}

zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# Descriptive only, not inferential
classify_residual_behavior <- function(cor_mitotic, cor_inflammatory, beta_adj, p_adj,
                                       alpha = alpha_sig) {
  if (is.na(p_adj) || p_adj >= alpha) return("Other/Weak")
  if (!is.na(cor_mitotic) && abs(cor_mitotic) < 0.10 &&
      !is.na(cor_inflammatory) && cor_inflammatory > 0.15) {
    return("Mitotic-independent contextual")
  }
  if (!is.na(cor_mitotic) && cor_mitotic > 0.15 &&
      !is.na(cor_inflammatory) && cor_inflammatory > 0.15) {
    return("Hybrid/coexistence")
  }
  if (!is.na(cor_mitotic) && cor_mitotic > 0.15) {
    return("Mitotic-associated")
  }
  return("Other/Weak")
}

# -----------------------------------------------------------
# 7) Core screening function
# -----------------------------------------------------------
screen_top_genes_robust <- function(top_gene_df,
                                    residual_name,
                                    min_cells_per_group = 20) {
  
  residual_dir <- file.path(outdir, residual_name)
  dir.create(residual_dir, showWarnings = FALSE, recursive = TRUE)
  
  top_genes <- top_gene_df %>%
    filter(gene %in% rownames(tex_obj))
  
  if (nrow(top_genes) == 0) {
    stop("No genes from ", residual_name, " top list were found in tex_obj.")
  }
  
  expr_long <- FetchData(tex_obj, vars = top_genes$gene) %>%
    rownames_to_column("cell_id") %>%
    pivot_longer(
      cols = -cell_id,
      names_to = "gene",
      values_to = "gene_expr"
    )
  
  df_cell_long <- meta_exhausted %>%
    left_join(expr_long, by = "cell_id") %>%
    mutate(detected = gene_expr > 0)
  
  write_csv(
    df_cell_long,
    file.path(residual_dir, paste0(residual_name, "_top_gene_cell_level_long.csv"))
  )
  
  patient_state_long <- df_cell_long %>%
    group_by(gene, PatientID, State) %>%
    summarise(
      mean_expr = mean(gene_expr, na.rm = TRUE),
      median_expr = median(gene_expr, na.rm = TRUE),
      pct_expr = mean(detected, na.rm = TRUE),
      G1S_score = mean(G1S_score, na.rm = TRUE),
      Checkpoint_score = mean(Checkpoint_score, na.rm = TRUE),
      Mitotic_score = mean(Mitotic_score, na.rm = TRUE),
      Inflammatory_score = mean(Inflammatory_score, na.rm = TRUE),
      Metabolic_score = mean(Metabolic_score, na.rm = TRUE),
      Exhaustion_score = mean(Exhaustion_score, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= min_cells_per_group)
  
  patient_state_paired <- patient_state_long %>%
    group_by(gene, PatientID) %>%
    filter(n_distinct(State) == 2) %>%
    ungroup()
  
  write_csv(
    patient_state_paired,
    file.path(residual_dir, paste0(residual_name, "_top_gene_patient_state_means.csv"))
  )
  
  pairwise_diffs <- patient_state_paired %>%
    select(
      gene, PatientID, State,
      mean_expr, median_expr, pct_expr,
      G1S_score, Checkpoint_score, Mitotic_score,
      Inflammatory_score, Metabolic_score, Exhaustion_score
    ) %>%
    pivot_wider(
      names_from = State,
      values_from = c(
        mean_expr, median_expr, pct_expr,
        G1S_score, Checkpoint_score, Mitotic_score,
        Inflammatory_score, Metabolic_score, Exhaustion_score
      )
    ) %>%
    mutate(
      diff_mean_expr = mean_expr_Prolif - mean_expr_NonProlif,
      diff_median_expr = median_expr_Prolif - median_expr_NonProlif,
      diff_pct_expr = pct_expr_Prolif - pct_expr_NonProlif,
      diff_G1S = G1S_score_Prolif - G1S_score_NonProlif,
      diff_Checkpoint = Checkpoint_score_Prolif - Checkpoint_score_NonProlif,
      diff_Mitotic = Mitotic_score_Prolif - Mitotic_score_NonProlif,
      diff_Inflammatory = Inflammatory_score_Prolif - Inflammatory_score_NonProlif,
      diff_Metabolic = Metabolic_score_Prolif - Metabolic_score_NonProlif,
      diff_Exhaustion = Exhaustion_score_Prolif - Exhaustion_score_NonProlif
    )
  
  write_csv(
    pairwise_diffs,
    file.path(residual_dir, paste0(residual_name, "_top_gene_pairwise_differences.csv"))
  )
  
  gene_summary <- map_dfr(top_genes$gene, function(g) {
    
    top_row <- top_genes %>% filter(gene == g)
    
    df_ps <- patient_state_paired %>% filter(gene == g)
    df_diff <- pairwise_diffs %>% filter(gene == g)
    
    n_paired_patients <- n_distinct(df_ps$PatientID)
    median_diff <- median(df_diff$diff_mean_expr, na.rm = TRUE)
    paired_p_value <- safe_wilcox_paired(df_diff$diff_mean_expr)
    
    frac_expr_prolif <- mean(df_ps$pct_expr[df_ps$State == "Prolif"], na.rm = TRUE)
    frac_expr_nonprolif <- mean(df_ps$pct_expr[df_ps$State == "NonProlif"], na.rm = TRUE)
    frac_expr_diff <- frac_expr_prolif - frac_expr_nonprolif
    
    cor_G1S <- safe_spearman(df_ps$mean_expr, df_ps$G1S_score)
    cor_checkpoint <- safe_spearman(df_ps$mean_expr, df_ps$Checkpoint_score)
    cor_mitotic <- safe_spearman(df_ps$mean_expr, df_ps$Mitotic_score)
    cor_metabolic <- safe_spearman(df_ps$mean_expr, df_ps$Metabolic_score)
    cor_inflammatory <- safe_spearman(df_ps$mean_expr, df_ps$Inflammatory_score)
    cor_exhaustion <- safe_spearman(df_ps$mean_expr, df_ps$Exhaustion_score)
    
    beta_adj_mitotic <- NA_real_
    p_adj_mitotic <- NA_real_
    
    if (n_distinct(df_ps$PatientID) >= 3 && nrow(df_ps) >= 6) {
      fit <- tryCatch(
        lm(mean_expr ~ State + Mitotic_score + PatientID, data = df_ps),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        fit_tidy <- broom::tidy(fit)
        row_state <- fit_tidy %>% filter(term == "StateProlif")
        if (nrow(row_state) == 1) {
          beta_adj_mitotic <- row_state$estimate
          p_adj_mitotic <- row_state$p.value
        }
      }
    }
    
    heuristic_class <- classify_residual_behavior(
      cor_mitotic = cor_mitotic,
      cor_inflammatory = cor_inflammatory,
      beta_adj = beta_adj_mitotic,
      p_adj = p_adj_mitotic
    )
    
    tibble(
      gene = g,
      rank_FDR = top_row$rank_FDR,
      logFC = top_row$logFC,
      FDR = top_row$FDR,
      n_paired_patients = n_paired_patients,
      median_diff = median_diff,
      paired_p_value = paired_p_value,
      frac_expr_prolif = frac_expr_prolif,
      frac_expr_nonprolif = frac_expr_nonprolif,
      frac_expr_diff = frac_expr_diff,
      cor_G1S = cor_G1S,
      cor_checkpoint = cor_checkpoint,
      cor_mitotic = cor_mitotic,
      cor_metabolic = cor_metabolic,
      cor_inflammatory = cor_inflammatory,
      cor_exhaustion = cor_exhaustion,
      beta_adj_mitotic = beta_adj_mitotic,
      p_adj_mitotic = p_adj_mitotic,
      heuristic_class = heuristic_class
    )
  }) %>%
    arrange(rank_FDR)
  
  heatmap_metrics <- c(
    "median_diff",
    "frac_expr_diff",
    "cor_G1S",
    "cor_checkpoint",
    "cor_mitotic",
    "cor_metabolic",
    "cor_inflammatory",
    "cor_exhaustion",
    "beta_adj_mitotic"
  )
  
  heatmap_long <- gene_summary %>%
    select(gene, all_of(heatmap_metrics)) %>%
    pivot_longer(
      cols = all_of(heatmap_metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(metric) %>%
    mutate(value_z = zscore_safe(value)) %>%
    ungroup()
  
  write_csv(
    gene_summary,
    file.path(residual_dir, paste0(residual_name, "_top_gene_screening_summary.csv"))
  )
  
  write_csv(
    heatmap_long,
    file.path(residual_dir, paste0(residual_name, "_top_gene_screening_heatmap_long.csv"))
  )
  
  return(list(
    summary = gene_summary,
    heatmap_long = heatmap_long,
    patient_state = patient_state_paired,
    pairwise_diffs = pairwise_diffs,
    cell_long = df_cell_long
  ))
}

# -----------------------------------------------------------
# 8) Run Residual A characterization
# -----------------------------------------------------------
cat("\nRunning characterization for ResidualA ", gene_set_label, " genes...\n", sep = "")

resA <- screen_top_genes_robust(
  top_gene_df = topA,
  residual_name = "ResidualA",
  min_cells_per_group = min_cells_per_group
)

# -----------------------------------------------------------
# 9) Run Residual B characterization
# -----------------------------------------------------------
cat("\nRunning characterization for ResidualB ", gene_set_label, " genes...\n", sep = "")

resB <- screen_top_genes_robust(
  top_gene_df = topB,
  residual_name = "ResidualB",
  min_cells_per_group = min_cells_per_group
)

# -----------------------------------------------------------
# 10) Save exhausted object with module scores
# -----------------------------------------------------------
saveRDS(tex_obj, file = file.path(outdir, "tex_obj_exhausted_with_module_scores.rds"))

# -----------------------------------------------------------
# 11) Discharge analysis: correlation structure of stringent survivors
# -----------------------------------------------------------
make_survivor_correlation <- function(patient_df, summary_df, residual_name,
                                      paired_cut = 0.01, adj_cut = 0.05) {
  
  residual_dir <- file.path(outdir, residual_name)
  
  robust_df <- summary_df %>%
    filter(!is.na(paired_p_value), !is.na(p_adj_mitotic)) %>%
    filter(paired_p_value < paired_cut, p_adj_mitotic < adj_cut)
  
  write_csv(
    robust_df,
    file.path(residual_dir, paste0(residual_name, "_stringent_survivors.csv"))
  )
  
  if (nrow(robust_df) < 2) {
    message("Skipping correlation matrix for ", residual_name, ": fewer than 2 stringent survivors.")
    return(invisible(NULL))
  }
  
  cor_wide <- patient_df %>%
    filter(gene %in% robust_df$gene) %>%
    mutate(patient_state = paste(PatientID, State, sep = "___")) %>%
    select(patient_state, gene, mean_expr) %>%
    distinct() %>%
    pivot_wider(names_from = gene, values_from = mean_expr)
  
  cor_input <- cor_wide %>% select(-patient_state)
  
  if (ncol(cor_input) < 2) {
    message("Skipping correlation matrix for ", residual_name, ": insufficient survivor columns.")
    return(invisible(NULL))
  }
  
  cor_mat <- suppressWarnings(cor(cor_input, use = "pairwise.complete.obs"))
  
  write_csv(
    as.data.frame(cor_mat) %>% tibble::rownames_to_column("gene"),
    file.path(residual_dir, paste0(residual_name, "_stringent_survivor_correlation_matrix.csv"))
  )
  
  invisible(cor_mat)
}

make_survivor_correlation(
  patient_df = resA$patient_state,
  summary_df = resA$summary,
  residual_name = "ResidualA"
)

make_survivor_correlation(
  patient_df = resB$patient_state,
  summary_df = resB$summary,
  residual_name = "ResidualB"
)

cat("\nScript 3 complete.\n")
cat("Outputs written to:\n", normalizePath(outdir), "\n")
cat("\nInterpretation reminder: heuristic labels are descriptive only; the key question is whether residual genes form coherent structure.\n")