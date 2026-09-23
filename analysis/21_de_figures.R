# ===========================================================
# 21_de_figures.R
# ===========================================================
# Thesis Methods section: 4.1.10 and 4.1.11 (figures)
# Produces: the direct and residual contrast figures - Figures 9 and 10
#
# Reads : outputs of scripts 18, 19 and 20 under Residual_pipeline_outputs/
# Writes: figure files under Residual_pipeline_outputs/
#
# Builds volcano plots for both direct and both residual contrasts, residual
# heterogeneity maps, and retained-gene heatmaps and correlation panels. It
# derives no groups and runs no differential expression; every input is an
# output of the three preceding scripts.
#
# Paths are relative to the working directory. See README.md.
# ===========================================================

# ===========================================================
# SCRIPT 4 — Final balanced storyline figures for thesis / paper
#
# Goal:
#   Generate final curated figures from Scripts 1–3 using the
#   same survivor-level logic for both Residual A and Residual B.
#
# Narrative:
#   direct contrasts -> dominant programs
#   residual contrasts -> remaining signal after filtering
#   survivor maps -> position of robust residual genes in broader space
#   survivor heatmaps/correlations -> test whether survivors form a module
#   conclusion -> limited, heterogeneous, non-modular residual structure
#
# Main outputs:
#   1) Direct contrast A volcano
#   2) Direct contrast B volcano
#   3) Residual A volcano
#   4) Residual B volcano
#   5) Residual A heterogeneity map
#   6) Residual B heterogeneity map
#   7) Residual A top-gene heterogeneity heatmap
#   8) Residual B top-gene heterogeneity heatmap
#   9) Residual A stringent-survivor heatmap
#   10) Residual B stringent-survivor heatmap
#   11) Residual A survivor correlation matrix
#   12) Residual B survivor correlation matrix
#   13) Optional survivor example panels for both residuals
#   14) Combined summary table for captioning
#
# Key point:
#   No residual survivor genes are manually selected in this script.
#   Stringent survivors are read from Script 3 outputs:
#     ResidualA_stringent_survivors.csv
#     ResidualB_stringent_survivors.csv
#
# Update:
#   All figures are exported as both PDF and PNG.
# ===========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(forcats)
  library(patchwork)
  library(scales)
  library(tibble)
})

# -----------------------------------------------------------
# 0) USER SETTINGS
# -----------------------------------------------------------

pipeline_dir <- "Residual_pipeline_outputs"

# CHANGE THIS ONLY:
gene_set_label <- "top100"   # use "top20" or "top100"

fig_dir <- file.path(
  pipeline_dir,
  paste0("Script4_storyline_figures_balanced_", gene_set_label)
)

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

char_dir <- file.path(
  pipeline_dir,
  paste0("Script3_residual_characterization_", gene_set_label)
)

resA_dir <- file.path(char_dir, "ResidualA")
resB_dir <- file.path(char_dir, "ResidualB")

state_colors <- c(
  "NonProlif" = "#1f77b4",
  "Prolif" = "#d62728"
)

# -----------------------------------------------------------
# 1) LOAD INPUT FILES
# -----------------------------------------------------------

# Script 1 outputs
directA_path <- file.path(
  pipeline_dir,
  "DirectContrastA_ProlifTex_vs_ProPlus_nonExh_full_edgeR.csv"
)

directB_path <- file.path(
  pipeline_dir,
  "DirectContrastB_ProlifTex_vs_NonProlifTex_full_edgeR.csv"
)

# Script 2 outputs
residA_path <- file.path(
  pipeline_dir,
  "ResidualA_ProlifTex_vs_ProPlus_nonExh_NO_EXHAUSTION_GENES.csv"
)

residB_path <- file.path(
  pipeline_dir,
  "ResidualB_ProlifTex_vs_NonProlifTex_NO_STANDARD_CELLCYCLE_HRG.csv"
)

# Script 3 Residual A outputs
A_summary_path <- file.path(resA_dir, "ResidualA_top_gene_screening_summary.csv")
A_heatmap_path <- file.path(resA_dir, "ResidualA_top_gene_screening_heatmap_long.csv")
A_patient_path <- file.path(resA_dir, "ResidualA_top_gene_patient_state_means.csv")
A_cell_path <- file.path(resA_dir, "ResidualA_top_gene_cell_level_long.csv")
A_corr_path <- file.path(resA_dir, "ResidualA_stringent_survivor_correlation_matrix.csv")
A_survivor_path <- file.path(resA_dir, "ResidualA_stringent_survivors.csv")

# Script 3 Residual B outputs
B_summary_path <- file.path(resB_dir, "ResidualB_top_gene_screening_summary.csv")
B_heatmap_path <- file.path(resB_dir, "ResidualB_top_gene_screening_heatmap_long.csv")
B_patient_path <- file.path(resB_dir, "ResidualB_top_gene_patient_state_means.csv")
B_cell_path <- file.path(resB_dir, "ResidualB_top_gene_cell_level_long.csv")
B_corr_path <- file.path(resB_dir, "ResidualB_stringent_survivor_correlation_matrix.csv")
B_survivor_path <- file.path(resB_dir, "ResidualB_stringent_survivors.csv")

needed_files <- c(
  directA_path, directB_path,
  residA_path, residB_path,
  A_summary_path, A_heatmap_path, A_patient_path, A_cell_path,
  B_summary_path, B_heatmap_path, B_patient_path, B_cell_path
)

missing_files <- needed_files[!file.exists(needed_files)]

if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

directA <- read_csv(directA_path, show_col_types = FALSE)
directB <- read_csv(directB_path, show_col_types = FALSE)

residA <- read_csv(residA_path, show_col_types = FALSE)
residB <- read_csv(residB_path, show_col_types = FALSE)

A_summary <- read_csv(A_summary_path, show_col_types = FALSE) %>%
  arrange(rank_FDR)

A_heatmap <- read_csv(A_heatmap_path, show_col_types = FALSE)
A_patient <- read_csv(A_patient_path, show_col_types = FALSE)
A_cell <- read_csv(A_cell_path, show_col_types = FALSE)

B_summary <- read_csv(B_summary_path, show_col_types = FALSE) %>%
  arrange(rank_FDR)

B_heatmap <- read_csv(B_heatmap_path, show_col_types = FALSE)
B_patient <- read_csv(B_patient_path, show_col_types = FALSE)
B_cell <- read_csv(B_cell_path, show_col_types = FALSE)

A_survivors <- if (file.exists(A_survivor_path)) {
  read_csv(A_survivor_path, show_col_types = FALSE)
} else {
  tibble()
}

B_survivors <- if (file.exists(B_survivor_path)) {
  read_csv(B_survivor_path, show_col_types = FALSE)
} else {
  tibble()
}

# -----------------------------------------------------------
# 2) DEFINE STRINGENT SURVIVORS
# -----------------------------------------------------------

get_survivor_genes <- function(survivor_df, summary_df, residual_name) {
  
  genes <- if (nrow(survivor_df) > 0 && "gene" %in% colnames(survivor_df)) {
    unique(survivor_df$gene[survivor_df$gene %in% summary_df$gene])
  } else {
    character(0)
  }
  
  cat("\n", residual_name, " stringent survivor genes:\n", sep = "")
  print(genes)
  
  if (length(genes) == 0) {
    warning(
      residual_name,
      ": no stringent survivors found. Survivor map labels, survivor heatmap, ",
      "correlation matrix, and example panels will be skipped for this residual."
    )
  }
  
  genes
}

label_genes_A <- get_survivor_genes(A_survivors, A_summary, "Residual A")
label_genes_B <- get_survivor_genes(B_survivors, B_summary, "Residual B")

write_csv(
  A_summary %>%
    filter(gene %in% label_genes_A),
  file.path(fig_dir, paste0("ResidualA_stringent_survivors_used_for_figures_", gene_set_label, ".csv"))
)

write_csv(
  B_summary %>%
    filter(gene %in% label_genes_B),
  file.path(fig_dir, paste0("ResidualB_stringent_survivors_used_for_figures_", gene_set_label, ".csv"))
)

# -----------------------------------------------------------
# 3) HELPER FUNCTIONS
# -----------------------------------------------------------

nice_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  format(p, scientific = TRUE, digits = 2)
}

zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

save_pdf_png <- function(filename_stub, plot, width, height, dpi = 300) {
  
  pdf_path <- file.path(fig_dir, paste0(filename_stub, ".pdf"))
  png_path <- file.path(fig_dir, paste0(filename_stub, ".png"))
  
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height
  )
  
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

make_volcano <- function(tt, title, subtitle = NULL, lfc = 0.5, fdr = 0.05, top_n = 15) {
  
  df <- tt %>%
    mutate(
      neglog10FDR = -log10(pmax(FDR, 1e-300)),
      sig = case_when(
        FDR < fdr & logFC >  lfc ~ "Up",
        FDR < fdr & logFC < -lfc ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  labs_df <- df %>%
    filter(sig != "NS") %>%
    arrange(FDR) %>%
    slice_head(n = top_n)
  
  ggplot(df, aes(x = logFC, y = neglog10FDR, color = sig)) +
    geom_point(alpha = 0.6, size = 1.1) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = 2) +
    geom_hline(yintercept = -log10(fdr), linetype = 2) +
    ggrepel::geom_text_repel(
      data = labs_df,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3
    ) +
    scale_color_manual(
      values = c(
        Up = "#D55E00",
        Down = "#0072B2",
        NS = "grey80"
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "log2 fold-change",
      y = "-log10(FDR)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

make_patient_paired_plot <- function(df_patient, gene_name, residual_name) {
  
  df <- df_patient %>%
    filter(gene == gene_name) %>%
    mutate(State = factor(State, levels = c("NonProlif", "Prolif")))
  
  if (nrow(df) == 0) return(NULL)
  
  ggplot(df, aes(x = State, y = mean_expr, group = PatientID)) +
    geom_line(alpha = 0.35) +
    geom_point(aes(color = State), size = 2) +
    scale_color_manual(values = state_colors) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(residual_name, " survivor: ", gene_name),
      x = NULL,
      y = "Mean expression per patient-state"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

make_corr_plot <- function(df_cell, gene_name, yvar, ylab = NULL) {
  
  df <- df_cell %>%
    filter(gene == gene_name) %>%
    mutate(State = factor(State, levels = c("NonProlif", "Prolif")))
  
  if (nrow(df) == 0) return(NULL)
  if (is.null(ylab)) ylab <- yvar
  
  ggplot(df, aes(x = gene_expr, y = .data[[yvar]], color = State)) +
    geom_point(alpha = 0.20, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    scale_color_manual(values = state_colors) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(gene_name, " vs ", ylab),
      x = paste0(gene_name, " expression"),
      y = ylab
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

make_summary_label <- function(summary_df, gene_name) {
  
  row <- summary_df %>% filter(gene == gene_name)
  
  if (nrow(row) != 1) return(gene_name)
  
  paste0(
    gene_name, "\n",
    "paired p=", nice_p(row$paired_p_value), "\n",
    "beta_adj=", round(row$beta_adj_mitotic, 3),
    ", p=", nice_p(row$p_adj_mitotic)
  )
}

# -----------------------------------------------------------
# 4) DIRECT AND RESIDUAL VOLCANO FIGURES
# -----------------------------------------------------------

p_directA <- make_volcano(
  directA,
  "Direct Contrast A: Prolif-Tex vs Pro+ non-exhausted",
  "Dominant signal before exhaustion-gene filtering"
)

p_directB <- make_volcano(
  directB,
  "Direct Contrast B: Prolif-Tex vs NonProlif-Tex",
  "Dominant signal before cell-cycle-gene filtering"
)

p_residA <- make_volcano(
  residA,
  "Residual A: canonical exhaustion genes removed",
  "Tests whether signal remains beyond canonical exhaustion markers"
)

p_residB <- make_volcano(
  residB,
  "Residual B: standardized cycling genes removed",
  "Tests whether signal remains beyond canonical cycling markers"
)

save_pdf_png("Figure_DirectContrastA_volcano", p_directA, width = 7, height = 5)
save_pdf_png("Figure_DirectContrastB_volcano", p_directB, width = 7, height = 5)
save_pdf_png("Figure_ResidualA_volcano", p_residA, width = 7, height = 5)
save_pdf_png("Figure_ResidualB_volcano", p_residB, width = 7, height = 5)

# -----------------------------------------------------------
# 5) RESIDUAL HETEROGENEITY MAPS
# -----------------------------------------------------------

make_residual_map <- function(summary_df, survivor_genes, residual_name, file_stub) {
  
  plot_df <- summary_df %>%
    mutate(
      is_stringent_survivor = gene %in% survivor_genes,
      gene_label = ifelse(is_stringent_survivor, gene, NA_character_)
    )
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
    
    geom_point(
      data = plot_df %>% filter(!is_stringent_survivor),
      aes(
        x = cor_mitotic,
        y = beta_adj_mitotic,
        size = abs(cor_inflammatory)
      ),
      color = "grey70",
      alpha = 0.65
    ) +
    
    geom_point(
      data = plot_df %>% filter(is_stringent_survivor),
      aes(
        x = cor_mitotic,
        y = beta_adj_mitotic,
        size = abs(cor_inflammatory)
      ),
      shape = 21,
      fill = "#1b9e77",
      color = "black",
      stroke = 0.8,
      alpha = 0.95
    ) +
    
    ggrepel::geom_text_repel(
      data = plot_df %>% filter(is_stringent_survivor),
      aes(
        x = cor_mitotic,
        y = beta_adj_mitotic,
        label = gene_label
      ),
      size = 4,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.4,
      min.segment.length = 0,
      force = 2,
      force_pull = 0.5
    ) +
    
    scale_size_continuous(
      name = "|cor inflammatory|",
      range = c(2.5, 6)
    ) +
    
    theme_classic(base_size = 12) +
    labs(
      title = paste0(residual_name, " gene heterogeneity map"),
      subtitle = "Grey points: all top residual genes; labelled points: stringent survivors",
      x = "Spearman correlation with Mitotic score",
      y = "State association after Mitotic adjustment\n(beta for Prolif vs NonProlif)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  save_pdf_png(
    paste0("Figure_", file_stub, "_heterogeneity_map_", gene_set_label),
    p,
    width = 10.5,
    height = 6.5
  )
  
  p
}

pA_map <- make_residual_map(
  summary_df = A_summary,
  survivor_genes = label_genes_A,
  residual_name = "Residual A",
  file_stub = "ResidualA"
)

pB_map <- make_residual_map(
  summary_df = B_summary,
  survivor_genes = label_genes_B,
  residual_name = "Residual B",
  file_stub = "ResidualB"
)

# -----------------------------------------------------------
# 6) TOP-GENE HETEROGENEITY HEATMAPS
# -----------------------------------------------------------

heatmap_metric_labels <- c(
  median_diff = "Paired mean diff",
  frac_expr_diff = "Fraction expressing diff",
  cor_G1S = "cor G1/S",
  cor_checkpoint = "cor Checkpoint",
  cor_mitotic = "cor Mitotic",
  cor_metabolic = "cor Metabolic",
  cor_inflammatory = "cor Inflammatory",
  cor_exhaustion = "cor Exhaustion",
  beta_adj_mitotic = "beta adj Mitotic"
)

make_top_gene_heatmap <- function(heatmap_df, summary_df, residual_name, file_stub) {
  
  gene_levels <- rev(summary_df$gene)
  
  plot_df <- heatmap_df %>%
    mutate(
      gene = factor(gene, levels = gene_levels),
      metric = factor(
        metric,
        levels = names(heatmap_metric_labels),
        labels = heatmap_metric_labels
      )
    )
  
  p <- ggplot(plot_df, aes(x = metric, y = gene, fill = value_z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      name = "Z-score"
    ) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0(residual_name, " top-gene heterogeneity heatmap (", gene_set_label, ")"),
      subtitle = "Top residual genes show mixed associations across biological axes",
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  save_pdf_png(
    paste0("Figure_", file_stub, "_top_gene_heatmap_", gene_set_label),
    p,
    width = 10,
    height = 7
  )
  
  p
}

pA_heatmap <- make_top_gene_heatmap(
  heatmap_df = A_heatmap,
  summary_df = A_summary,
  residual_name = "Residual A",
  file_stub = "ResidualA"
)

pB_heatmap <- make_top_gene_heatmap(
  heatmap_df = B_heatmap,
  summary_df = B_summary,
  residual_name = "Residual B",
  file_stub = "ResidualB"
)

# -----------------------------------------------------------
# 7) STRINGENT-SURVIVOR HEATMAPS
# -----------------------------------------------------------

make_survivor_heatmap <- function(summary_df, survivor_genes, residual_name, file_stub) {
  
  if (length(survivor_genes) == 0) {
    message("Skipping ", residual_name, " survivor heatmap: no stringent survivors.")
    return(NULL)
  }
  
  reduced_metrics <- c(
    "median_diff",
    "frac_expr_diff",
    "cor_mitotic",
    "cor_inflammatory",
    "cor_exhaustion",
    "beta_adj_mitotic"
  )
  
  reduced_labels <- c(
    median_diff = "Patient-paired\nexpression difference",
    frac_expr_diff = "Fraction-expressing\ndifference",
    cor_mitotic = "Mitotic\ncorrelation",
    cor_inflammatory = "Inflammatory\ncorrelation",
    cor_exhaustion = "Exhaustion\ncorrelation",
    beta_adj_mitotic = "Mitotic-adjusted\nstate effect"
  )
  
  subtitle_text <- case_when(
    file_stub == "ResidualA" ~
      "Genes retained after exhaustion-gene filtering and predefined residual-survivor criteria",
    file_stub == "ResidualB" ~
      "Genes retained after cycling-gene filtering and predefined residual-survivor criteria",
    TRUE ~
      "Genes retained after predefined residual-survivor criteria"
  )
  
  plot_df <- summary_df %>%
    filter(gene %in% survivor_genes) %>%
    select(gene, all_of(reduced_metrics)) %>%
    pivot_longer(
      cols = all_of(reduced_metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    group_by(metric) %>%
    mutate(value_z = zscore_safe(value)) %>%
    ungroup() %>%
    mutate(
      gene = factor(gene, levels = rev(survivor_genes)),
      metric = factor(
        metric,
        levels = names(reduced_labels),
        labels = reduced_labels
      )
    )
  
  p <- ggplot(plot_df, aes(x = metric, y = gene, fill = value_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      name = "Z-score"
    ) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0(residual_name, ": stringent residual gene profiles"),
      subtitle = subtitle_text,
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10.5),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  save_pdf_png(
    paste0("Figure_", file_stub, "_stringent_survivor_heatmap_", gene_set_label),
    p,
    width = 8.5,
    height = max(4.8, 0.48 * length(survivor_genes) + 2.8)
  )
  
  p
}

pA_survivor_heatmap <- make_survivor_heatmap(
  summary_df = A_summary,
  survivor_genes = label_genes_A,
  residual_name = "Residual A",
  file_stub = "ResidualA"
)

pB_survivor_heatmap <- make_survivor_heatmap(
  summary_df = B_summary,
  survivor_genes = label_genes_B,
  residual_name = "Residual B",
  file_stub = "ResidualB"
)

# -----------------------------------------------------------
# 8) SURVIVOR CORRELATION MATRICES
# -----------------------------------------------------------

make_survivor_corr_plot <- function(corr_path, residual_name, file_stub) {
  
  if (!file.exists(corr_path)) {
    message("Skipping ", residual_name, " correlation matrix: file not found.")
    return(NULL)
  }
  
  corr_df <- read_csv(corr_path, show_col_types = FALSE)
  
  if (!"gene" %in% colnames(corr_df) || ncol(corr_df) < 3) {
    message("Skipping ", residual_name, " correlation matrix: insufficient survivor genes.")
    return(NULL)
  }
  
  corr_long <- corr_df %>%
    pivot_longer(
      cols = -gene,
      names_to = "Gene2",
      values_to = "Correlation"
    ) %>%
    rename(Gene1 = gene)
  
  gene_order <- colnames(corr_df)[-1]
  
  corr_long <- corr_long %>%
    mutate(
      Gene1 = factor(Gene1, levels = gene_order),
      Gene2 = factor(Gene2, levels = gene_order)
    )
  
  p <- ggplot(corr_long, aes(Gene1, Gene2, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0(residual_name, " survivor correlation structure (", gene_set_label, ")"),
      subtitle = "Weak or mixed co-variation argues against a coherent residual module",
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  save_pdf_png(
    paste0("Figure_", file_stub, "_survivor_correlation_", gene_set_label),
    p,
    width = 6,
    height = 5
  )
  
  p
}

pA_corr <- make_survivor_corr_plot(
  corr_path = A_corr_path,
  residual_name = "Residual A",
  file_stub = "ResidualA"
)

pB_corr <- make_survivor_corr_plot(
  corr_path = B_corr_path,
  residual_name = "Residual B",
  file_stub = "ResidualB"
)

# -----------------------------------------------------------
# 9) OPTIONAL SURVIVOR EXAMPLE PANELS
# -----------------------------------------------------------

make_survivor_examples <- function(patient_df, cell_df, summary_df, survivor_genes,
                                   residual_name, file_stub) {
  
  if (length(survivor_genes) == 0) {
    message("Skipping ", residual_name, " example panels: no stringent survivors.")
    return(invisible(NULL))
  }
  
  for (gene in survivor_genes) {
    
    if (!(gene %in% summary_df$gene)) next
    
    p_pair <- make_patient_paired_plot(
      df_patient = patient_df,
      gene_name = gene,
      residual_name = residual_name
    )
    
    p_mit <- make_corr_plot(cell_df, gene, "Mitotic_score", "Mitotic score")
    p_inf <- make_corr_plot(cell_df, gene, "Inflammatory_score", "Inflammatory score")
    p_exh <- make_corr_plot(cell_df, gene, "Exhaustion_score", "Exhaustion score")
    
    label_text <- make_summary_label(summary_df, gene)
    
    p_title <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = label_text,
        hjust = 0,
        vjust = 1,
        size = 4.2
      ) +
      theme_void()
    
    plots_exist <- list(p_pair, p_mit, p_inf, p_exh)
    
    if (any(vapply(plots_exist, is.null, logical(1)))) next
    
    fig_gene <- (p_pair | p_title) / (p_mit | p_inf) / p_exh
    
    save_pdf_png(
      paste0(
        "Supplementary_",
        file_stub,
        "_",
        gene,
        "_stringent_survivor_example_",
        gene_set_label
      ),
      fig_gene,
      width = 12,
      height = 12
    )
  }
}

make_survivor_examples(
  patient_df = A_patient,
  cell_df = A_cell,
  summary_df = A_summary,
  survivor_genes = label_genes_A,
  residual_name = "Residual A",
  file_stub = "ResidualA"
)

make_survivor_examples(
  patient_df = B_patient,
  cell_df = B_cell,
  summary_df = B_summary,
  survivor_genes = label_genes_B,
  residual_name = "Residual B",
  file_stub = "ResidualB"
)

# -----------------------------------------------------------
# 10) COMBINED SUMMARY TABLE FOR CAPTIONING / WRITING
# -----------------------------------------------------------

summary_cols <- c(
  "source", "stringent_survivor",
  "gene", "rank_FDR", "logFC", "FDR",
  "median_diff", "paired_p_value",
  "frac_expr_diff",
  "cor_G1S", "cor_checkpoint", "cor_mitotic",
  "cor_metabolic", "cor_inflammatory", "cor_exhaustion",
  "beta_adj_mitotic", "p_adj_mitotic", "heuristic_class"
)

concept_summary_A <- A_summary %>%
  mutate(
    source = "ResidualA",
    stringent_survivor = gene %in% label_genes_A
  ) %>%
  select(any_of(summary_cols))

concept_summary_B <- B_summary %>%
  mutate(
    source = "ResidualB",
    stringent_survivor = gene %in% label_genes_B
  ) %>%
  select(any_of(summary_cols))

combined_summary <- bind_rows(concept_summary_A, concept_summary_B)

write_csv(
  combined_summary,
  file.path(
    fig_dir,
    paste0("Combined_residual_gene_summary_for_captioning_", gene_set_label, ".csv")
  )
)

# -----------------------------------------------------------
# 11) FINAL REPORT
# -----------------------------------------------------------

cat("\nScript 4 complete.\n")
cat("Figures written to:\n", normalizePath(fig_dir), "\n")

cat("\nResidual A stringent survivors used for labelled figures:\n")
print(label_genes_A)

cat("\nResidual B stringent survivors used for labelled figures:\n")
print(label_genes_B)