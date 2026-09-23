# ===========================================================
# 18_pseudobulk_direct_contrasts.R
# ===========================================================
# Thesis Methods section: 4.1.10 (pseudobulk differential expression)
#
# Requires: CD8_TumourOnly_scVI in the R session (script 09).
# Writes  : Residual_pipeline_outputs/ - State_thresholds.csv, group counts,
#           DirectContrastA_*_full_edgeR.csv, DirectContrastB_*_full_edgeR.csv,
#           volcano plots, ranked .rnk files, and Residual_pipeline_obj.rds
#           carrying the CMP_A / CMP_B labels used by scripts 19-21.
#
# Contrasts:
#   A: Prolif-Tex vs proliferating non-exhausted (Exh-low, Prolif-high)
#   B: Prolif-Tex vs NonProlif-Tex
#
# Gene panels: these contrasts use the extended canonical panels specified in
# Methods 4.1.10 - fourteen exhaustion genes and twelve proliferation genes -
# rather than the nine-gene panels used for configuration assignment in
# Methods 4.1.7. Scores are z-normalised within dataset and the same
# percentile thresholds are applied. Both panel sets are in panels/.
#
# Note: AddModuleScore draws control genes at random and no seed is set here,
# so a re-run shifts a small number of borderline cells between groups.
# See README.
#
# Paths are relative to the working directory. Input data are not included in
# this repository; see README.md.
# ===========================================================

# ===========================================================
# SCRIPT 1 — Direct contrasts in tumour CD8 T cells
#
# Goal:
#   Define Prolif-Tex-related comparison groups and quantify
#   transcriptomic differences using paired patient-level
#   pseudobulk differential expression.
#
# Biological purpose:
#   This script establishes the dominant signal in the direct
#   contrasts before any residual filtering is applied.
#
# Direct contrasts:
#   Contrast A: Prolif-Tex vs Pro+ non-exhausted
#   Contrast B: Prolif-Tex vs NonProlif-Tex
#
# Key outputs:
#   1) State thresholds and group counts
#   2) Full edgeR DE tables for both direct contrasts
#   3) Direct volcano plots
#   4) Ranked .rnk files for downstream direct GSEA
#   5) Saved Seurat object with CMP_A / CMP_B labels
#
# Downstream logic:
#   Script 1 = direct contrasts
#   Script 2 = residual contrasts after removing dominant programs
#
# Input object:
#   CD8_TumourOnly_scVI
#
# Required metadata:
#   PatientID, Dataset
#
# Will compute if missing:
#   Exh_z, Prolif_z
# ===========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(edgeR)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

# -----------------------------------------------------------
# OUTPUT DIRECTORY
# -----------------------------------------------------------
outdir <- "Residual_pipeline_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------
# 0) INPUT OBJECT
# -----------------------------------------------------------
obj <- CD8_TumourOnly_scVI
DefaultAssay(obj) <- "RNA"

meta <- obj@meta.data
stopifnot(all(c("PatientID", "Dataset") %in% colnames(meta)))

cts <- GetAssayData(obj, slot = "counts")
stopifnot(inherits(cts, "dgCMatrix"))

# -----------------------------------------------------------
# 1) Define core signatures and compute Exh_z / Prolif_z if missing
# -----------------------------------------------------------
if (!all(c("Exh_z", "Prolif_z") %in% colnames(meta))) {
  message("Exh_z / Prolif_z not found. Computing from module scores...")
  
  exhaustion_genes <- c(
    "PDCD1","LAYN","HAVCR2","CTLA4","TOX",
    "TIGIT","LAG3","CXCL13","ENTPD1",
    "NR4A1","NR4A2","NR4A3","EOMES","BATF"
  )
  
  proliferation_genes <- c(
    "MKI67","TOP2A","PCNA","HMGB2","CDC20",
    "UBE2C","BIRC5","AURKB","CCNB1",
    "MCM2","MCM5","TYMS"
  )
  
  exhaustion_genes    <- intersect(exhaustion_genes, rownames(obj))
  proliferation_genes <- intersect(proliferation_genes, rownames(obj))
  
  stopifnot(length(exhaustion_genes) >= 5)
  stopifnot(length(proliferation_genes) >= 5)
  
  if (!any(grepl("^RAW_Exhaustion", colnames(obj@meta.data)))) {
    obj <- AddModuleScore(
      obj,
      features = list(exhaustion_genes),
      name = "RAW_Exhaustion",
      assay = "RNA"
    )
  }
  
  if (!any(grepl("^RAW_Proliferation", colnames(obj@meta.data)))) {
    obj <- AddModuleScore(
      obj,
      features = list(proliferation_genes),
      name = "RAW_Proliferation",
      assay = "RNA"
    )
  }
  
  meta <- obj@meta.data
  raw_exh_col  <- grep("^RAW_Exhaustion", colnames(meta), value = TRUE)[1]
  raw_prol_col <- grep("^RAW_Proliferation", colnames(meta), value = TRUE)[1]
  
  z_within <- function(x, g) {
    sp <- split(x, g)
    as.numeric(unsplit(lapply(sp, function(v) {
      m <- mean(v, na.rm = TRUE)
      s <- stats::sd(v, na.rm = TRUE)
      if (is.na(s) || s == 0) rep(0, length(v)) else (v - m) / s
    }), g))
  }
  
  meta$Exh_z    <- z_within(meta[[raw_exh_col]],  meta$Dataset)
  meta$Prolif_z <- z_within(meta[[raw_prol_col]], meta$Dataset)
  
  obj@meta.data <- meta
}

meta <- obj@meta.data
stopifnot(all(c("Exh_z", "Prolif_z") %in% colnames(meta)))

# -----------------------------------------------------------
# 2) Define quartile-based states
# -----------------------------------------------------------
q75_exh <- quantile(meta$Exh_z,    0.75, na.rm = TRUE)
q25_exh <- quantile(meta$Exh_z,    0.25, na.rm = TRUE)
q75_pro <- quantile(meta$Prolif_z, 0.75, na.rm = TRUE)
q25_pro <- quantile(meta$Prolif_z, 0.25, na.rm = TRUE)

is_exh_pos <- meta$Exh_z    >= q75_exh
is_exh_neg <- meta$Exh_z    <= q25_exh
is_pro_pos <- meta$Prolif_z >= q75_pro
is_pro_neg <- meta$Prolif_z <= q25_pro

# Contrast A:
# Prolif-Tex (Exh+ Pro+) vs Pro+ non-exhausted (Exh- Pro+)
meta$CMP_A <- NA_character_
meta$CMP_A[is_exh_pos & is_pro_pos] <- "Prolif-Tex"
meta$CMP_A[is_exh_neg & is_pro_pos] <- "ProPlus_nonExh"

# Contrast B:
# Prolif-Tex (Exh+ Pro+) vs NonProlif-Tex (Exh+ Pro-)
meta$CMP_B <- NA_character_
meta$CMP_B[is_exh_pos & is_pro_pos] <- "Prolif-Tex"
meta$CMP_B[is_exh_pos & is_pro_neg] <- "NonProlif-Tex"

obj@meta.data <- meta

cat("CMP_A counts:\n")
print(table(obj$CMP_A, useNA = "ifany"))
cat("\nCMP_B counts:\n")
print(table(obj$CMP_B, useNA = "ifany"))

state_summary <- tibble(
  q75_exh = q75_exh,
  q25_exh = q25_exh,
  q75_pro = q75_pro,
  q25_pro = q25_pro
)

write_csv(state_summary, file.path(outdir, "State_thresholds.csv"))

cmpA_counts <- as.data.frame(table(obj$CMP_A, useNA = "ifany"))
cmpB_counts <- as.data.frame(table(obj$CMP_B, useNA = "ifany"))

write_csv(cmpA_counts, file.path(outdir, "CMP_A_counts.csv"))
write_csv(cmpB_counts, file.path(outdir, "CMP_B_counts.csv"))

# -----------------------------------------------------------
# 3) Helper: pseudobulk by patient × group
# -----------------------------------------------------------
make_pseudobulk <- function(obj, counts, group_col, patient_col = "PatientID",
                            ref_level = NULL, min_cpm = 1, min_samples = 3) {
  md <- obj@meta.data
  stopifnot(all(c(patient_col, group_col) %in% names(md)))
  
  ok <- !is.na(md[[group_col]]) & !is.na(md[[patient_col]])
  if (!any(ok)) stop("No cells available for ", group_col)
  
  g_all   <- droplevels(factor(md[[group_col]][ok]))
  pat_all <- droplevels(factor(md[[patient_col]][ok]))
  
  tab_pg   <- table(pat_all, g_all)
  keep_pat <- rownames(tab_pg)[rowSums(tab_pg > 0) == ncol(tab_pg)]
  
  keep <- ok
  keep[ok] <- md[[patient_col]][ok] %in% keep_pat
  if (!any(keep)) stop("No paired patients for ", group_col)
  
  g   <- droplevels(factor(md[[group_col]][keep]))
  pat <- droplevels(factor(md[[patient_col]][keep]))
  
  SAFE_SEP <- "|||"
  sample_id <- factor(paste(pat, g, sep = SAFE_SEP))
  
  J <- Matrix::sparseMatrix(
    i = seq_along(sample_id),
    j = as.integer(sample_id),
    x = 1,
    dims = c(length(sample_id), nlevels(sample_id))
  )
  
  counts_pb <- counts[, keep, drop = FALSE] %*% J
  colnames(counts_pb) <- levels(sample_id)
  
  parts <- do.call(rbind, strsplit(colnames(counts_pb), SAFE_SEP, fixed = TRUE))
  samp_df <- tibble(
    sample  = colnames(counts_pb),
    patient = factor(parts[, 1]),
    group   = factor(parts[, 2])
  )
  
  if (nlevels(samp_df$group) != 2) {
    print(table(samp_df$group))
    stop(group_col, " requires exactly two groups.")
  }
  
  if (!is.null(ref_level)) {
    samp_df$group <- relevel(samp_df$group, ref = ref_level)
  }
  
  dge <- DGEList(counts = counts_pb)
  keep_genes <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  
  list(dge = dge, samp = samp_df)
}

# -----------------------------------------------------------
# 4) Run edgeR DE for both direct contrasts
# -----------------------------------------------------------
run_edgeR_pb <- function(obj, counts, group_col, ref_level, out_prefix) {
  pb <- make_pseudobulk(
    obj, counts,
    group_col = group_col,
    patient_col = "PatientID",
    ref_level = ref_level
  )
  
  message(out_prefix, " paired patients: ", length(unique(pb$samp$patient)))
  print(table(pb$samp$group))
  
  design <- model.matrix(~ patient + group, data = pb$samp)
  pb$dge <- estimateDisp(pb$dge, design)
  fit <- glmQLFit(pb$dge, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = ncol(design))
  tt <- topTags(qlf, n = Inf)$table %>% rownames_to_column("gene")
  
  write_csv(tt, file.path(outdir, paste0(out_prefix, "_full_edgeR.csv")))
  
  list(pb = pb, design = design, fit = fit, qlf = qlf, tt = tt)
}

# Direct Contrast A
resA <- run_edgeR_pb(
  obj, cts,
  group_col = "CMP_A",
  ref_level = "ProPlus_nonExh",
  out_prefix = "DirectContrastA_ProlifTex_vs_ProPlus_nonExh"
)

# Direct Contrast B
resB <- run_edgeR_pb(
  obj, cts,
  group_col = "CMP_B",
  ref_level = "NonProlif-Tex",
  out_prefix = "DirectContrastB_ProlifTex_vs_NonProlifTex"
)

ttA <- resA$tt
ttB <- resB$tt

# Interpretation:
# A logFC > 0 => higher in Prolif-Tex
# B logFC > 0 => higher in Prolif-Tex

# -----------------------------------------------------------
# 5) Volcano plots for direct contrasts
# -----------------------------------------------------------
make_volcano <- function(tt, title, lfc = 0.5, fdr = 0.05, top_n = 15) {
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
    scale_color_manual(values = c(Up = "#D55E00", Down = "#0072B2", NS = "grey80")) +
    labs(title = title, x = "log2 fold-change", y = "-log10(FDR)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

pA_direct <- make_volcano(
  ttA,
  "Direct Contrast A: Prolif-Tex vs Pro+ non-exhausted"
)

pB_direct <- make_volcano(
  ttB,
  "Direct Contrast B: Prolif-Tex vs NonProlif-Tex"
)

print(pA_direct)
print(pB_direct)

ggsave(file.path(outdir, "DirectContrastA_volcano.png"), pA_direct, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "DirectContrastB_volcano.png"), pB_direct, width = 7, height = 5, dpi = 300)

# -----------------------------------------------------------
# 6) Ranked lists for downstream direct GSEA
# -----------------------------------------------------------
mk_rank_table <- function(tt) {
  tt %>%
    transmute(
      gene = toupper(str_replace(gene, "\\.\\d+$", "")),
      score = sign(logFC) * -log10(pmax(PValue, 1e-300))
    ) %>%
    filter(!is.na(gene), is.finite(score)) %>%
    group_by(gene) %>%
    summarise(score = max(score), .groups = "drop") %>%
    arrange(desc(score))
}

rankA_direct <- mk_rank_table(ttA)
rankB_direct <- mk_rank_table(ttB)

write.table(
  rankA_direct,
  file = file.path(outdir, "DirectContrastA_ranked_for_GSEA.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  rankB_direct,
  file = file.path(outdir, "DirectContrastB_ranked_for_GSEA.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# -----------------------------------------------------------
# 7) Export significant genes for quick review
# -----------------------------------------------------------
export_ranked <- function(tt, prefix, fdr_cut = 0.05, lfc_cut = 0.25) {
  up   <- tt %>% filter(FDR < fdr_cut, logFC >  lfc_cut) %>% arrange(FDR)
  down <- tt %>% filter(FDR < fdr_cut, logFC < -lfc_cut) %>% arrange(FDR)
  
  write_csv(up[, c("gene", "logFC", "FDR")],   file.path(outdir, paste0(prefix, "_UP_sig.csv")))
  write_csv(down[, c("gene", "logFC", "FDR")], file.path(outdir, paste0(prefix, "_DOWN_sig.csv")))
}

export_ranked(ttA, "DirectContrastA")
export_ranked(ttB, "DirectContrastB")

# -----------------------------------------------------------
# 8) Quick summaries
# -----------------------------------------------------------
message("Direct Contrast A: top UP genes in Prolif-Tex")
print(
  ttA %>%
    arrange(FDR) %>%
    filter(logFC > 0) %>%
    select(gene, logFC, FDR) %>%
    head(20)
)

message("Direct Contrast B: top UP genes in Prolif-Tex")
print(
  ttB %>%
    arrange(FDR) %>%
    filter(logFC > 0) %>%
    select(gene, logFC, FDR) %>%
    head(20)
)

# -----------------------------------------------------------
# 9) Save labeled object for downstream scripts
# -----------------------------------------------------------
saveRDS(obj, file = file.path(outdir, "Residual_pipeline_obj.rds"))

message("\nScript 1 complete: direct contrasts established.")
message("All outputs written to: ", normalizePath(outdir))
message("Next step: Script 2 performs residual filtering of these direct contrasts.")