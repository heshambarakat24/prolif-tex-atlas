# ===========================================================
# 19_residual_contrasts.R
# ===========================================================
# Thesis Methods section: 4.1.11 (residual differential expression)
#
# Reads : Residual_pipeline_outputs/DirectContrastA_*_full_edgeR.csv
#         Residual_pipeline_outputs/DirectContrastB_*_full_edgeR.csv
# Writes: residual DE tables, ranked .rnk files, and the two exclusion lists
#         actually used - Canonical_exhaustion_exclusion_genes.txt and
#         Standardized_cellcycle_exclusion_union_Hallmark_Reactome_GOBP.txt
#
# This script filters genes out of the direct-contrast tables produced by
# script 18. It does not re-run differential expression and does not
# re-derive cell groups.
#
#   Residual A: contrast A with canonical exhaustion genes removed
#   Residual B: contrast B with a standardised cycling programme removed
#               (Hallmark, Reactome and GO:BP union, resolved via msigdbr)
#
# The cycling exclusion set is built at run time from MSigDB rather than
# hard-coded. The resolved list is written to the output folder and should be
# deposited alongside the results, since it cannot be reconstructed from this
# file alone.
#
# Paths are relative to the working directory. See README.md.
# ===========================================================

# ===========================================================
# SCRIPT 2 — Residual contrasts after dominant-program filtering
#
# Goal:
#   Starting from the full direct-contrast edgeR outputs produced
#   in Script 1, remove dominant canonical programs downstream and
#   test whether coherent residual structure remains.
#
# Residual contrasts:
#   Residual A:
#     Direct Contrast A after removal of canonical exhaustion genes
#
#   Residual B:
#     Direct Contrast B after removal of a standardized canonical
#     cycling program (Hallmark + Reactome + GO:BP)
#
# Important interpretation:
#   This is downstream gene filtering, not latent-variable regression.
#   The goal is to test whether coherent structure remains beyond
#   dominant canonical markers.
#
# Key outputs:
#   1) Residual DE tables
#   2) Residual volcano plots
#   3) Ranked .rnk files for residual GSEA
#   4) Standardized cycling exclusion membership files
#   5) Top20 / Top100 residual up-gene lists
#
# Inputs (from Script 1):
#   Residual_pipeline_outputs/
#     DirectContrastA_ProlifTex_vs_ProPlus_nonExh_full_edgeR.csv
#     DirectContrastB_ProlifTex_vs_NonProlifTex_full_edgeR.csv
#     Residual_pipeline_obj.rds
# ===========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(msigdbr)
  library(tibble)
})

# -----------------------------------------------------------
# 0) SETTINGS
# -----------------------------------------------------------
outdir <- "Residual_pipeline_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

A_path <- file.path(outdir, "DirectContrastA_ProlifTex_vs_ProPlus_nonExh_full_edgeR.csv")
B_path <- file.path(outdir, "DirectContrastB_ProlifTex_vs_NonProlifTex_full_edgeR.csv")
obj_path <- file.path(outdir, "Residual_pipeline_obj.rds")

if (!file.exists(A_path)) stop("Missing input file: ", A_path)
if (!file.exists(B_path)) stop("Missing input file: ", B_path)
if (!file.exists(obj_path)) stop("Missing input file: ", obj_path)

ttA <- read_csv(A_path, show_col_types = FALSE)
ttB <- read_csv(B_path, show_col_types = FALSE)

# -----------------------------------------------------------
# 1) Define dominant programs to remove downstream
# -----------------------------------------------------------

# Residual A: canonical exhaustion markers
canonical_exhaustion <- unique(toupper(c(
  "PDCD1","HAVCR2","CTLA4","TIGIT","LAG3","TOX","LAYN","ENTPD1","CXCL13",
  "NR4A1","NR4A2","NR4A3","EOMES","BATF","CD244","CD160","LGALS3"
)))

write.table(
  sort(canonical_exhaustion),
  file = file.path(outdir, "Canonical_exhaustion_exclusion_genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Residual B: standardized canonical cycling genes
hallmark_cc <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  filter(gs_name %in% c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS")) %>%
  transmute(
    source = "Hallmark",
    pathway = gs_name,
    gene = toupper(gene_symbol)
  ) %>%
  distinct()

hallmark_cc_genes <- unique(hallmark_cc$gene)

reactome_all <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
)

reactome_keep_patterns <- c(
  "CELL_CYCLE",
  "DNA_REPLICATION",
  "S_PHASE",
  "G2_M_CHECKPOINT",
  "MITOTIC",
  "M_PHASE",
  "CHROMOSOME_SEGREGATION",
  "SISTER_CHROMATID",
  "SPINDLE",
  "METAPHASE",
  "ANAPHASE",
  "CYTOKINESIS"
)

reactome_regex <- stringr::str_c(reactome_keep_patterns, collapse = "|")

reactome_cc <- reactome_all %>%
  filter(stringr::str_detect(gs_name, reactome_regex)) %>%
  transmute(
    source = "Reactome",
    pathway = gs_name,
    gene = toupper(gene_symbol)
  ) %>%
  distinct()

reactome_cc_genes <- unique(reactome_cc$gene)

gobp_all <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:BP"
)

gobp_keep_patterns <- c(
  "MITOTIC_CELL_CYCLE",
  "CELL_DIVISION",
  "NUCLEAR_DIVISION",
  "DNA_REPLICATION",
  "CHROMOSOME_SEGREGATION",
  "SISTER_CHROMATID_SEGREGATION",
  "SPINDLE_ORGANIZATION",
  "MITOTIC_SPINDLE",
  "CHROMOSOME_ORGANIZATION",
  "MITOTIC_NUCLEAR_DIVISION"
)

gobp_regex <- stringr::str_c(gobp_keep_patterns, collapse = "|")

gobp_cc <- gobp_all %>%
  filter(stringr::str_detect(gs_name, gobp_regex)) %>%
  transmute(
    source = "GO_BP",
    pathway = gs_name,
    gene = toupper(gene_symbol)
  ) %>%
  distinct()

gobp_cc_genes <- unique(gobp_cc$gene)

standard_cc_genes <- unique(c(
  hallmark_cc_genes,
  reactome_cc_genes,
  gobp_cc_genes
))

cat("\nCycling exclusion sizes:\n")
cat("Hallmark genes:", length(hallmark_cc_genes), "\n")
cat("Reactome genes:", length(reactome_cc_genes), "\n")
cat("GO:BP genes:", length(gobp_cc_genes), "\n")
cat("Standardized union:", length(standard_cc_genes), "\n")

standard_membership <- bind_rows(
  hallmark_cc,
  reactome_cc,
  gobp_cc
) %>%
  distinct() %>%
  arrange(source, pathway, gene)

write.table(
  sort(standard_cc_genes),
  file = file.path(outdir, "Standardized_cellcycle_exclusion_union_Hallmark_Reactome_GOBP.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write_csv(
  standard_membership,
  file.path(outdir, "Standardized_cellcycle_exclusion_membership_Hallmark_Reactome_GOBP.csv")
)

# -----------------------------------------------------------
# 2) Residual filtering
# -----------------------------------------------------------
ttA_resid <- ttA %>%
  filter(!toupper(gene) %in% canonical_exhaustion)

ttB_resid <- ttB %>%
  filter(!toupper(gene) %in% standard_cc_genes)

write_csv(
  ttA_resid,
  file.path(outdir, "ResidualA_ProlifTex_vs_ProPlus_nonExh_NO_EXHAUSTION_GENES.csv")
)

write_csv(
  ttB_resid,
  file.path(outdir, "ResidualB_ProlifTex_vs_NonProlifTex_NO_STANDARD_CELLCYCLE_HRG.csv")
)

# -----------------------------------------------------------
# 3) Residual volcano plots
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

pA <- make_volcano(
  ttA_resid,
  "Residual A: Prolif-Tex vs Pro+ non-exhausted\n(canonical exhaustion genes removed)"
)

pB <- make_volcano(
  ttB_resid,
  "Residual B: Prolif-Tex vs NonProlif-Tex\n(standardized cycling genes removed)"
)

print(pA)
print(pB)

ggsave(file.path(outdir, "ResidualA_volcano.png"), pA, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "ResidualB_volcano.png"), pB, width = 7, height = 5, dpi = 300)

# -----------------------------------------------------------
# 4) Ranked lists for residual GSEA
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

rankA <- mk_rank_table(ttA_resid)
rankB <- mk_rank_table(ttB_resid)

write.table(
  rankA,
  file = file.path(outdir, "ResidualA_ranked_for_GSEA.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  rankB,
  file = file.path(outdir, "ResidualB_ranked_for_GSEA.rnk"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# -----------------------------------------------------------
# 5) Export significant residual genes
# -----------------------------------------------------------
export_ranked <- function(tt, prefix, fdr_cut = 0.05, lfc_cut = 0.25) {
  up   <- tt %>% filter(FDR < fdr_cut, logFC >  lfc_cut) %>% arrange(FDR)
  down <- tt %>% filter(FDR < fdr_cut, logFC < -lfc_cut) %>% arrange(FDR)
  
  write_csv(up[, c("gene", "logFC", "FDR")],   file.path(outdir, paste0(prefix, "_UP_sig.csv")))
  write_csv(down[, c("gene", "logFC", "FDR")], file.path(outdir, paste0(prefix, "_DOWN_sig.csv")))
}

export_ranked(ttA_resid, "ResidualA")
export_ranked(ttB_resid, "ResidualB")

# -----------------------------------------------------------
# 6) Export Top20 and Top100 residual up-genes
# -----------------------------------------------------------
export_top_ranked <- function(tt, prefix, n_top = 20) {
  tt %>%
    arrange(FDR) %>%
    filter(logFC > 0) %>%
    slice_head(n = n_top) %>%
    mutate(rank_FDR = row_number()) %>%
    select(gene, rank_FDR, logFC, FDR) %>%
    write_csv(file.path(outdir, paste0(prefix, "_top", n_top, "_up_genes.csv")))
}

export_top_ranked(ttA_resid, "ResidualA", 20)
export_top_ranked(ttA_resid, "ResidualA", 100)
export_top_ranked(ttB_resid, "ResidualB", 20)
export_top_ranked(ttB_resid, "ResidualB", 100)

# -----------------------------------------------------------
# 7) Quick summaries
# -----------------------------------------------------------
message("Residual A: top residual UP genes in Prolif-Tex")
print(
  ttA_resid %>%
    arrange(FDR) %>%
    filter(logFC > 0) %>%
    select(gene, logFC, FDR) %>%
    head(20)
)

message("Residual B: top residual UP genes in Prolif-Tex")
print(
  ttB_resid %>%
    arrange(FDR) %>%
    filter(logFC > 0) %>%
    select(gene, logFC, FDR) %>%
    head(20)
)

message("\nScript 2 complete: residual contrasts established.")
message("All outputs written to: ", normalizePath(outdir))
message("Next step: Script 3 characterizes residual survivor genes within exhausted cells.")