# ============================================================
# 12_patient_prevalence_compartments.R
# ============================================================
# Thesis Methods section: 4.1.8 (patient-level prevalence across compartments)
#
# Requires: CD8_Full_scVI in the R session with TexClass attached, i.e. run
#           08 then 11 first. This script does not load anything from disk.
# Produces: patient_tex, the per-patient/per-compartment prevalence table used
#           by scripts 13 and 14.
#
# Patient, not cell, is the unit of analysis: prevalence is computed per
# patient and tissue compartment to avoid pseudo-replication.
# ============================================================

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# CHUNK 1 — Extract metadata from Seurat object
# Purpose: Create a working dataframe containing patient,
#          tissue compartment, and Tex classification.
# ------------------------------------------------------------

md <- CD8_Full_scVI@meta.data

df <- md %>%
  select(PatientID, Tissue, TexClass)


# ------------------------------------------------------------
# CHUNK 2 — Compute per-patient CD8 counts
# Purpose: Determine total CD8 cells per patient and tissue
#          compartment.
# ------------------------------------------------------------

patient_totals <- df %>%
  group_by(PatientID, Tissue) %>%
  summarise(
    total_cd8 = n(),
    .groups = "drop"
  )


# ------------------------------------------------------------
# CHUNK 3 — Compute per-patient counts of Tex populations
# Purpose: Count Prolif-Tex and Non-Prolif-Tex cells for each
#          patient in each tissue.
# ------------------------------------------------------------

patient_tex <- df %>%
  filter(TexClass %in% c("Prolif-Tex", "Non-Prolif-Tex")) %>%
  group_by(PatientID, Tissue, TexClass) %>%
  summarise(
    tex_cells = n(),
    .groups = "drop"
  )


# ------------------------------------------------------------
# CHUNK 4 — Merge totals with Tex counts
# Purpose: Calculate fraction of Tex cells among all CD8 cells
#          for each patient.
# ------------------------------------------------------------

patient_tex <- patient_tex %>%
  left_join(patient_totals, by = c("PatientID", "Tissue")) %>%
  mutate(
    prevalence = tex_cells / total_cd8
  )


# ------------------------------------------------------------
# CHUNK 5 — Inspect resulting dataset
# Purpose: Quick sanity check of resulting table.
# ------------------------------------------------------------

head(patient_tex)

table(patient_tex$Tissue)


# ------------------------------------------------------------
# CHUNK 6 — Boxplot of patient-level prevalence
# Purpose: Visualise distribution of Tex prevalence across
#          tissue compartments.
# ------------------------------------------------------------

ggplot(patient_tex,
       aes(x = Tissue,
           y = prevalence,
           fill = TexClass)) +
  
  geom_boxplot(outlier.shape = NA,
               alpha = 0.7) +
  
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.8) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  
  labs(
    title = "Patient-level prevalence of Prolif-Tex and Non-Prolif-Tex",
    subtitle = "Each point represents one patient",
    x = "Tissue compartment",
    y = "Fraction of CD8 cells"
  ) +
  
  theme_classic()


# ------------------------------------------------------------
# CHUNK 7 — Optional: combine Tex populations
# Purpose: Evaluate total Tex prevalence regardless of
#          proliferative status.
# ------------------------------------------------------------

patient_tex_total <- df %>%
  filter(TexClass != "Unlabelled") %>%
  group_by(PatientID, Tissue) %>%
  summarise(
    tex_cells = n(),
    .groups = "drop"
  ) %>%
  left_join(patient_totals, by = c("PatientID", "Tissue")) %>%
  mutate(
    prevalence = tex_cells / total_cd8
  )


ggplot(patient_tex_total,
       aes(x = Tissue, y = prevalence)) +
  
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.8) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1)
  ) +
  
  labs(
    title = "Patient-level prevalence of exhausted CD8 T cells",
    subtitle = "Combined Prolif-Tex and Non-Prolif-Tex",
    x = "Tissue compartment",
    y = "Fraction of CD8 cells"
  ) +
  
  theme_classic()