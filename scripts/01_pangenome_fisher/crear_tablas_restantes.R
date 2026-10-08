#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Crear Table 2 (ganados y perdidos), Table S4 y Table S5

setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")

# ---------- 1. Table 2: Genes with extreme association (Gained and Lost) ----------
cat("Generando Table 2...\n")

# Leer resultados anotados de genes extremos (113 genes)
sig_ann <- read.delim("R_analysis/significant_extreme_annotated_R.tsv",
                      sep = "\t", stringsAsFactors = FALSE)

# Seleccionar columnas de interés
cols_t2 <- c("gene_family", "Description", "Preferred_name",
             "present_AHPND_pos", "absent_AHPND_pos",
             "present_AHPND_neg", "absent_AHPND_neg",
             "odds_ratio", "p_value", "FDR")

# --- Genes Gained (odds ratio Inf, presentes en AHPND+, ausentes en AHPND-) ---
gained <- sig_ann[is.infinite(sig_ann$odds_ratio) & sig_ann$FDR < 0.05, ]
gained <- gained[, cols_t2]
gained$Direction <- "Gained in AHPND+"
# Ordenar por p
gained <- gained[order(gained$p_value), ]

# --- Genes Lost (odds ratio 0, ausentes en AHPND+, presentes en AHPND-) ---
lost <- sig_ann[sig_ann$odds_ratio == 0 & sig_ann$p_value < 0.05, ]
# Seleccionar los 5 con menor p
lost <- lost[order(lost$p_value), ][1:min(5, nrow(lost)), ]
lost <- lost[, cols_t2]
lost$Direction <- "Lost in AHPND+"

# Combinar
table2 <- rbind(gained, lost)

# Renombrar columnas
colnames(table2) <- c("Gene", "Product", "Preferred_name",
                      "AHPND+_present", "AHPND+_absent",
                      "AHPND−_present", "AHPND−_absent",
                      "Odds_ratio", "p_value", "FDR", "Direction")

# Escribir
write.table(table2, "Table_2_summary.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Table 2 guardada: Table_2_summary.tsv\n")
cat("Contiene", nrow(gained), "genes ganados y", nrow(lost), "genes perdidos.\n")

# ---------- 2. Table S4: Validation of absences ----------
cat("Generando Table S4...\n")
tabla_s4 <- data.frame(
  Gene = c("vpadF", "hcp-2"),
  Strain = "15_CESAIBC",
  BLASTp_hit = "No",
  TBLASTN_hit = "No",
  Synteny_break = "Yes",
  Read_mapping = "Not done",
  PCR = "Not done",
  Conclusion = "Not detected"
)
write.table(tabla_s4, "Table_S4_validation_absences.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Table S4 guardada: Table_S4_validation_absences.tsv\n")

# ---------- 3. Table S5: Replicate-level mortality ----------
cat("Generando plantilla Table S5...\n")
plantilla_s5 <- data.frame(
  Strain = c("15_CESAIBC", "15_CESAIBC", "15_CESAIBC", "15_CESAIBC", "15_CESAIBC",
             "11_VM", "11_VM", "11_VM", "11_VM", "11_VM",
             "6_VM", "6_VM", "6_VM", "6_VM", "6_VM",
             "Control", "Control", "Control", "Control", "Control"),
  Replicate = rep(1:5, times = 4),
  Initial_n = 18,
  Deaths = NA_integer_,
  Mortality = NA_real_
)
write.table(plantilla_s5, "Table_S5_mortality_template.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Plantilla Table S5 guardada: Table_S5_mortality_template.tsv\n")
cat("Completa las columnas Deaths y Mortality con tus datos reales.\n")
