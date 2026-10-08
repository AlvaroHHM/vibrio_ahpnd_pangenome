#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Añadir anotaciones a los 15 genes FDR<0.05 del análisis de sensibilidad (excluyendo R13)

# ---------- Configuración ----------
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")

# Archivos de entrada
sens_file <- "sensitivity_excluding_R13.tsv"
pirate_tsv <- "PIRATE.gene_families.ordered.tsv"
egg_annot_file <- "Galaxy9_113_annotations.tabular"

# Leer resultados de sensibilidad
sens <- read.delim(sens_file, sep = "\t", stringsAsFactors = FALSE)

# Filtrar FDR < 0.05
sig <- sens[sens$FDR < 0.05, ]

# ---------- Obtener anotación de PIRATE ----------
cat("Leyendo PIRATE...\n")
pirate <- read.delim(pirate_tsv, sep = "\t", quote = "", comment.char = "",
                     stringsAsFactors = FALSE, check.names = FALSE)
# Columnas típicas: gene_family (col 2), consensus_product (col 4)
pirate_ann <- data.frame(
  gene_family = pirate$gene_family,
  PIRATE_product = pirate$consensus_product,
  stringsAsFactors = FALSE
)

# ---------- Obtener anotación de eggNOG ----------
cat("Leyendo eggNOG...\n")
egg <- read.delim(egg_annot_file, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE, check.names = FALSE)
# Extraer gene_id de #query (antes del '|')
egg$gene_id <- sub("\\|.*", "", egg[["#query"]])
egg_ann <- egg[, c("gene_id", "Description", "Preferred_name", "COG_category", "PFAMs")]

# ---------- Fusionar ----------
sig_ann <- merge(sig, pirate_ann, by = "gene_family", all.x = TRUE)
sig_ann <- merge(sig_ann, egg_ann, by.x = "gene_family", by.y = "gene_id", all.x = TRUE)

# Reordenar columnas
col_order <- c("gene_family", "present_AHPND_pos", "absent_AHPND_pos",
               "present_AHPND_neg", "absent_AHPND_neg", "odds_ratio",
               "p_value", "FDR", "PIRATE_product", "Description", 
               "Preferred_name", "COG_category", "PFAMs")
sig_ann <- sig_ann[, col_order]

# Guardar
write.table(sig_ann, "sensitivity_excluding_R13_FDR005_annotated.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Tabla anotada guardada: sensitivity_excluding_R13_FDR005_annotated.tsv\n")
print(sig_ann)
