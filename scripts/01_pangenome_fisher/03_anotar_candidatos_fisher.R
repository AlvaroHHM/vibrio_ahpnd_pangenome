#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Integración de anotaciones eggNOG con genes candidatos de Fisher (19 cepas)
# Entrada: PIRATE_results_v3/R_analysis/significant_extreme_R.tsv
#          Galaxy9_113_annotations.tabular
#          Galaxy8_113_seed_orthologs.tabular
# Salida: R_analysis/significant_extreme_annotated_R.tsv

# ---------- Configuración ----------
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")
output_dir <- "R_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ---------- 1. Leer genes significativos extremos ----------
sig <- read.delim(file.path(output_dir, "significant_extreme_R.tsv"),
                  sep = "\t", stringsAsFactors = FALSE)
genes <- sig$gene_family

# ---------- 2. Leer anotaciones eggNOG ----------
cat("Leyendo anotaciones eggNOG...\n")
ann <- read.delim("Galaxy9_113_annotations.tabular", sep = "\t", quote = "",
                  comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
ortho <- read.delim("Galaxy8_113_seed_orthologs.tabular", sep = "\t", quote = "",
                    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

# Extraer ID del gen (antes del '|')
ann$gene <- sub("\\|.*", "", ann[["#query"]])
if ("#qseqid" %in% colnames(ortho)) {
  ortho$gene <- sub("\\|.*", "", ortho[["#qseqid"]])
} else if ("qseqid" %in% colnames(ortho)) {
  ortho$gene <- sub("\\|.*", "", ortho[["qseqid"]])
}

# Seleccionar columnas relevantes de anotación
cols_ann <- c("gene", "Description", "Preferred_name", "COG_category",
              "GOs", "EC", "KEGG_ko", "PFAMs")
cols_ann <- cols_ann[cols_ann %in% colnames(ann)]
ann_sel <- ann[, cols_ann, drop = FALSE]

# Ortólogo: quedarse con el mejor hit por gen (asumimos que el archivo ya está ordenado)
cols_ortho <- c("gene", "sseqid", "pident", "qcov", "evalue", "bitscore")
cols_ortho <- cols_ortho[cols_ortho %in% colnames(ortho)]
ortho_sel <- ortho[, cols_ortho, drop = FALSE]
# Para cada gen, tomar la primera fila
ortho_first <- ortho_sel[!duplicated(ortho_sel$gene), ]

# ---------- 3. Fusionar ----------
merged <- merge(sig, ann_sel, by.x = "gene_family", by.y = "gene", all.x = TRUE)
merged <- merge(merged, ortho_first, by.x = "gene_family", by.y = "gene", all.x = TRUE)

# Reordenar columnas: primero gene, conteos, OR, p, FDR, luego anotaciones
# (ya está bien, solo guardamos)

# ---------- 4. Guardar ----------
write.table(merged, file.path(output_dir, "significant_extreme_annotated_R.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Tabla anotada guardada en", file.path(output_dir, "significant_extreme_annotated_R.tsv"), "\n")
cat("Total genes:", nrow(merged), "\n")
cat("Genes con descripción:", sum(!is.na(merged$Description) & merged$Description != "-"), "\n")
