#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Análisis de sensibilidad: Fisher excluyendo L2181

setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")

# Leer fenotipo
pheno <- read.csv("phenotype.csv", stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")

# Excluir L2181
pheno_sens <- pheno[pheno$strain != "L2181", ]

# Leer PIRATE TSV
pirate <- read.delim("PIRATE.gene_families.ordered.tsv", sep = "\t",
                     quote = "", comment.char = "", stringsAsFactors = FALSE,
                     check.names = FALSE)

# Identificar columnas de cepas en la lista filtrada
cepas_sens <- pheno_sens$strain
cepas_cols <- intersect(cepas_sens, colnames(pirate))
if (length(cepas_cols) != length(cepas_sens)) {
  warning("No todas las cepas están en PIRATE")
}

# Construir matriz de presencia/ausencia
presence <- pirate[, cepas_cols, drop = FALSE]
presence_bin <- ifelse(presence == "", 0, 1)
rownames(presence_bin) <- pirate$gene_family

# Preparar grupos
AHPND_pos <- pheno_sens$AHPND == 1
AHPND_neg <- pheno_sens$AHPND == 0

# Test exacto de Fisher para cada gen
n_genes <- nrow(presence_bin)
resultados_sens <- data.frame(
  gene_family = rownames(presence_bin),
  present_AHPND_pos = integer(n_genes),
  absent_AHPND_pos = integer(n_genes),
  present_AHPND_neg = integer(n_genes),
  absent_AHPND_neg = integer(n_genes),
  odds_ratio = numeric(n_genes),
  p_value = numeric(n_genes),
  stringsAsFactors = FALSE
)

for (i in seq_len(n_genes)) {
  vec <- as.numeric(presence_bin[i, ])
  a <- sum(vec[AHPND_pos] == 1)
  b <- sum(vec[AHPND_pos] == 0)
  c <- sum(vec[AHPND_neg] == 1)
  d <- sum(vec[AHPND_neg] == 0)

  resultados_sens$present_AHPND_pos[i] <- a
  resultados_sens$absent_AHPND_pos[i] <- b
  resultados_sens$present_AHPND_neg[i] <- c
  resultados_sens$absent_AHPND_neg[i] <- d

  if ((a+b) > 0 && (c+d) > 0 && (a+c) > 0 && (b+d) > 0) {
    ft <- fisher.test(matrix(c(a,b,c,d), nrow = 2))
    resultados_sens$odds_ratio[i] <- ft$estimate
    resultados_sens$p_value[i] <- ft$p.value
  } else {
    resultados_sens$odds_ratio[i] <- ifelse(a+b > 0 && c+d > 0, Inf, 0)
    resultados_sens$p_value[i] <- 1
  }
}

# FDR
resultados_sens$FDR <- p.adjust(resultados_sens$p_value, method = "BH")

# Guardar resultados completos
write.table(resultados_sens, "sensitivity_excluding_L2181.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Filtrar genes con FDR < 0.05
sig_fdr_sens <- resultados_sens[resultados_sens$FDR < 0.05, ]
write.table(sig_fdr_sens, "sensitivity_excluding_L2181_FDR005.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("Análisis de sensibilidad completado.\n")
cat("Genes con FDR < 0.05 excluyendo L2181:", nrow(sig_fdr_sens), "\n")
print(sig_fdr_sens[, c("gene_family", "odds_ratio", "p_value", "FDR")])
