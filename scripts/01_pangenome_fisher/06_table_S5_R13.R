#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Sensibilidad excluyendo R13 (cepa AHPND- con pirB)

setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")

# Leer fenotipo
pheno <- read.csv("phenotype.csv", stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")
pheno_sens <- pheno[pheno$strain != "R13", ]

# Leer PIRATE TSV
pirate <- read.delim("PIRATE.gene_families.ordered.tsv", sep = "\t",
                     quote = "", comment.char = "", stringsAsFactors = FALSE,
                     check.names = FALSE)

# Seleccionar columnas de cepas
cepas_sens <- pheno_sens$strain
cepas_cols <- intersect(cepas_sens, colnames(pirate))

# Construir matriz binaria
presence <- pirate[, cepas_cols, drop = FALSE]
presence_bin <- ifelse(presence == "", 0, 1)
rownames(presence_bin) <- pirate$gene_family

# Grupos
AHPND_pos <- pheno_sens$AHPND == 1
AHPND_neg <- pheno_sens$AHPND == 0

# Test de Fisher para cada gen
n_genes <- nrow(presence_bin)
resultados <- data.frame(
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

  resultados$present_AHPND_pos[i] <- a
  resultados$absent_AHPND_pos[i] <- b
  resultados$present_AHPND_neg[i] <- c
  resultados$absent_AHPND_neg[i] <- d

  if ((a+b) > 0 && (c+d) > 0 && (a+c) > 0 && (b+d) > 0) {
    ft <- fisher.test(matrix(c(a,b,c,d), nrow = 2))
    resultados$odds_ratio[i] <- ft$estimate
    resultados$p_value[i] <- ft$p.value
  } else {
    resultados$odds_ratio[i] <- ifelse(a+b > 0 && c+d > 0, Inf, 0)
    resultados$p_value[i] <- 1
  }
}

# FDR
resultados$FDR <- p.adjust(resultados$p_value, method = "BH")

# Guardar
write.table(resultados, "sensitivity_excluding_R13.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Genes con FDR < 0.05
sig <- resultados[resultados$FDR < 0.05, ]
write.table(sig, "sensitivity_excluding_R13_FDR005.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("Análisis de sensibilidad excluyendo R13 completado.\n")
cat("Genes con FDR < 0.05:", nrow(sig), "\n")
print(sig[, c("gene_family", "odds_ratio", "p_value", "FDR")])
