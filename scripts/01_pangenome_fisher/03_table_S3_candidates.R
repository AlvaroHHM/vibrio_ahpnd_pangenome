#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Añadir IC95% solo a los 113 genes con p nominal < 0.05 y odds ratio extremo

setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/R_analysis")

# Leer resultados completos
fisher <- read.delim("fisher_results_R.tsv", sep = "\t", stringsAsFactors = FALSE)

# Filtrar genes con odds ratio extremo Y p < 0.05
candidatos <- fisher[is.finite(fisher$odds_ratio) & (is.infinite(fisher$odds_ratio) | fisher$odds_ratio == 0) & fisher$p_value < 0.05, ]

# Si hay odds ratio Inf, se maneja con is.infinite
candidatos <- fisher[(is.infinite(fisher$odds_ratio) | fisher$odds_ratio == 0) & fisher$p_value < 0.05, ]

# Calcular IC95%
ic_lower <- numeric(nrow(candidatos))
ic_upper <- numeric(nrow(candidatos))
for (i in seq_len(nrow(candidatos))) {
  a <- candidatos$present_AHPND_pos[i]
  b <- candidatos$absent_AHPND_pos[i]
  c <- candidatos$present_AHPND_neg[i]
  d <- candidatos$absent_AHPND_neg[i]
  if ((a+b) > 0 && (c+d) > 0 && (a+c) > 0 && (b+d) > 0) {
    ft <- fisher.test(matrix(c(a,b,c,d), nrow = 2))
    ic_lower[i] <- ft$conf.int[1]
    ic_upper[i] <- ft$conf.int[2]
  } else {
    ic_lower[i] <- NA
    ic_upper[i] <- NA
  }
}

candidatos$IC95_lower <- ic_lower
candidatos$IC95_upper <- ic_upper

# Seleccionar columnas finales
tabla_s2 <- candidatos[, c("gene_family", "present_AHPND_pos", "absent_AHPND_pos",
                           "present_AHPND_neg", "absent_AHPND_neg",
                           "odds_ratio", "IC95_lower", "IC95_upper",
                           "p_value", "FDR")]

# Guardar
write.table(tabla_s2, "Table_S2_candidatos.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("Tabla S2 corregida guardada con", nrow(tabla_s2), "genes candidatos.\n")
cat("Genes con FDR < 0.05:", sum(tabla_s2$FDR < 0.05), "\n")
