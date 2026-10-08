#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Curvas de acumulación del pangenoma con aleatorización

# Cargar matriz de presencia/ausencia
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")
pres <- read.delim("PIRATE.gene_families.ordered.tsv", sep = "\t", quote = "", comment.char = "",
                   stringsAsFactors = FALSE, check.names = FALSE)
pheno <- read.csv("phenotype.csv", stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")
cepas <- pheno$strain

# Matriz binaria (genes x cepas)
mat <- pres[, cepas]
mat_bin <- ifelse(mat == "", 0, 1)
rownames(mat_bin) <- pres$gene_family

# Función para calcular pangenoma/core dado un orden
calc_curves <- function(orden) {
  pg <- numeric(length(orden))
  core <- numeric(length(orden))
  genes_acumulados <- c()
  genes_core <- rep(TRUE, nrow(mat_bin))
  for (i in seq_along(orden)) {
    genes_presentes <- rownames(mat_bin)[mat_bin[, orden[i]] == 1]
    genes_acumulados <- union(genes_acumulados, genes_presentes)
    pg[i] <- length(genes_acumulados)
    genes_core <- genes_core & (mat_bin[, orden[i]] == 1)
    core[i] <- sum(genes_core)
  }
  list(pg = pg, core = core)
}

# Orden original (según fenotipo)
orden_orig <- cepas
res_orig <- calc_curves(orden_orig)

# Aleatorizaciones
n_perm <- 100
pg_matrix <- matrix(NA, nrow = length(cepas), ncol = n_perm)
core_matrix <- matrix(NA, nrow = length(cepas), ncol = n_perm)
for (i in 1:n_perm) {
  orden_rand <- sample(cepas)
  res <- calc_curves(orden_rand)
  pg_matrix[, i] <- res$pg
  core_matrix[, i] <- res$core
}

# Calcular medianas y cuantiles
pg_med <- apply(pg_matrix, 1, median)
pg_li <- apply(pg_matrix, 1, quantile, 0.025)
pg_ls <- apply(pg_matrix, 1, quantile, 0.975)
core_med <- apply(core_matrix, 1, median)
core_li <- apply(core_matrix, 1, quantile, 0.025)
core_ls <- apply(core_matrix, 1, quantile, 0.975)

# Graficar
png("curvas_acumulacion_aleatorizadas.png", width = 3000, height = 2000, res = 300)
plot(1:length(cepas), res_orig$pg, type = "l", col = "blue", ylim = c(0, max(pg_ls)),
     xlab = "Número de genomas", ylab = "Número de familias génicas", lwd = 2)
lines(1:length(cepas), pg_med, col = "gray", lwd = 2)
lines(1:length(cepas), pg_li, col = "gray", lty = 2)
lines(1:length(cepas), pg_ls, col = "gray", lty = 2)
lines(1:length(cepas), res_orig$core, col = "red", lwd = 2)
lines(1:length(cepas), core_med, col = "gray", lwd = 2)
lines(1:length(cepas), core_li, col = "gray", lty = 2)
lines(1:length(cepas), core_ls, col = "gray", lty = 2)
legend("topleft", legend = c("Pangenoma (orden original)", "Pangenoma (mediana permutada)", "Core (orden original)", "Core (mediana permutada)"),
       col = c("blue", "gray", "red", "gray"), lty = c(1,1,1,2), lwd = 2)
dev.off()

cat("Figura guardada: curvas_acumulacion_aleatorizadas.png\n")
