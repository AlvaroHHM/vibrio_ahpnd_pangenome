#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Calcular alpha de Heaps para el pangenoma

setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")

# Leer fenotipo y matriz de presencia/ausencia
pheno <- read.csv("phenotype.csv", stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")

pirate <- read.delim("PIRATE.gene_families.ordered.tsv", sep = "\t",
                     quote = "", comment.char = "", stringsAsFactors = FALSE,
                     check.names = FALSE)

cepas <- pheno$strain
cepas_cols <- intersect(cepas, colnames(pirate))

# Matriz binaria genes x cepas
presence <- pirate[, cepas_cols, drop = FALSE]
presence_bin <- ifelse(presence == "", 0, 1)
rownames(presence_bin) <- pirate$gene_family

# Función para calcular tamaño del pangenoma dado un orden
pangenome_size <- function(orden) {
  genes_acumulados <- c()
  tamanos <- numeric(length(orden))
  for (i in seq_along(orden)) {
    genes_presentes <- rownames(presence_bin)[presence_bin[, orden[i]] == 1]
    genes_acumulados <- union(genes_acumulados, genes_presentes)
    tamanos[i] <- length(genes_acumulados)
  }
  return(tamanos)
}

# Orden original (según fenotipo)
orden_original <- cepas_cols
tam_original <- pangenome_size(orden_original)

# Ajuste de Heaps: log(y) ~ log(k) + alpha * log(x)
x <- seq_along(tam_original)
y <- tam_original
# Evitar ceros
validos <- y > 0 & x > 0
modelo <- lm(log(y[validos]) ~ log(x[validos]))
alpha <- coef(modelo)[2]
alpha_ci <- confint(modelo, "log(x[validos])", level = 0.95)
k <- exp(coef(modelo)[1])

# Guardar resultado
resultado <- data.frame(
  alpha = alpha,
  alpha_lower = alpha_ci[1],
  alpha_upper = alpha_ci[2],
  k = k
)
write.table(resultado, "alpha_heaps.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("Alpha de Heaps estimado:", round(alpha, 3), "\n")
cat("IC 95% para alpha:", round(alpha_ci[1], 3), "-", round(alpha_ci[2], 3), "\n")
cat("Guardado en alpha_heaps.tsv\n")
