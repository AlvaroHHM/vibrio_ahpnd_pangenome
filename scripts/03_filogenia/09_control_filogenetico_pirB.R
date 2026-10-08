#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Control filogenético solo para pirB (g05919)

library(ape)
library(phangorn)
if (!requireNamespace("phylolm", quietly = TRUE)) {
  install.packages("phylolm", repos = "https://cloud.r-project.org")
}
library(phylolm)

# Rutas
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")
tree_file <- "core_snp_phylogeny/iqtree_out/core_snp.treefile"
pirate_tsv <- "PIRATE.gene_families.ordered.tsv"
pheno_file <- "phenotype.csv"

# Leer árbol y enraizar
tree <- read.tree(tree_file)
tree <- midpoint(tree)
tree <- makeNodeLabel(tree, prefix = "node")

# Leer fenotipos
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")
pheno$AHPND <- as.numeric(pheno$AHPND)

# Leer cabecera PIRATE y localizar columnas de cepas
pirate_header <- read.delim(pirate_tsv, sep = "\t", nrows = 1, header = FALSE,
                            stringsAsFactors = FALSE, check.names = FALSE)
cepas_arbol <- tree$tip.label
cepas_cols <- match(cepas_arbol, pirate_header)

# Fenotipo en el orden del árbol
fen <- pheno[match(cepas_arbol, pheno$strain), "AHPND"]
names(fen) <- cepas_arbol

# Obtener presencia/ausencia de pirB (g05919)
fam <- "g05919"
cmd <- sprintf("awk -F'\\t' '$2==\"%s\"' %s", fam, pirate_tsv)
linea <- system(cmd, intern = TRUE)
if (length(linea) == 0) stop("No se encontró g05919 en PIRATE")

partes <- strsplit(linea[1], "\t")[[1]]

estado <- numeric(length(cepas_arbol))
names(estado) <- cepas_arbol
for (i in seq_along(cepas_arbol)) {
  val <- partes[cepas_cols[i]]
  estado[i] <- ifelse(is.na(val) | val == "", 0, 1)
}
estado <- estado[tree$tip.label]
estado <- as.numeric(estado)
names(estado) <- tree$tip.label

# Data frame
df <- data.frame(AHPND = fen, gene = estado)
rownames(df) <- names(fen)

cat("Ajustando modelo phylolm para pirB...\n")
fit <- tryCatch(
  phylolm(AHPND ~ gene, data = df, phy = tree, model = "logistic_MPLE"),
  error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    return(NULL)
  }
)

if (is.null(fit)) {
  cat("No se pudo ajustar el modelo. Posible separación completa o falta de variación.\n")
  writeLines("familia\tcoef\tp_valor\tnota\ng05919\tNA\tNA\tError/separación perfecta", "control_filogenetico_pirB_R.tsv")
} else {
  coef <- summary(fit)$coefficients["gene", "Estimate"]
  pval <- summary(fit)$coefficients["gene", "p.value"]
  cat("Coeficiente:", coef, "\n")
  cat("p-valor:", pval, "\n")
  writeLines(sprintf("familia\tcoef\tp_valor\tnota\ng05919\t%.4f\t%.4e\t", coef, pval), "control_filogenetico_pirB_R.tsv")
}
