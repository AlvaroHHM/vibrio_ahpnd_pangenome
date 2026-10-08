#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Visualización del árbol core-SNP con ape + phangorn

# Cargar librerías
library(ape)
if (!requireNamespace("phangorn", quietly = TRUE)) {
  install.packages("phangorn", repos = "https://cloud.r-project.org")
}
library(phangorn)

# Ruta al árbol
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/core_snp_phylogeny")
tree_file <- "iqtree_out/core_snp.treefile"

# Leer árbol
tree <- read.tree(tree_file)

# Enraizar en el punto medio
tree <- midpoint(tree)

# Guardar PNG
png("core_snp_tree_ape.png", width = 3000, height = 2000, res = 300)
plot(tree, show.tip.label = TRUE, cex = 0.8, no.margin = TRUE)
add.scale.bar(length = 0.01)
dev.off()

# Guardar SVG
if (requireNamespace("svglite", quietly = TRUE)) {
  library(svglite)
  svglite("core_snp_tree_ape.svg", width = 10, height = 8)
  plot(tree, show.tip.label = TRUE, cex = 0.8, no.margin = TRUE)
  add.scale.bar(length = 0.01)
  dev.off()
}

cat("Árbol guardado como core_snp_tree_ape.png\n")
