#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Panel 2x2 de reconstrucción ancestral para los 4 genes con asociación perfecta a AHPND+

library(ape)
library(phangorn)

# Rutas
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")
pirate_tsv <- "PIRATE.gene_families.ordered.tsv"
tree_file <- "core_snp_phylogeny/iqtree_out/core_snp.treefile"

# Leer árbol y enraizar
tree <- read.tree(tree_file)
tree <- midpoint(tree)
tree <- makeNodeLabel(tree, prefix = "node")

cepas_arbol <- tree$tip.label

# Leer cabecera del TSV
pirate_header <- read.delim(pirate_tsv, sep = "\t", nrows = 1, header = FALSE,
                            stringsAsFactors = FALSE, check.names = FALSE)
cepas_cols <- match(cepas_arbol, pirate_header)

# Genes de interés (family IDs de PIRATE ya mapeados)
genes <- list(
  pirA       = "g07679",
  `g07720`   = "g07720",
  `g06662`   = "g06662",
  `g07221`   = "g07221"
)

# Preparar lista de estados ancestrales y datos
anc_list <- list()
estados_list <- list()

for (gen in names(genes)) {
  fam <- genes[[gen]]
  cat("Procesando", gen, "(", fam, ")...\n")

  # Buscar fila con awk (segunda columna == familia)
  cmd <- sprintf("awk -F'\\t' '$2==\"%s\"' %s", fam, pirate_tsv)
  linea <- system(cmd, intern = TRUE)
  if (length(linea) == 0) {
    cat("  No se encontró fila para", fam, "\n")
    next
  }
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

  fit <- ace(estado, tree, type = "discrete", method = "ML", model = "ARD")
  anc <- fit$lik.anc
  anc_list[[gen]] <- anc
  estados_list[[gen]] <- estado
}

# Generar panel 2x2 en SVG y PNG
svg("panel_ancestral_4genes.svg", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))

for (gen in names(anc_list)) {
  plot(tree, show.tip.label = TRUE, cex = 0.7, no.margin = TRUE)
  nodelabels(pie = anc_list[[gen]], piecol = c("white", "black"), cex = 0.5)
  tiplabels(pch = 21, bg = ifelse(estados_list[[gen]], "black", "white"), cex = 1.0)
  title(gen, line = 0.5)
}
dev.off()

png("panel_ancestral_4genes.png", width = 3600, height = 3000, res = 300)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
for (gen in names(anc_list)) {
  plot(tree, show.tip.label = TRUE, cex = 0.7, no.margin = TRUE)
  nodelabels(pie = anc_list[[gen]], piecol = c("white", "black"), cex = 0.5)
  tiplabels(pch = 21, bg = ifelse(estados_list[[gen]], "black", "white"), cex = 1.0)
  title(gen, line = 0.5)
}
dev.off()

cat("Panel guardado: panel_ancestral_4genes.svg / .png\n")
