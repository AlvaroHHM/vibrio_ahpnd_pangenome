#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Figure 2: Recombination-filtered core-genome SNP phylogeny
# Node labels: SH-aLRT / UFBoot values (already in the .treefile)

library(ape)
library(phangorn)

# ---------- Rutas relativas ----------
tree_file <- "outputs/filogenia/core_snp.treefile"
out_dir   <- "outputs/figures"
out_pdf   <- file.path(out_dir, "Figure_2_core_snp_phylogeny.pdf")
out_png   <- file.path(out_dir, "Figure_2_core_snp_phylogeny.png")

if (!file.exists(tree_file)) {
  stop("No se encontró el árbol: ", tree_file)
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Leer y enraizar ----------
tree <- read.tree(tree_file)
tree <- midpoint(tree)

# ---------- Verificar labels ----------
cat("Nodos internos:", tree$Nnode, "\n")
cat("Labels en el árbol:", length(tree$node.label), "\n")
if (length(tree$node.label) > 0) {
  cat("Primeros 5 labels:", paste(head(tree$node.label, 5), collapse = " | "), "\n")
} else {
  cat("⚠ El árbol no tiene labels de soporte\n")
}

# ---------- Función de ploteo ----------
plot_tree <- function() {
  par(mar = c(2, 2, 3, 8), xpd = NA)
  plot(tree,
       show.tip.label = TRUE,
       cex = 0.9,
       no.margin = FALSE,
       align.tip.label = TRUE,
       main = "Core-genome SNP phylogeny (19 V. parahaemolyticus strains)")

  # Añadir labels de soporte en los nodos internos
  if (length(tree$node.label) > 0) {
    nodelabels(tree$node.label,
               frame = "none",
               cex = 0.7,
               col = "black",
               adj = c(1.2, 0.5))
  }

  add.scale.bar(length = 0.01, cex = 0.8)
}

# ---------- Guardar ----------
pdf(out_pdf, width = 8, height = 10)
plot_tree()
dev.off()

png(out_png, width = 2400, height = 3000, res = 300)
plot_tree()
dev.off()

cat("\n=========================================\n")
cat("Figure 2 generated\n")
cat("=========================================\n")
cat("PDF :", out_pdf, "\n")
cat("PNG :", out_png, "\n")
