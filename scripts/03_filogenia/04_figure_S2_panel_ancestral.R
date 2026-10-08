#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Figure S2: Ancestral state reconstruction panel (2x2)
# 4 genes with perfect association to AHPND+ (pirA, g07720, g06662, g07221)

library(ape)
library(phangorn)

# ---------- Rutas relativas ----------
tree_file  <- "outputs/filogenia/core_snp.treefile"
pirate_tsv <- "input_data/pirate/PIRATE.gene_families.ordered.tsv"
out_dir    <- "outputs/supplementary_figures"
out_pdf    <- file.path(out_dir, "Figure_S2_ancestral_reconstruction.pdf")
out_png    <- file.path(out_dir, "Figure_S2_ancestral_reconstruction.png")
out_stats  <- file.path(out_dir, "Figure_S2_ancestral_stats.tsv")

if (!file.exists(tree_file))  stop("Falta: ", tree_file)
if (!file.exists(pirate_tsv)) stop("Falta: ", pirate_tsv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Leer árbol ----------
tree <- read.tree(tree_file)
tree <- midpoint(tree)
tree <- makeNodeLabel(tree, prefix = "node")

cepas_arbol <- tree$tip.label

# ---------- Leer cabecera del TSV ----------
pirate_header <- read.delim(pirate_tsv, sep = "\t", nrows = 1, header = FALSE,
                            stringsAsFactors = FALSE, check.names = FALSE)
cepas_cols <- match(cepas_arbol, pirate_header)

# ---------- Genes de interés ----------
genes <- list(
  pirA     = "g07679",
  g07720   = "g07720",
  g06662   = "g06662",
  g07221   = "g07221"
)

# ---------- Calcular reconstrucciones ancestrales ----------
anc_list     <- list()
estados_list <- list()

for (gen in names(genes)) {
  fam <- genes[[gen]]
  cat("Procesando", gen, "(", fam, ")...\n")

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
  anc_list[[gen]]     <- fit$lik.anc
  estados_list[[gen]] <- estado
  cat("  OK\n")
}

# ---------- Cálculo de probabilidades ancestrales ----------
cat("\n=========================================\n")
cat("RESUMEN NUMÉRICO DE ESTADOS ANCESTRALES\n")
cat("=========================================\n")

stats_df <- data.frame(
  gene                  = character(),
  n_internal_nodes      = integer(),
  n_with_prob_present   = integer(),
  n_with_prob_absent    = integer(),
  max_prob_present      = numeric(),
  mean_prob_present     = numeric(),
  min_prob_present      = numeric(),
  stringsAsFactors      = FALSE
)

for (gen in names(anc_list)) {
  cat("\n", gen, ":\n", sep = "")
  # Columna 2 = probabilidad de "presente" (estado 1)
  probs <- anc_list[[gen]][, 2]
  n_total    <- length(probs)
  n_present  <- sum(probs > 0.5)
  n_absent   <- sum(probs <= 0.5)
  max_p      <- max(probs)
  mean_p     <- mean(probs)
  min_p      <- min(probs)

  cat("  Nodos internos totales:              ", n_total, "\n")
  cat("  Nodos con prob(presente) > 0.5:      ", n_present, "\n")
  cat("  Nodos con prob(presente) <= 0.5:     ", n_absent, "\n")
  cat("  Max prob(presente):                  ", round(max_p, 3), "\n")
  cat("  Media prob(presente):                ", round(mean_p, 3), "\n")
  cat("  Min prob(presente):                  ", round(min_p, 3), "\n")

  stats_df <- rbind(stats_df, data.frame(
    gene                = gen,
    n_internal_nodes    = n_total,
    n_with_prob_present = n_present,
    n_with_prob_absent  = n_absent,
    max_prob_present    = round(max_p, 4),
    mean_prob_present   = round(mean_p, 4),
    min_prob_present    = round(min_p, 4),
    stringsAsFactors    = FALSE
  ))
}

# Guardar tabla de estadísticas
write.table(stats_df, out_stats, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nTabla de estadísticas guardada en:", out_stats, "\n")

# ---------- Función de ploteo del panel 2x2 ----------
plot_panel <- function() {
  par(mfrow = c(2, 2), mar = c(2, 2, 3.5, 2), oma = c(0, 0, 2, 0))

  # Etiquetas A), B), C), D) en orden
  letters_labels <- c("A)", "B)", "C)", "D)")

  for (i in seq_along(names(anc_list))) {
    gen <- names(anc_list)[i]
    probs <- anc_list[[gen]][, 2]

    plot(tree, show.tip.label = TRUE, cex = 0.6, no.margin = FALSE)
    nodelabels(pie = anc_list[[gen]], piecol = c("white", "black"), cex = 0.45)
    tiplabels(pch = 21, bg = ifelse(estados_list[[gen]], "black", "white"), cex = 0.9)

    # Título con nombre del gen + media prob(presente)
    title(sprintf("%s  [mean P(present) = %.3f]",
                  gen, mean(probs)),
          line = 0.8, cex.main = 1.1)

    # Añadir etiqueta A), B), C), D)
    mtext(letters_labels[i], side = 3, adj = 0, line = 2.3,
          font = 2, cex = 1.6)
  }
}

# ---------- Guardar PDF ----------
pdf(out_pdf, width = 14, height = 12)
plot_panel()
dev.off()

# ---------- Preview PNG ----------
png(out_png, width = 4200, height = 3600, res = 300)
plot_panel()
dev.off()

cat("\n=========================================\n")
cat("Figure S2 generated\n")
cat("=========================================\n")
cat("PDF :", out_pdf, "\n")
cat("PNG :", out_png, "\n")
cat("Stats:", out_stats, "\n")
cat("Genes:", paste(names(anc_list), collapse = ", "), "\n")
