#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Control de estructura poblacional en asociación gen-fenotipo
# Regresión filogenética (PGLS) con ape + nlme + phangorn

# Cargar librerías
library(ape)
library(nlme)
if (!requireNamespace("phangorn", quietly = TRUE)) {
  install.packages("phangorn", repos = "https://cloud.r-project.org")
}
library(phangorn)

# Rutas
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3")
tree_file <- "core_snp_phylogeny/iqtree_out/core_snp.treefile"
pirate_tsv <- "PIRATE.gene_families.ordered.tsv"
pheno_file <- "phenotype.csv"

# Leer árbol y enraizar
tree <- read.tree(tree_file)
tree <- midpoint(tree)                # phangorn::midpoint
tree <- makeNodeLabel(tree, prefix = "node")   # opcional para evitar conflictos

# Leer fenotipos
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")
pheno$AHPND <- as.numeric(pheno$AHPND)

# Leer cabecera PIRATE y localizar columnas de cepas
pirate_header <- read.delim(pirate_tsv, sep = "\t", nrows = 1, header = FALSE,
                            stringsAsFactors = FALSE, check.names = FALSE)
cepas_arbol <- tree$tip.label
cepas_cols <- match(cepas_arbol, pirate_header)

# Preparar fenotipo en el orden del árbol
fen <- pheno[match(cepas_arbol, pheno$strain), "AHPND"]
names(fen) <- cepas_arbol

# Genes de interés
genes_top <- c("g07679", "g07720", "g06662", "g07221", "g05919")
resultados <- data.frame()

for (fam in genes_top) {
  # Buscar fila con awk (segunda columna == familia)
  cmd <- sprintf("awk -F'\\t' '$2==\"%s\"' %s", fam, pirate_tsv)
  linea <- system(cmd, intern = TRUE)
  if (length(linea) == 0) next

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
  df <- data.frame(AHPND = fen, gene = factor(estado))
  rownames(df) <- names(fen)

  # PGLS con correlación de Pagel
  fit <- gls(AHPND ~ gene,
             correlation = corPagel(0.5, phy = tree, fixed = FALSE),
             data = df,
             method = "ML")

  pval <- summary(fit)$tTable["gene1", "p-value"]
  coef <- summary(fit)$tTable["gene1", "Value"]

  resultados <- rbind(resultados,
                      data.frame(familia = fam, coef = coef, p_valor = pval))
}

# Guardar
write.table(resultados, "control_filogenetico_R.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(resultados)
