#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Análisis de pangenoma y asociación fenotipo-gen en R
# Entrada: PIRATE_results_v3
# Salida: carpeta R_analysis/ con resultados

# ---------- Configuración de rutas ----------
base_dir <- "~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3"
setwd(base_dir)

pheno_file <- "phenotype.csv"
pirate_tsv <- "PIRATE.gene_families.ordered.tsv"

# Crear carpeta de resultados
output_dir <- file.path(base_dir, "R_analysis")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Carpeta creada:", output_dir, "\n")
} else {
  cat("Carpeta ya existe:", output_dir, "\n")
}

# ---------- 1. Leer fenotipo ----------
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")
pheno$AHPND <- ifelse(pheno$AHPND == 1, 1, 0)

# ---------- 2. Leer la tabla de PIRATE ----------
cat("Leyendo PIRATE.gene_families.ordered.tsv ...\n")
pirate <- read.delim(pirate_tsv, header = TRUE, sep = "\t", 
                     quote = "", comment.char = "", 
                     stringsAsFactors = FALSE, check.names = FALSE)

# Identificar las columnas de cepas
strain_cols <- intersect(colnames(pirate), pheno$strain)
if (length(strain_cols) != nrow(pheno)) {
  warning("No se encontraron todas las cepas del fenotipo en la tabla PIRATE")
}

# Crear matriz de presencia/ausencia
cat("Construyendo matriz de presencia/ausencia ...\n")
presence <- pirate[, strain_cols, drop = FALSE]
presence_bin <- ifelse(presence == "", 0, 1)
rownames(presence_bin) <- pirate$gene_family

# ---------- 3. Estadísticas de pangenoma ----------
n_strains <- ncol(presence_bin)
n_genes <- nrow(presence_bin)
gene_freq <- rowSums(presence_bin)

core_thresh <- 0.99 * n_strains
softcore_thresh <- 0.95 * n_strains
shell_thresh <- 0.15 * n_strains

core_genes <- sum(gene_freq >= core_thresh)
softcore_genes <- sum(gene_freq >= softcore_thresh & gene_freq < core_thresh)
shell_genes <- sum(gene_freq >= shell_thresh & gene_freq < softcore_thresh)
cloud_genes <- sum(gene_freq < shell_thresh)

cat("\nResumen del pangenoma:\n")
cat("  Total de familias:", n_genes, "\n")
cat("  Core (≥99%):", core_genes, "\n")
cat("  Soft-core (95-99%):", softcore_genes, "\n")
cat("  Shell (15-95%):", shell_genes, "\n")
cat("  Cloud (<15%):", cloud_genes, "\n")
cat("  Accesorio (shell+cloud):", shell_genes + cloud_genes, "\n")

# ---------- 4. Fisher exact test por familia ----------
cat("\nRealizando test exacto de Fisher ...\n")

pheno_ordered <- pheno[match(strain_cols, pheno$strain), ]
AHPND_pos <- pheno_ordered$AHPND == 1
AHPND_neg <- pheno_ordered$AHPND == 0

results <- data.frame(
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

  results$present_AHPND_pos[i] <- a
  results$absent_AHPND_pos[i] <- b
  results$present_AHPND_neg[i] <- c
  results$absent_AHPND_neg[i] <- d

  if ((a+b) > 0 && (c+d) > 0 && (a+c) > 0 && (b+d) > 0) {
    ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
    results$odds_ratio[i] <- ft$estimate
    results$p_value[i] <- ft$p.value
  } else {
    results$odds_ratio[i] <- ifelse(a+b > 0 && c+d > 0, Inf, 0)
    results$p_value[i] <- 1
  }
}

# Añadir FDR
results$FDR <- p.adjust(results$p_value, method = "BH")

# ---------- 5. Guardar resultados en R_analysis ----------
write.table(results, file.path(output_dir, "fisher_results_R.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

sig_nominal <- results[results$p_value < 0.05, ]
sig_extreme <- sig_nominal[is.infinite(sig_nominal$odds_ratio) | sig_nominal$odds_ratio == 0, ]

cat("\nGenes con p < 0.05 nominal:", nrow(sig_nominal), "\n")
cat("Genes con p < 0.05 y odds ratio extremo:", nrow(sig_extreme), "\n")
cat("Genes con FDR < 0.05:", sum(results$FDR < 0.05), "\n")

write.table(sig_nominal, file.path(output_dir, "significant_fisher_R.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sig_extreme, file.path(output_dir, "significant_extreme_R.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Guardar también un resumen de pangenoma
pangenome_summary <- data.frame(
  Statistic = c("Total families", "Core (>=99%)", "Soft-core (95-99%)",
                "Shell (15-95%)", "Cloud (<15%)", "Accessory (shell+cloud)"),
  Count = c(n_genes, core_genes, softcore_genes, shell_genes, cloud_genes,
            shell_genes + cloud_genes)
)
write.table(pangenome_summary, file.path(output_dir, "pangenome_summary_R.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResultados guardados en:", output_dir, "\n")
cat("  fisher_results_R.tsv\n")
cat("  significant_fisher_R.tsv\n")
cat("  significant_extreme_R.tsv\n")
cat("  pangenome_summary_R.tsv\n")
