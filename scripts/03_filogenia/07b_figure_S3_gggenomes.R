#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ============================================================
# Figure S3 with gggenomes (Hackl et al. 2024)
# ============================================================

suppressPackageStartupMessages({
  library(gggenomes)
  library(ggplot2)
})

IN_GENES <- "input_data/strain_metadata/gggenomes_genes.tsv"
IN_SEQS  <- "input_data/strain_metadata/gggenomes_seqs.tsv"
OUT_DIR  <- "outputs/supplementary_figures"
OUT_PDF  <- file.path(OUT_DIR, "Figure_S3_colocalization.pdf")
OUT_PNG  <- file.path(OUT_DIR, "Figure_S3_colocalization.png")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

genes <- read.delim(IN_GENES, sep = "\t", stringsAsFactors = FALSE)
seqs  <- read.delim(IN_SEQS,  sep = "\t", stringsAsFactors = FALSE)

strain_order <- unique(genes$bin_id)
seqs$bin_id  <- factor(seqs$bin_id,  levels = strain_order)
genes$bin_id <- factor(genes$bin_id, levels = strain_order)
genes$gene_id <- factor(genes$gene_id,
                        levels = c("pirA","pirB","g07720","g06662","g07221"))

gene_colors <- c(
  pirA   = "#EE0000",
  pirB   = "#8B0000",
  g07720 = "#3B4992",
  g06662 = "#5F559B",
  g07221 = "#008B45"
)

p <- gggenomes(genes = genes, seqs = seqs) +
  geom_seq() +
  geom_gene(aes(fill = gene_id), size = 4, color = "black") +
  geom_seq_label(aes(label = seq_id), size = 3.5) +
  scale_fill_manual(values = gene_colors, name = NULL) +
  theme_gggenomes_clean(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(OUT_PDF, p, width = 14, height = 9, dpi = 300)
ggsave(OUT_PNG, p, width = 14, height = 9, dpi = 300)

cat("\n=========================================\n")
cat("Figure S3 (gggenomes) generated\n")
cat("=========================================\n")
cat("PDF :", OUT_PDF, "\n")
cat("PNG :", OUT_PNG, "\n")
