#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ============================================================
# Figure S4: Top-50 filtered gene fingerprint (Roary pattern)
# ============================================================
#
# NOTE: This script requires feature_sequences/ from PIRATE_results_v3,
# which is not included in this repository (47 MB, 10,424 files).
# The script performs BLASTp mapping of Roary gene probes against
# PIRATE feature sequences to build the presence/absence matrix.
#
# External input required:
#   ~/AHPND_Vp_study/02_pangenome/PIRATE_results_v3/feature_sequences/
#
# Inputs (from repo):
#   - input_data/pangenome/exclusivos_15.tsv
#   - input_data/pangenome/ausentes_15_presentes_11_6.tsv
#   - input_data/pirate/PIRATE.gene_families.ordered.tsv
#   - input_data/strain_metadata/phenotype.csv
#
# Outputs:
#   - outputs/supplementary_figures/Figure_S4_top50_heatmap.pdf
#   - outputs/supplementary_figures/Figure_S4_top50_heatmap.png
#   - outputs/huella/huella_pirate_filtrada_R.tsv
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
})

# ---------- Paths (relative to repo root) ----------
pirate_tsv  <- "input_data/pirate/PIRATE.gene_families.ordered.tsv"
pheno_file  <- "input_data/strain_metadata/phenotype.csv"
excl_faa    <- "input_data/pangenome/exclusivos_15.faa"
aus_faa     <- "input_data/pangenome/ausentes_15_presentes_11_6.faa"

# External (not in repo) — see header
feature_dir <- path.expand("~/AHPND_Vp_study/02_pangenome/PIRATE_results_v3/feature_sequences")

# Output dirs
output_dir <- "outputs/huella"
figure_dir <- "outputs/supplementary_figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# Sanity checks
for (f in c(pirate_tsv, pheno_file)) {
  if (!file.exists(f)) stop("Missing: ", f)
}
if (!dir.exists(feature_dir)) {
  stop("Missing external input: ", feature_dir,
       "\n  This script requires feature_sequences/ from PIRATE_results_v3.")
}

# ---------- 1. Build combined PIRATE FASTA ----------
combined_fasta <- file.path(output_dir, "pirate_families.faa")
if (!file.exists(combined_fasta)) {
  cat("Concatenating PIRATE feature sequences...\n")
  archivos <- list.files(feature_dir, pattern = "\\.aa\\.fasta$", full.names = TRUE)
  con <- file(combined_fasta, "w")
  for (f in archivos) {
    fam <- sub("\\.aa\\.fasta$", "", basename(f))
    lines <- readLines(f)
    headers <- grep("^>", lines)
    for (i in seq_along(headers)) {
      header <- sub("^>", "", lines[headers[i]])
      seq_start <- headers[i] + 1
      seq_end <- ifelse(i < length(headers), headers[i+1] - 1, length(lines))
      seq <- paste(lines[seq_start:seq_end], collapse = "")
      writeLines(paste0(">", fam, "|", header), con)
      writeLines(seq, con)
    }
  }
  close(con)
}

# ---------- 2. BLAST database ----------
db_name <- file.path(output_dir, "pirate_fam_db")
if (!file.exists(paste0(db_name, ".phr"))) {
  cat("Building BLAST database...\n")
  system(sprintf("makeblastdb -in %s -dbtype prot -out %s",
                 shQuote(combined_fasta), shQuote(db_name)))
}

# ---------- 3. Read Roary probes ----------
leer_sondas <- function(faa_file, label) {
  lines <- readLines(faa_file)
  headers <- grep("^>", lines)
  sondas <- list()
  for (i in seq_along(headers)) {
    h <- sub("^>", "", lines[headers[i]])
    gene_id <- sub("\\|.*", "", h)
    gene_id <- sub(" .*", "", gene_id)
    seq_start <- headers[i] + 1
    seq_end <- ifelse(i < length(headers), headers[i+1] - 1, length(lines))
    seq <- paste(lines[seq_start:seq_end], collapse = "")
    sondas[[gene_id]] <- list(seq = seq, etiqueta = label,
                              annotation = sub(".*annotation=", "", h))
  }
  sondas
}

sondas_excl <- leer_sondas(excl_faa, "Exclusive_15")
sondas_aus  <- leer_sondas(aus_faa,  "Absent_15")
sondas <- c(sondas_excl, sondas_aus)
cat("Probes loaded:", length(sondas), "\n")

# ---------- 4. BLAST each probe against PIRATE DB ----------
blast_sonda <- function(seq, id) {
  qfile <- tempfile(pattern = "q_", fileext = ".faa")
  writeLines(paste0(">", id, "\n", seq), qfile)
  cmd <- sprintf(
    "blastp -query %s -db %s -outfmt '6 qseqid sseqid pident qcovs evalue bitscore' -evalue 1e-5 -max_target_seqs 1",
    shQuote(qfile), shQuote(db_name))
  out <- suppressWarnings(system(cmd, intern = TRUE))
  unlink(qfile)
  if (length(out) == 0) return(NULL)
  r <- strsplit(out[1], "\t")[[1]]
  if (as.numeric(r[3]) >= 50 && as.numeric(r[4]) >= 50) {
    return(sub("\\|.*", "", r[2]))
  }
  NULL
}

cat("Mapping probes to PIRATE families...\n")
mapeo <- list()
for (gene in names(sondas)) {
  fam <- blast_sonda(sondas[[gene]]$seq, gene)
  if (!is.null(fam)) mapeo[[gene]] <- fam
}
cat("Mapped probes:", length(mapeo), "\n")

# ---------- 5. Presence/absence matrix ----------
pirate <- read.delim(pirate_tsv, header = TRUE, sep = "\t",
                     quote = "", comment.char = "",
                     stringsAsFactors = FALSE, check.names = FALSE)
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")
pheno$strain <- gsub("-", "_", pheno$strain)
cepas <- pheno$strain
cepas_presentes <- intersect(cepas, colnames(pirate))

matriz <- matrix(0, nrow = length(mapeo), ncol = length(cepas_presentes),
                 dimnames = list(names(mapeo), cepas_presentes))
for (gene in names(mapeo)) {
  fam <- mapeo[[gene]]
  idx <- which(pirate$gene_family == fam)
  if (length(idx) == 1) {
    for (cepa in cepas_presentes) {
      val <- pirate[idx, cepa]
      if (!is.na(val) && val != "") matriz[gene, cepa] <- 1
    }
  }
}

# ---------- 6. Filter by Roary pattern ----------
df <- as.data.frame(matriz)
df$Set <- sapply(names(mapeo), function(g) sondas[[g]]$etiqueta)
df$annotation <- sapply(names(mapeo), function(g) sondas[[g]]$annotation)

mask_excl <- (df$Set == "Exclusive_15") & (df[["15_CESAIBC"]] == 1) &
             (df[["11_VM"]] == 0) & (df[["6_VM"]] == 0)
mask_aus  <- (df$Set == "Absent_15")   & (df[["15_CESAIBC"]] == 0) &
             (df[["11_VM"]] == 1) & (df[["6_VM"]] == 1)
df_filtered <- df[mask_excl | mask_aus, ]

# ---------- 7. Score and top 50 ----------
AHPND_pos <- intersect(pheno$strain[pheno$AHPND == 1], colnames(matriz))
AHPND_neg <- intersect(pheno$strain[pheno$AHPND == 0], colnames(matriz))
cat("AHPND+ strains in matrix:", length(AHPND_pos), "\n")
cat("AHPND- strains in matrix:", length(AHPND_neg), "\n")
score_pos <- rowMeans(df_filtered[, AHPND_pos, drop = FALSE]) * 100
score_neg <- rowMeans(df_filtered[, AHPND_neg, drop = FALSE]) * 100
df_filtered$score <- score_pos - score_neg
df_filtered <- df_filtered[order(-abs(df_filtered$score)), ]

genes_seleccionados <- head(rownames(df_filtered), 50)
for (g in c("pirA", "pirB")) {
  if (g %in% rownames(df_filtered) && !(g %in% genes_seleccionados)) {
    genes_seleccionados <- c(genes_seleccionados, g)
  }
}
top_genes <- df_filtered[genes_seleccionados, , drop = FALSE]

# ---------- 8. Heatmap data ----------
ordered_cols <- c(AHPND_pos, AHPND_neg)
heat_data <- as.matrix(top_genes[, ordered_cols, drop = FALSE])
heat_data <- apply(heat_data, 2, as.numeric)
rownames(heat_data) <- rownames(top_genes)

ann_col <- data.frame(AHPND = ifelse(ordered_cols %in% AHPND_pos,
                                     "Positive", "Negative"))
rownames(ann_col) <- ordered_cols
ann_colors <- list(AHPND = c(Positive = "#EE0000", Negative = "#BBBBBB"))

ann_row <- data.frame(Set = top_genes$Set)
rownames(ann_row) <- rownames(top_genes)
ann_row_colors <- list(Set = c(Exclusive_15 = "#3B4992", Absent_15 = "#EE0000"))

etiquetas_fila <- ifelse(nchar(top_genes$annotation) > 40,
                         paste0(substr(top_genes$annotation, 1, 37), "..."),
                         top_genes$annotation)
etiquetas_fila <- paste0(etiquetas_fila, " (", rownames(top_genes), ")")

# ---------- 9. Plot: PDF + PNG (no title) ----------
plot_heatmap <- function() {
  pheatmap(heat_data,
           color = colorRampPalette(c("white", "#3B4992"))(50),
           annotation_col = ann_col,
           annotation_row = ann_row,
           annotation_colors = c(ann_colors, ann_row_colors),
           labels_row = etiquetas_fila,
           fontsize_row = 6,
           fontsize_col = 8,
           angle_col = 45,
           main = "",
           legend = FALSE,
           silent = FALSE)
}

pdf(file.path(figure_dir, "Figure_S4_top50_heatmap.pdf"),
    width = 12, height = 8)
plot_heatmap()
dev.off()

png(file.path(figure_dir, "Figure_S4_top50_heatmap.png"),
    width = 12 * 200, height = 8 * 200, res = 200)
plot_heatmap()
dev.off()

# ---------- 10. Save filtered matrix ----------
write.table(top_genes, file.path(output_dir, "huella_pirate_filtrada_R.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("\n=========================================\n")
cat("Figure S4 generated\n")
cat("=========================================\n")
cat("PDF :", file.path(figure_dir, "Figure_S4_top50_heatmap.pdf"), "\n")
cat("PNG :", file.path(figure_dir, "Figure_S4_top50_heatmap.png"), "\n")
cat("TSV :", file.path(output_dir, "huella_pirate_filtrada_R.tsv"), "\n")
