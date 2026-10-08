#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Huella digital filtrada (Roary -> PIRATE) en R
# Entrada: 
#   - ~/bacterial-genomics-tutorial/pangenome/exclusivos_15.faa
#   - ~/bacterial-genomics-tutorial/pangenome/ausentes_15_presentes_11_6.faa
#   - ~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/feature_sequences/
#   - ~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/PIRATE.gene_families.ordered.tsv
#   - ~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3/phenotype.csv
# Salida: carpeta R_analysis_huella/

# ---------- Configuración ----------
setwd("~/bacterial-genomics-tutorial")
pirate_dir <- "reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3"
pangenome_dir <- "pangenome"

output_dir <- file.path(pangenome_dir, "R_analysis_huella")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------- 1. Preparar base de datos PIRATE (todas las secuencias de feature_sequences) ----------
feature_dir <- file.path(pirate_dir, "feature_sequences")
combined_fasta <- file.path(output_dir, "pirate_families.faa")

cat("Concatenando secuencias de feature_sequences ...\n")
archivos <- list.files(feature_dir, pattern = "\\.aa\\.fasta$", full.names = TRUE)
if (length(archivos) == 0) stop("No se encontraron archivos .aa.fasta en feature_sequences")

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

# Crear base de datos BLASTp
db_name <- file.path(output_dir, "pirate_fam_db")
if (!file.exists(paste0(db_name, ".phr"))) {
  cat("Creando base de datos BLASTp ...\n")
  system(sprintf("makeblastdb -in %s -dbtype prot -out %s", shQuote(combined_fasta), shQuote(db_name)))
} else {
  cat("Base de datos BLAST ya existe.\n")
}

# ---------- 2. Cargar secuencias sonda de Roary ----------
leer_sondas <- function(archivo, etiqueta) {
  lines <- readLines(archivo)
  headers <- grep("^>", lines)
  sondas <- list()
  for (i in seq_along(headers)) {
    h <- sub("^>", "", lines[headers[i]])
    if (grepl("\\|", h)) {
      parts <- strsplit(h, "\\|", fixed = TRUE)[[1]]
      gene_id <- trimws(parts[1])
      annotation <- trimws(paste(parts[-1], collapse = "|"))
    } else {
      gene_id <- trimws(h)
      annotation <- ""
    }
    seq_start <- headers[i] + 1
    seq_end <- ifelse(i < length(headers), headers[i+1] - 1, length(lines))
    seq <- paste(lines[seq_start:seq_end], collapse = "")
    sondas[[gene_id]] <- list(etiqueta = etiqueta, annotation = annotation, seq = seq)
  }
  return(sondas)
}

cat("Cargando secuencias sonda ...\n")
sondas_excl <- leer_sondas(file.path(pangenome_dir, "exclusivos_15.faa"), "Exclusive_15")
sondas_aus  <- leer_sondas(file.path(pangenome_dir, "ausentes_15_presentes_11_6.faa"), "Absent_15")
sondas <- c(sondas_excl, sondas_aus)
cat("Total sondas:", length(sondas), " (Exclusivas:", length(sondas_excl), ", Ausentes:", length(sondas_aus), ")\n")

# ---------- 3. BLASTp de cada sonda contra PIRATE ----------
blast_sonda <- function(query_seq, query_id) {
  qfile <- tempfile(pattern = "query_", fileext = ".faa")
  writeLines(paste0(">", query_id, "\n", query_seq), qfile)
  cmd <- sprintf("blastp -query %s -db %s -outfmt \"6 qseqid sseqid pident qcovs evalue bitscore\" -evalue 1e-5 -max_target_seqs 1",
                 shQuote(qfile), shQuote(db_name))
  out <- system(cmd, intern = TRUE)
  unlink(qfile)
  if (length(out) == 0) return(NULL)
  res <- strsplit(out[1], "\t")[[1]]
  pident <- as.numeric(res[3])
  qcovs <- as.numeric(res[4])
  evalue <- as.numeric(res[5])
  if (pident >= 50 && qcovs >= 50) {
    family_id <- sub("\\|.*", "", res[2])
    return(family_id)
  }
  return(NULL)
}

cat("Mapeando sondas a familias PIRATE ...\n")
mapeo <- list()
for (gene in names(sondas)) {
  fam <- blast_sonda(sondas[[gene]]$seq, gene)
  if (!is.null(fam)) {
    mapeo[[gene]] <- fam
  }
}
cat("Sondas mapeadas:", length(mapeo), "\n")

# ---------- 4. Leer presencia/ausencia desde PIRATE TSV ----------
cat("Leyendo PIRATE.gene_families.ordered.tsv ...\n")
pirate <- read.delim(file.path(pirate_dir, "PIRATE.gene_families.ordered.tsv"),
                     header = TRUE, sep = "\t", quote = "", comment.char = "",
                     stringsAsFactors = FALSE, check.names = FALSE)

pheno <- read.csv(file.path(pirate_dir, "phenotype.csv"), stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")
cepas <- pheno$strain
cepas_presentes <- intersect(cepas, colnames(pirate))
cat("Cepas encontradas:", length(cepas_presentes), "\n")

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

# ---------- 5. Filtrar por patrón Roary ----------
df <- as.data.frame(matriz)
df$Set <- sapply(names(mapeo), function(g) sondas[[g]]$etiqueta)
df$annotation <- sapply(names(mapeo), function(g) sondas[[g]]$annotation)

mask_excl <- (df$Set == "Exclusive_15") & (df[["15_CESAIBC"]] == 1) & (df[["11_VM"]] == 0) & (df[["6_VM"]] == 0)
mask_aus  <- (df$Set == "Absent_15")   & (df[["15_CESAIBC"]] == 0) & (df[["11_VM"]] == 1) & (df[["6_VM"]] == 1)
df_filtered <- df[mask_excl | mask_aus, ]

cat("Filas después del filtro Roary:", nrow(df_filtered), "\n")

# ---------- 6. Calcular score ----------
AHPND_pos <- pheno$strain[pheno$AHPND == 1]
AHPND_neg <- pheno$strain[pheno$AHPND == 0]

score_pos <- rowMeans(df_filtered[, AHPND_pos, drop = FALSE]) * 100
score_neg <- rowMeans(df_filtered[, AHPND_neg, drop = FALSE]) * 100
df_filtered$score <- score_pos - score_neg
df_filtered <- df_filtered[order(-abs(df_filtered$score)), ]

# ---------- 7. Seleccionar top 50 y añadir pirA/pirB ----------
top_n <- 50
genes_importantes <- c("pirA", "pirB")
genes_seleccionados <- head(rownames(df_filtered), top_n)
for (g in genes_importantes) {
  if (g %in% rownames(df_filtered) && !(g %in% genes_seleccionados)) {
    genes_seleccionados <- c(genes_seleccionados, g)
  }
}
top_genes <- df_filtered[genes_seleccionados, , drop = FALSE]

# ---------- 8. Heatmap ----------
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}
library(pheatmap)

ordered_cols <- c(AHPND_pos, AHPND_neg)
heat_data <- top_genes[, ordered_cols, drop = FALSE]
heat_data <- apply(heat_data, 2, as.numeric)
rownames(heat_data) <- rownames(top_genes)

ann_col <- data.frame(AHPND = ifelse(ordered_cols %in% AHPND_pos, "Positive", "Negative"))
rownames(ann_col) <- ordered_cols
ann_colors <- list(AHPND = c(Positive = "#EE0000", Negative = "#BBBBBB"))

ann_row <- data.frame(Set = top_genes$Set)
rownames(ann_row) <- rownames(top_genes)
ann_row_colors <- list(Set = c(Exclusive_15 = "#3B4992", Absent_15 = "#EE0000"))

etiquetas_fila <- ifelse(nchar(top_genes$annotation) > 40,
                         paste0(substr(top_genes$annotation, 1, 37), "..."),
                         top_genes$annotation)
etiquetas_fila <- paste0(etiquetas_fila, " (", rownames(top_genes), ")")

png(file.path(output_dir, "heatmap_top50_pirate_R.png"), width = 12, height = 8, units = "in", res = 300)
pheatmap(heat_data,
         color = colorRampPalette(c("white", "#3B4992"))(50),
         annotation_col = ann_col,
         annotation_row = ann_row,
         annotation_colors = c(ann_colors, ann_row_colors),
         labels_row = etiquetas_fila,
         fontsize_row = 6,
         fontsize_col = 8,
         angle_col = 45,
         main = "Top 50 genes filtrados (patrón Roary en trío mexicano)",
         legend = FALSE)
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  library(svglite)
  svglite::svglite(file.path(output_dir, "heatmap_top50_pirate_R.svg"), width = 12, height = 8)
  pheatmap(heat_data,
           color = colorRampPalette(c("white", "#3B4992"))(50),
           annotation_col = ann_col,
           annotation_row = ann_row,
           annotation_colors = c(ann_colors, ann_row_colors),
           labels_row = etiquetas_fila,
           fontsize_row = 6,
           fontsize_col = 8,
           angle_col = 45,
           main = "Top 50 genes filtrados (patrón Roary en trío mexicano)",
           legend = FALSE)
  dev.off()
}

# Guardar matriz filtrada
write.table(top_genes, file.path(output_dir, "huella_pirate_filtrada_R.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("Resultados guardados en", output_dir, "\n")
