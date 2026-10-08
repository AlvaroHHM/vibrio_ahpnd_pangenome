#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Reconstrucción ancestral automática de genes candidatos (SVG)
# Mapea secuencias Roary -> familias PIRATE usando BLASTp

library(ape)
library(phangorn)

# ---------- Rutas ----------
setwd("~/bacterial-genomics-tutorial")
pirate_dir <- "reference_genomes_pipeline/pirate_input_pirAB/PIRATE_results_v3"
pangenome_dir <- "pangenome"
output_dir <- file.path(pirate_dir, "ancestral_svg")
dir.create(output_dir, showWarnings = FALSE)

# Árbol core-SNP
tree_file <- file.path(pirate_dir, "core_snp_phylogeny/iqtree_out/core_snp.treefile")
tree <- read.tree(tree_file)
tree <- midpoint(tree)
tree <- makeNodeLabel(tree, prefix = "node")

# Archivos de secuencias Roary
excl_faa <- file.path(pangenome_dir, "exclusivos_15.faa")
aus_faa  <- file.path(pangenome_dir, "ausentes_15_presentes_11_6.faa")

# Base de datos PIRATE
feature_dir <- file.path(pirate_dir, "feature_sequences")
combined_fasta <- file.path(output_dir, "pirate_families.faa")
db_name <- file.path(output_dir, "pirate_fam_db")

# Crear FASTA combinado de PIRATE si no existe
if (!file.exists(combined_fasta)) {
  cat("Creando FASTA combinado de PIRATE...\n")
  archivos <- list.files(feature_dir, pattern = "\\.aa\\.fasta$", full.names = TRUE)
  con <- file(combined_fasta, "w")
  for (f in archivos) {
    fam <- sub("\\.aa\\.fasta$", "", basename(f))
    lines <- readLines(f)
    headers <- grep("^>", lines)
    for (i in seq_along(headers)) {
      h <- sub("^>", "", lines[headers[i]])
      seq_start <- headers[i] + 1
      seq_end <- ifelse(i < length(headers), headers[i+1] - 1, length(lines))
      seq <- paste(lines[seq_start:seq_end], collapse = "")
      writeLines(paste0(">", fam, "|", h), con)
      writeLines(seq, con)
    }
  }
  close(con)
}

# Crear BLASTdb si no existe
if (!file.exists(paste0(db_name, ".phr"))) {
  cat("Creando base de datos BLASTp...\n")
  system(sprintf("makeblastdb -in %s -dbtype prot -out %s", shQuote(combined_fasta), shQuote(db_name)))
}

# ---------- Función para extraer secuencia desde FASTA Roary ----------
leer_secuencia <- function(archivo, patron) {
  lines <- readLines(archivo)
  headers <- grep("^>", lines)
  for (i in seq_along(headers)) {
    h <- sub("^>", "", lines[headers[i]])
    if (grepl(patron, h, ignore.case = TRUE)) {
      seq_start <- headers[i] + 1
      seq_end <- ifelse(i < length(headers), headers[i+1] - 1, length(lines))
      seq <- paste(lines[seq_start:seq_end], collapse = "")
      return(seq)
    }
  }
  return(NULL)
}

# ---------- Función para mapear secuencia a familia PIRATE ----------
blast_a_familia <- function(query_seq, query_id) {
  qfile <- tempfile(pattern = "query_", fileext = ".faa")
  writeLines(paste0(">", query_id, "\n", query_seq), qfile)
  cmd <- sprintf(
    "blastp -query %s -db %s -outfmt \"6 qseqid sseqid pident qcovs evalue bitscore\" -evalue 1e-5 -max_target_seqs 1",
    shQuote(qfile), shQuote(db_name)
  )
  out <- system(cmd, intern = TRUE)
  unlink(qfile)
  if (length(out) == 0) return(NULL)
  res <- strsplit(out[1], "\t")[[1]]
  pident <- as.numeric(res[3])
  qcovs <- as.numeric(res[4])
  if (pident >= 50 && qcovs >= 50) {
    return(sub("\\|.*", "", res[2]))  # family ID
  }
  return(NULL)
}

# ---------- Genes de interés (patrón en header Roary) ----------
genes_interes <- c(
  "pirA" = "pirA\\|",
  "pirB" = "pirB\\|",
  "group_1485" = "group_1485\\|",
  "group_1486" = "group_1486\\|",
  "group_1487" = "group_1487\\|"
)

# ---------- Leer cabecera PIRATE ----------
pirate_tsv <- file.path(pirate_dir, "PIRATE.gene_families.ordered.tsv")
pirate_header <- read.delim(pirate_tsv, sep = "\t", nrows = 1, header = FALSE,
                            stringsAsFactors = FALSE, check.names = FALSE)
cepas_arbol <- tree$tip.label
cepas_cols <- match(cepas_arbol, pirate_header)

# ---------- Procesar cada gen ----------
for (gen in names(genes_interes)) {
  patron <- genes_interes[gen]
  cat("Procesando", gen, "...\n")

  # Extraer secuencia de Roary (buscar en exclusivos y ausentes)
  seq <- leer_secuencia(excl_faa, patron)
  if (is.null(seq)) seq <- leer_secuencia(aus_faa, patron)
  if (is.null(seq)) {
    cat("  No se encontró secuencia para", gen, "\n")
    next
  }

  # Mapear a familia PIRATE
  fam <- blast_a_familia(seq, gen)
  if (is.null(fam)) {
    cat("  No se pudo mapear a PIRATE\n")
    next
  }
  cat("  Familia PIRATE:", fam, "\n")

  # Extraer presencia/ausencia desde PIRATE TSV
  cmd <- sprintf("awk -F'\\t' '$2==\"%s\"' %s", fam, pirate_tsv)
  linea <- system(cmd, intern = TRUE)
  if (length(linea) == 0) {
    cat("  No se encontró fila de la familia en PIRATE\n")
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

  # Reconstrucción ancestral
  fit <- ace(estado, tree, type = "discrete", method = "ML", model = "ARD")
  anc <- fit$lik.anc

  # Guardar SVG
  svg_file <- file.path(output_dir, paste0("ancestral_", gen, ".svg"))
  svg(svg_file, width = 10, height = 8)
  plot(tree, show.tip.label = TRUE, cex = 0.8, no.margin = TRUE)
  nodelabels(pie = anc, piecol = c("white", "black"), cex = 0.6)
  tiplabels(pch = 21, bg = ifelse(estado, "black", "white"), cex = 1.5)
  title(paste("Ancestral reconstruction of", gen))
  dev.off()
  cat("  SVG guardado:", svg_file, "\n")
}
