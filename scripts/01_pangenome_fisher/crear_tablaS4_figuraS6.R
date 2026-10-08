#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Crear Table S4 y Figure S6

# ---------- Configuración ----------
setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB")

pirate_tsv <- "PIRATE_results_v3/PIRATE.gene_families.ordered.tsv"
base_genomas <- "nuevos_genomas"

# Cepas AHPND+
cepas_ahpnd_pos <- c("15_CESAIBC", "T13P", "CGVP22", "L2181",
                     "13028_A3", "CGVP3", "CGVP8", "L2171")

# Genes de interés
genes_familias <- c(
  pirA   = "g07679",
  pirB   = "g05919",
  g07720 = "g07720",
  g06662 = "g06662",
  g07221 = "g07221"
)

# Función para encontrar archivo GFF
encontrar_gff <- function(cepa) {
  carpeta <- file.path(base_genomas, cepa)
  if (dir.exists(carpeta)) {
    archivos <- list.files(carpeta, pattern = "\\.gff$", full.names = TRUE)
    if (length(archivos) > 0) return(archivos[1])
  }
  carpeta2 <- file.path(base_genomas, paste0(cepa, "_prokka"))
  if (dir.exists(carpeta2)) {
    archivos <- list.files(carpeta2, pattern = "\\.gff$", full.names = TRUE)
    if (length(archivos) > 0) return(archivos[1])
  }
  return(NA)
}

# Función para extraer contig y posición de un locus_tag (número)
extraer_info_gff <- function(gff_file, locus_pirate) {
  if (is.na(gff_file)) return(NULL)
  partes <- strsplit(locus_pirate, "_")[[1]]
  numero <- tail(partes, n = 1)
  gff <- read.delim(gff_file, sep = "\t", header = FALSE,
                    comment.char = "#", stringsAsFactors = FALSE, quote = "")
  colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  patron <- paste0("locus_tag=", ".*", numero)
  idx <- grep(patron, gff$attributes)
  if (length(idx) == 0) return(NULL)
  return(list(
    contig = gff$seqid[idx[1]],
    start = gff$start[idx[1]],
    end = gff$end[idx[1]],
    strand = gff$strand[idx[1]]
  ))
}

# Leer cabecera PIRATE
pirate_header <- read.delim(pirate_tsv, sep = "\t", nrows = 1, header = FALSE,
                            stringsAsFactors = FALSE, check.names = FALSE)
col_cepas <- setNames(match(cepas_ahpnd_pos, pirate_header), cepas_ahpnd_pos)

# ---------- Table S4: validación de ausencias ----------
cat("Creando Table S4: validación de ausencias...\n")
tabla_s4 <- data.frame(
  Gene = c("vpadF", "hcp-2"),
  Cepa = "15_CESAIBC",
  BLASTp_hit = "No",
  TBLASTN_hit = "No",
  Synteny = "Not consistent",
  Read_mapping = "Not done",
  PCR = "Not done",
  Conclusion = "Not detected",
  stringsAsFactors = FALSE
)
write.table(tabla_s4, "Table_S4_validation_absences.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Table S4 guardada: Table_S4_validation_absences.tsv\n")

# ---------- Figure S6: co-localización ----------
cat("Creando Figure S6: co-localización...\n")
# Preparar datos para la figura
plot_data <- data.frame()
for (cepa in cepas_ahpnd_pos) {
  gff_file <- encontrar_gff(cepa)
  cat("  Procesando", cepa, ":", ifelse(is.na(gff_file), "GFF no encontrado", gff_file), "\n")
  if (is.na(gff_file)) next

  col_idx <- col_cepas[cepa]
  for (gen in names(genes_familias)) {
    fam <- genes_familias[gen]
    cmd <- sprintf("awk -F'\\t' '$2==\"%s\"' %s", fam, pirate_tsv)
    linea <- system(cmd, intern = TRUE)
    if (length(linea) == 0) next
    partes <- strsplit(linea[1], "\t")[[1]]
    locus_pirate <- partes[col_idx]
    if (is.na(locus_pirate) || locus_pirate == "") next

    info <- extraer_info_gff(gff_file, locus_pirate)
    if (!is.null(info)) {
      plot_data <- rbind(plot_data,
                         data.frame(cepa = cepa, gen = gen,
                                    contig = info$contig,
                                    start = as.numeric(info$start),
                                    end = as.numeric(info$end),
                                    stringsAsFactors = FALSE))
    }
  }
}

# Guardar datos
write.table(plot_data, "Figure_S6_colocalization_data.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

# Gráfico esquemático
if (nrow(plot_data) > 0) {
  library(ggplot2)
  # Simplificar: mostrar como puntos por cepa/contig
  p <- ggplot(plot_data, aes(x = start, y = cepa, colour = contig, shape = gen)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(x = "Posición en el contig", y = "Cepa",
         colour = "Contig", shape = "Gen",
         title = "Co-localización de pirA, pirB y genes asociados")
  ggsave("Figure_S6_colocalization.png", p, width = 10, height = 6, dpi = 300)
  ggsave("Figure_S6_colocalization.svg", p, width = 10, height = 6)
  cat("Figure S6 guardada: Figure_S6_colocalization.png / .svg\n")
} else {
  cat("No se generaron datos para Figure S6.\n")
}
