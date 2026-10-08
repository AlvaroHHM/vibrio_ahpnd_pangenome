#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Localizar contigs reales de genes de interés
# Versión corregida: busca por el número final del locus_tag (ej. 04690)

setwd("~/bacterial-genomics-tutorial/reference_genomes_pipeline/pirate_input_pirAB")

pirate_tsv <- "PIRATE_results_v3/PIRATE.gene_families.ordered.tsv"
cepas <- c("15_CESAIBC", "T13P")

gff_paths <- list(
  "15_CESAIBC" = "nuevos_genomas/15_CESAIBC/15_CESAIBC.gff",
  "T13P"       = "nuevos_genomas/T13P_prokka/T13P.gff"
)

genes_familias <- c(
  pirA   = "g07679",
  pirB   = "g05919",
  g07720 = "g07720",
  g06662 = "g06662",
  g07221 = "g07221"
)

# Función para extraer el contig desde GFF buscando por el número final
extraer_contig_por_numero <- function(gff_file, locus_pirate) {
  if (!file.exists(gff_file)) return(NA)

  # Extraer el número final del locus (ej. "15_CESAIBC_04690" -> "04690")
  partes <- strsplit(locus_pirate, "_")[[1]]
  numero <- tail(partes, n = 1)
  # Buscar en atributos "locus_tag=...<numero>"
  gff <- read.delim(gff_file, sep = "\t", header = FALSE,
                    comment.char = "#", stringsAsFactors = FALSE, quote = "")
  colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

  patron <- paste0("locus_tag=", ".*", numero)
  idx <- grep(patron, gff$attributes)
  if (length(idx) == 0) return(NA)
  return(gff$seqid[idx[1]])
}

# Obtener locus_tags desde PIRATE
pirate_header <- read.delim(pirate_tsv, sep = "\t", nrows = 1, header = FALSE,
                            stringsAsFactors = FALSE, check.names = FALSE)
col_cepas <- setNames(match(cepas, pirate_header), cepas)

resultados <- data.frame()

for (cepa in cepas) {
  cat("\n=== Cepa:", cepa, "===\n")
  col_idx <- col_cepas[cepa]
  gff_file <- gff_paths[[cepa]]

  for (gen in names(genes_familias)) {
    fam <- genes_familias[gen]
    cmd <- sprintf("awk -F'\\t' '$2==\"%s\"' %s", fam, pirate_tsv)
    linea <- system(cmd, intern = TRUE)
    if (length(linea) == 0) {
      cat(sprintf("  %-12s (%-8s): familia no encontrada\n", gen, fam))
      next
    }
    partes <- strsplit(linea[1], "\t")[[1]]
    locus_pirate <- partes[col_idx]
    if (is.na(locus_pirate) || locus_pirate == "") {
      cat(sprintf("  %-12s (%-8s): locus ausente en PIRATE\n", gen, fam))
      next
    }

    contig <- extraer_contig_por_numero(gff_file, locus_pirate)
    cat(sprintf("  %-12s (%-8s): locus %-20s -> contig %s\n", gen, fam, locus_pirate, contig))

    resultados <- rbind(resultados,
                        data.frame(cepa = cepa, gen = gen, familia = fam,
                                   locus_pirate = locus_pirate, contig = contig,
                                   stringsAsFactors = FALSE))
  }
}

write.table(resultados, "contigs_pirAB_genes_reales_v3.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nTabla guardada: contigs_pirAB_genes_reales_v3.tsv\n")

# Co-localización
cat("\n--- Co-localización ---\n")
for (cepa in cepas) {
  sub <- resultados[resultados$cepa == cepa, ]
  contigs <- unique(sub$contig[!is.na(sub$contig)])
  cat(sprintf("%s: %d contig(s) únicos -> %s\n", cepa, length(contigs),
              paste(contigs, collapse = ", ")))
}
