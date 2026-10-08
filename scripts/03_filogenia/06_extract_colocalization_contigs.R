#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ============================================================
# Extract real genomic coordinates for the 5 pVA1-borne genes.
# Tries multiple GFF directories: returns the first GFF that
# actually CONTAINS the requested locus_tag.
# ============================================================

IN_TSV  <- "input_data/strain_metadata/colocalizacion_AHPND_pos_v2.tsv"
OUT_TSV <- "input_data/strain_metadata/colocalization_with_coordinates.tsv"

GFF_DIRS <- c(
  path.expand("~/AHPND_Vp_study/01_genomes/annotated_prokka"),
  path.expand("~/AHPND_Vp_study/02_pangenome/PIRATE_results_v3/modified_gffs"),
  path.expand("~/AHPND_Vp_study/01_genomes/nuevos_genomas")
)

list_gff_candidates <- function(strain) {
  candidates <- character(0)
  for (gdir in GFF_DIRS) {
    candidates <- c(candidates,
      file.path(gdir, strain, paste0(strain, ".gff")),
      file.path(gdir, paste0(strain, ".gff"))
    )
  }
  candidates[file.exists(candidates)]
}

read_gff_cds <- function(gff_file) {
  gff <- read.delim(gff_file, sep = "\t", header = FALSE,
                    comment.char = "#", stringsAsFactors = FALSE, quote = "")
  colnames(gff) <- c("seqid","source","type","start","end",
                     "score","strand","phase","attributes")
  gff[gff$type == "CDS", ]
}

find_coords <- function(strain, locus_pirate) {
  num <- sub(".*_", "", locus_pirate)
  pat <- paste0("locus_tag=[^;]*", num, "($|;)")
  gffs <- list_gff_candidates(strain)
  for (gff_file in gffs) {
    cds <- read_gff_cds(gff_file)
    idx <- grep(pat, cds$attributes)
    if (length(idx) > 0) {
      g <- cds[idx[1], ]
      return(data.frame(contig_gff = g$seqid,
                        start      = as.numeric(g$start),
                        end        = as.numeric(g$end),
                        strand     = g$strand,
                        gff_source = basename(dirname(gff_file)),
                        stringsAsFactors = FALSE))
    }
  }
  return(NULL)
}

df <- read.delim(IN_TSV, sep = "\t", stringsAsFactors = FALSE)

rows <- list()
for (i in seq_len(nrow(df))) {
  s <- df$cepa[i]
  locus <- df$locus_pirate[i]
  coords <- find_coords(s, locus)
  if (is.null(coords)) {
    coords <- data.frame(contig_gff = NA, start = NA, end = NA,
                         strand = NA, gff_source = NA,
                         stringsAsFactors = FALSE)
  }
  rows[[i]] <- cbind(df[i, ], coords)
}

out <- do.call(rbind, rows)
write.table(out, OUT_TSV, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=========================================\n")
cat("Coordinates extracted\n")
cat("=========================================\n")
cat("Output:", OUT_TSV, "\n")
cat("Rows:  ", nrow(out), "\n")
cat("Missing coords:", sum(is.na(out$start)), "\n")

# Show which GFF was used per strain
cat("\nGFF source used per strain:\n")
for (s in unique(out$cepa)) {
  src <- unique(out$gff_source[out$cepa == s & !is.na(out$gff_source)])
  cat(" ", s, ":", ifelse(length(src) == 0, "(none)", paste(src, collapse=", ")), "\n")
}
