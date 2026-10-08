#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ============================================================
# Prepare gggenomes input: genes + seqs (WINDOWED around cluster)
# ============================================================

IN_TSV    <- "input_data/strain_metadata/colocalization_with_coordinates.tsv"
OUT_DIR   <- "input_data/strain_metadata"
OUT_GENES <- file.path(OUT_DIR, "gggenomes_genes.tsv")
OUT_SEQS  <- file.path(OUT_DIR, "gggenomes_seqs.tsv")

PAD <- 2000  # bp padding on each side of the gene cluster

df <- read.delim(IN_TSV, sep = "\t", stringsAsFactors = FALSE)

clean_contig <- function(cg) {
  cg <- sub("^gnl\\|Prokka\\|", "", cg)
  cg <- sub("^NZ_", "", cg)
  cg <- sub("\\.1$", "", cg)
  cg
}
df$contig_display <- clean_contig(df$contig_gff)
df$strand_num <- ifelse(df$strand == "+", 1,
                 ifelse(df$strand == "-", -1, NA))
df$seq_id <- paste(df$cepa, df$contig_display, sep = "__")

df <- df[!is.na(df$start), ]

# ---- Compute window per seq_id (loop avoids aggregate() issues) ----
seq_ids <- unique(df$seq_id)
windows <- data.frame(
  seq_id   = seq_ids,
  bin_id   = sapply(seq_ids, function(s) df$cepa[df$seq_id == s][1]),
  gene_min = sapply(seq_ids, function(s) min(df$start[df$seq_id == s])),
  gene_max = sapply(seq_ids, function(s) max(df$end[df$seq_id == s])),
  stringsAsFactors = FALSE
)

windows$offset    <- pmax(0, windows$gene_min - PAD)
windows$seq_start <- 0
windows$seq_end   <- windows$gene_max - windows$offset + PAD

# ---- Shift gene coordinates ----
df$window_offset <- windows$offset[match(df$seq_id, windows$seq_id)]
df$start_shift   <- df$start - df$window_offset
df$end_shift     <- df$end   - df$window_offset

# ---- Build GENES table ----
genes <- data.frame(
  seq_id  = df$seq_id,
  start   = df$start_shift,
  end     = df$end_shift,
  strand  = df$strand_num,
  gene_id = df$gen,
  bin_id  = df$cepa,
  stringsAsFactors = FALSE
)

# ---- Build SEQS table ----
seqs <- data.frame(
  seq_id = windows$seq_id,
  bin_id = windows$bin_id,
  start  = windows$seq_start,
  end    = windows$seq_end,
  stringsAsFactors = FALSE
)

# ---- Order ----
strain_order <- unique(df$cepa)
if ("15_CESAIBC" %in% strain_order) {
  strain_order <- c("15_CESAIBC", sort(setdiff(strain_order, "15_CESAIBC")))
}
seqs$bin_id  <- factor(seqs$bin_id,  levels = strain_order)
genes$bin_id <- factor(genes$bin_id, levels = strain_order)
seqs  <- seqs[order(seqs$bin_id, seqs$seq_id), ]
genes <- genes[order(genes$bin_id, genes$seq_id, genes$start), ]

seqs$bin_id  <- as.character(seqs$bin_id)
genes$bin_id <- as.character(genes$bin_id)

write.table(genes, OUT_GENES, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(seqs,  OUT_SEQS,  sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=========================================\n")
cat("gggenomes input tables (windowed)\n")
cat("=========================================\n")
cat("Genes :", OUT_GENES, " (", nrow(genes), "rows )\n")
cat("Seqs  :", OUT_SEQS,  " (", nrow(seqs),  "rows )\n\n")
cat("Window sizes per contig (bp):\n")
print(seqs)
