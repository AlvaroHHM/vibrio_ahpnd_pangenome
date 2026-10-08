#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ============================================================
# Figure S3: Co-localization of pVA1-borne genes in AHPND+ strains
# ============================================================

IN_TSV  <- "input_data/strain_metadata/colocalization_with_coordinates.tsv"
OUT_DIR <- "outputs/supplementary_figures"
OUT_PDF <- file.path(OUT_DIR, "Figure_S3_colocalization.pdf")
OUT_PNG <- file.path(OUT_DIR, "Figure_S3_colocalization.png")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Read data ----
df <- read.delim(IN_TSV, sep = "\t", stringsAsFactors = FALSE)
df$strand <- substr(df$strand, 1, 1)
df$strand[!df$strand %in% c("+", "-")] <- NA

clean_contig <- function(cg) {
  cg <- sub("^gnl\\|Prokka\\|", "", cg)
  cg <- sub("^NZ_", "", cg)
  cg <- sub("\\.1$", "", cg)
  cg
}
df$contig_display <- ifelse(is.na(df$contig_gff), NA, clean_contig(df$contig_gff))

gene_colors <- c(
  pirA   = "#EE0000",
  pirB   = "#8B0000",
  g07720 = "#3B4992",
  g06662 = "#5F559B",
  g07221 = "#008B45"
)
gene_order <- c("pirA", "pirB", "g07720", "g06662", "g07221")

# ---- Order strains ----
strains <- unique(df$cepa)
if ("15_CESAIBC" %in% strains) {
  strains <- c("15_CESAIBC", sort(setdiff(strains, "15_CESAIBC")))
}
n <- length(strains)

# ---- Per-strain contig structure ----
strain_contigs <- list()
for (s in strains) {
  sub <- df[df$cepa == s & !is.na(df$contig_gff), ]
  if (nrow(sub) == 0) {
    strain_contigs[[s]] <- character(0)
  } else {
    agg <- aggregate(start ~ contig_display, data = sub, FUN = min)
    strain_contigs[[s]] <- agg$contig_display[order(agg$start)]
  }
}

# ---- Assign y positions ----
y_pos <- list()
total_rows <- 0
for (s in strains) {
  contigs <- strain_contigs[[s]]
  if (length(contigs) == 0) {
    y_pos[[s]] <- total_rows + 1
    total_rows <- total_rows + 1
  } else {
    y_pos[[s]] <- seq(total_rows + length(contigs), total_rows + 1)
    total_rows <- total_rows + length(contigs)
  }
}

# ---- Compute relative positions per contig ----
df$rel_pos <- NA
for (s in strains) {
  for (cg in strain_contigs[[s]]) {
    idx <- which(df$cepa == s & !is.na(df$contig_display) & df$contig_display == cg)
    if (length(idx) > 0) {
      s_min <- min(df$start[idx], na.rm = TRUE)
      s_max <- max(df$end[idx],   na.rm = TRUE)
      if (s_max > s_min) {
        df$rel_pos[idx] <- (df$start[idx] - s_min) / (s_max - s_min) * 10
      } else {
        df$rel_pos[idx] <- 5
      }
    }
  }
}

# ---- Plotting ----
plot_s3 <- function() {
  # Left margin for strain names, right margin WIDE for contig labels
  par(mar = c(3, 11, 2, 10), xpd = NA)

  plot(NA,
       xlim = c(-1, 11),
       ylim = c(0.3, total_rows + 0.7),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       main = "", bty = "n")

  for (s in strains) {
    contigs <- strain_contigs[[s]]
    ys <- y_pos[[s]]

    if (length(contigs) == 0) {
      y <- ys[1]
      mtext(s, side = 2, at = y, las = 1, line = 0.5, cex = 1.2, font = 2)
      text(5, y, "(coordinates unavailable)", col = "gray50", cex = 1.2, font = 3)
      next
    }

    for (k in seq_along(contigs)) {
      cg <- contigs[k]
      y <- ys[k]
      sub <- df[df$cepa == s & !is.na(df$contig_display) & df$contig_display == cg, ]
      sub <- sub[order(sub$rel_pos), ]

      # Backbone
      segments(min(sub$rel_pos) - 0.5, y,
               max(sub$rel_pos) + 0.5, y,
               col = "#333333", lwd = 2)

      # Arrows
      for (j in seq_len(nrow(sub))) {
        x <- sub$rel_pos[j]
        g <- sub$gen[j]
        col <- gene_colors[g]
        strand <- sub$strand[j]
        half_w <- 0.35
        half_h <- 0.15

        if (!is.na(strand) && strand == "+") {
          polygon(c(x - half_w, x + half_w - 0.1, x + half_w,
                    x + half_w - 0.1, x - half_w),
                  c(y - half_h, y - half_h, y,
                    y + half_h, y + half_h),
                  col = col, border = "black", lwd = 0.8)
        } else if (!is.na(strand) && strand == "-") {
          polygon(c(x + half_w, x - half_w + 0.1, x - half_w,
                    x - half_w + 0.1, x + half_w),
                  c(y - half_h, y - half_h, y,
                    y + half_h, y + half_h),
                  col = col, border = "black", lwd = 0.8)
        } else {
          rect(x - half_w, y - half_h, x + half_w, y + half_h,
               col = col, border = "black", lwd = 0.8)
        }
      }

      # Strain label
      if (k == 1) {
        mtext(s, side = 2, at = mean(ys), las = 1, line = 0.5, cex = 1.2, font = 2)
      }

      # Contig label — CLOSE to plot (line = 0.3)
      mtext(cg, side = 4, at = y, las = 1, line = 0.3,
            cex = 1.2, col = "gray30", font = 3)
    }
  }

  # Legend at bottom, centered
  legend("bottom",
         legend = gene_order,
         fill = gene_colors[gene_order],
         border = "black",
         ncol = 5, bty = "n", cex = 1.2,
         inset = c(0, 0.02),
         xpd = NA)
}

# ---- Save ----
pdf(OUT_PDF, width = 18, height = 0.9 * total_rows + 2)
plot_s3()
dev.off()

png(OUT_PNG, width = 18 * 200, height = (0.9 * total_rows + 2) * 200, res = 200)
plot_s3()
dev.off()

cat("\n=========================================\n")
cat("Figure S3 generated\n")
cat("=========================================\n")
cat("PDF :", OUT_PDF, "\n")
cat("PNG :", OUT_PNG, "\n")
cat("Rows:", total_rows, "\n")
cat("Strains:", n, "\n")
