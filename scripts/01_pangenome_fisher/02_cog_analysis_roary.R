#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Análisis de categorías COG en genes exclusivos y ausentes de 15_CESAIBC (Roary)
# Entrada: archivos en ~/bacterial-genomics-tutorial/pangenome
# Salida: carpeta R_analysis_roary/

# ---------- Configuración ----------
setwd("~/bacterial-genomics-tutorial/pangenome")

# Crear carpeta de resultados
output_dir <- "R_analysis_roary"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ---------- 1. Leer listas de genes ----------
excl <- read.delim("exclusivos_15.tsv", sep = "\t", stringsAsFactors = FALSE)
aus  <- read.delim("ausentes_15_presentes_11_6.tsv", sep = "\t", stringsAsFactors = FALSE)

# Verificar columnas
if (!"gene" %in% colnames(excl) | !"gene" %in% colnames(aus)) {
  stop("Las tablas no tienen columna 'gene'")
}

genes_excl <- excl$gene
genes_aus <- aus$gene

# ---------- 2. Leer anotaciones eggNOG ----------
cat("Leyendo anotaciones eggNOG...\n")

leer_eggnog <- function(archivo) {
  dat <- read.delim(archivo, sep = "\t", quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
  # Extraer el ID del gen: columna #query o query, antes del '|'
  if ("#query" %in% colnames(dat)) {
    q <- dat[["#query"]]
  } else if ("query" %in% colnames(dat)) {
    q <- dat[["query"]]
  } else {
    stop("No se encontró columna #query")
  }
  gen_id <- sub("\\|.*", "", q)
  dat$gene <- gen_id
  # Seleccionar columna COG_category si existe
  if ("COG_category" %in% colnames(dat)) {
    dat <- dat[, c("gene", "COG_category")]
  } else {
    dat$COG_category <- NA
  }
  return(dat)
}

ann_excl <- leer_eggnog("eggNOG_exclusivos_15_annotations.tabular")
ann_aus  <- leer_eggnog("eggNOG_ausentes_15_annotations.tabular")

# Fusionar por gen (para cada conjunto)
excl_ann <- merge(data.frame(gene = genes_excl), ann_excl, by = "gene", all.x = TRUE)
aus_ann  <- merge(data.frame(gene = genes_aus),  ann_aus,  by = "gene", all.x = TRUE)

# ---------- 3. Preparar datos para COG ----------
# COG_category puede contener varias letras (ej. "KL"), las separamos
library(stringr)  # Si no está instalado, install.packages("stringr")

split_cog <- function(cog_vec) {
  # Devuelve lista de vectores de letras
  lapply(cog_vec, function(x) {
    if (is.na(x) | x == "-" | x == "") return(character(0))
    strsplit(x, "")[[1]]
  })
}

cog_excl <- unlist(split_cog(excl_ann$COG_category))
cog_aus  <- unlist(split_cog(aus_ann$COG_category))

# También contar genes sin COG
n_excl_total <- nrow(excl_ann)
n_aus_total  <- nrow(aus_ann)
n_excl_sin_cog <- sum(is.na(excl_ann$COG_category) | excl_ann$COG_category == "-")
n_aus_sin_cog  <- sum(is.na(aus_ann$COG_category)  | aus_ann$COG_category == "-")

# ---------- 4. Test exacto de Fisher por categoría COG ----------
letras <- LETTERS  # A-Z

resultados <- data.frame(
  COG = letras,
  Presente_Exclusivos = integer(length(letras)),
  Total_Exclusivos = n_excl_total,
  Presente_Ausentes = integer(length(letras)),
  Total_Ausentes = n_aus_total,
  odds_ratio = numeric(length(letras)),
  p_value = numeric(length(letras)),
  stringsAsFactors = FALSE
)

for (i in seq_along(letras)) {
  letra <- letras[i]
  n_excl_con <- sum(sapply(split_cog(excl_ann$COG_category), function(x) letra %in% x))
  n_aus_con  <- sum(sapply(split_cog(aus_ann$COG_category),  function(x) letra %in% x))
  
  # Tabla 2x2: [con letra, sin letra; con letra, sin letra]
  a <- n_excl_con
  b <- n_excl_total - n_excl_con
  c <- n_aus_con
  d <- n_aus_total - n_aus_con
  
  resultados$Presente_Exclusivos[i] <- a
  resultados$Presente_Ausentes[i]   <- c
  
  if (a+b > 0 & c+d > 0 & a+c > 0 & b+d > 0) {
    ft <- fisher.test(matrix(c(a,b,c,d), nrow=2))
    resultados$odds_ratio[i] <- ft$estimate
    resultados$p_value[i]    <- ft$p.value
  } else {
    resultados$odds_ratio[i] <- ifelse(a+b > 0 & c+d > 0, Inf, 0)
    resultados$p_value[i]    <- 1
  }
}

# Ajuste BH
resultados$FDR <- p.adjust(resultados$p_value, method = "BH")

# ---------- 5. Guardar resultados ----------
write.table(resultados, file.path(output_dir, "enriquecimiento_COG_roary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResumen de enriquecimiento COG (top significativos):\n")
print(head(resultados[order(resultados$p_value), ], 10))

cat("\nCategorías con FDR < 0.05:", sum(resultados$FDR < 0.05), "\n")

# ---------- 6. Gráfico de barras de proporciones ----------
# Calcular proporciones para graficar
resultados$Prop_Exclusivos <- resultados$Presente_Exclusivos / n_excl_total
resultados$Prop_Ausentes   <- resultados$Presente_Ausentes / n_aus_total

# Preparar datos para ggplot2 (si está instalado)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  df_plot <- data.frame(
    COG = rep(letras, 2),
    Conjunto = rep(c("Exclusive_15", "Absent_15"), each = length(letras)),
    Proporcion = c(resultados$Prop_Exclusivos, resultados$Prop_Ausentes)
  )
  
  p <- ggplot(df_plot, aes(x = COG, y = Proporcion, fill = Conjunto)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Categoría COG", y = "Proporción de genes con categoría",
         title = "Distribución de categorías COG en genes exclusivos y ausentes de 15_CESAIBC") +
    theme_minimal() +
    theme(legend.position = "top")
  
  ggsave(file.path(output_dir, "COG_barplot_roary.png"), plot = p,
         width = 12, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "COG_barplot_roary.svg"), plot = p,
         width = 12, height = 6)
  cat("Gráfico guardado en", output_dir, "\n")
} else {
  warning("ggplot2 no está instalado; no se generó el gráfico.")
}
