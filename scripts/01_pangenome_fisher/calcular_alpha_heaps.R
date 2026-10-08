#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Calcular alpha de Heaps para el pangenoma con randomización + bootstrap

# Rutas relativas al repo
# Ejecutar desde la raíz del repo: Rscript scripts/01_pangenome_fisher/calcular_alpha_heaps.R

set.seed(42)  # reproducibilidad

# Rutas relativas
input_pheno <- "input_data/strain_metadata/phenotype.csv"
input_pirate <- "input_data/pirate/PIRATE.gene_families.ordered.tsv"
output_file <- "outputs/pangenome/alpha_heaps.tsv"

# Verificar que existen
if (!file.exists(input_pheno)) stop("Falta: ", input_pheno)
if (!file.exists(input_pirate)) stop("Falta: ", input_pirate)

# Leer datos
pheno <- read.csv(input_pheno, stringsAsFactors = FALSE)
colnames(pheno) <- c("strain", "AHPND")

pirate <- read.delim(input_pirate, sep = "\t",
                     quote = "", comment.char = "", stringsAsFactors = FALSE,
                     check.names = FALSE)

cepas <- pheno$strain
cepas_cols <- intersect(cepas, colnames(pirate))

# Matriz binaria
presence <- pirate[, cepas_cols, drop = FALSE]
presence_bin <- ifelse(presence == "", 0, 1)
rownames(presence_bin) <- pirate$gene_family

# Función: tamaño del pangenoma dado un orden
pangenome_size <- function(orden) {
  genes_acumulados <- c()
  tamanos <- numeric(length(orden))
  for (i in seq_along(orden)) {
    genes_presentes <- rownames(presence_bin)[presence_bin[, orden[i]] == 1]
    genes_acumulados <- union(genes_acumulados, genes_presentes)
    tamanos[i] <- length(genes_acumulados)
  }
  return(tamanos)
}

# === RANDOMIZACIÓN ===
n_perm <- 1000
n_gen  <- length(cepas_cols)

# Matriz para guardar las 1000 curvas
curvas <- matrix(NA, nrow = n_gen, ncol = n_perm)

cat("Randomizando", n_perm, "veces...\n")
pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
for (i in 1:n_perm) {
  orden <- sample(cepas_cols)
  curvas[, i] <- pangenome_size(orden)
  setTxtProgressBar(pb, i)
}
close(pb)

# Curva promedio
y_mean <- rowMeans(curvas)
x <- seq_along(y_mean)

# === AJUSTE DE HEAPS ===
validos <- y_mean > 0 & x > 0
modelo <- lm(log(y_mean[validos]) ~ log(x[validos]))
alpha <- coef(modelo)[2]
k <- exp(coef(modelo)[1])

# === IC DE ALPHA DESDE LAS CURVAS RANDOMIZADAS ===
# Cada permutación ya representa un remuestreo del orden de genomas.
# Ajustamos Heaps a cada una de las 1000 curvas randomizadas para obtener
# la distribución empírica de alpha bajo distintos órdenes de adición.
cat("\nAjustando Heaps a las", n_perm, "curvas randomizadas...\n")
alpha_boot <- numeric(n_perm)
for (i in 1:n_perm) {
  y_i <- curvas[, i]
  x_i <- seq_along(y_i)
  validos_i <- y_i > 0 & x_i > 0
  fit_i <- lm(log(y_i[validos_i]) ~ log(x_i[validos_i]))
  alpha_boot[i] <- coef(fit_i)[2]
}
alpha_ci <- quantile(alpha_boot, c(0.025, 0.975), na.rm = TRUE)

# === ÚLTIMO INCREMENTO ===
last_inc <- (y_mean[n_gen] - y_mean[n_gen - 1]) / y_mean[n_gen - 1]
plateau <- last_inc < 0.05

# === GUARDAR ===
resultado <- data.frame(
  alpha = alpha,
  alpha_lower = alpha_ci[1],
  alpha_upper = alpha_ci[2],
  k = k,
  last_increment = last_inc,
  plateau_reached = plateau
)

write.table(resultado, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n========================================\n")
cat("Alpha de Heaps estimado:", round(alpha, 4), "\n")
cat("IC 95% (bootstrap percentil):", round(alpha_ci[1], 4), "-", round(alpha_ci[2], 4), "\n")
cat("Último incremento:", round(last_inc * 100, 2), "%\n")
cat("¿Meseta alcanzada (<5%)?:", plateau, "\n")
cat("Guardado en:", output_file, "\n")
