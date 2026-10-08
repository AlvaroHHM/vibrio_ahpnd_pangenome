#!/usr/bin/env Rscript
# =============================================================================
# 05_kaplan_meier.R
# Kaplan-Meier survival curves for the AHPND pathogenicity bioassay
# Run from repo root: Rscript scripts/04_bioensayo/05_kaplan_meier.R
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(ggsci)
})

# --- 1. Paths ----------------------------------------------------------------
input_file  <- file.path("input_data", "bioensayo", "bioensayo_datos.csv")
output_dir  <- file.path("outputs", "bioensayo")
figure_dir  <- file.path("outputs", "figures")

if (!file.exists(input_file)) {
  stop("❌ File not found: ", input_file,
       "\n   Run this script from the repo root")
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

# --- 2. Read data ------------------------------------------------------------
cat("📥 Reading data from:", input_file, "\n")

datos <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Treatment = recode(Sample,
                       "vp15VM" = "15_CESAIBC",
                       "vp11VM" = "11_VM",
                       "vp06VM" = "6_VM",
                       "C-"     = "Control"),
    Treatment = factor(Treatment,
                       levels = c("15_CESAIBC", "11_VM", "6_VM", "Control")),
    status = ifelse(evento == 1, 2, 1)
  )

cat("   Total individuals :", nrow(datos), "\n")
cat("   Treatments        :", paste(levels(datos$Treatment), collapse = ", "), "\n")

# --- 3. Survival object ------------------------------------------------------
surv_obj <- Surv(time = datos$HPI, event = datos$status)
km_fit   <- survfit(surv_obj ~ Treatment, data = datos, conf.type = "log-log")

# --- 4. Log-rank test --------------------------------------------------------
logrank      <- survdiff(surv_obj ~ Treatment, data = datos)
pval_logrank <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
df_logrank   <- length(logrank$n) - 1

p_annotation <- if (pval_logrank < 1e-16) {
  sprintf("Log-rank: chi2 = %.1f, df = %d, p < 2.2e-16",
          logrank$chisq, df_logrank)
} else if (pval_logrank < 0.001) {
  sprintf("Log-rank: chi2 = %.1f, df = %d, p < 0.001",
          logrank$chisq, df_logrank)
} else {
  sprintf("Log-rank: chi2 = %.1f, df = %d, p = %.3f",
          logrank$chisq, df_logrank, pval_logrank)
}

cat("📊", p_annotation, "\n")

# --- 5. Median survival (LT50) -----------------------------------------------
km_median <- surv_median(km_fit)
cat("\n📋 Median survival (LT50):\n")
print(km_median)

lt50_15 <- km_median$median[km_median$strata == "Treatment=15_CESAIBC"]
cat("\n   LT50 for 15_CESAIBC:", lt50_15, "HPI\n")

# --- 6. Colors ---------------------------------------------------------------
colores <- c(
  "15_CESAIBC" = "#EE0000",
  "11_VM"      = "#3B4992",
  "6_VM"       = "#008B45",
  "Control"    = "#7F7F7F"
)

# --- 7. Kaplan-Meier plot ----------------------------------------------------
km_plot <- ggsurvplot(
  km_fit,
  data              = datos,
  palette           = unname(colores),
  linewidth         = 1.0,
  censor.size       = 2.8,
  conf.int          = TRUE,
  conf.int.alpha    = 0.15,
  conf.int.style    = "ribbon",
  pval              = FALSE,
  risk.table        = "percentage",
  risk.table.height = 0.22,
  risk.table.col    = "strata",
  risk.table.y.text = FALSE,
  break.time.by     = 20,
  xlim              = c(0, 122),
  xlab              = "Hours post-inoculation (HPI)",
  ylab              = "Survival probability",
  legend.title      = "Treatment",
  legend.labs       = c("15_CESAIBC (AHPND+)", "11_VM", "6_VM", "Control"),
  ggtheme           = theme_cowplot(font_size = 12),
  font.main         = c(12, "bold"),
  tables.theme      = theme_cleantable() +
                      theme(
                        axis.text.y  = element_text(size = 8),
                        axis.text.x  = element_text(size = 8),
                        plot.title   = element_text(size = 9),
                        axis.title.x = element_text(size = 9)
                      )
)

# --- 8. LT50 lines + p-value annotation --------------------------------------
km_plot$plot <- km_plot$plot +
  annotate("segment", x = lt50_15, xend = lt50_15,
           y = 0, yend = 0.5,
           linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("segment", x = 0, xend = lt50_15,
           y = 0.5, yend = 0.5,
           linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  annotate("text", x = lt50_15 + 2.5, y = 0.53,
           label = sprintf("LT50 = %.1f h", lt50_15),
           hjust = 0, size = 3.6, colour = "grey20") +
  annotate("text", x = 65, y = 0.30,
           label = p_annotation,
           hjust = 0, size = 3.8, colour = "black") +
  theme(
    legend.position      = c(100 / 122, 0.50),
    legend.justification = c("center", "center"),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.key           = element_rect(fill = "white", colour = NA),
    legend.key.size      = unit(0.9, "lines"),
    legend.text          = element_text(size = 10)
  )

# --- 9. Export PDF + PNG -----------------------------------------------------
cat("\n💾 Saving figure to:", figure_dir, "\n")

cairo_pdf(file.path(figure_dir, "Figure_3_kaplan_meier.pdf"),
          width = 7, height = 6.5)
print(km_plot)
dev.off()

png(file.path(figure_dir, "Figure_3_kaplan_meier.png"),
    width = 7 * 200, height = 6.5 * 200, res = 200)
print(km_plot)
dev.off()

# --- 10. Summary tables ------------------------------------------------------
medianas <- km_median %>%
  mutate(
    Treatment = gsub("Treatment=", "", strata),
    median    = ifelse(is.na(median), "NR", as.character(median)),
    lower     = ifelse(is.na(lower),  "NR", as.character(lower)),
    upper     = ifelse(is.na(upper),  "NR", as.character(upper))
  ) %>%
  select(Treatment, median, lower, upper)

write_csv(medianas, file.path(output_dir, "kaplan_meier_summary.csv"))

surv_122 <- summary(km_fit, times = 122)
surv_tabla <- data.frame(
  Treatment     = gsub("Treatment=", "", surv_122$strata),
  HPI           = surv_122$time,
  n_risk        = surv_122$n.risk,
  n_event       = surv_122$n.event,
  Survival      = round(surv_122$surv,  4),
  Lower_95      = round(surv_122$lower, 4),
  Upper_95      = round(surv_122$upper, 4),
  Survival_pct  = round(surv_122$surv  * 100, 2),
  Lower_95_pct  = round(surv_122$lower * 100, 2),
  Upper_95_pct  = round(surv_122$upper * 100, 2)
)
write_csv(surv_tabla, file.path(output_dir, "survival_at_122HPI.csv"))

surv_full <- as.data.frame(summary(km_fit)$table)
write_csv(surv_full, file.path(output_dir, "survival_summary_full.csv"))

cat("\n📋 Median survival (NR = not reached):\n")
print(medianas)
cat("\n📋 Survival at 122 HPI (with 95% CI):\n")
print(surv_tabla)

cat("\n✅ Done. Output in:", output_dir, "\n")
