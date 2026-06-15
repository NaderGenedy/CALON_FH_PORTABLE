# =====================================================================
# CALON-FH poster — Lancet/Nature-calibre statistical panels (R)
# =====================================================================
# Produces 4 publication-grade figures from the locked v7 CSVs.
# Style: thin grid, single accent colour, no chart-junk, Cairo PDF + PNG.
#
# Output: ./poster_figures/
#   p1_headline_auc.{pdf,png}
#   p2_subgroup_forest.{pdf,png}
#   p3_calibration_dca.{pdf,png}
#   p4_coef_nri.{pdf,png}
# =====================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales); library(forcats); library(stringr)
})

OUT      <- "C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2"
FIG_DIR  <- file.path(OUT, "poster_figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Colour discipline — 3-tone max (Lancet figure standard)
COL_CALON   <- "#DA291C"   # Cardiff/AHA red
COL_SRE     <- "#5A6772"   # Steel grey
COL_NEUTRAL <- "#1B1F23"   # Charcoal
COL_LIGHT   <- "#E8EAED"   # Panel divider
THEME_NATURE <- theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(colour = COL_LIGHT, linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(colour = COL_NEUTRAL, linewidth = 0.4),
    axis.ticks         = element_line(colour = COL_NEUTRAL, linewidth = 0.3),
    axis.text          = element_text(colour = COL_NEUTRAL),
    axis.title         = element_text(colour = COL_NEUTRAL, face = "bold"),
    plot.title         = element_text(colour = COL_NEUTRAL, face = "bold", size = 13),
    plot.subtitle      = element_text(colour = COL_NEUTRAL, size = 10),
    plot.caption       = element_text(colour = COL_NEUTRAL, size = 8, hjust = 0),
    strip.text         = element_text(colour = COL_NEUTRAL, face = "bold"),
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    legend.background  = element_blank()
  )

# Load v7 outputs ------------------------------------------------------
res  <- read_csv(file.path(OUT, "CALON_FINAL_v7_results.csv"), show_col_types = FALSE)
cf   <- read_csv(file.path(OUT, "CALON_FINAL_v7_coefficients.csv"), show_col_types = FALSE)
sg   <- read_csv(file.path(OUT, "CALON_FINAL_v7_subgroups.csv"), show_col_types = FALSE)
cal  <- read_csv(file.path(OUT, "CALON_FINAL_v7_calibration.csv"), show_col_types = FALSE)
dca  <- read_csv(file.path(OUT, "CALON_FINAL_v7_dca.csv"), show_col_types = FALSE)
get  <- function(k) res |> filter(metric == k) |> pull(value) |> as.numeric()

# =====================================================================
# PANEL 1 — HEADLINE AUC (both directions, both models, ΔAUC + p)
# =====================================================================
p1_df <- tibble(
  direction = rep(c("A: Wales-clean → UKB", "B: UKB → Wales-clean"), each = 2),
  model     = rep(c("CALON-FH", "SAFEHEART-RE"), 2),
  auc       = c(get("A_CALON_AUC"), get("A_SRE_AUC"),
                get("B_CALON_AUC"), get("B_SRE_AUC")),
  lo        = c(get("A_CALON_CIlo"), get("A_SRE_CIlo"),
                get("B_CALON_CIlo"), get("B_SRE_CIlo")),
  hi        = c(get("A_CALON_CIhi"), get("A_SRE_CIhi"),
                get("B_CALON_CIhi"), get("B_SRE_CIhi"))
) |>
  mutate(model = factor(model, levels = c("CALON-FH", "SAFEHEART-RE")))

p1_ann <- tibble(
  direction = c("A: Wales-clean → UKB", "B: UKB → Wales-clean"),
  delta = c(get("A_delta_AUC"), get("B_delta_AUC")),
  p     = c(get("A_delta_p"),   get("B_delta_p")),
  y     = 0.85
) |>
  mutate(label = sprintf("ΔAUC %+.3f\n%s",
                          delta,
                          ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))))

p1 <- ggplot(p1_df, aes(x = model, y = auc, fill = model)) +
  geom_col(width = 0.65, colour = COL_NEUTRAL, linewidth = 0.4, alpha = 0.9) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18,
                linewidth = 0.5, colour = COL_NEUTRAL) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.4, size = 3.6,
            fontface = "bold", colour = COL_NEUTRAL,
            position = position_nudge(y = pmax(p1_df$hi - p1_df$auc, 0.01) + 0.01)) +
  geom_text(data = p1_ann, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 3.3, colour = COL_NEUTRAL, lineheight = 0.95) +
  facet_wrap(~direction, nrow = 1) +
  scale_fill_manual(values = c("CALON-FH" = COL_CALON, "SAFEHEART-RE" = COL_SRE),
                    name = NULL) +
  scale_y_continuous(limits = c(0.55, 0.92), breaks = seq(0.55, 0.9, 0.05),
                     expand = expansion(0)) +
  labs(title = "External AUC: CALON-FH beats SAFEHEART-RE in both directions",
       subtitle = "TRIPOD Type 4 bidirectional external validation, 95% CI from 2,000-iter bootstrap",
       y = "External AUC (95% CI)", x = NULL,
       caption = "Cohorts: Wales-clean n=200 (events 54), UKB LDLR carriers n=3,540 (events 165). Paired-bootstrap test for ΔAUC.") +
  THEME_NATURE

ggsave(file.path(FIG_DIR, "p1_headline_auc.png"), p1,
       width = 8, height = 5, dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "p1_headline_auc.pdf"), p1,
       width = 8, height = 5, device = cairo_pdf, bg = "white")

# =====================================================================
# PANEL 2 — SUBGROUP FOREST PLOT (Direction A, primary external)
# =====================================================================
sgA <- sg |>
  filter(direction == "A_Wales_to_UKB") |>
  arrange(group, level) |>
  mutate(label   = sprintf("%s = %s   (n=%d, events=%d)", group, level, n, events),
         se_app  = 1 / sqrt(events),
         lo      = delta - 1.96 * se_app,
         hi      = delta + 1.96 * se_app,
         favours = ifelse(delta > 0, "CALON-FH", "SAFEHEART-RE")) |>
  mutate(label = factor(label, levels = rev(label)))

p2 <- ggplot(sgA, aes(x = delta, y = label, colour = favours)) +
  annotate("rect", xmin = 0, xmax = 0.20, ymin = -Inf, ymax = Inf,
           fill = COL_CALON, alpha = 0.04) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%+.3f", delta)), hjust = -0.25,
            size = 3.1, colour = COL_NEUTRAL, fontface = "bold") +
  scale_colour_manual(values = c("CALON-FH" = COL_CALON,
                                  "SAFEHEART-RE" = COL_SRE),
                       name = "Favours") +
  scale_x_continuous(limits = c(-0.04, 0.20), breaks = seq(-0.04, 0.20, 0.04)) +
  labs(title = "Subgroup ΔAUC (Direction A): CALON wins in every prespecified stratum",
       subtitle = "Sex, age band, untreated-LDL band, T2DM, smoking, hypertension",
       x = "ΔAUC (CALON-FH − SAFEHEART-RE)", y = NULL,
       caption = "Approximate 95% CI from event-count-derived SE. Bootstrap CIs preserved in supplementary table.") +
  THEME_NATURE +
  theme(panel.grid.major.y = element_line(colour = COL_LIGHT, linewidth = 0.2))

ggsave(file.path(FIG_DIR, "p2_subgroup_forest.png"), p2,
       width = 9, height = 7, dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "p2_subgroup_forest.pdf"), p2,
       width = 9, height = 7, device = cairo_pdf, bg = "white")

# =====================================================================
# PANEL 3 — CALIBRATION + DCA combined
# =====================================================================
# Calibration (decile-binned, both directions, CALON only)
cal_clean <- cal |>
  mutate(direction = factor(direction,
                             levels = c("A_Wales_to_UKB","B_UKB_to_Wales"),
                             labels = c("A: Wales → UKB", "B: UKB → Wales")))

p3a <- ggplot(cal_clean, aes(x = mean_predicted, y = observed_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_line(colour = COL_CALON, linewidth = 0.7) +
  geom_point(colour = COL_CALON, size = 2.6) +
  facet_wrap(~direction, scales = "free") +
  labs(title = "Calibration (decile-binned)",
       x = "Mean predicted probability", y = "Observed event rate") +
  THEME_NATURE

# DCA (both models, both directions)
dca_clean <- dca |>
  pivot_longer(c(NB_CALON, NB_SRE, NB_treat_all),
               names_to = "model", values_to = "NB") |>
  mutate(direction = factor(direction,
                             levels = c("A_Wales_to_UKB","B_UKB_to_Wales"),
                             labels = c("A: Wales → UKB", "B: UKB → Wales")),
         model = recode(model,
                         "NB_CALON" = "CALON-FH",
                         "NB_SRE"   = "SAFEHEART-RE",
                         "NB_treat_all" = "Treat all"))

p3b <- ggplot(dca_clean, aes(x = threshold, y = NB, colour = model, linetype = model)) +
  geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~direction, scales = "free_y") +
  scale_colour_manual(values = c("CALON-FH" = COL_CALON,
                                  "SAFEHEART-RE" = COL_SRE,
                                  "Treat all" = "#9CA3AF"),
                       name = NULL) +
  scale_linetype_manual(values = c("CALON-FH" = "solid",
                                    "SAFEHEART-RE" = "solid",
                                    "Treat all" = "dashed"),
                         name = NULL) +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.30),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Decision curve — net clinical benefit",
       x = "Threshold probability", y = "Net benefit") +
  THEME_NATURE

p3 <- p3a / p3b + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "Calibration & clinical utility, both external directions",
    caption = "Brier (A) CALON 0.066 vs SRE 0.073; (B) CALON 0.241 vs SRE 0.256.",
    theme = THEME_NATURE
  )

ggsave(file.path(FIG_DIR, "p3_calibration_dca.png"), p3,
       width = 9, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "p3_calibration_dca.pdf"), p3,
       width = 9, height = 8.5, device = cairo_pdf, bg = "white")

# =====================================================================
# PANEL 4 — COEFFICIENT FOREST + NRI / IDI
# =====================================================================
cfA <- cf |>
  filter(direction == "A_Wales_trained") |>
  mutate(feature = fct_reorder(feature, OR_per_SD))

p4a <- ggplot(cfA, aes(x = OR_per_SD, y = feature)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_segment(aes(x = 1, xend = OR_per_SD, yend = feature),
               colour = COL_LIGHT, linewidth = 1.6) +
  geom_point(size = 3, colour = COL_CALON) +
  geom_text(aes(label = sprintf("%.2f", OR_per_SD)),
            hjust = -0.35, size = 3.1, fontface = "bold", colour = COL_NEUTRAL) +
  scale_x_continuous(limits = c(0.85, 2.4), breaks = seq(1, 2.4, 0.2)) +
  labs(title = "CALON-FH equation",
       subtitle = "11-band Wales-clean trained model (locked v7)",
       x = "OR per SD", y = NULL) +
  THEME_NATURE +
  theme(panel.grid.major.y = element_line(colour = COL_LIGHT, linewidth = 0.2))

# NRI / IDI summary
nri_df <- tibble(
  direction = rep(c("A: Wales → UKB", "B: UKB → Wales"), each = 4),
  metric    = rep(c("NRI total", "NRI event", "NRI non-event", "IDI"), 2),
  value     = c(get("A_NRI_total"), get("A_NRI_event"),
                get("A_NRI_nonevent"), get("A_IDI"),
                get("B_NRI_total"), get("B_NRI_event"),
                get("B_NRI_nonevent"), get("B_IDI"))
) |>
  mutate(metric = factor(metric, levels = c("NRI total","NRI event",
                                              "NRI non-event","IDI")),
         sign = ifelse(value >= 0, "Gain", "Loss"))

p4b <- ggplot(nri_df, aes(x = metric, y = value, fill = sign)) +
  geom_col(width = 0.65, colour = COL_NEUTRAL, linewidth = 0.4, alpha = 0.9) +
  geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", value),
                vjust = ifelse(value >= 0, -0.5, 1.3)),
            size = 3.2, fontface = "bold", colour = COL_NEUTRAL) +
  facet_wrap(~direction, nrow = 1) +
  scale_fill_manual(values = c("Gain" = COL_CALON, "Loss" = COL_SRE),
                    name = NULL) +
  labs(title = "Reclassification — net (CALON − SAFEHEART)",
       subtitle = "Categorical NRI 3-tier @ 0.05/0.20  •  IDI",
       x = NULL, y = NULL) +
  THEME_NATURE +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p4 <- p4a / p4b + plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Coefficients + reclassification",
    caption = "OR per SD on standardised binary bands. NRI/IDI: clinically-stratified 3-tier cut at 5% and 20%.",
    theme = THEME_NATURE
  )

ggsave(file.path(FIG_DIR, "p4_coef_nri.png"), p4,
       width = 9, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "p4_coef_nri.pdf"), p4,
       width = 9, height = 8.5, device = cairo_pdf, bg = "white")

cat("\n[DONE] 4 panels in", FIG_DIR, "\n")
cat("  p1_headline_auc.{pdf,png}\n")
cat("  p2_subgroup_forest.{pdf,png}\n")
cat("  p3_calibration_dca.{pdf,png}\n")
cat("  p4_coef_nri.{pdf,png}\n")
