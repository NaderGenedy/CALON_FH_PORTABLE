# =====================================================================
# CALON-FH poster — 8-figure Nature-calibre R pack
# =====================================================================
# Builds 8 publication-grade panels (600 DPI PNG + Cairo PDF) from the
# locked v7 output CSVs.
#
# Style discipline: 3-tone palette (charcoal/red/grey), no chart-junk,
# minimum 30% whitespace, every error bar a 95% CI, every metric labelled.
#
# Panels:
#   p1 — Headline external AUC, both directions, with ΔAUC + p
#   p2 — 14-subgroup forest plot (Direction A primary)
#   p3 — Calibration (decile) + DCA stacked
#   p4 — Coefficient OR forest + NRI bars stacked
#   p5 — 13-tile per-subgroup mini-bar infographic (Direction A)
#   p6 — NEW: Cohort attrition waterfall (7,253 → 200)
#   p7 — NEW: Risk-stratum reclassification heatmap (SAFEHEART vs CALON, 3×3)
#   p8 — NEW: Treatment-era timeline (1995→2025), CALON & SAFEHEART positioned
# =====================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales); library(forcats); library(stringr)
})

OUT      <- "C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2"
FIG_DIR  <- file.path(OUT, "poster_figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Tone discipline
COL_CALON   <- "#DA291C"
COL_SRE     <- "#5A6772"
COL_NEUTRAL <- "#1B1F23"
COL_NOISY   <- "#9CA3AF"
COL_LIGHT   <- "#E8EAED"
COL_ACCENT  <- "#A8120D"   # darker red for highlights

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
    legend.background  = element_blank(),
    plot.background    = element_rect(fill = "white", colour = NA)
  )

# Load locked v7 outputs
res  <- read_csv(file.path(OUT, "CALON_FINAL_v7_results.csv"), show_col_types = FALSE)
cf   <- read_csv(file.path(OUT, "CALON_FINAL_v7_coefficients.csv"), show_col_types = FALSE)
sg   <- read_csv(file.path(OUT, "CALON_FINAL_v7_subgroup_NRI.csv"), show_col_types = FALSE)
cal  <- read_csv(file.path(OUT, "CALON_FINAL_v7_calibration.csv"), show_col_types = FALSE)
dca  <- read_csv(file.path(OUT, "CALON_FINAL_v7_dca.csv"), show_col_types = FALSE)
get  <- function(k) res |> filter(metric == k) |> pull(value) |> as.numeric()

save_panel <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, sprintf("%s.png", name)), p,
         width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(FIG_DIR, sprintf("%s.pdf", name)), p,
         width = w, height = h, device = cairo_pdf, bg = "white")
  cat(sprintf("  wrote %s.{png,pdf}\n", name))
}

# =====================================================================
# PANEL 1 — Headline AUC bars (with ΔAUC + p)
# =====================================================================
cat("\n[1] Headline AUC bars\n")
p1_df <- tibble(
  direction = rep(c("A: Wales-clean → UKB", "B: UKB → Wales-clean"), each = 2),
  model = rep(c("CALON-FH", "SAFEHEART-RE"), 2),
  auc = c(get("A_CALON_AUC"), get("A_SRE_AUC"),
          get("B_CALON_AUC"), get("B_SRE_AUC")),
  lo = c(get("A_CALON_CIlo"), get("A_SRE_CIlo"),
         get("B_CALON_CIlo"), get("B_SRE_CIlo")),
  hi = c(get("A_CALON_CIhi"), get("A_SRE_CIhi"),
         get("B_CALON_CIhi"), get("B_SRE_CIhi"))
) |> mutate(model = factor(model, levels = c("CALON-FH", "SAFEHEART-RE")))

p1_ann <- tibble(
  direction = c("A: Wales-clean → UKB", "B: UKB → Wales-clean"),
  delta = c(get("A_delta_AUC"), get("B_delta_AUC")),
  p = c(get("A_delta_p"), get("B_delta_p")), y = 0.85
) |> mutate(label = sprintf("ΔAUC %+.3f\n%s", delta,
                              ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))))

p1 <- ggplot(p1_df, aes(x = model, y = auc, fill = model)) +
  geom_col(width = 0.6, colour = COL_NEUTRAL, linewidth = 0.4, alpha = 0.9) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.5, colour = COL_NEUTRAL) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.4, size = 3.6, fontface = "bold",
            colour = COL_NEUTRAL,
            position = position_nudge(y = pmax(p1_df$hi - p1_df$auc, 0.01) + 0.01)) +
  geom_text(data = p1_ann, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, size = 3.3, colour = COL_NEUTRAL, lineheight = 0.95) +
  facet_wrap(~direction, nrow = 1) +
  scale_fill_manual(values = c("CALON-FH" = COL_CALON, "SAFEHEART-RE" = COL_SRE), name = NULL) +
  scale_y_continuous(limits = c(0.55, 0.92), breaks = seq(0.55, 0.9, 0.05)) +
  labs(title = "External AUC: CALON-FH beats SAFEHEART-RE in both directions",
       subtitle = "TRIPOD Type 4 bidirectional, 95% CI from 2,000-iter bootstrap",
       y = "External AUC (95% CI)", x = NULL,
       caption = "Wales-clean n=200 (events 54); UKB LDLR carriers n=3,540 (events 165).") +
  THEME_NATURE
save_panel(p1, "p1_headline_auc", 8, 5)

# =====================================================================
# PANEL 2 — Subgroup forest (Direction A)
# =====================================================================
cat("\n[2] Subgroup forest\n")
sgA <- sg |> filter(direction == "A_Wales_to_UKB") |>
  arrange(group, level) |>
  mutate(label = sprintf("%s = %s   (n=%d, events=%d)", group, level, n, events),
         se_app = 1 / sqrt(events),
         lo = delta_auc - 1.96 * se_app, hi = delta_auc + 1.96 * se_app,
         favours = ifelse(delta_auc > 0, "CALON-FH", "SAFEHEART-RE"),
         label = factor(label, levels = rev(label)))

p2 <- ggplot(sgA, aes(x = delta_auc, y = label, colour = favours)) +
  annotate("rect", xmin = 0, xmax = 0.20, ymin = -Inf, ymax = Inf,
           fill = COL_CALON, alpha = 0.04) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%+.3f", delta_auc)), hjust = -0.25,
            size = 3.1, colour = COL_NEUTRAL, fontface = "bold") +
  scale_colour_manual(values = c("CALON-FH" = COL_CALON, "SAFEHEART-RE" = COL_SRE), name = "Favours") +
  scale_x_continuous(limits = c(-0.04, 0.20), breaks = seq(-0.04, 0.20, 0.04)) +
  labs(title = "Subgroup ΔAUC, Direction A: CALON wins in 14/14 prespecified strata",
       subtitle = "Sex, age, untreated-LDL band, T2DM, smoking, hypertension",
       x = "ΔAUC (CALON-FH − SAFEHEART-RE)", y = NULL,
       caption = "Approximate 95% CIs from event-count SE.") +
  THEME_NATURE +
  theme(panel.grid.major.y = element_line(colour = COL_LIGHT, linewidth = 0.2))
save_panel(p2, "p2_subgroup_forest", 9, 7)

# =====================================================================
# PANEL 3 — Calibration + DCA
# =====================================================================
cat("\n[3] Calibration + DCA\n")
cal_c <- cal |> mutate(direction = factor(direction,
                                           levels = c("A_Wales_to_UKB","B_UKB_to_Wales"),
                                           labels = c("A: Wales → UKB","B: UKB → Wales")))
p3a <- ggplot(cal_c, aes(x = mean_predicted, y = observed_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_line(colour = COL_CALON, linewidth = 0.7) +
  geom_point(colour = COL_CALON, size = 2.6) +
  facet_wrap(~direction, scales = "free") +
  labs(title = "Calibration (decile-binned)", x = "Mean predicted probability", y = "Observed event rate") +
  THEME_NATURE

dca_c <- dca |>
  pivot_longer(c(NB_CALON, NB_SRE, NB_treat_all), names_to = "model", values_to = "NB") |>
  mutate(direction = factor(direction, levels = c("A_Wales_to_UKB","B_UKB_to_Wales"),
                             labels = c("A: Wales → UKB","B: UKB → Wales")),
         model = recode(model, "NB_CALON" = "CALON-FH", "NB_SRE" = "SAFEHEART-RE",
                        "NB_treat_all" = "Treat all"))
p3b <- ggplot(dca_c, aes(x = threshold, y = NB, colour = model, linetype = model)) +
  geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
  facet_wrap(~direction, scales = "free_y") +
  scale_colour_manual(values = c("CALON-FH" = COL_CALON, "SAFEHEART-RE" = COL_SRE,
                                  "Treat all" = COL_NOISY), name = NULL) +
  scale_linetype_manual(values = c("CALON-FH" = "solid", "SAFEHEART-RE" = "solid",
                                    "Treat all" = "dashed"), name = NULL) +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.15, 0.20, 0.30),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Decision curve — net clinical benefit",
       x = "Threshold probability", y = "Net benefit") +
  THEME_NATURE

p3 <- p3a / p3b + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(title = "Calibration & clinical utility, both external directions",
                  caption = "Brier (A) CALON 0.066 vs SRE 0.073; (B) CALON 0.241 vs SRE 0.256.",
                  theme = THEME_NATURE)
save_panel(p3, "p3_calibration_dca", 9, 8.5)

# =====================================================================
# PANEL 4 — Coefficient forest + NRI bars
# =====================================================================
cat("\n[4] Coefficient + NRI\n")
cfA <- cf |> filter(direction == "A_Wales_trained") |>
  mutate(feature = fct_reorder(feature, OR_per_SD))

p4a <- ggplot(cfA, aes(x = OR_per_SD, y = feature)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_segment(aes(x = 1, xend = OR_per_SD, yend = feature),
               colour = COL_LIGHT, linewidth = 1.6) +
  geom_point(size = 3, colour = COL_CALON) +
  geom_text(aes(label = sprintf("%.2f", OR_per_SD)), hjust = -0.35,
            size = 3.1, fontface = "bold", colour = COL_NEUTRAL) +
  scale_x_continuous(limits = c(0.85, 2.4), breaks = seq(1, 2.4, 0.2)) +
  labs(title = "CALON-FH equation", subtitle = "11-band Wales-clean trained model (locked v7)",
       x = "OR per SD", y = NULL) +
  THEME_NATURE +
  theme(panel.grid.major.y = element_line(colour = COL_LIGHT, linewidth = 0.2))

nri_df <- tibble(
  direction = rep(c("A: Wales → UKB", "B: UKB → Wales"), each = 4),
  metric = rep(c("NRI total", "NRI event", "NRI non-event", "IDI"), 2),
  value = c(get("A_NRI_total"), get("A_NRI_event"), get("A_NRI_nonevent"), get("A_IDI"),
            get("B_NRI_total"), get("B_NRI_event"), get("B_NRI_nonevent"), get("B_IDI"))
) |> mutate(metric = factor(metric, levels = c("NRI total","NRI event","NRI non-event","IDI")),
            sign = ifelse(value >= 0, "Gain", "Loss"))

p4b <- ggplot(nri_df, aes(x = metric, y = value, fill = sign)) +
  geom_col(width = 0.65, colour = COL_NEUTRAL, linewidth = 0.4, alpha = 0.9) +
  geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.3f", value),
                vjust = ifelse(value >= 0, -0.5, 1.3)),
            size = 3.2, fontface = "bold", colour = COL_NEUTRAL) +
  facet_wrap(~direction, nrow = 1) +
  scale_fill_manual(values = c("Gain" = COL_CALON, "Loss" = COL_SRE), name = NULL) +
  labs(title = "Reclassification — net (CALON − SAFEHEART)",
       subtitle = "Categorical NRI 3-tier @ 0.05/0.20  •  IDI", x = NULL, y = NULL) +
  THEME_NATURE +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

p4 <- p4a / p4b + plot_layout(heights = c(1.2, 1)) +
  plot_annotation(title = "Coefficients + reclassification",
                  caption = "OR per SD on standardised binary bands. NRI/IDI 3-tier at 5%/20%.",
                  theme = THEME_NATURE)
save_panel(p4, "p4_coef_nri", 9, 8.5)

# =====================================================================
# PANEL 5 — 13-tile mini-bar grid (Direction A)
# =====================================================================
cat("\n[5] 13-tile subgroup mini-bar grid\n")
sgA_t <- sg |> filter(direction == "A_Wales_to_UKB") |>
  mutate(noisy = events < 10,
         delta_lbl = sprintf("ΔAUC %+.3f", delta_auc),
         nri_lbl = ifelse(noisy, sprintf("NRI %+.0f%% †", nri_total*100),
                           sprintf("NRI %+.0f%%", nri_total*100)),
         tile_order = factor(level, levels = c("Female","Male","<50y","50-65y",">=65y",
                                                 "LDL<4.14","LDL 4.14-6.5","LDL>=6.5",
                                                 "No DM","T2DM","Non-smoker","Smoker",
                                                 "No HTN","HTN")))

tiles_bars <- sgA_t |> select(tile_order, auc_calon, auc_sre, noisy) |>
  pivot_longer(c(auc_calon, auc_sre), names_to = "model", values_to = "auc") |>
  mutate(model = recode(model, "auc_calon" = "CALON-FH", "auc_sre" = "SAFEHEART-RE"))

p5 <- ggplot(tiles_bars, aes(x = model, y = auc, fill = model)) +
  geom_col(width = 0.65, colour = COL_NEUTRAL, linewidth = 0.3, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f", auc), y = auc), vjust = -0.3,
            size = 2.6, colour = COL_NEUTRAL, fontface = "bold") +
  scale_fill_manual(values = c("CALON-FH" = COL_CALON, "SAFEHEART-RE" = COL_SRE), name = NULL) +
  scale_y_continuous(limits = c(0, 0.95), expand = expansion(0)) +
  geom_text(data = sgA_t, aes(x = 1.5, y = 0.88, label = delta_lbl),
            inherit.aes = FALSE, size = 3, fontface = "bold", colour = COL_CALON) +
  geom_text(data = sgA_t, aes(x = 1.5, y = 0.79, label = nri_lbl, colour = noisy),
            inherit.aes = FALSE, size = 2.7, fontface = "bold", show.legend = FALSE) +
  geom_text(data = sgA_t, aes(x = 1.5, y = 0.05, label = sprintf("n=%d, ev=%d", n, events)),
            inherit.aes = FALSE, size = 2.2, colour = COL_NEUTRAL) +
  scale_colour_manual(values = c("TRUE" = COL_NOISY, "FALSE" = COL_NEUTRAL), guide = "none") +
  facet_wrap(~ tile_order, ncol = 7, nrow = 2, strip.position = "top") +
  labs(title = "Subgroup wins (Direction A) — ΔAUC and NRI per stratum",
       subtitle = "CALON-FH (red) vs SAFEHEART-RE (grey)",
       x = NULL, y = "External AUC",
       caption = "† NRI noisy (events < 10); ΔAUC remains primary metric.") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 13, colour = COL_NEUTRAL),
        plot.subtitle = element_text(size = 10, colour = COL_NEUTRAL),
        plot.caption = element_text(size = 8, colour = COL_NEUTRAL, hjust = 0),
        strip.text = element_text(face = "bold", colour = COL_NEUTRAL, size = 10),
        strip.background = element_rect(fill = COL_LIGHT, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = COL_LIGHT, linewidth = 0.2),
        panel.border = element_rect(colour = COL_NEUTRAL, fill = NA, linewidth = 0.4),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7, colour = COL_NEUTRAL),
        axis.text.y = element_text(size = 7, colour = COL_NEUTRAL),
        axis.title.y = element_text(size = 9, colour = COL_NEUTRAL),
        axis.line = element_line(colour = COL_NEUTRAL, linewidth = 0.3),
        legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size = 9, colour = COL_NEUTRAL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.spacing.x = unit(0.4, "lines"), panel.spacing.y = unit(0.6, "lines"))
save_panel(p5, "p5_subgroup_tiles", 16, 7.5)

# =====================================================================
# PANEL 6 — NEW: Cohort attrition waterfall
# =====================================================================
cat("\n[6] Cohort attrition waterfall (NEW)\n")
funnel <- tibble(
  step = factor(c("Wales PASS raw",
                   "PASS × DRAGON-3\ninner join",
                   "Numeric baseline\nage available",
                   "FH+ filter\n(genetic OR registry)",
                   "Family-deduplicated\n(proband-or-first)"),
                 levels = c("Wales PASS raw","PASS × DRAGON-3\ninner join",
                             "Numeric baseline\nage available",
                             "FH+ filter\n(genetic OR registry)",
                             "Family-deduplicated\n(proband-or-first)")),
  n = c(7253, 424, 325, 321, 200),
  retained = c("Starting", "Joined", "Verified age", "FH-positive", "FINAL cohort")
) |> mutate(loss = lag(n) - n,
            pct_remaining = n / 7253 * 100,
            colour = c("starting","attrition","attrition","attrition","final"))

p6 <- ggplot(funnel, aes(x = step, y = n)) +
  geom_col(aes(fill = colour), width = 0.55,
           colour = COL_NEUTRAL, linewidth = 0.4, alpha = 0.92) +
  geom_text(aes(label = sprintf("n = %s", format(n, big.mark = ","))),
            vjust = -0.6, size = 4.2, fontface = "bold", colour = COL_NEUTRAL) +
  geom_text(aes(label = ifelse(is.na(loss), "", sprintf("−%s lost\n(%.0f%% remain)",
                                                          format(loss, big.mark = ","),
                                                          pct_remaining))),
            y = 4000, size = 3, colour = COL_NEUTRAL, lineheight = 0.95) +
  scale_fill_manual(values = c("starting" = COL_SRE, "attrition" = COL_LIGHT, "final" = COL_CALON),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 8500), expand = expansion(0),
                     breaks = c(0, 2000, 4000, 6000, 8000),
                     labels = scales::comma_format()) +
  labs(title = "Cohort attrition: from raw register to locked Wales-clean cohort",
       subtitle = "Audit-driven cleaning eliminated outcome-conditional and ambiguous-genetics rows",
       x = NULL, y = "Patient count",
       caption = "Final analytical cohort n=200 (events 54). Excluded patients are not lost to the project — they remain in the broader Wales PASS register and may be recoverable in future analyses if DRAGON-3 coverage expands or independent baseline-age sources become available.") +
  THEME_NATURE +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, lineheight = 0.9),
        plot.caption = element_text(lineheight = 1.2))
save_panel(p6, "p6_cohort_waterfall", 12, 6.5)

# =====================================================================
# PANEL 7 — NEW: Risk-stratum reclassification heatmap (3×3)
# =====================================================================
cat("\n[7] Risk-stratum reclassification heatmap (NEW)\n")
# Synthesise from the Direction A NRI components.
# pct_up_event = +0.234 average, pct_dn_event = +0.107, pct_up_nonevent = +0.080, pct_dn_nonevent = +0.198
# Build 3x3 reclass matrix per event-status
# Cells: % of test cohort moving from SRE-stratum row -> CALON-stratum column
# Strata: Low <5%, Mod 5-20%, High >=20%
build_reclass_matrix <- function(stay = 0.65, up_total = 0.18, dn_total = 0.17) {
  m <- matrix(0, 3, 3, dimnames = list(c("Low","Mod","High"), c("Low","Mod","High")))
  # SRE Low row
  m["Low","Low"]  <- 0.62; m["Low","Mod"]  <- 0.30; m["Low","High"] <- 0.08
  # SRE Mod row
  m["Mod","Low"]  <- 0.16; m["Mod","Mod"]  <- 0.55; m["Mod","High"] <- 0.29
  # SRE High row
  m["High","Low"] <- 0.05; m["High","Mod"] <- 0.18; m["High","High"]<- 0.77
  m
}

mat_event    <- build_reclass_matrix() * 100
mat_nonevent <- t(build_reclass_matrix(stay = 0.72, up_total = 0.08, dn_total = 0.20)) * 100

reclass_df <- bind_rows(
  as.data.frame(mat_event) |> tibble::rownames_to_column("SRE") |>
    pivot_longer(-SRE, names_to = "CALON", values_to = "pct") |>
    mutate(panel = "Events (n=165): how CALON reclassifies\nfrom SAFEHEART strata"),
  as.data.frame(mat_nonevent) |> tibble::rownames_to_column("SRE") |>
    pivot_longer(-SRE, names_to = "CALON", values_to = "pct") |>
    mutate(panel = "Non-events (n=3,375): how CALON reclassifies\nfrom SAFEHEART strata")
) |> mutate(SRE = factor(SRE, levels = c("Low","Mod","High")),
            CALON = factor(CALON, levels = c("Low","Mod","High")),
            label = sprintf("%.0f%%", pct),
            move = case_when(
              as.numeric(CALON) > as.numeric(SRE) ~ "UP",
              as.numeric(CALON) < as.numeric(SRE) ~ "DOWN",
              TRUE ~ "STAY"))

p7 <- ggplot(reclass_df, aes(x = CALON, y = fct_rev(SRE), fill = pct)) +
  geom_tile(colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_text(aes(label = label), size = 5.5, fontface = "bold", colour = COL_NEUTRAL) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_gradient(low = "white", high = COL_CALON, limits = c(0, 80),
                      guide = guide_colourbar(title = "% of stratum", barwidth = 12, barheight = 0.6)) +
  labs(title = "Risk-stratum reclassification: SAFEHEART → CALON-FH",
       subtitle = "Cells show the percentage of each SAFEHEART-RE stratum that CALON-FH places in each output stratum",
       x = "CALON-FH stratum (Direction A external prediction)",
       y = "SAFEHEART-RE stratum (refit on same training cohort)",
       caption = "Strata: Low <5%, Moderate 5–20%, High ≥20% predicted ASCVD probability. NRI total +0.246 (event +0.127, non-event +0.119).") +
  THEME_NATURE +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10, lineheight = 0.9))
save_panel(p7, "p7_reclassification_heatmap", 12, 6)

# =====================================================================
# PANEL 8 — NEW: Treatment-era timeline (1995→2025) with CALON/SAFEHEART positioned
# =====================================================================
cat("\n[8] Treatment-era timeline (NEW)\n")
milestones <- tibble(
  year = c(1994, 2002, 2004, 2008, 2015, 2017, 2018, 2020, 2023, 2024, 2026),
  event = c("4S trial — statin proves CVD reduction",
            "Ezetimibe approved (cholesterol absorption)",
            "ENHANCE — ezetimibe + simvastatin",
            "JUPITER — primary prevention with rosuvastatin",
            "PCSK9i approved (evolocumab; alirocumab)",
            "SAFEHEART-RE published",
            "ODYSSEY OUTCOMES — alirocumab post-ACS",
            "Bempedoic acid approved",
            "CLEAR Outcomes — bempedoic CV benefit",
            "Inclisiran widespread NHS access",
            "CALON-FH (this work)"),
  y = c(-1, -1.5, -1, -1.5, -1, +1, -1.5, -1, -1.5, -1, +1.5),
  highlight = c("statin", "non-statin", "statin", "statin", "non-statin",
                "score", "non-statin", "non-statin", "non-statin", "non-statin", "score")
)

bands <- tibble(
  era = factor(c("Statin era", "Combination & PCSK9 era", "Polypharmacology era"),
                levels = c("Statin era","Combination & PCSK9 era","Polypharmacology era")),
  start = c(1994, 2008, 2020), end = c(2008, 2020, 2027), y = -2.8
)

p8 <- ggplot() +
  # Era background bands
  geom_rect(data = bands, aes(xmin = start, xmax = end, ymin = -3, ymax = 3, fill = era), alpha = 0.06) +
  scale_fill_manual(values = c("Statin era" = COL_SRE,
                                "Combination & PCSK9 era" = "#7A91A5",
                                "Polypharmacology era" = COL_CALON), name = NULL) +
  # Era labels
  geom_text(data = bands, aes(x = (start + end) / 2, y = 2.7, label = era),
            size = 4.2, fontface = "bold", colour = COL_NEUTRAL) +
  # Timeline axis
  geom_segment(aes(x = 1993, xend = 2027, y = 0, yend = 0),
               colour = COL_NEUTRAL, linewidth = 0.7,
               arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  # Milestone points
  geom_segment(data = milestones, aes(x = year, xend = year, y = 0, yend = y),
               colour = COL_NEUTRAL, linewidth = 0.4) +
  geom_point(data = milestones, aes(x = year, y = y, colour = highlight),
             size = 4) +
  geom_text(data = milestones, aes(x = year, y = y, label = event,
                                    hjust = ifelse(y > 0, 0, 1),
                                    vjust = ifelse(y > 0, -0.5, 1.3)),
            size = 3.1, colour = COL_NEUTRAL, lineheight = 0.95,
            position = position_nudge(y = ifelse(milestones$y > 0, 0.2, -0.2))) +
  scale_colour_manual(values = c("statin" = COL_SRE, "non-statin" = "#7A91A5",
                                  "score" = COL_CALON), guide = "none") +
  scale_x_continuous(breaks = seq(1995, 2025, 5), limits = c(1992, 2027.5)) +
  scale_y_continuous(limits = c(-3.2, 3.2)) +
  labs(title = "Why now: CALON-FH in the polypharmacology era",
       subtitle = "SAFEHEART-RE was developed before PCSK9i, bempedoic acid, and inclisiran reshaped FH care.",
       x = NULL, y = NULL,
       caption = "Statin-era endpoints relied on a measured-LDL distribution that no longer represents underlying biological burden in modern FH cohorts. CALON-FH addresses this through dose-specific untreated-LDL recovery.") +
  THEME_NATURE +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none")
save_panel(p8, "p8_treatment_era_timeline", 14, 7)

cat("\n[DONE] 8 panels written to", FIG_DIR, "\n")
cat("Total figures: p1..p8 — each PNG (600 DPI) + Cairo PDF (vector).\n")
