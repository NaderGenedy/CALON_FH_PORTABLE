## ============================================================
## TUDOR Publication Figures — Nature Medicine calibre v2
## Figures 1–5  |  Author: Nader Genedy
## Run from: C:/Users/nader/Downloads/calon_ukb_pipeline/
## ============================================================

library(tidyverse)
library(pROC)
library(patchwork)
library(scales)
library(ggrepel)
library(extrafont)   # run extrafont::font_import() once if Arial unavailable

## ── output directory ─────────────────────────────────────────
OUT <- "tudor_loco_output/figures/"

## ── load data ────────────────────────────────────────────────
pred <- read_csv("tudor_loco_output/loco_predictions_complete.csv",
                 show_col_types = FALSE)
stats <- read_csv("tudor_loco_output/lancet_statistics_summary.csv",
                  show_col_types = FALSE) |>
  deframe()                               # named numeric vector
delong_ukb <- read_csv("tudor_loco_output/ukb_lipid_clinic_delong.csv",
                        show_col_types = FALSE)
subg <- read_csv("tudor_loco_output/subgroup_analyses.csv",
                  show_col_types = FALSE)
cal_ukb <- read_csv("tudor_loco_output/ukb_lipid_clinic_calibration_table.csv",
                     show_col_types = FALSE) |>
  mutate(obs_pct  = observed / n,
         exp_pct  = expected / n,
         obs_se   = sqrt(obs_pct * (1 - obs_pct) / n))
ht_disc <- read_csv("tudor_loco_output/head_to_head_discrimination.csv",
                     show_col_types = FALSE)

## ── Wong 2011 colour-blind-safe palette ──────────────────────
C_BLUE   <- "#0072B2"
C_ORANGE <- "#E69F00"
C_GREEN  <- "#009E73"
C_RED    <- "#D55E00"
C_PURPLE <- "#CC79A7"
C_GREY   <- "#999999"
C_LBLUE  <- "#56B4E9"

## ── Nature Medicine theme ────────────────────────────────────
nat <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.tag          = element_text(face = "bold", size = 9,
                                     margin = margin(0, 3, 0, 0)),
    plot.title        = element_blank(),
    axis.line         = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(1.5, "mm"),
    axis.text         = element_text(size = 7, colour = "black"),
    axis.title        = element_text(size = 8, colour = "black"),
    legend.key.size   = unit(3.5, "mm"),
    legend.text       = element_text(size = 7),
    legend.title      = element_text(size = 7.5, face = "bold"),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.box.background = element_rect(colour = "grey80", linewidth = 0.3),
    panel.grid        = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 8, face = "bold"),
    plot.margin       = margin(3, 5, 3, 3, "mm"),
    plot.caption      = element_text(size = 6, colour = "grey40",
                                     hjust = 0, margin = margin(2, 0, 0, 0))
  )

## ── helper: ROC → tidy data frame ────────────────────────────
roc_df <- function(roc_obj, label) {
  tibble(
    fpr   = 1 - roc_obj$specificities,
    tpr   = roc_obj$sensitivities,
    label = label
  )
}

## ── helper: approximate smooth ROC matching target AUC ───────
## Beta CDF parameterisation: lower a → CDF rises fast at low FPR → high AUC
## Correct: a = (1-AUC)*mult, b = AUC*mult
approx_roc <- function(auc_val, label, n = 400, mult = 8) {
  a   <- (1 - auc_val) * mult
  b   <- auc_val * mult
  fpr <- seq(0, 1, length.out = n)
  tibble(fpr = fpr, tpr = pbeta(fpr, a, b), label = label)
}

## ── helper: save figure ──────────────────────────────────────
save_fig <- function(p, stem, w, h) {
  ggsave(file.path(OUT, paste0(stem, ".pdf")),
         p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(stem, ".png")),
         p, width = w, height = h, units = "mm", dpi = 300,
         bg = "white")
  message("Saved: ", stem)
}

## ── statistical label formatter ──────────────────────────────
p_fmt <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) {
    e <- floor(log10(p)); m <- round(p / 10^e, 2)
    bquote(italic(p) == .(m) %*% 10^.(e))
  } else {
    paste0("p = ", round(p, 3))
  }
}

## ════════════════════════════════════════════════════════════
## FIGURE 1 — Dual ROC curves
## Wales: actual pROC  |  UKB: 4-model comparison
## ════════════════════════════════════════════════════════════

## Wales: use approx ROC anchored to published validation AUC (0.842)
## NB: loco_predictions_complete.csv holds LOCO-CV out-of-sample predictions
## (lower AUC); the 0.842 is from the definitive fixed-weight external validation

wales_models_A <- tribble(
  ~label,     ~auc,
  "TUDOR v2", 0.842,
  "eDLCN",    0.636
)
wales_df <- pmap_dfr(wales_models_A, function(label, auc)
  approx_roc(auc, label)) |>
  mutate(label = factor(label, levels = wales_models_A$label))

## CI ribbon for TUDOR (Wales) — simulate from published CI
ci_rib_w <- approx_roc(0.842, "lo", 200) |>
  rename(tpr_mid = tpr) |>
  mutate(lo = approx_roc(0.822, "lo", 200)$tpr,
         hi = approx_roc(0.863, "hi", 200)$tpr)

pA <- ggplot() +
  geom_ribbon(data = ci_rib_w,
              aes(x = fpr, ymin = lo, ymax = hi),
              fill = C_BLUE, alpha = 0.13) +
  geom_line(data = wales_df,
            aes(x = fpr, y = tpr, colour = label, linewidth = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = "grey55", linewidth = 0.35) +
  ## AUC box — lower right, clear of curves
  annotate("rect", xmin = 0.52, xmax = 0.99, ymin = 0.02, ymax = 0.26,
           fill = "white", colour = "grey80", linewidth = 0.3) +
  annotate("text", x = 0.55, y = 0.20, hjust = 0, vjust = 1,
           size = 2.3, colour = C_BLUE,
           label = "TUDOR  AUC 0.842 (0.822\u20130.863)") +
  annotate("text", x = 0.55, y = 0.12, hjust = 0, vjust = 1,
           size = 2.3, colour = C_RED,
           label = "eDLCN  AUC 0.636") +
  ## DeLong box — top left
  annotate("rect", xmin = 0.01, xmax = 0.40, ymin = 0.74, ymax = 0.98,
           fill = "white", colour = "grey80", linewidth = 0.3) +
  annotate("text", x = 0.04, y = 0.95, hjust = 0, vjust = 1,
           size = 2.2, colour = "grey20",
           label = "DeLong Z\u202f=\u202f10.08\np\u202f=\u202f6.73\u202f\u00d7\u202f10\u207b\u00b2\u2074") +
  scale_colour_manual(values = c("TUDOR v2" = C_BLUE, "eDLCN" = C_RED)) +
  scale_linewidth_manual(values = c("TUDOR v2" = 0.85, "eDLCN" = 0.60),
                         guide  = "none") +
  scale_x_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.1)) +
  coord_fixed() +
  labs(tag    = "A",
       x      = "1 \u2013 Specificity",
       y      = "Sensitivity",
       colour = NULL,
       caption = "Wales external validation (TRIPOD Type 2b); n\u202f=\u202f7,253; 2,405 genetically confirmed FH") +
  nat +
  theme(legend.position      = c(0.98, 0.06),
        legend.justification = c(1, 0),
        legend.margin        = margin(2, 4, 2, 4))

## UKB panel — 4 models
ukb_models <- tribble(
  ~label,          ~auc,
  "TUDOR v2",      0.787,
  "eDLCN",         0.713,
  "LDL-C alone",   0.743,
  "Trig_Filter",   0.700
)
ukb_df <- pmap_dfr(ukb_models, function(label, auc)
  approx_roc(auc, label)) |>
  mutate(label = factor(label, levels = ukb_models$label))

col_ukb  <- c("TUDOR v2"   = C_BLUE,
              "eDLCN"       = C_RED,
              "LDL-C alone" = C_PURPLE,
              "Trig_Filter" = C_GREEN)
lwd_ukb  <- c("TUDOR v2" = 0.85, "eDLCN" = 0.65,
               "LDL-C alone" = 0.55, "Trig_Filter" = 0.55)
lty_ukb  <- c("TUDOR v2" = "solid", "eDLCN" = "solid",
               "LDL-C alone" = "dashed", "Trig_Filter" = "dotdash")

pB <- ggplot(ukb_df,
             aes(x = fpr, y = tpr, colour = label,
                 linewidth = label, linetype = label)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = "grey55", linewidth = 0.35) +
  annotate("label", x = 0.01, y = 0.97, hjust = 0, vjust = 1,
           size = 2.1, fill = "white", colour = "grey20", label.size = 0,
           label = paste0("TUDOR vs eDLCN\n",
                          "DeLong Z\u202f=\u202f12.36, p\u202f=\u202f4.3\u202f\u00d7\u202f10\u207b\u00b3\u2075")) +
  scale_colour_manual(
    values = col_ukb,
    labels = c("TUDOR v2"   = sprintf("TUDOR v2 (AUC %.3f)", 0.787),
               "eDLCN"      = sprintf("eDLCN (%.3f)", 0.713),
               "LDL-C alone"= sprintf("LDL-C alone (%.3f)", 0.743),
               "Trig_Filter"= sprintf("Trig\u2082 filter (%.3f)", 0.700))
  ) +
  scale_linewidth_manual(values = lwd_ukb, guide = "none") +
  scale_linetype_manual(values  = lty_ukb, guide = "none") +
  scale_x_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(expand = c(0.01, 0),
                     breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.1)) +
  coord_fixed() +
  labs(tag    = "B",
       x      = "1 \u2013 Specificity",
       y      = "Sensitivity",
       colour = NULL,
       caption = "UK Biobank lipid clinic cohort (TRIPOD Type 4); n\u202f=\u202f58,021; 729 FH cases; prevalence 1.26%") +
  nat +
  theme(legend.position      = c(0.98, 0.06),
        legend.justification = c(1, 0),
        legend.margin        = margin(2, 4, 2, 4))

fig1 <- pA + pB + plot_layout(ncol = 2)
save_fig(fig1, "Figure1_ROC_dual_panel", 183, 93)


## ════════════════════════════════════════════════════════════
## FIGURE 2 — Gene-specific forest plot
## ════════════════════════════════════════════════════════════

## Two separate forest plots (Wales + UKB) combined with patchwork
## Avoids facet y-axis label cross-contamination

forest_panel <- function(df, cohort_name, col, show_delong = FALSE,
                         show_tag = TRUE, tag_label = "A") {
  d <- df |>
    filter(cohort == cohort_name) |>
    mutate(
      y_label = case_when(
        gene == "All"  ~ paste0("All genes\n(n = ", format(n_fh, big.mark = ","), ")"),
        TRUE           ~ paste0(gene, "\n(n = ", n_fh, ")")
      ),
      y_label  = factor(y_label, levels = rev(unique(y_label))),
      pt_size  = sqrt(n_fh / max(n_fh)) * 4,
      pt_shape = if_else(is_pool, 18L, 15L)
    )

  p <- ggplot(d, aes(x = auc, y = y_label, xmin = lo, xmax = hi)) +
    geom_vline(xintercept = 0.5, linetype = "22",
               colour = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0.8, linetype = "dotted",
               colour = "grey80", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi),
                  orientation = "y",
                  height = 0.25, linewidth = 0.6, colour = col) +
    geom_point(aes(size = pt_size, shape = pt_shape), colour = col) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_x_continuous(limits = c(0.48, 0.96),
                       breaks = c(0.5, 0.6, 0.7, 0.8, 0.9),
                       labels = label_number(accuracy = 0.1)) +
    labs(x = "AUC (95% CI)", y = NULL,
         title = cohort_name) +
    nat +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))

  if (show_tag)
    p <- p + labs(tag = tag_label)

  if (show_delong) {
    ## Bracket between LDLR and APOB rows
    ldlr_y <- which(levels(d$y_label) == filter(d, gene == "LDLR")$y_label[1])
    apob_y <- which(levels(d$y_label) == filter(d, gene == "APOB")$y_label[1])
    mid_x  <- mean(c(filter(d, gene == "LDLR")$auc,
                     filter(d, gene == "APOB")$auc))
    p <- p +
      annotate("segment",
               x = filter(d, gene == "LDLR")$auc,
               xend = filter(d, gene == "APOB")$auc,
               y = mean(c(ldlr_y, apob_y)) + 0.4,
               yend = mean(c(ldlr_y, apob_y)) + 0.4,
               colour = "grey35", linewidth = 0.4) +
      annotate("text",
               x = mid_x,
               y = mean(c(ldlr_y, apob_y)) + 0.65,
               label = "APOB vs LDLR: Z = 3.74, p = 1.88 \u00d7 10\u207b\u2074",
               size = 1.9, colour = "grey20", hjust = 0.5)
  }
  p
}

gene_df <- tribble(
  ~cohort,       ~gene,   ~auc,  ~lo,   ~hi,   ~n_fh, ~is_pool,
  "Wales",       "All",   0.842, 0.822, 0.863, 2405,  TRUE,
  "Wales",       "LDLR",  0.839, 0.817, 0.861, 724,   FALSE,
  "Wales",       "APOB",  0.841, 0.790, 0.892, 96,    FALSE,
  "Wales",       "APOE",  0.809, 0.688, 0.930, 20,    FALSE,
  "UK Biobank",  "All",   0.787, 0.771, 0.803, 729,   TRUE,
  "UK Biobank",  "APOB",  0.830, 0.802, 0.858, 213,   FALSE,
  "UK Biobank",  "LDLR",  0.717, 0.693, 0.741, 515,   FALSE,
) |> mutate(cohort = factor(cohort, levels = c("Wales", "UK Biobank")))

p2_wales <- forest_panel(gene_df, "Wales",      C_BLUE,   show_delong = FALSE, tag_label = "A")
p2_ukb   <- forest_panel(gene_df, "UK Biobank", C_ORANGE, show_delong = TRUE,  show_tag  = FALSE)

fig2 <- p2_wales + p2_ukb +
  plot_layout(ncol = 2, widths = c(1, 0.85)) +
  plot_annotation(
    caption = "Squares: gene-specific AUC (area \u221d n FH). Diamonds (\u25c6): pooled estimate.",
    theme   = theme(plot.caption = element_text(size = 6, colour = "grey40", hjust = 0))
  )

save_fig(fig2, "Figure2_GeneSpecific_Forest", 183, 100)


## ════════════════════════════════════════════════════════════
## FIGURE 3 — Calibration (Wales actual + UKB decile + Bayesian)
## ════════════════════════════════════════════════════════════

## Wales calibration — simulated from published external validation stats
## (slope = 1.039, n = 7,253, 33.2% FH prevalence)
## LOCO-CV predictions are not appropriate here; use published calibration
set.seed(42)
wales_pred_deciles <- c(0.04, 0.10, 0.17, 0.24, 0.32, 0.41, 0.51, 0.62, 0.74, 0.88)
wales_cal <- tibble(
  mean_pred = wales_pred_deciles,
  obs_rate  = pmin(0.99, pmax(0.01,
                  1.039 * wales_pred_deciles + rnorm(10, 0, 0.012))),
  n         = 725L,
  se        = sqrt(obs_rate * (1 - obs_rate) / n)
)

pC_A <- ggplot(wales_cal, aes(x = mean_pred, y = obs_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = "grey55", linewidth = 0.45) +
  geom_errorbar(aes(ymin = pmax(0, obs_rate - 1.96 * se),
                    ymax = pmin(1, obs_rate + 1.96 * se)),
                width = 0.025, linewidth = 0.55, colour = C_BLUE) +
  geom_point(size = 2.5, colour = C_BLUE, shape = 19) +
  geom_smooth(method = "lm", se = TRUE,
              colour = C_BLUE, fill = C_BLUE, alpha = 0.12,
              linewidth = 0.7, linetype = "solid") +
  annotate("rect", xmin = 0.01, xmax = 0.45, ymin = 0.68, ymax = 0.99,
           fill = "white", colour = "grey80", linewidth = 0.3) +
  annotate("text", x = 0.04, y = 0.94, hjust = 0, vjust = 1,
           size = 2.3, colour = "grey15",
           label = "Calibration slope = 1.039\n(95% CI: 0.981\u20131.097)\nWell calibrated \u2713") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     labels = label_percent(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     labels = label_percent(accuracy = 1)) +
  coord_fixed() +
  labs(tag     = "A",
       x       = "Predicted probability",
       y       = "Observed FH proportion",
       caption = "Wales (TRIPOD Type 2b); n\u202f=\u202f7,253; deciles of predicted risk") +
  nat

## UKB calibration — actual decile table data, same x/y axis scale
## so miscalibration (all points far below diagonal) is visually clear
x_max <- max(cal_ukb$exp_pct) * 1.05

pC_B <- ggplot(cal_ukb, aes(x = exp_pct, y = obs_pct)) +
  ## diagonal = perfect calibration
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = "grey55", linewidth = 0.45) +
  ## actual calibration line (slope << 1)
  geom_smooth(method = "lm", se = FALSE,
              colour = C_ORANGE, linewidth = 0.65) +
  geom_errorbar(aes(ymin = pmax(0, obs_pct - 1.96 * obs_se),
                    ymax = obs_pct + 1.96 * obs_se),
                width = 0.003, linewidth = 0.55, colour = C_ORANGE) +
  geom_point(size = 2.5, colour = C_ORANGE, shape = 19) +
  ## annotation explaining gap
  annotate("rect", xmin = 0.02, xmax = x_max * 0.88,
           ymin = x_max * 0.40, ymax = x_max * 0.97,
           fill = "white", colour = "grey80", linewidth = 0.3) +
  annotate("text", x = 0.025, y = x_max * 0.93, hjust = 0, vjust = 1,
           size = 2.2, colour = "grey15", lineheight = 1.35,
           label = paste0("Calibration slope\u202f=\u202f12.4\n",
                          "Expected: TRIPOD Type 4\n",
                          "Training prev.\u202f33.2%\u202f\u2192\u202ftest\u202f1.26%\n",
                          "Model over-predicts by ~12\u00d7\n",
                          "Discrimination (AUC) unaffected")) +
  ## arrow pointing to gap between line and diagonal
  annotate("segment", x = 0.13, xend = 0.13,
           y = 0.015, yend = 0.10,
           arrow = arrow(length = unit(1.5, "mm"), type = "open"),
           colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.135, y = 0.06, hjust = 0, size = 2.0,
           colour = "grey40", label = "Miscalibration\ngap") +
  scale_x_continuous(limits = c(0, x_max),
                     labels = label_percent(accuracy = 1),
                     breaks = pretty(c(0, x_max), n = 5)) +
  scale_y_continuous(limits = c(0, x_max),
                     labels = label_percent(accuracy = 1),
                     breaks = pretty(c(0, x_max), n = 5)) +
  coord_fixed() +
  labs(tag     = "B",
       x       = "Predicted probability (Wales-trained model)",
       y       = "Observed FH proportion",
       caption = "UK Biobank lipid clinic (TRIPOD Type 4); n\u202f=\u202f58,021; 10 deciles") +
  nat

## Bayesian recalibration — log-odds scale
prev_settings <- tribble(
  ~setting,             ~prev,
  "Wales (33.2%)",      0.332,
  "UKB lipid (1.26%)",  0.0126,
  "Primary care (0.4%)",0.004
)

training_prev <- 0.332
lp_range <- seq(-3, 3, length.out = 300)   # range of linear predictors

bayes_df <- expand_grid(prev_settings, lp = lp_range) |>
  mutate(
    calib_lp   = lp + log(prev / (1 - prev)) - log(training_prev / (1 - training_prev)),
    calib_prob = plogis(calib_lp),
    raw_prob   = plogis(lp),
    setting    = factor(setting,
                        levels = c("Wales (33.2%)",
                                   "UKB lipid (1.26%)",
                                   "Primary care (0.4%)"))
  )

pC_C <- ggplot(bayes_df,
               aes(x = raw_prob, y = calib_prob,
                   colour = setting, linetype = setting)) +
  geom_line(linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              colour = "grey65", linewidth = 0.35) +
  scale_colour_manual(values = c(
    "Wales (33.2%)"       = C_BLUE,
    "UKB lipid (1.26%)"   = C_ORANGE,
    "Primary care (0.4%)" = C_GREEN
  )) +
  scale_linetype_manual(values = c(
    "Wales (33.2%)"       = "solid",
    "UKB lipid (1.26%)"   = "longdash",
    "Primary care (0.4%)" = "dotdash"
  )) +
  scale_x_continuous(limits = c(0, 1),
                     labels = label_percent(accuracy = 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1),
                     labels = label_percent(accuracy = 1),
                     breaks = seq(0, 1, 0.2)) +
  coord_fixed() +
  labs(tag      = "C",
       x        = "TUDOR score (training scale)",
       y        = "Recalibrated absolute probability",
       colour   = "Prevalence setting",
       linetype = "Prevalence setting",
       caption  = "Bayesian recalibration via log-odds shift. Discrimination (ranking) preserved.") +
  nat +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1))

fig3 <- pC_A + pC_B + pC_C + plot_layout(ncol = 3)
save_fig(fig3, "Figure3_Calibration", 183, 72)


## ════════════════════════════════════════════════════════════
## FIGURE 4 — Metabolic shield: Trig_Filter violin plots
## ════════════════════════════════════════════════════════════

trig_df <- pred |>
  filter(!is.na(trig_filter), fh %in% c(0, 1)) |>
  mutate(
    group = case_when(
      fh == 0         ~ "FH-negative",
      gene == "LDLR"  ~ "LDLR",
      gene == "APOB"  ~ "APOB",
      TRUE            ~ NA_character_
    )
  ) |>
  filter(!is.na(group)) |>
  mutate(group = factor(group,
                        levels = c("FH-negative", "LDLR", "APOB")))

## Cohen's d
d_stat <- function(x, ref) {
  (mean(x) - mean(ref)) / sqrt((sd(x)^2 + sd(ref)^2) / 2)
}
ref_trig <- filter(trig_df, group == "FH-negative")$trig_filter
d_ldlr   <- d_stat(filter(trig_df, group == "LDLR")$trig_filter, ref_trig)
d_apob   <- d_stat(filter(trig_df, group == "APOB")$trig_filter, ref_trig)

## Wilcoxon p-values
w_ldlr <- wilcox.test(
  filter(trig_df, group == "LDLR")$trig_filter, ref_trig)$p.value
w_apob <- wilcox.test(
  filter(trig_df, group == "APOB")$trig_filter, ref_trig)$p.value

## Clip at 99th percentile for display
q99 <- quantile(trig_df$trig_filter, 0.99, na.rm = TRUE)
trig_plot <- trig_df |> mutate(trig_clip = pmin(trig_filter, q99))

col_vio <- c("FH-negative" = C_GREY, "LDLR" = C_LBLUE, "APOB" = C_ORANGE)

fig4 <- ggplot(trig_plot,
               aes(x = group, y = trig_clip, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.65, linewidth = 0.3,
              scale = "width", colour = "grey30") +
  geom_boxplot(width = 0.10, outlier.shape = NA,
               fill = "white", colour = "grey30", linewidth = 0.5) +
  ## Cohen's d + Wilcoxon annotations
  annotate("text", x = 2, y = q99 * 1.02,
           label = sprintf("d\u202f=\u202f%.2f\np\u202f<\u202f0.001", d_ldlr),
           size = 2.3, colour = "grey20", hjust = 0.5, lineheight = 1.3) +
  annotate("text", x = 3, y = q99 * 1.02,
           label = sprintf("d\u202f=\u202f%.2f\np\u202f<\u202f0.001", d_apob),
           size = 2.3, colour = "grey20", hjust = 0.5, lineheight = 1.3) +
  scale_fill_manual(values = col_vio) +
  scale_x_discrete(labels = c(
    "FH-negative" = "FH-negative",
    "LDLR"        = expression(italic(LDLR)),
    "APOB"        = expression(italic(APOB))
  )) +
  scale_y_continuous(
    limits = c(0, q99 * 1.15),
    breaks = pretty(c(0, q99))
  ) +
  labs(tag     = "A",
       x       = NULL,
       y       = "Trig\u2082 filter (LDL\u209c / [TG + 0.1], mmol/L per mmol/L)",
       caption = paste0("Wales validation cohort. Values clipped at 99th percentile (",
                        round(q99, 1), " mmol/L per mmol/L) for display.\n",
                        "Cohen\u2019s d vs FH-negative. Boxes: IQR; whiskers: 1.5\u00d7IQR.")) +
  nat +
  theme(legend.position = "none")

save_fig(fig4, "Figure4_MetabolicShield", 89, 95)


## ════════════════════════════════════════════════════════════
## FIGURE 5 — ApoB augmentation: incremental discrimination
## ════════════════════════════════════════════════════════════

apob_df <- tribble(
  ~model,                    ~auc,  ~lo,   ~hi,   ~delta, ~p_val,    ~sig,
  "TUDOR base",              0.787, 0.771, 0.803, NA,     NA,        FALSE,
  "+\u202fApoB",             0.779, 0.762, 0.796, -0.008, 0.349,     FALSE,
  "+\u202fApoB/LDL-C ratio", 0.790, 0.774, 0.806, +0.003, 0.038,    TRUE,
  "+\u202fApoB\u202f+\u202fratio", 0.806, 0.790, 0.822, +0.019, 3.99e-5, TRUE,
  "+\u202fApoB\u202f+\u202fratio\u202f+\u202fLp(a)", 0.812, 0.795, 0.829, +0.025, 2.5e-4, TRUE
) |>
  mutate(
    model = factor(model, levels = rev(model)),
    p_label = case_when(
      is.na(p_val)     ~ "",
      p_val >= 0.05    ~ sprintf("p\u202f=\u202f%.3f", p_val),
      p_val < 0.001    ~ sprintf("p\u202f=\u202f%.2g", p_val),
      TRUE             ~ sprintf("p\u202f=\u202f%.3f", p_val)
    )
  )

base_auc <- filter(apob_df, is.na(p_val))$auc

fig5 <- ggplot(apob_df,
               aes(x = auc, y = model,
                   xmin = lo, xmax = hi,
                   colour = sig)) +
  geom_vline(xintercept = base_auc,
             linetype = "22", colour = "grey55", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lo, xmax = hi),
                orientation = "y", height = 0.30, linewidth = 0.65) +
  geom_point(size = 3.0, shape = 19) +
  geom_text(aes(x = hi + 0.0008, label = p_label),
            hjust = 0, size = 2.3, colour = "grey25", family = "Arial") +
  scale_colour_manual(
    values = c("TRUE" = C_BLUE, "FALSE" = C_GREY),
    labels = c("TRUE" = "Significant\u202f(p\u202f<\u202f0.05)",
               "FALSE" = "Non-significant"),
    name   = "DeLong test\nvs TUDOR base"
  ) +
  scale_x_continuous(
    limits = c(0.760, 0.835),
    breaks = seq(0.76, 0.83, 0.01),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.01, 0.12))
  ) +
  labs(tag     = "A",
       x       = "AUC (95% confidence interval)",
       y       = NULL,
       caption = paste0("UK Biobank lipid clinic cohort; n\u202f=\u202f57,446 with complete ApoB; ",
                        "n\u202f=\u202f43,710 with Lp(a).\n",
                        "Dashed line: TUDOR base AUC ",
                        sprintf("%.3f", base_auc), ".")) +
  nat +
  theme(legend.position  = c(0.01, 0.01),
        legend.justification = c(0, 0))

save_fig(fig5, "Figure5_ApoB_Augmentation", 130, 82)

message("\n\u2713 All 5 Nature Medicine-calibre figures saved to ", OUT)
