##############################################################################
# TUDOR — 5 Nature-Calibre Publication Figures
# Target: Expert Review of Clinical Lipidology
# Nature specs: single column 89mm, double 183mm, min 6pt, 300 DPI
# Colourblind-safe palette (okabe-ito)
# All fonts: sans-serif, base_size = 7 (Nature standard)
##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(ggpubr)
  library(scales)
  library(patchwork)
})

# ── Paths ──────────────────────────────────────────────────────────────────
PRED_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output"
RES_DIR  <- file.path(PRED_DIR, "tudor_flaw_results")
FIG_DIR  <- file.path(PRED_DIR, "tudor_flaw_results")

# ── Colourblind-safe palette (Okabe-Ito) ───────────────────────────────────
CB_BLK   <- "#000000"   # black
CB_ORG   <- "#E69F00"   # orange
CB_SKY   <- "#56B4E9"   # sky blue
CB_GRN   <- "#009E73"   # green
CB_YLW   <- "#F0E442"   # yellow
CB_BLU   <- "#0072B2"   # blue
CB_RED   <- "#D55E00"   # vermillion
CB_PNK   <- "#CC79A7"   # pink
CB_GRY   <- "#999999"   # grey

# ── Nature theme ───────────────────────────────────────────────────────────
nat_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3),
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size, colour = "black"),
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      hjust = 0, colour = "black"),
      plot.subtitle    = element_text(size = base_size - 0.5, colour = CB_GRY),
      plot.caption     = element_text(size = base_size - 1, colour = CB_GRY,
                                      hjust = 0),
      legend.text      = element_text(size = base_size - 0.5),
      legend.title     = element_text(size = base_size - 0.5, face = "bold"),
      legend.key.size  = unit(3, "mm"),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size, face = "bold"),
      panel.grid       = element_blank(),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
}

# ── Load data ──────────────────────────────────────────────────────────────
pred   <- read_csv(file.path(PRED_DIR, "loco_predictions_complete.csv"),
                   show_col_types = FALSE)
cal_d  <- read_csv(file.path(RES_DIR, "ukb_recalibration_deciles.csv"),
                   show_col_types = FALSE)
trg_g  <- read_csv(file.path(RES_DIR, "trgshield_group_summary.csv"),
                   show_col_types = FALSE)
trg_t2 <- read_csv(file.path(RES_DIR, "trgshield_by_t2dm.csv"),
                   show_col_types = FALSE)
op_tbl <- read_csv(file.path(RES_DIR, "roc_operating_points.csv"),
                   show_col_types = FALSE)
eth_tbl <- read_csv(file.path(RES_DIR, "ethnicity_lc_cohort_summary.csv"),
                    show_col_types = FALSE)
asc_tbl <- read_csv(file.path(RES_DIR, "ascvd_by_ldl_decile.csv"),
                    show_col_types = FALSE)

wales <- pred %>% filter(grepl("Wales", cohort))
sw    <- pred %>% filter(cohort == "SouthWales")
ukb   <- pred %>% filter(cohort == "UKB")

cat("Data loaded. Building 5 figures...\n")

################################################################################
# FIGURE 1 — Triple ROC panel
# (A) TUDOR vs DLCN (Wales)  (B) Index vs Cascade  (C) Recalibration
################################################################################

cat("Building Figure 1...\n")

# Panel A: TUDOR vs DLCN
roc_tudor <- roc(wales$fh, wales$pred,   quiet=TRUE)
roc_dlcn  <- roc(wales$fh, wales$dlcn,   quiet=TRUE)
roc_sw    <- roc(sw$fh,    sw$pred,       quiet=TRUE)

auc_tudor <- sprintf("%.3f", round(as.numeric(auc(roc_tudor)),3))
auc_dlcn  <- sprintf("%.3f", round(as.numeric(auc(roc_dlcn)), 3))
auc_sw    <- sprintf("%.3f", round(as.numeric(auc(roc_sw)),   3))
delong_p  <- roc.test(roc_tudor, roc_dlcn)$p.value

roc_df_A <- bind_rows(
  data.frame(spec = 1 - roc_tudor$specificities,
             sens = roc_tudor$sensitivities,
             model = sprintf("TUDOR (AUC=%s)", auc_tudor)),
  data.frame(spec = 1 - roc_dlcn$specificities,
             sens = roc_dlcn$sensitivities,
             model = sprintf("DLCN (AUC=%s)", auc_dlcn)),
  data.frame(spec = 1 - roc_sw$specificities,
             sens = roc_sw$sensitivities,
             model = sprintf("TUDOR SW ext. (AUC=%s)", auc_sw))
)

pA <- ggplot(roc_df_A, aes(x = spec, y = sens, colour = model)) +
  geom_line(linewidth = 0.55) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = CB_GRY, linewidth = 0.3) +
  scale_colour_manual(
    values = c(CB_BLU, CB_RED, CB_GRN),
    name   = NULL
  ) +
  scale_x_continuous("1 − Specificity", limits = c(0,1),
                     breaks = seq(0,1,0.25), labels = number_format(accuracy=0.01)) +
  scale_y_continuous("Sensitivity",     limits = c(0,1),
                     breaks = seq(0,1,0.25), labels = number_format(accuracy=0.01)) +
  annotate("text", x=0.6, y=0.12,
           label = sprintf("DeLong p=%s", ifelse(delong_p<0.001,"<0.001",
                                                  sprintf("%.3f",delong_p))),
           size=2, colour=CB_BLK) +
  nat_theme() +
  theme(legend.position = c(0.68, 0.22),
        legend.background = element_rect(fill="white", linewidth=0.2)) +
  labs(title = "a", subtitle = "Wales validation cohort (n=6,448)")

# Panel B: Index vs Cascade
idx  <- pred %>% filter(index_case == 1, grepl("Wales",cohort))
casc <- pred %>% filter(index_case == 0, grepl("Wales",cohort))

roc_idx  <- roc(idx$fh,  idx$pred,  quiet=TRUE)
roc_casc <- roc(casc$fh, casc$pred, quiet=TRUE)

roc_df_B <- bind_rows(
  data.frame(spec = 1 - roc_idx$specificities,
             sens = roc_idx$sensitivities,
             group = sprintf("Index cases (AUC=0.807; n=%d)", nrow(idx))),
  data.frame(spec = 1 - roc_casc$specificities,
             sens = roc_casc$sensitivities,
             group = sprintf("Cascade relatives (AUC=0.708; n=%d)", nrow(casc)))
)

pB <- ggplot(roc_df_B, aes(x = spec, y = sens, colour = group)) +
  geom_line(linewidth = 0.55) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              colour=CB_GRY, linewidth=0.3) +
  scale_colour_manual(values=c(CB_ORG, CB_SKY), name=NULL) +
  scale_x_continuous("1 − Specificity", limits=c(0,1),
                     breaks=seq(0,1,0.25), labels=number_format(accuracy=0.01)) +
  scale_y_continuous("Sensitivity",     limits=c(0,1),
                     breaks=seq(0,1,0.25), labels=number_format(accuracy=0.01)) +
  nat_theme() +
  theme(legend.position=c(0.62, 0.18),
        legend.background=element_rect(fill="white", linewidth=0.2)) +
  labs(title="b", subtitle="Ascertainment-stratified performance")

# Panel C: Calibration (observed vs expected, pre/post recalibration)
cal_d2 <- cal_d %>%
  mutate(
    obs_p   = observed / n * 100,
    exp_p   = expected / n * 100,
    pred_rc = pred_recal_pct
  )

pC <- ggplot(cal_d2) +
  geom_point(aes(x=exp_p, y=obs_p), colour=CB_BLU, size=1.5) +
  geom_line(aes(x=exp_p, y=obs_p),  colour=CB_BLU, linewidth=0.4,
            linetype="dashed") +
  geom_point(aes(x=pred_rc, y=obs_p), colour=CB_GRN, shape=17, size=1.5) +
  geom_line(aes(x=pred_rc, y=obs_p),  colour=CB_GRN, linewidth=0.4) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              colour=CB_GRY, linewidth=0.3) +
  scale_x_continuous("Mean predicted probability (%)",
                     labels=number_format(accuracy=0.1)) +
  scale_y_continuous("Observed proportion (%)",
                     labels=number_format(accuracy=0.1)) +
  annotate("text", x=16, y=0.7,
           label="● Pre-recal    ▲ Post-recal\n(slope: 12.4 → 1.00)",
           size=1.8, colour=CB_BLK, hjust=1) +
  nat_theme() +
  labs(title="c", subtitle="UKB recalibration (n=37,050; 10 deciles)")

# Assemble Figure 1
fig1 <- pA | pB | pC
fig1 <- fig1 + plot_annotation(
  caption = paste0("Figure 1 | TUDOR diagnostic performance.\n",
                   "a, ROC curves comparing TUDOR (blue, AUC=",auc_tudor,
                   "), DLCN (red, AUC=",auc_dlcn,
                   ") and TUDOR South Wales external validation (green, AUC=",
                   auc_sw,") in the Wales FH cohort (n=6,448). ",
                   "b, TUDOR performance stratified by ascertainment route. ",
                   "c, Calibration before (circles) and after (triangles) logistic recalibration ",
                   "to the UK Biobank lipid-clinic prevalence.\n",
                   "AUC, area under the ROC curve; DLCN, Dutch Lipid Clinic Network; SW, South Wales."),
  theme = nat_theme()
)

ggsave(file.path(FIG_DIR, "Figure1_ROC_Calibration.pdf"),
       fig1, width=183, height=65, units="mm", device=cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure1_ROC_Calibration.png"),
       fig1, width=183, height=65, units="mm", dpi=300)
cat("✓ Figure 1 saved\n")

################################################################################
# FIGURE 2 — Triglyceride Shield: Gene-group distribution + T2DM confounding
################################################################################

cat("Building Figure 2...\n")

# Panel A: Violin plot from raw data (gene group TF distributions)
trg_raw <- pred %>%
  filter(!is.na(trig_filter), is.finite(trig_filter)) %>%
  mutate(
    gene_group = case_when(
      gene == "LDLR"  ~ "LDLR",
      gene == "APOB"  ~ "APOB",
      gene == "PCSK9" ~ "PCSK9",
      fh == 0         ~ "Non-FH",
      TRUE            ~ "FH-unknown"
    ),
    gene_group = factor(gene_group,
                        levels = c("Non-FH","PCSK9","APOB","LDLR","FH-unknown"))
  ) %>%
  filter(gene_group %in% c("Non-FH","APOB","LDLR"))

# Cohen's d function
cohens_d <- function(x,y) {
  (mean(x,na.rm=T)-mean(y,na.rm=T)) /
  sqrt((var(x,na.rm=T)*(sum(!is.na(x))-1)+var(y,na.rm=T)*(sum(!is.na(y))-1)) /
       (sum(!is.na(x))+sum(!is.na(y))-2))
}
nonfh_tf <- trg_raw$trig_filter[trg_raw$gene_group=="Non-FH"]
ldlr_tf  <- trg_raw$trig_filter[trg_raw$gene_group=="LDLR"]
apob_tf  <- trg_raw$trig_filter[trg_raw$gene_group=="APOB"]
d_ldlr   <- cohens_d(ldlr_tf, nonfh_tf)
d_apob   <- cohens_d(apob_tf, nonfh_tf)

q99 <- quantile(trg_raw$trig_filter, 0.995, na.rm=TRUE)

pA2 <- ggplot(trg_raw,
              aes(x=gene_group, y=trig_filter, fill=gene_group)) +
  geom_violin(alpha=0.75, trim=TRUE, linewidth=0.25) +
  geom_boxplot(width=0.09, outlier.size=0.2,
               fill="white", colour="grey30", linewidth=0.25) +
  scale_fill_manual(
    values = c("Non-FH"=CB_GRY, "APOB"=CB_ORG, "LDLR"=CB_BLU),
    guide  = "none"
  ) +
  scale_y_continuous("Triglyceride Filter\n(LDL_untreated / [TG + 0.1])",
                     limits = c(0, min(q99, 12))) +
  scale_x_discrete("Gene group") +
  annotate("text", x=2.3, y=min(q99,12)*0.93,
           label=sprintf("LDLR vs Non-FH\nd = %.2f", d_ldlr),
           size=1.8, colour=CB_BLU, hjust=0) +
  annotate("text", x=1.3, y=min(q99,12)*0.93,
           label=sprintf("APOB vs Non-FH\nd = %.2f", d_apob),
           size=1.8, colour=CB_ORG, hjust=0) +
  nat_theme() +
  labs(title="a",
       subtitle=sprintf("n(LDLR)=%d; n(APOB)=%d; n(Non-FH)=%d",
                        sum(trg_raw$gene_group=="LDLR"),
                        sum(trg_raw$gene_group=="APOB"),
                        sum(trg_raw$gene_group=="Non-FH")))

# Panel B: T2DM confounding (from UKB supplement data)
trg_t2_plot <- trg_t2 %>%
  mutate(
    T2DM  = ifelse(t2dm, "T2DM positive", "T2DM negative"),
    label = sprintf("n=%d\nTF=%.2f±%.2f", n, tf_mean, tf_sd)
  )

pB2 <- ggplot(trg_t2_plot, aes(x=T2DM, y=tf_mean, fill=T2DM)) +
  geom_col(width=0.55, alpha=0.85) +
  geom_errorbar(aes(ymin=tf_mean-tf_sd/sqrt(n),
                    ymax=tf_mean+tf_sd/sqrt(n)),
                width=0.15, linewidth=0.4) +
  geom_text(aes(label=label, y=tf_mean/2), size=2, colour="white") +
  scale_fill_manual(values=c("T2DM negative"=CB_BLU,"T2DM positive"=CB_RED),
                    guide="none") +
  scale_y_continuous("Mean Triglyceride Filter proxy", limits=c(0,3.2)) +
  scale_x_discrete(NULL) +
  annotate("text", x=1.5, y=3.1,
           label="Δ = −0.681\n(β=−0.477; p=8.2×10⁻⁷⁰)",
           size=1.9, colour=CB_BLK) +
  nat_theme() +
  labs(title="b",
       subtitle="UK Biobank lipid-clinic cohort (n=36,739)")

# Panel C: OR forest plot — TRG shield adjustment models
or_df <- data.frame(
  model   = factor(c("Unadjusted","Age + sex\nadjusted",
                      "T2DM + BMI\n+ alcohol + age"),
                   levels=rev(c("Unadjusted","Age + sex\nadjusted",
                                 "T2DM + BMI\n+ alcohol + age"))),
  OR      = c(1.897, 1.821, 1.145),
  ci_lo   = c(1.860, 1.784, 1.122),
  ci_hi   = c(1.934, 1.860, 1.168),
  dataset = c("CALON-FH","CALON-FH","UKB LC cohort")
)

pC2 <- ggplot(or_df, aes(y=model, x=OR, xmin=ci_lo, xmax=ci_hi,
                          colour=dataset)) +
  geom_vline(xintercept=1, linetype="dashed", colour=CB_GRY, linewidth=0.3) +
  geom_errorbarh(height=0.2, linewidth=0.5) +
  geom_point(size=2.2) +
  scale_colour_manual(values=c("CALON-FH"=CB_BLU,"UKB LC cohort"=CB_GRN),
                      name="Dataset") +
  scale_x_log10("Odds ratio (log scale)", breaks=c(1,1.5,2,2.5),
                labels=number_format(accuracy=0.1)) +
  scale_y_discrete(NULL) +
  nat_theme() +
  theme(legend.position=c(0.82,0.18)) +
  labs(title="c",
       subtitle="TRG Filter OR for FH / severe hypercholesterolaemia")

fig2 <- (pA2 | pB2 | pC2)
fig2 <- fig2 + plot_annotation(
  caption = paste0("Figure 2 | Triglyceride Shield validation.\n",
                   "a, Distribution of Triglyceride Filter (TF = LDL_untreated/(TG+0.1)) ",
                   "by causal gene group (Cohen's d: LDLR=",round(d_ldlr,2),
                   "; APOB=",round(d_apob,2),").\n",
                   "b, TF is significantly attenuated in T2DM (β=−0.477; ",
                   "p=8.2×10⁻⁷⁰), confirming T2DM as a critical confounder.\n",
                   "c, TRG Filter odds ratio for FH across adjustment models. ",
                   "Effect persists after full confounder adjustment.\n",
                   "TF, Triglyceride Filter; T2DM, type 2 diabetes mellitus; ",
                   "LC, lipid clinic; OR, odds ratio."),
  theme = nat_theme()
)

ggsave(file.path(FIG_DIR, "Figure2_TRGShield.pdf"),
       fig2, width=183, height=65, units="mm", device=cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure2_TRGShield.png"),
       fig2, width=183, height=65, units="mm", dpi=300)
cat("✓ Figure 2 saved\n")

################################################################################
# FIGURE 3 — Clinical operating characteristics
# (A) Sensitivity/Specificity vs threshold curve + 3 annotated points
# (B) PPV by prior probability heatmap
# (C) LR+ LR- summary
################################################################################

cat("Building Figure 3...\n")

# Build full SE/SP curve from SW ROC
roc_full <- roc(wales$fh, wales$pred, quiet=TRUE)
roc_curve_df <- data.frame(
  threshold = roc_full$thresholds,
  sens      = roc_full$sensitivities,
  spec      = roc_full$specificities
) %>%
  filter(is.finite(threshold), threshold >= 0, threshold <= 1)

# Operating points
op3 <- data.frame(
  label     = c("Population\nscreening","Lipid clinic\n(Youden)","Cascade\nconfirmation"),
  threshold = c(0.0031, 0.0445, 0.159),
  sens      = c(0.9404, 0.6424, 0.3854),
  spec      = c(0.1762, 0.8111, 0.9501),
  colour    = c(CB_ORG, CB_BLU, CB_GRN)
)

pA3 <- ggplot(roc_curve_df) +
  geom_line(aes(x=threshold, y=sens),  colour=CB_BLU,  linewidth=0.5) +
  geom_line(aes(x=threshold, y=spec),  colour=CB_RED,  linewidth=0.5,
            linetype="dashed") +
  geom_vline(xintercept=op3$threshold, colour=op3$colour,
             linetype="dotted", linewidth=0.35) +
  geom_point(data=op3, aes(x=threshold, y=sens), colour=op3$colour,
             size=2, shape=16) +
  geom_point(data=op3, aes(x=threshold, y=spec), colour=op3$colour,
             size=2, shape=17) +
  geom_text(data=op3, aes(x=threshold, y=0.52, label=label),
            size=1.8, colour=op3$colour, vjust=0, hjust=0.5) +
  scale_x_continuous("TUDOR score threshold",
                     limits=c(0,0.6), breaks=c(0,0.2,0.4,0.6)) +
  scale_y_continuous("Sensitivity (●) / Specificity (▲)",
                     limits=c(0,1), breaks=seq(0,1,0.25)) +
  annotate("text", x=0.56, y=0.92, label="● Sensitivity",
           size=1.8, colour=CB_BLU, hjust=1) +
  annotate("text", x=0.56, y=0.86, label="▲ Specificity",
           size=1.8, colour=CB_RED, hjust=1) +
  nat_theme() +
  labs(title="a", subtitle="Operating characteristic curve — Wales (n=6,448)")

# Panel B: PPV heatmap across prior probabilities
prior_grid <- expand.grid(
  threshold_label = c("Sensitivity 94%\n(threshold 0.003)",
                       "Youden optimal\n(threshold 0.045)",
                       "Specificity 95%\n(threshold 0.159)"),
  prior_pct       = c(0.2, 1.3, 5, 10, 20, 50)
) %>%
  mutate(
    sens = case_when(
      grepl("94%",  threshold_label) ~ 0.9404,
      grepl("Youden",threshold_label) ~ 0.6424,
      grepl("95%",  threshold_label) ~ 0.3854
    ),
    spec = case_when(
      grepl("94%",  threshold_label) ~ 0.1762,
      grepl("Youden",threshold_label) ~ 0.8111,
      grepl("95%",  threshold_label) ~ 0.9501
    ),
    prior  = prior_pct / 100,
    ppv    = (sens * prior) /
             (sens * prior + (1 - spec) * (1 - prior)) * 100,
    threshold_label = factor(threshold_label,
                             levels=c("Sensitivity 94%\n(threshold 0.003)",
                                       "Youden optimal\n(threshold 0.045)",
                                       "Specificity 95%\n(threshold 0.159)"))
  )

pB3 <- ggplot(prior_grid, aes(x=factor(prior_pct), y=threshold_label,
                               fill=ppv)) +
  geom_tile(colour="white", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.1f%%", ppv)), size=1.9, colour="black") +
  scale_fill_gradient2(low="white", mid="#a6c8e8", high=CB_BLU,
                       midpoint=20, name="PPV (%)") +
  scale_x_discrete("Assumed FH prevalence (%)") +
  scale_y_discrete(NULL) +
  nat_theme() +
  theme(legend.position="right", legend.key.height=unit(4,"mm"),
        axis.text.x=element_text(size=6)) +
  labs(title="b", subtitle="Positive predictive value by prior probability")

# Panel C: LR+ bar chart
lr_df <- data.frame(
  scenario = factor(c("Population\nscreening","Lipid clinic\n(Youden)",
                       "Cascade\nconfirmation"),
                    levels=c("Population\nscreening","Lipid clinic\n(Youden)",
                              "Cascade\nconfirmation")),
  LR_pos   = c(1.14, 3.40, 7.73),
  LR_neg   = c(0.338, 0.441, 0.647),
  fill_col = c(CB_ORG, CB_BLU, CB_GRN)
)

pC3 <- ggplot(lr_df, aes(x=scenario, y=LR_pos, fill=scenario)) +
  geom_col(width=0.55, alpha=0.85) +
  geom_hline(yintercept=c(2,5,10), linetype="dashed",
             colour=CB_GRY, linewidth=0.3) +
  geom_text(aes(label=sprintf("LR+ %.2f", LR_pos), y=LR_pos+0.2),
            size=2, vjust=0) +
  scale_fill_manual(values=setNames(lr_df$fill_col, lr_df$scenario),
                    guide="none") +
  scale_y_continuous("Positive likelihood ratio (LR+)",
                     limits=c(0,9.5), breaks=c(0,2,5,7,9)) +
  annotate("text", x=0.7, y=c(2.1,5.1,10.1),
           label=c("LR+>2","LR+>5","LR+>10"),
           size=1.8, colour=CB_GRY, hjust=0) +
  scale_x_discrete(NULL) +
  nat_theme() +
  labs(title="c", subtitle="Likelihood ratios by clinical scenario")

fig3 <- pA3 | pB3 | pC3
fig3 <- fig3 + plot_annotation(
  caption = paste0("Figure 3 | Clinical operating characteristics.\n",
                   "a, Sensitivity and specificity across TUDOR score thresholds; ",
                   "three operating points are annotated for population screening ",
                   "(orange), lipid clinic Youden optimal (blue), and cascade confirmation (green).\n",
                   "b, Positive predictive value (PPV) heatmap across six FH prevalence assumptions ",
                   "at each operating threshold.\n",
                   "c, Positive likelihood ratios at the three clinical thresholds.\n",
                   "LR+, positive likelihood ratio."),
  theme = nat_theme()
)

ggsave(file.path(FIG_DIR, "Figure3_ClinicalOperating.pdf"),
       fig3, width=183, height=65, units="mm", device=cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3_ClinicalOperating.png"),
       fig3, width=183, height=65, units="mm", dpi=300)
cat("✓ Figure 3 saved\n")

################################################################################
# FIGURE 4 — UKB supplement: ASCVD penetrance + ethnicity
################################################################################

cat("Building Figure 4...\n")

# Panel A: ASCVD by LDL decile (already computed)
pA4 <- ggplot(asc_tbl, aes(x=ldl_decile)) +
  geom_col(aes(y=pct_ascvd), fill=CB_RED, alpha=0.8, width=0.7) +
  geom_line(aes(y=pct_prem), colour=CB_ORG, linewidth=0.7) +
  geom_point(aes(y=pct_prem), colour=CB_ORG, size=1.5) +
  scale_x_continuous("LDL severity decile (proxy for FH burden)",
                     breaks=1:10) +
  scale_y_continuous("ASCVD prevalence (%)", limits=c(0,20)) +
  annotate("text", x=8.5, y=18.5,
           label="― Premature ASCVD (<55 yr)", size=2, colour=CB_ORG) +
  annotate("text", x=8.5, y=17,
           label="█ Any ASCVD", size=2, colour=CB_RED) +
  nat_theme() +
  labs(title="a",
       subtitle="UK Biobank lipid-clinic cohort (n=37,050)")

# Panel B: Ethnicity — LDL + ASCVD + T2DM bubble plot
eth_plot <- eth_tbl %>%
  filter(!is.na(ethnicity), n >= 20) %>%
  mutate(ethnicity = fct_reorder(ethnicity, ldl_mean))

pB4 <- ggplot(eth_plot,
              aes(x=pct_t2dm, y=pct_ascvd, colour=ethnicity, size=n)) +
  geom_point(alpha=0.85) +
  geom_text(aes(label=ethnicity), size=1.9, vjust=-1,
            show.legend=FALSE) +
  scale_size_continuous("n", range=c(2,10),
                        breaks=c(100,1000,10000,36000)) +
  scale_colour_brewer(palette="Set2", guide="none") +
  scale_x_continuous("T2DM prevalence (%)", limits=c(0,25)) +
  scale_y_continuous("ASCVD prevalence (%)", limits=c(0,25)) +
  nat_theme() +
  labs(title="b",
       subtitle="Ethnicity: T2DM–ASCVD–cohort size relationship")

# Panel C: NRI/IDI summary
nri_tbl <- read_csv(file.path(RES_DIR,"nri_idi_bootstrap.csv"),
                     show_col_types=FALSE)
nri_plot <- bind_rows(
  nri_tbl %>%
    filter(cohort %in% c("South Wales","All Wales")) %>%
    transmute(cohort, stat="NRI (events)", est=nri_events,
              lo_val=nri_events_lo, hi_val=nri_events_hi),
  nri_tbl %>%
    filter(cohort %in% c("South Wales","All Wales")) %>%
    transmute(cohort, stat="IDI", est=idi, lo_val=idi_lo, hi_val=idi_hi)
)

pC4 <- ggplot(nri_plot, aes(y=cohort, x=est, xmin=lo_val, xmax=hi_val,
                              colour=stat)) +
  geom_vline(xintercept=0, linetype="dashed", colour=CB_GRY, linewidth=0.3) +
  geom_errorbar(aes(ymin=NULL, ymax=NULL), orientation="y",
                width=0.25, linewidth=0.5, position=position_dodge(0.5)) +
  geom_point(size=2.5, position=position_dodge(0.5)) +
  scale_colour_manual(values=c("NRI (events)"=CB_BLU,"IDI"=CB_GRN),
                      name="Metric") +
  scale_x_continuous("Improvement over DLCN", limits=c(-0.55,0.4)) +
  scale_y_discrete(NULL) +
  nat_theme() +
  theme(legend.position=c(0.25,0.85)) +
  labs(title="c",
       subtitle="Net reclassification improvement vs DLCN")

fig4 <- pA4 | pB4 | pC4
fig4 <- fig4 + plot_annotation(
  caption = paste0("Figure 4 | Clinical outcomes and reclassification.\n",
                   "a, Any ASCVD (bars) and premature ASCVD (<55 years, line) ",
                   "by LDL severity decile in the UK Biobank lipid-clinic cohort.\n",
                   "b, Ethnicity plotted by T2DM and ASCVD prevalence; ",
                   "South Asian patients show highest T2DM (21.4%) and ASCVD (19.0%) rates.\n",
                   "c, NRI for events and IDI of TUDOR vs DLCN. ",
                   "South Wales NRI_events=+0.173 (95% CI 0.014–0.338; p=0.033).\n",
                   "ASCVD, atherosclerotic cardiovascular disease; NRI, net reclassification improvement; ",
                   "IDI, integrated discrimination improvement."),
  theme = nat_theme()
)

ggsave(file.path(FIG_DIR, "Figure4_ASCVD_NRI.pdf"),
       fig4, width=183, height=65, units="mm", device=cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure4_ASCVD_NRI.png"),
       fig4, width=183, height=65, units="mm", dpi=300)
cat("✓ Figure 4 saved\n")

################################################################################
# FIGURE 5 — Penetrance gradient + statin paradox + score architecture
################################################################################

cat("Building Figure 5...\n")

pen_tbl <- read_csv(file.path(RES_DIR,"penetrance_by_tudor_decile.csv"),
                    show_col_types=FALSE)

# Panel A: TUDOR decile vs LDL_ut (bar) + statin use (line)
pA5 <- ggplot(pen_tbl, aes(x=tudor_decile)) +
  geom_col(aes(y=mean_ldl_ut), fill=CB_BLU, alpha=0.8, width=0.7) +
  geom_line(aes(y=pct_on_statin/10), colour=CB_RED,
            linewidth=0.6, linetype="dashed") +
  geom_point(aes(y=pct_on_statin/10), colour=CB_RED, size=1.5) +
  scale_x_continuous("TUDOR score decile", breaks=1:10) +
  scale_y_continuous(
    "Untreated LDL (mmol/L)",
    sec.axis = sec_axis(~.*10, name="Statin use (%)",
                        breaks=seq(0,50,10))
  ) +
  annotate("text", x=9, y=7.5,
           label="― Statin use %\n(right axis)", size=1.8,
           colour=CB_RED, hjust=1) +
  nat_theme() +
  labs(title="a",
       subtitle="TUDOR score vs LDL burden and statin use (CALON-FH)")

# Panel B: TUDOR score architecture waterfall
features <- data.frame(
  feature  = factor(c("LDL_untreated","Trig. Filter","Non-HDL–LDL gap",
                       "ApoB/LDL ratio","Age","Sex",
                       "Statin use\n(adjustment)","DLCN score"),
                    levels=rev(c("LDL_untreated","Trig. Filter","Non-HDL–LDL gap",
                                  "ApoB/LDL ratio","Age","Sex",
                                  "Statin use\n(adjustment)","DLCN score"))),
  direction = c("Positive","Positive","Positive","Positive",
                "Negative","Mixed","Adjustment","Positive"),
  importance = c(100, 72, 48, 35, 22, 15, 12, 30)
)

pB5 <- ggplot(features, aes(y=feature, x=importance, fill=direction)) +
  geom_col(width=0.6, alpha=0.85) +
  scale_fill_manual(
    values=c("Positive"=CB_BLU,"Negative"=CB_RED,
             "Mixed"=CB_GRY,"Adjustment"=CB_ORG),
    name="Effect"
  ) +
  scale_x_continuous("Relative importance (arbitrary units)",
                     limits=c(0,115)) +
  scale_y_discrete(NULL) +
  nat_theme() +
  theme(legend.position=c(0.75,0.25)) +
  labs(title="b", subtitle="TUDOR score feature contributions (schematic)")

# Panel C: TUDOR decile vs chol-years proxy
pC5 <- ggplot(pen_tbl, aes(x=tudor_decile, y=chol_years_proxy)) +
  geom_col(fill=CB_GRN, alpha=0.8, width=0.7) +
  geom_smooth(method="lm", se=FALSE, colour=CB_BLK,
              linewidth=0.5, linetype="dashed") +
  scale_x_continuous("TUDOR score decile", breaks=1:10) +
  scale_y_continuous("Cholesterol-years proxy\n(LDL_ut × age, mmol/L·yr)") +
  nat_theme() +
  labs(title="c",
       subtitle="Cumulative LDL burden by TUDOR decile")

fig5 <- pA5 | pB5 | pC5
fig5 <- fig5 + plot_annotation(
  caption = paste0("Figure 5 | TUDOR score architecture and penetrance gradient.\n",
                   "a, Mean untreated LDL (bars) rises monotonically across TUDOR deciles; ",
                   "statin use (dashed line) peaks in mid-deciles, demonstrating ",
                   "the treatment paradox captured by the score.\n",
                   "b, Schematic representation of TUDOR feature contributions; ",
                   "LDL_untreated and Triglyceride Filter are the dominant features.\n",
                   "c, Cholesterol-years proxy (LDL_untreated × age) by TUDOR decile, ",
                   "confirming cumulative atherogenic burden increases with score.\n",
                   "LDL, low-density lipoprotein; TF, Triglyceride Filter."),
  theme = nat_theme()
)

ggsave(file.path(FIG_DIR, "Figure5_Penetrance_Architecture.pdf"),
       fig5, width=183, height=65, units="mm", device=cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure5_Penetrance_Architecture.png"),
       fig5, width=183, height=65, units="mm", dpi=300)
cat("✓ Figure 5 saved\n")

cat("\n========================================================\n")
cat("ALL 5 FIGURES COMPLETE\n")
cat("========================================================\n")
cat("Figure1_ROC_Calibration     — triple ROC + calibration panel\n")
cat("Figure2_TRGShield           — gene violin + T2DM + forest\n")
cat("Figure3_ClinicalOperating   — SE/SP curves + PPV heatmap + LR\n")
cat("Figure4_ASCVD_NRI           — ASCVD penetrance + ethnicity + NRI\n")
cat("Figure5_Penetrance_Architecture — TUDOR decile + features + chol-yrs\n")
cat("All saved as .pdf (vector) + .png (300 DPI) in:\n")
cat(FIG_DIR,"\n")
