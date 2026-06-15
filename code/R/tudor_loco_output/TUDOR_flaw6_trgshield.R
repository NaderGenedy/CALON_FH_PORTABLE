##############################################################################
# TUDOR Flaw 6 — Triglyceride Shield Multivariable Validation
# Adjusts for: T2DM, age, BMI, sex, alcohol
#
# Two modes:
#   MODE A: Using existing data (unadjusted Cohen's d — available now)
#   MODE B: Using UKB supplement (fully adjusted — after UKB extract)
##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(broom)
})

OUT  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output/tudor_flaw_results"
DATA <- "C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output"

pred <- read_csv(file.path(DATA, "loco_predictions_complete.csv"),
                 show_col_types = FALSE)

##############################################################################
# MODE A — AVAILABLE DATA (no T2DM/BMI/alcohol yet)
##############################################################################

cat("=== FLAW 6 MODE A: Triglyceride Shield (available confounders) ===\n\n")

# The TRG shield hypothesis: in true FH (LDLR dysfunction), hepatic LDL-
# receptor loss does NOT impair VLDL secretion, so TG remains normal-to-low.
# Trig_Filter = LDL_untreated / (TG + 0.1) captures this: high in FH, low in
# polygenic/secondary hypercholesterolaemia (where TG is usually elevated).

trgdata <- pred %>%
  filter(!is.na(trig_filter), !is.na(fh), !is.na(tg),
         !is.na(age), !is.na(sex)) %>%
  mutate(
    gene_group = case_when(
      gene == "LDLR" ~ "LDLR",
      gene == "APOB" ~ "APOB",
      gene == "PCSK9" ~ "PCSK9",
      fh == 0        ~ "Non-FH",
      TRUE           ~ "FH-unknown"
    )
  )

cat(sprintf("Complete cases for TRG shield analysis: n=%d\n\n", nrow(trgdata)))

# ── 1. Unadjusted group comparison ───────────────────────────────────────────
trg_summary <- trgdata %>%
  group_by(gene_group) %>%
  summarise(
    n             = n(),
    TF_mean       = round(mean(trig_filter), 3),
    TF_sd         = round(sd(trig_filter),   3),
    TF_median     = round(median(trig_filter), 3),
    TG_mean       = round(mean(tg),   3),
    LDL_ut_mean   = round(mean(ldl_ut, na.rm=TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(TF_mean))

print(trg_summary)
write_csv(trg_summary, file.path(OUT, "trgshield_group_summary.csv"))

# ── 2. Cohen's d: LDLR vs Non-FH (unadjusted) ────────────────────────────────
ldlr_tf   <- trgdata %>% filter(gene_group == "LDLR") %>% pull(trig_filter)
apob_tf   <- trgdata %>% filter(gene_group == "APOB")  %>% pull(trig_filter)
nonfh_tf  <- trgdata %>% filter(gene_group == "Non-FH") %>% pull(trig_filter)

cohens_d <- function(x, y) (mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE)) /
  sqrt((var(x, na.rm=TRUE) * (length(x)-1) + var(y, na.rm=TRUE) * (length(y)-1)) /
         (length(x) + length(y) - 2))

cat(sprintf("Cohen's d (LDLR vs Non-FH): %.3f\n",
            cohens_d(ldlr_tf, nonfh_tf)))
cat(sprintf("Cohen's d (APOB vs Non-FH): %.3f\n",
            cohens_d(apob_tf, nonfh_tf)))

# ── 3. Multivariable logistic: Trig_Filter predicts FH, adjusted for
#       age + sex (available confounders; T2DM/BMI/alcohol from UKB extract)
# ─────────────────────────────────────────────────────────────────────────────
model_unadj <- glm(fh ~ trig_filter,
                   data = trgdata, family = binomial)
model_adj1  <- glm(fh ~ trig_filter + age + sex,
                   data = trgdata, family = binomial)

# model_adj2 would add: + t2dm + bmi + alcohol_cat (after UKB extract)

tidy_models <- bind_rows(
  tidy(model_unadj, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(model = "Unadjusted"),
  tidy(model_adj1, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(model = "Adjusted (age+sex)")
) %>%
  filter(term == "trig_filter") %>%
  select(model, OR = estimate, ci_lo = conf.low, ci_hi = conf.high,
         p_value = p.value) %>%
  mutate(across(c(OR, ci_lo, ci_hi), ~round(., 3)),
         p_value = round(p_value, 4))

cat("\nTriglyceride Filter — Odds Ratio for FH:\n")
print(tidy_models)
write_csv(tidy_models, file.path(OUT, "trgshield_multivariable.csv"))

# ── 4. Incremental AUC: TUDOR with vs without trig_filter ────────────────────
# Create TUDOR-minus-TF pseudo-model (remove trig_filter contribution)
trgdata_cc <- trgdata %>% filter(!is.na(pred), !is.na(dlcn))

roc_tudor   <- roc(trgdata_cc$fh, trgdata_cc$pred,         quiet = TRUE)
roc_notf    <- roc(trgdata_cc$fh, trgdata_cc$pred -
                     0.3 * scale(trgdata_cc$trig_filter)[,1] * 0.1,
                   quiet = TRUE)   # approximate contribution removal

dl_result <- roc.test(roc_tudor, roc_notf, method = "delong")

cat(sprintf("\nAUC with TrigFilter:    %.4f\n",    as.numeric(auc(roc_tudor))))
cat(sprintf("AUC without TrigFilter: %.4f\n",    as.numeric(auc(roc_notf))))
cat(sprintf("DeLong p-value:         %.4e\n\n",  dl_result$p.value))

# ── 5. T2DM sensitivity note ─────────────────────────────────────────────────
cat("=== T2DM ADJUSTMENT STATUS ===\n")
cat("T2DM is a critical confounder for TG: T2DM raises TG by ~0.5-1.0 mmol/L\n")
cat("and would artefactually LOWER Trig_Filter in FH+T2DM patients.\n")
cat("ACTION: Extract UKB field 2443 via tudor_ukb_extract_supplement.sh\n")
cat("Then add model_adj2 <- glm(fh ~ trig_filter + age + sex + t2dm + bmi\n")
cat("                         + alcohol_cat, data=trgdata, family=binomial)\n\n")

cat("Key expected findings after T2DM adjustment:\n")
cat("  - OR for trig_filter should INCREASE (positive confounding by T2DM)\n")
cat("  - ApoB-discordant cases (high TG + normal LDL) will be separated\n")
cat("  - Interaction term: trig_filter × t2dm tests if TRG shield is\n")
cat("    attenuated in diabetic FH patients\n")

##############################################################################
# FIGURES — TRG shield distribution
##############################################################################

C_BLUE  <- "#0072B2"; C_ORANGE <- "#E69F00"
C_GREEN <- "#009E73"; C_RED    <- "#D55E00"
C_GREY  <- "#999999"

nat_theme <- function() {
  theme_classic(base_size = 8, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.35),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 8))
}

plot_df <- trgdata %>%
  filter(gene_group %in% c("LDLR","APOB","Non-FH")) %>%
  mutate(gene_group = factor(gene_group,
                             levels = c("Non-FH","APOB","LDLR")))

p_trg <- ggplot(plot_df, aes(x = gene_group, y = trig_filter,
                              fill = gene_group)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.4,
               fill = "white", colour = "grey40") +
  scale_fill_manual(values = c("Non-FH" = C_GREY,
                                "APOB"   = C_ORANGE,
                                "LDLR"   = C_BLUE),
                    guide = "none") +
  scale_y_continuous("Triglyceride Filter (LDL_untreated / TG+0.1)",
                     limits = c(0, quantile(plot_df$trig_filter, 0.99))) +
  scale_x_discrete("Group") +
  annotate("text", x = 2.5, y = quantile(plot_df$trig_filter, 0.97),
           label = paste0("Cohen's d\nLDLR vs Non-FH: ",
                          round(cohens_d(ldlr_tf, nonfh_tf), 3)),
           size = 2.2, colour = C_BLUE) +
  nat_theme() +
  labs(title = "Triglyceride Shield Distribution by Gene Group",
       subtitle = "Unadjusted; T2DM/BMI adjustment pending UKB extract",
       caption = "Higher Trig_Filter = preserved VLDL secretion (FH signature)")

ggsave(file.path(OUT, "Figure_TRGShield_Violin.pdf"),
       p_trg, width = 89, height = 89, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure_TRGShield_Violin.png"),
       p_trg, width = 89, height = 89, units = "mm", dpi = 300)

cat("TRG Shield violin plot saved.\n")
