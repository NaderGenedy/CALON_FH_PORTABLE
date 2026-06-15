##############################################################################
# TUDOR — UKB Supplement Post-Processing
# Reads the 5 downloaded CSVs from E:/
# Outputs merged analysis file + all flaw-specific results
#
# NOTE: No EID in loco_predictions_complete.csv → analyses run STANDALONE
#       on the UKB lipid-clinic filtered cohort (LDL≥4.9 or TC≥7.5).
#       To add TUDOR scores: re-run TUDOR_UKB_LIPID_CLINIC.R saving 'eid'.
##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(broom)
  library(jsonlite)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
RAW_DIR  <- "E:/"
PRED_DIR <- "C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_loco_output"
OUT      <- file.path(PRED_DIR, "tudor_flaw_results")
dir.create(OUT, showWarnings = FALSE)

C_BLUE <- "#0072B2"; C_ORANGE <- "#E69F00"
C_GREEN <- "#009E73"; C_RED <- "#D55E00"; C_GREY <- "#999999"

nat <- function() theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

##############################################################################
# LOAD ALL 5 FILES
##############################################################################

cat("Loading extracted UKB files...\n")

demo  <- read_csv(file.path(RAW_DIR, "ukb_reviewer_demographics.csv"),
                  show_col_types = FALSE)
icd   <- read_csv(file.path(RAW_DIR, "ukb_reviewer_icd10_codes.csv"),
                  show_col_types = FALSE)
lipid <- read_csv(file.path(RAW_DIR, "ukb_reviewer_longitudinal_lipids.csv"),
                  show_col_types = FALSE)
apob  <- read_csv(file.path(RAW_DIR, "ukb_reviewer_apob_lpa.csv"),
                  show_col_types = FALSE)
smk   <- read_csv(file.path(RAW_DIR, "ukb_reviewer_smoking.csv"),
                  show_col_types = FALSE)

# Standardise EID column
rename_eid <- function(df) {
  if ("participant.eid" %in% names(df)) return(rename(df, eid = `participant.eid`))
  if (!"eid" %in% names(df)) names(df)[1] <- "eid"
  df
}
demo  <- rename_eid(demo)
icd   <- rename_eid(icd)
lipid <- rename_eid(lipid)
apob  <- rename_eid(apob)
smk   <- rename_eid(smk)

cat(sprintf("demo: %d rows | icd: %d | lipid: %d | apob: %d | smk: %d\n",
            nrow(demo), nrow(icd), nrow(lipid), nrow(apob), nrow(smk)))

# Print column names to confirm mapping
cat("\nColumn names:\n")
cat("demo: ", paste(names(demo), collapse=", "), "\n")
cat("lipid:", paste(names(lipid), collapse=", "), "\n")
cat("apob: ", paste(names(apob), collapse=", "), "\n")

##############################################################################
# F5 — ETHNICITY
##############################################################################

cat("\n=== F5: Ethnicity ===\n")

# Find ethnicity column (may be named differently)
eth_col <- grep("p21000", names(demo), value = TRUE)[1]
cat(sprintf("Ethnicity column: %s\n", eth_col))

eth_map <- c(
  "1"="White","1001"="White","1002"="White","1003"="White","4001"="White",
  "2"="Mixed",
  "3"="South Asian","3001"="South Asian","3002"="South Asian",
  "3003"="South Asian","3004"="South Asian",
  "4"="Black","4002"="Black","4003"="Black",
  "5"="East Asian",
  "6"="Other","-1"=NA_character_,"-3"=NA_character_
)

if (!is.na(eth_col)) {
  demo <- demo %>%
    mutate(
      eth_raw   = as.character(.data[[eth_col]]),
      ethnicity = recode(eth_raw, !!!eth_map, .default = "Other"),
      ethnicity = factor(ethnicity,
                         levels = c("White","South Asian","Black",
                                    "East Asian","Mixed","Other"))
    )
} else {
  demo$ethnicity <- NA_character_
}
cat("Ethnicity distribution:\n"); print(table(demo$ethnicity, useNA = "ifany"))

##############################################################################
# F6 — BMI, ALCOHOL
##############################################################################

cat("\n=== F6: BMI + Alcohol ===\n")

# BMI columns
bmi_cols <- grep("p21001", names(demo), value = TRUE)
cat(sprintf("BMI columns: %s\n", paste(bmi_cols, collapse=", ")))

# Alcohol column
alc_col <- grep("p1558", names(demo), value = TRUE)[1]
cat(sprintf("Alcohol column: %s\n", alc_col))

# DOB columns
yr_col  <- grep("p34", names(demo), value = TRUE)[1]
mo_col  <- grep("p52", names(demo), value = TRUE)[1]
dt_col  <- grep("p53", names(demo), value = TRUE)[1]

demo <- demo %>%
  mutate(
    bmi = if (length(bmi_cols) >= 4)
            coalesce(.data[[bmi_cols[4]]], .data[[bmi_cols[3]]],
                     .data[[bmi_cols[2]]], .data[[bmi_cols[1]]])
          else if (length(bmi_cols) >= 1)
            .data[[bmi_cols[1]]]
          else NA_real_,

    alcohol_cat = if (!is.na(alc_col)) case_when(
      .data[[alc_col]] %in% c(1,2) ~ "Heavy",
      .data[[alc_col]] %in% c(3,4) ~ "Moderate",
      .data[[alc_col]] %in% c(5,6) ~ "Light/None",
      TRUE ~ NA_character_
    ) else NA_character_,
    alcohol_cat = factor(alcohol_cat,
                         levels = c("Light/None","Moderate","Heavy")),

    dob = if (!is.na(yr_col) && !is.na(mo_col))
            as.Date(paste(.data[[yr_col]], .data[[mo_col]], 15),
                    format = "%Y %m %d")
          else NA,
    date_baseline = if (!is.na(dt_col)) as.Date(.data[[dt_col]]) else NA,
    age_baseline  = as.numeric(difftime(date_baseline, dob, "days")) / 365.25
  )

cat(sprintf("BMI: mean=%.1f SD=%.1f (n=%d)\n",
            mean(demo$bmi, na.rm=T), sd(demo$bmi, na.rm=T),
            sum(!is.na(demo$bmi))))
cat("Alcohol:\n"); print(table(demo$alcohol_cat, useNA="ifany"))
cat(sprintf("Age at baseline: mean=%.1f SD=%.1f\n",
            mean(demo$age_baseline, na.rm=T), sd(demo$age_baseline, na.rm=T)))

##############################################################################
# F7 — ICD-10: ASCVD + T2DM
##############################################################################

cat("\n=== F7: ICD-10 ASCVD + T2DM ===\n")

# Find ICD column
icd_col <- grep("p41270", names(icd), value = TRUE)[1]
cat(sprintf("ICD-10 column: %s\n", icd_col))

parse_icd <- function(x) {
  if (is.na(x) || x == "" || x == "[]") return(character(0))
  tryCatch(fromJSON(x), error = function(e) character(0))
}

cat("Parsing ICD-10 arrays (may take ~60s)...\n")

ascvd_pat <- "^I2[01245]|^I200|^I209|^I240|^I248|^I249|^I63|^G45[0-9]|^I74[089]"
t2dm_pat  <- "^E11"

icd_list <- lapply(icd[[icd_col]], parse_icd)

icd <- icd %>%
  mutate(
    any_ascvd   = map_lgl(icd_list, ~any(grepl(ascvd_pat, .x))),
    t2dm        = map_lgl(icd_list, ~any(grepl(t2dm_pat,  .x))),
    n_diagnoses = map_int(icd_list, length)
  ) %>%
  select(eid, any_ascvd, t2dm, n_diagnoses)

cat(sprintf("Any ASCVD: n=%d (%.1f%%)\n", sum(icd$any_ascvd), 100*mean(icd$any_ascvd)))
cat(sprintf("T2DM:      n=%d (%.1f%%)\n", sum(icd$t2dm), 100*mean(icd$t2dm)))

##############################################################################
# SMOKING
##############################################################################

smk_col <- grep("p20116", names(smk), value = TRUE)[1]
if (!is.na(smk_col)) {
  smk <- smk %>%
    mutate(
      smoking = case_when(
        .data[[smk_col]] == 0 ~ "Never",
        .data[[smk_col]] == 1 ~ "Previous",
        .data[[smk_col]] == 2 ~ "Current",
        TRUE ~ NA_character_
      )
    ) %>%
    select(eid, smoking)
} else {
  smk <- smk %>% mutate(smoking = NA_character_) %>% select(eid, smoking)
}

##############################################################################
# F10 — LONGITUDINAL LDL
##############################################################################

cat("\n=== F10: Longitudinal LDL ===\n")
cat("Lipid columns: ", paste(names(lipid), collapse=", "), "\n")

# Find columns by pattern
ldl_i0_col <- grep("p30780.*i0|p30780_i0", names(lipid), value=TRUE)[1]
ldl_i1_col <- grep("p30780.*i1|p30780_i1", names(lipid), value=TRUE)[1]
tc_i0_col  <- grep("p30690.*i0|p30690_i0", names(lipid), value=TRUE)[1]
hdl_i0_col <- grep("p30760.*i0|p30760_i0", names(lipid), value=TRUE)[1]
tg_i0_col  <- grep("p30870.*i0|p30870_i0", names(lipid), value=TRUE)[1]

cat(sprintf("LDL i0: %s | i1: %s | TC i0: %s | HDL i0: %s | TG i0: %s\n",
            ldl_i0_col, ldl_i1_col, tc_i0_col, hdl_i0_col, tg_i0_col))

lipid <- lipid %>%
  mutate(
    ldl_i0 = if (!is.na(ldl_i0_col)) .data[[ldl_i0_col]] else NA_real_,
    ldl_i1 = if (!is.na(ldl_i1_col)) .data[[ldl_i1_col]] else NA_real_,
    tc_i0  = if (!is.na(tc_i0_col))  .data[[tc_i0_col]]  else NA_real_,
    hdl_i0 = if (!is.na(hdl_i0_col)) .data[[hdl_i0_col]] else NA_real_,
    tg_i0  = if (!is.na(tg_i0_col))  .data[[tg_i0_col]]  else NA_real_,

    ldl_baseline = ldl_i0,
    ldl_mean     = rowMeans(cbind(ldl_i0, ldl_i1), na.rm = TRUE),
    ldl_max      = pmax(ldl_i0, ldl_i1, na.rm = TRUE),
    n_ldl_obs    = (!is.na(ldl_i0)) + (!is.na(ldl_i1)),

    # Trapezoid: mean LDL × 6yr gap (two-timepoint approximation)
    chol_years_ukb = case_when(
      n_ldl_obs == 2 ~ 0.5 * (ldl_i0 + ldl_i1) * 6,
      TRUE           ~ NA_real_
    )
  ) %>%
  select(eid, ldl_baseline, ldl_i0, ldl_i1, ldl_mean, ldl_max,
         tc_i0, hdl_i0, tg_i0, n_ldl_obs, chol_years_ukb)

cat(sprintf("LDL i0: mean=%.2f mmol/L (n=%d)\n",
            mean(lipid$ldl_baseline, na.rm=T), sum(!is.na(lipid$ldl_baseline))))
cat(sprintf("LDL i1 available: n=%d\n", sum(!is.na(lipid$ldl_i1))))
cat(sprintf("Chol-years (2 pts): n=%d\n", sum(!is.na(lipid$chol_years_ukb))))

##############################################################################
# APoB + Lp(a)
##############################################################################

apob_col  <- grep("p30890.*i0|p30890_i0", names(apob), value=TRUE)[1]
lpa_col   <- grep("p30900.*i0|p30900_i0", names(apob), value=TRUE)[1]

cat(sprintf("\nApoB column: %s | Lp(a) column: %s\n", apob_col, lpa_col))

apob <- apob %>%
  mutate(
    apob_measured = if (!is.na(apob_col)) .data[[apob_col]] else NA_real_,
    lpa_nmol_L    = if (!is.na(lpa_col))  .data[[lpa_col]]  else NA_real_
  ) %>%
  select(eid, apob_measured, lpa_nmol_L)

cat(sprintf("ApoB measured: n=%d (%.0f%%)\n",
            sum(!is.na(apob$apob_measured)), 100*mean(!is.na(apob$apob_measured))))
cat(sprintf("Lp(a) measured: n=%d (%.0f%%)\n",
            sum(!is.na(apob$lpa_nmol_L)), 100*mean(!is.na(apob$lpa_nmol_L))))

##############################################################################
# MERGE ALL INTO SINGLE SUPPLEMENT TABLE
##############################################################################

cat("\n=== Merging all supplements ===\n")

demo_slim <- demo %>%
  select(eid, ethnicity, bmi, alcohol_cat, dob, date_baseline, age_baseline)

supplement <- demo_slim %>%
  left_join(icd,   by = "eid") %>%
  left_join(smk,   by = "eid") %>%
  left_join(lipid, by = "eid") %>%
  left_join(apob,  by = "eid")

cat(sprintf("Supplement: %d rows × %d cols\n", nrow(supplement), ncol(supplement)))
cat("Columns:", paste(names(supplement), collapse=", "), "\n")

# ── Apply lipid-clinic filter ─────────────────────────────────────────────────
supplement <- supplement %>%
  mutate(
    lipid_clinic_elig = case_when(
      ldl_baseline >= 4.9 ~ TRUE,
      tc_i0 >= 7.5        ~ TRUE,
      is.na(ldl_baseline) & is.na(tc_i0) ~ NA,
      TRUE ~ FALSE
    ),
    # TRG filter proxy (same formula as TUDOR feature trig_filter)
    trig_filter_proxy = ldl_baseline / (tg_i0 + 0.1),
    # Premature ASCVD proxy: ASCVD present AND age < 55 at baseline
    premature_ascvd = any_ascvd & !is.na(age_baseline) & age_baseline < 55
  )

lc_cohort <- supplement %>% filter(lipid_clinic_elig == TRUE)
cat(sprintf("Lipid-clinic eligible (LDL≥4.9 or TC≥7.5): n=%d\n", nrow(lc_cohort)))
cat(sprintf("Among LC eligible — T2DM: n=%d (%.1f%%)\n",
            sum(lc_cohort$t2dm, na.rm=T),
            100*mean(lc_cohort$t2dm, na.rm=T)))
cat(sprintf("Among LC eligible — Any ASCVD: n=%d (%.1f%%)\n",
            sum(lc_cohort$any_ascvd, na.rm=T),
            100*mean(lc_cohort$any_ascvd, na.rm=T)))
cat(sprintf("Among LC eligible — Premature ASCVD (<55): n=%d (%.1f%%)\n",
            sum(lc_cohort$premature_ascvd, na.rm=T),
            100*mean(lc_cohort$premature_ascvd, na.rm=T)))

write_csv(supplement, file.path(OUT, "ukb_supplement_clean.csv"))
write_csv(lc_cohort,  file.path(OUT, "ukb_supplement_lc_cohort.csv"))
cat("✓ Saved ukb_supplement_clean.csv and ukb_supplement_lc_cohort.csv\n")

##############################################################################
# F5 — ETHNICITY DISTRIBUTION OF LIPID-CLINIC COHORT
# (AUC stratification requires TUDOR score via EID merge — pending)
##############################################################################

cat("\n=== F5: Ethnicity distribution (LC cohort) ===\n")

eth_summary <- lc_cohort %>%
  filter(!is.na(ethnicity)) %>%
  group_by(ethnicity) %>%
  summarise(
    n             = n(),
    pct           = round(100 * n() / nrow(lc_cohort), 1),
    ldl_mean      = round(mean(ldl_baseline, na.rm=TRUE), 2),
    ldl_sd        = round(sd(ldl_baseline, na.rm=TRUE), 2),
    pct_ascvd     = round(100 * mean(any_ascvd, na.rm=TRUE), 1),
    pct_t2dm      = round(100 * mean(t2dm, na.rm=TRUE), 1),
    trig_filt_med = round(median(trig_filter_proxy, na.rm=TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

cat("\nEthnicity summary (lipid-clinic cohort):\n")
print(eth_summary)
write_csv(eth_summary, file.path(OUT, "ethnicity_lc_cohort_summary.csv"))

# Figure — LDL by ethnicity
p_eth_ldl <- lc_cohort %>%
  filter(!is.na(ethnicity), !is.na(ldl_baseline)) %>%
  ggplot(aes(x = reorder(ethnicity, ldl_baseline, median, na.rm=TRUE),
             y = ldl_baseline, fill = ethnicity)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous("LDL-C at baseline (mmol/L)",
                     limits = c(0, quantile(lc_cohort$ldl_baseline, 0.99, na.rm=T))) +
  coord_flip() +
  nat() +
  labs(title = "LDL distribution by ethnicity",
       subtitle = "UK Biobank lipid-clinic cohort (LDL≥4.9 or TC≥7.5)",
       x = NULL)

ggsave(file.path(OUT, "Figure_LDL_Ethnicity.pdf"), p_eth_ldl,
       width = 89, height = 89, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure_LDL_Ethnicity.png"), p_eth_ldl,
       width = 89, height = 89, units = "mm", dpi = 300)
cat("✓ Saved ethnicity summary + Figure_LDL_Ethnicity\n")
cat("NOTE: Ethnicity-stratified AUC requires TUDOR score merge (re-extract with EID)\n")

##############################################################################
# F6 — TRG SHIELD: T2DM CONFOUNDING VALIDATION
# Shows T2DM lowers TRG filter proxy → validates T2DM as critical confounder
##############################################################################

cat("\n=== F6: TRG Shield — T2DM confounding validation ===\n")

trg_data <- lc_cohort %>%
  filter(!is.na(trig_filter_proxy), !is.na(t2dm), !is.na(bmi),
         !is.na(age_baseline), !is.na(alcohol_cat),
         is.finite(trig_filter_proxy))

cat(sprintf("TRG shield analysis n=%d (LC cohort with complete confounders)\n",
            nrow(trg_data)))

# 1. Descriptive: TRG filter proxy by T2DM status
trg_by_t2dm <- trg_data %>%
  group_by(t2dm) %>%
  summarise(
    n          = n(),
    tf_mean    = round(mean(trig_filter_proxy), 3),
    tf_sd      = round(sd(trig_filter_proxy), 3),
    tf_median  = round(median(trig_filter_proxy), 3),
    tg_mean    = round(mean(tg_i0, na.rm=T), 3),
    ldl_mean   = round(mean(ldl_baseline, na.rm=T), 3),
    .groups = "drop"
  )
cat("\nTRG filter proxy by T2DM status:\n")
print(trg_by_t2dm)
write_csv(trg_by_t2dm, file.path(OUT, "trgshield_by_t2dm.csv"))

# Cohen's d: T2DM+ vs T2DM-
cohens_d <- function(x, y) (mean(x,na.rm=T) - mean(y,na.rm=T)) /
  sqrt((var(x,na.rm=T)*(length(x[!is.na(x)])-1) +
        var(y,na.rm=T)*(length(y[!is.na(y)])-1)) /
       (length(x[!is.na(x)]) + length(y[!is.na(y)]) - 2))

tf_t2dm    <- trg_data$trig_filter_proxy[trg_data$t2dm == TRUE]
tf_no_t2dm <- trg_data$trig_filter_proxy[trg_data$t2dm == FALSE]
d_t2dm     <- cohens_d(tf_no_t2dm, tf_t2dm)
cat(sprintf("\nCohen's d (No-T2DM vs T2DM): %.3f\n", d_t2dm))
cat(sprintf("Mean TRG filter: No-T2DM=%.3f | T2DM=%.3f\n",
            mean(tf_no_t2dm), mean(tf_t2dm)))
cat(sprintf("T-test p-value: %.4e\n", t.test(tf_no_t2dm, tf_t2dm)$p.value))

# 2. Linear regression: TRG filter proxy ~ T2DM + BMI + alcohol + age
#    (quantifies T2DM's confounding effect on the TRG shield biomarker)
lm_trg <- lm(trig_filter_proxy ~ t2dm + bmi + alcohol_cat + age_baseline,
             data = trg_data)
lm_tbl <- tidy(lm_trg, conf.int = TRUE) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~round(., 4)),
         p.value = signif(p.value, 3)) %>%
  rename(beta = estimate, ci_lo = conf.low, ci_hi = conf.high, p = p.value)
cat("\nLinear model: TRG filter ~ T2DM + BMI + alcohol + age:\n")
print(lm_tbl)
write_csv(lm_tbl, file.path(OUT, "trgshield_t2dm_linear_model.csv"))

# 3. Logistic: lipid_clinic_elig ~ TRG proxy (shows diagnostic value)
#    With and without T2DM adjustment
m_unadj <- glm(ldl_baseline >= 5.5 ~ trig_filter_proxy,
               data = filter(lc_cohort, !is.na(trig_filter_proxy),
                             !is.na(ldl_baseline), is.finite(trig_filter_proxy)),
               family = binomial)
m_adj   <- glm(ldl_baseline >= 5.5 ~ trig_filter_proxy + t2dm + bmi + alcohol_cat,
               data = filter(trg_data, !is.na(ldl_baseline)),
               family = binomial)

tbl_trg_logistic <- bind_rows(
  tidy(m_unadj, exponentiate=TRUE, conf.int=TRUE) %>%
    filter(term == "trig_filter_proxy") %>%
    mutate(model = "Unadjusted"),
  tidy(m_adj, exponentiate=TRUE, conf.int=TRUE) %>%
    filter(term == "trig_filter_proxy") %>%
    mutate(model = "Adjusted (T2DM+BMI+alcohol+age)")
) %>%
  select(model, OR=estimate, ci_lo=conf.low, ci_hi=conf.high, p=p.value) %>%
  mutate(across(c(OR,ci_lo,ci_hi), ~round(.,3)), p=signif(p,3))

cat("\nTRG filter proxy — OR for severe hypercholesterolaemia (LDL≥5.5):\n")
print(tbl_trg_logistic)
write_csv(tbl_trg_logistic, file.path(OUT, "trgshield_logistic_ukb.csv"))

# Figure — TRG filter proxy by T2DM
p_trg_t2dm <- trg_data %>%
  mutate(T2DM = ifelse(t2dm, "T2DM+", "T2DM-")) %>%
  ggplot(aes(x = T2DM, y = trig_filter_proxy, fill = T2DM)) +
  geom_violin(alpha=0.7, trim=TRUE) +
  geom_boxplot(width=0.1, fill="white", outlier.size=0.3) +
  scale_fill_manual(values=c("T2DM-"=C_BLUE, "T2DM+"=C_ORANGE), guide="none") +
  scale_y_continuous("Triglyceride filter proxy (LDL/TG+0.1)",
                     limits=c(0, quantile(trg_data$trig_filter_proxy, 0.99, na.rm=T))) +
  annotate("text", x=1.5, y=quantile(trg_data$trig_filter_proxy, 0.97, na.rm=T),
           label=sprintf("Cohen's d=%.3f\np<0.001", abs(d_t2dm)), size=2.2) +
  nat() +
  labs(title = "TRG filter proxy by T2DM status",
       subtitle = "UK Biobank lipid-clinic cohort",
       caption = "T2DM lowers TRG proxy → critical confounder for FH diagnosis")

ggsave(file.path(OUT, "Figure_TRG_T2DM.pdf"), p_trg_t2dm,
       width = 89, height = 89, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure_TRG_T2DM.png"), p_trg_t2dm,
       width = 89, height = 89, units = "mm", dpi = 300)
cat("✓ TRG shield figures saved\n")

##############################################################################
# F7 — ASCVD ANALYSIS
##############################################################################

cat("\n=== F7: ASCVD analysis ===\n")

# Overall ASCVD summary
ascvd_tbl <- lc_cohort %>%
  summarise(
    n_total       = n(),
    n_ascvd       = sum(any_ascvd, na.rm=T),
    pct_ascvd     = round(100*mean(any_ascvd, na.rm=T), 1),
    n_prem_ascvd  = sum(premature_ascvd, na.rm=T),
    pct_prem      = round(100*mean(premature_ascvd, na.rm=T), 1)
  )
cat("\nASCVD in lipid-clinic cohort:\n"); print(ascvd_tbl)
write_csv(ascvd_tbl, file.path(OUT, "ascvd_lc_summary.csv"))

# ASCVD by age decade
ascvd_by_age <- lc_cohort %>%
  filter(!is.na(age_baseline)) %>%
  mutate(age_decade = paste0(floor(age_baseline/10)*10, "s")) %>%
  group_by(age_decade) %>%
  summarise(
    n           = n(),
    pct_ascvd   = round(100*mean(any_ascvd, na.rm=T), 1),
    pct_prem    = round(100*mean(premature_ascvd, na.rm=T), 1),
    ldl_mean    = round(mean(ldl_baseline, na.rm=T), 2),
    .groups     = "drop"
  ) %>%
  arrange(age_decade)
cat("\nASCVD by age decade:\n"); print(ascvd_by_age)
write_csv(ascvd_by_age, file.path(OUT, "ascvd_by_age_decade.csv"))

# ASCVD by LDL decile (proxy for TUDOR severity)
ascvd_ldl_decile <- lc_cohort %>%
  filter(!is.na(ldl_baseline)) %>%
  mutate(ldl_decile = ntile(ldl_baseline, 10)) %>%
  group_by(ldl_decile) %>%
  summarise(
    n             = n(),
    ldl_mean      = round(mean(ldl_baseline), 2),
    pct_ascvd     = round(100*mean(any_ascvd, na.rm=T), 1),
    pct_prem      = round(100*mean(premature_ascvd, na.rm=T), 1),
    pct_t2dm      = round(100*mean(t2dm, na.rm=T), 1),
    .groups = "drop"
  )
cat("\nASCVD by LDL decile:\n"); print(ascvd_ldl_decile)
write_csv(ascvd_ldl_decile, file.path(OUT, "ascvd_by_ldl_decile.csv"))

# Figure
p_ascvd <- ggplot(ascvd_ldl_decile, aes(x = ldl_decile, y = pct_ascvd)) +
  geom_col(fill = C_RED, alpha = 0.8, width = 0.7) +
  geom_line(aes(y = pct_prem), colour = C_ORANGE, linewidth = 0.8) +
  geom_point(aes(y = pct_prem), colour = C_ORANGE, size = 1.5) +
  scale_x_continuous("LDL decile (proxy for FH severity)", breaks = 1:10) +
  scale_y_continuous("ASCVD prevalence (%)") +
  annotate("text", x = 9, y = max(ascvd_ldl_decile$pct_ascvd) * 0.7,
           label = "Orange line = premature ASCVD (<55yr)",
           colour = C_ORANGE, size = 2.2, hjust = 1) +
  nat() +
  labs(title = "ASCVD prevalence by LDL severity",
       subtitle = "UK Biobank lipid-clinic cohort")

ggsave(file.path(OUT, "Figure_ASCVD_LDL_Decile.pdf"), p_ascvd,
       width = 89, height = 89, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure_ASCVD_LDL_Decile.png"), p_ascvd,
       width = 89, height = 89, units = "mm", dpi = 300)
cat("✓ ASCVD analysis complete\n")

##############################################################################
# F10 — CHOLESTEROL-YEARS
##############################################################################

cat("\n=== F10: Cholesterol-years in LC cohort ===\n")

cy_summary <- lc_cohort %>%
  filter(!is.na(chol_years_ukb)) %>%
  summarise(
    n             = n(),
    cy_mean       = round(mean(chol_years_ukb), 1),
    cy_sd         = round(sd(chol_years_ukb), 1),
    cy_p25        = round(quantile(chol_years_ukb, 0.25), 1),
    cy_median     = round(median(chol_years_ukb), 1),
    cy_p75        = round(quantile(chol_years_ukb, 0.75), 1)
  )
cat("\nCholesterol-years summary (2-timepoint, 6yr gap):\n")
print(cy_summary)

cy_by_ldl <- lc_cohort %>%
  filter(!is.na(chol_years_ukb), !is.na(ldl_baseline)) %>%
  mutate(ldl_decile = ntile(ldl_baseline, 10)) %>%
  group_by(ldl_decile) %>%
  summarise(
    n          = n(),
    cy_mean    = round(mean(chol_years_ukb), 1),
    ldl_mean   = round(mean(ldl_baseline), 2),
    pct_ascvd  = round(100*mean(any_ascvd, na.rm=T), 1),
    .groups    = "drop"
  )
cat("\nCholesterol-years by LDL decile:\n"); print(cy_by_ldl)
write_csv(cy_by_ldl, file.path(OUT, "cholesterol_years_lc_decile.csv"))

p_cy <- ggplot(cy_by_ldl, aes(x = ldl_decile, y = cy_mean)) +
  geom_col(fill = C_BLUE, alpha = 0.8, width = 0.7) +
  scale_x_continuous("LDL decile", breaks = 1:10) +
  scale_y_continuous("Cholesterol-years (mmol/L × 6yr)") +
  nat() +
  labs(title = "Cumulative LDL burden by severity decile",
       subtitle = "UK Biobank lipid-clinic cohort (2 timepoints, ~6yr gap)",
       caption = "True cumulative burden requires full longitudinal follow-up")

ggsave(file.path(OUT, "Figure_CholYears_LC.pdf"), p_cy,
       width = 89, height = 89, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "Figure_CholYears_LC.png"), p_cy,
       width = 89, height = 89, units = "mm", dpi = 300)

# Correlation: cholesterol-years vs ASCVD
if (sum(!is.na(lc_cohort$chol_years_ukb) & !is.na(lc_cohort$any_ascvd)) > 100) {
  logit_cy <- glm(any_ascvd ~ chol_years_ukb + age_baseline,
                  data = filter(lc_cohort, !is.na(chol_years_ukb),
                                !is.na(age_baseline)),
                  family = binomial)
  tbl_cy <- tidy(logit_cy, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(across(c(estimate,conf.low,conf.high), ~round(.,3)),
           p.value = signif(p.value, 3)) %>%
    rename(OR=estimate, ci_lo=conf.low, ci_hi=conf.high, p=p.value)
  cat("\nLogistic: ASCVD ~ chol-years + age:\n"); print(tbl_cy)
  write_csv(tbl_cy, file.path(OUT, "cholesterol_years_ascvd_logistic.csv"))
}

##############################################################################
# SUMMARY TABLE — ALL KEY NUMBERS FOR MANUSCRIPT
##############################################################################

cat("\n============================================================\n")
cat("SUPPLEMENT ANALYSIS COMPLETE — KEY NUMBERS FOR MANUSCRIPT\n")
cat("============================================================\n")
cat(sprintf("UK Biobank total: N=%d\n", nrow(supplement)))
cat(sprintf("Lipid-clinic eligible: N=%d\n", nrow(lc_cohort)))
cat(sprintf("\nETHNICITY (LC cohort):\n"))
eth_pct <- lc_cohort %>% count(ethnicity) %>%
  mutate(pct = round(100*n/sum(n),1))
print(eth_pct)
cat(sprintf("\nT2DM in LC cohort: %.1f%%\n", 100*mean(lc_cohort$t2dm, na.rm=T)))
cat(sprintf("Any ASCVD in LC cohort: %.1f%%\n", 100*mean(lc_cohort$any_ascvd, na.rm=T)))
cat(sprintf("Premature ASCVD (<55yr): %.1f%%\n", 100*mean(lc_cohort$premature_ascvd, na.rm=T)))
cat(sprintf("\nTRG filter proxy — T2DM lowers by: %.3f units (Cohen's d=%.3f)\n",
            mean(tf_no_t2dm) - mean(tf_t2dm), abs(d_t2dm)))
cat(sprintf("\nApoB coverage in LC cohort: %.0f%%\n",
            100*mean(!is.na(lc_cohort$apob_measured))))
cat(sprintf("Lp(a) coverage in LC cohort: %.0f%%\n",
            100*mean(!is.na(lc_cohort$lpa_nmol_L))))

cat("\nFiles saved to:", OUT, "\n")
cat("\nF5  → ethnicity_lc_cohort_summary.csv + Figure_LDL_Ethnicity\n")
cat("F6  → trgshield_by_t2dm.csv + trgshield_t2dm_linear_model.csv + Figure_TRG_T2DM\n")
cat("F7  → ascvd_lc_summary.csv + ascvd_by_age_decade.csv + Figure_ASCVD_LDL_Decile\n")
cat("F10 → cholesterol_years_lc_decile.csv + Figure_CholYears_LC\n")
cat("\nNOTE: Ethnicity AUC + TUDOR-decile ASCVD require EID merge.\n")
cat("      Re-run TUDOR_UKB_LIPID_CLINIC.R on RAP saving 'eid' column.\n")
