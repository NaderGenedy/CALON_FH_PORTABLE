################################################################################
#                                                                              #
#  CALON-2 v2: PUBLICATION-READY SUMMARY TABLES                              #
#  Consolidate results from scripts 08 and 09                                 #
#                                                                              #
#  Generates CSV tables suitable for paper / supplementary materials           #
#                                                                              #
#  IMPORTANT: Run this script AFTER running scripts 08 and 09                 #
#             and pasting the hardcoded results below.                         #
#                                                                              #
#  Author:  Dr Nader Genedy                                                    #
#  Date:    March 2026                                                         #
#                                                                              #
################################################################################

rm(list = ls())

# =============================================================================
# CONFIGURATION
# =============================================================================

CLOUD_PATH  <- "/cloud/project/"
LOCAL_PATH  <- "C:/Users/nader/Downloads/calon_ukb_pipeline/"

if (dir.exists(CLOUD_PATH) && file.exists(paste0(CLOUD_PATH, "project.Rproj"))) {
    BASE_DIR <- CLOUD_PATH
    cat("  Environment: Posit Cloud\n")
} else if (dir.exists(LOCAL_PATH)) {
    BASE_DIR <- LOCAL_PATH
    cat("  Environment: Local\n")
} else {
    stop("No valid base directory found")
}

TAB_DIR <- paste0(BASE_DIR, "output/tables/")
if (!dir.exists(TAB_DIR)) dir.create(TAB_DIR, recursive = TRUE)

cat("\n")
cat("======================================================================\n")
cat("  CALON-2 v2: PUBLICATION-READY TABLES                              \n")
cat("======================================================================\n\n")

# =============================================================================
# TABLE 1: COHORT BASELINE CHARACTERISTICS
# =============================================================================

cat("===================================================================\n")
cat("TABLE 1: COHORT BASELINE CHARACTERISTICS\n")
cat("===================================================================\n\n")

# Hardcoded from data exploration (update with exact values after final runs)
table1 <- data.frame(
    Variable = c("N", "ASCVD events, n (%)",
                 "Age, years, mean (SD)", "Male, n (%)",
                 "LDL-C, mmol/L, mean (SD)", "HDL-C, mmol/L, mean (SD)",
                 "Total cholesterol, mmol/L, mean (SD)",
                 "Triglycerides, mmol/L, mean (SD)",
                 "Current/ex-smoker, n (%)", "Diabetes, n (%)",
                 "Hypertension, n (%)", "BMI, kg/m2, mean (SD)",
                 "APOB mutation, n (%)"),
    UKB_Development = c("1,619", "397 (24.5%)",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE"),
    Wales_External1 = c("7,253 (3,038 CC)", "826 (27.2%, CC)",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE", "UPDATE",
                        "UPDATE"),
    DRAGON3_External2 = c("424", "62 (14.6%)",
                          "UPDATE", "UPDATE",
                          "UPDATE", "UPDATE",
                          "UPDATE", "UPDATE",
                          "UPDATE", "UPDATE",
                          "UPDATE", "UPDATE",
                          "UPDATE"),
    stringsAsFactors = FALSE
)

cat("  Table 1 template generated (fill UPDATE values from data).\n")
cat("  NOTE: Run the fill_table1() function below after loading actual data.\n\n")

# Function to fill Table 1 from actual data
fill_table1 <- function(ukb, wales, dragon3) {
    fmt_mean_sd <- function(x) {
        x <- x[!is.na(x)]
        sprintf("%.1f (%.1f)", mean(x), sd(x))
    }
    fmt_n_pct <- function(x) {
        x <- x[!is.na(x)]
        sprintf("%d (%.1f%%)", sum(x == 1), 100 * mean(x == 1))
    }

    # UKB column
    ukb_vals <- c(
        as.character(nrow(ukb)),
        sprintf("%d (%.1f%%)", sum(ukb$ascvd_combined), 100*mean(ukb$ascvd_combined)),
        fmt_mean_sd(ukb$age), fmt_n_pct(ukb$sex),
        fmt_mean_sd(ukb$re_ldl), fmt_mean_sd(ukb$hdl),
        fmt_mean_sd(ukb$tc), fmt_mean_sd(ukb$trig),
        fmt_n_pct(ukb$smoking_binary), fmt_n_pct(ukb$diabetes),
        fmt_n_pct(ukb$hypertension), fmt_mean_sd(ukb$bmi),
        fmt_n_pct(ukb$gene_apob))

    # Wales column
    wales_vals <- c(
        as.character(nrow(wales)),
        sprintf("%d (%.1f%%)", sum(wales$ascvd_combined), 100*mean(wales$ascvd_combined)),
        fmt_mean_sd(wales$age), fmt_n_pct(wales$sex),
        fmt_mean_sd(wales$re_ldl), fmt_mean_sd(wales$hdl),
        fmt_mean_sd(wales$tc), fmt_mean_sd(wales$trig),
        fmt_n_pct(wales$smoking_binary), fmt_n_pct(wales$diabetes),
        fmt_n_pct(wales$hypertension), fmt_mean_sd(wales$bmi),
        fmt_n_pct(wales$gene_apob))

    # DRAGON3 column
    d3_vals <- c(
        as.character(nrow(dragon3)),
        sprintf("%d (%.1f%%)", sum(dragon3$ascvd_combined), 100*mean(dragon3$ascvd_combined)),
        fmt_mean_sd(dragon3$age), fmt_n_pct(dragon3$sex),
        fmt_mean_sd(dragon3$re_ldl), fmt_mean_sd(dragon3$hdl),
        fmt_mean_sd(dragon3$tc), fmt_mean_sd(dragon3$trig),
        fmt_n_pct(dragon3$smoking_binary), fmt_n_pct(dragon3$diabetes),
        fmt_n_pct(dragon3$hypertension),
        ifelse(sum(!is.na(dragon3$bmi)) > 10, fmt_mean_sd(dragon3$bmi), "N/A"),
        fmt_n_pct(dragon3$gene_apob))

    data.frame(Variable = table1$Variable,
               UKB_Development = ukb_vals,
               Wales_External1 = wales_vals,
               DRAGON3_External2 = d3_vals,
               stringsAsFactors = FALSE)
}

# =============================================================================
# TABLE 2: MODEL DISCRIMINATION — ALL TIERS ACROSS COHORTS
# =============================================================================

cat("===================================================================\n")
cat("TABLE 2: MODEL DISCRIMINATION ACROSS COHORTS\n")
cat("===================================================================\n\n")

# Results from script 08 (Posit Cloud output)
table2 <- data.frame(
    Model = c("EN v1-A (original, 10 vars)",
              "EN v2-A (lean clean, 8 vars)",
              "EN v2-A (recalibrated)",
              "EN v2-B (enriched, 11 vars)",
              "EN v2-C (LDL-focused, 12 vars)",
              "EN v2-D (full, 15 vars)",
              "XGB v2-A (lean clean)",
              "SAFEHEART-RE (published)"),
    Active_Vars = c(10, 7, 7, 7, 7, 9, 8, 6),
    Max_VIF = c("13.49*", "2.30", "2.30", "2.37", "72.00*", "230.28*", "-", "-"),
    UKB_CV_AUC = c("0.7250", "0.7247", "0.7247", "0.7287", "0.7265", "0.7287",
                   "0.7162", "0.7095"),
    Wales_AUC = c("0.7522", "0.7551", "0.7551", "0.7479", "0.7411", "0.7327",
                  "0.7389", "0.7211"),
    Wales_95CI = c("0.7341-0.7703", "0.7369-0.7733", "0.7369-0.7733",
                   "0.7293-0.7664", "0.7221-0.7602", "0.7134-0.7520",
                   "0.7199-0.7580", "0.6954-0.7468"),
    DRAGON3_AUC = c("RUN 09", "RUN 09", "RUN 09", "RUN 09", "-", "-",
                    "-", "RUN 09"),
    DRAGON3_95CI = c("RUN 09", "RUN 09", "RUN 09", "RUN 09", "-", "-",
                     "-", "RUN 09"),
    stringsAsFactors = FALSE
)

cat("  * = collinear (VIF > 5)\n")
cat("  v2-A AUC unchanged after recalibration (monotonic intercept shift)\n\n")

write.csv(table2, paste0(TAB_DIR, "Table2_discrimination.csv"), row.names = FALSE)
cat(sprintf("  Saved: Table2_discrimination.csv\n"))

# Print formatted
cat("\n")
cat(sprintf("  %-35s %4s %8s %8s %8s %-18s %8s\n",
            "Model", "Vars", "VIF", "UKB", "Wales", "Wales 95%CI", "DRAGON3"))
cat(paste(rep("-", 110), collapse = ""), "\n")
for (i in 1:nrow(table2)) {
    cat(sprintf("  %-35s %4d %8s %8s %8s %-18s %8s\n",
                table2$Model[i], table2$Active_Vars[i], table2$Max_VIF[i],
                table2$UKB_CV_AUC[i], table2$Wales_AUC[i], table2$Wales_95CI[i],
                table2$DRAGON3_AUC[i]))
}

# =============================================================================
# TABLE 3: CALIBRATION METRICS
# =============================================================================

cat("\n===================================================================\n")
cat("TABLE 3: CALIBRATION METRICS (WALES)\n")
cat("===================================================================\n\n")

table3 <- data.frame(
    Model = c("EN v1-A (original)", "EN v2-A (lean clean)",
              "EN v2-A (recalibrated)", "EN v2-B (enriched)",
              "SAFEHEART-RE"),
    Wales_Cal_Intercept = c("+0.0318", "-0.3655", "UPDATE", "-0.0879", "UPDATE"),
    Wales_Cal_Slope = c("0.8987", "0.9219", "UPDATE", "0.9007", "UPDATE"),
    Brier_Score = c("UPDATE", "UPDATE", "UPDATE", "UPDATE", "UPDATE"),
    DRAGON3_Cal_Intercept = c("RUN 09", "RUN 09", "RUN 09", "RUN 09", "RUN 09"),
    DRAGON3_Cal_Slope = c("RUN 09", "RUN 09", "RUN 09", "RUN 09", "RUN 09"),
    stringsAsFactors = FALSE
)

write.csv(table3, paste0(TAB_DIR, "Table3_calibration.csv"), row.names = FALSE)
cat(sprintf("  Saved: Table3_calibration.csv\n"))

cat("\n")
cat(sprintf("  %-28s %12s %10s %10s %12s %10s\n",
            "Model", "W.Intercept", "W.Slope", "Brier", "D3.Intercept", "D3.Slope"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in 1:nrow(table3)) {
    cat(sprintf("  %-28s %12s %10s %10s %12s %10s\n",
                table3$Model[i], table3$Wales_Cal_Intercept[i],
                table3$Wales_Cal_Slope[i], table3$Brier_Score[i],
                table3$DRAGON3_Cal_Intercept[i], table3$DRAGON3_Cal_Slope[i]))
}

# =============================================================================
# TABLE 4: HEAD-TO-HEAD — v2-A vs SAFEHEART-RE
# =============================================================================

cat("\n===================================================================\n")
cat("TABLE 4: HEAD-TO-HEAD — CALON-2 v2-A vs SAFEHEART-RE\n")
cat("===================================================================\n\n")

table4 <- data.frame(
    Metric = c("AUC (CALON-2 v2-A)", "AUC (SAFEHEART-RE)",
               "Delta AUC", "DeLong P-value",
               "Bootstrap Delta 95% CI (2000 rep)",
               "Category-free NRI", "NRI 95% CI",
               "IDI", "IDI P-value", "IDI 95% CI",
               "Category NRI (10%/20% cuts)", "NRI P-value"),
    Wales = c("0.7551", "0.7211",
              "+0.0340", "0.0052",
              "(+0.0066, +0.0361)",
              "+0.0476", "(-0.0176, +0.1128)",
              "+0.0333", "<0.0001", "(+0.0221, +0.0445)",
              "-0.0100", "UPDATE"),
    DRAGON3 = c("RUN 09", "RUN 09",
                "RUN 09", "RUN 09",
                "RUN 09",
                "RUN 09", "RUN 09",
                "RUN 09", "RUN 09", "RUN 09",
                "RUN 09", "RUN 09"),
    stringsAsFactors = FALSE
)

write.csv(table4, paste0(TAB_DIR, "Table4_headtohead.csv"), row.names = FALSE)
cat(sprintf("  Saved: Table4_headtohead.csv\n"))

cat("\n")
cat(sprintf("  %-40s %-24s %-20s\n", "Metric", "Wales", "DRAGON3"))
cat(paste(rep("-", 85), collapse = ""), "\n")
for (i in 1:nrow(table4)) {
    cat(sprintf("  %-40s %-24s %-20s\n",
                table4$Metric[i], table4$Wales[i], table4$DRAGON3[i]))
}

# =============================================================================
# TABLE 5: SUBGROUP ANALYSIS — delta AUC by subgroup
# =============================================================================

cat("\n===================================================================\n")
cat("TABLE 5: SUBGROUP ANALYSIS (v2-A vs SAFEHEART-RE, Wales)\n")
cat("===================================================================\n\n")

table5 <- data.frame(
    Subgroup = c("HDL < 1.0 mmol/L", "Female", "Male",
                 "Age < 50", "Age 50-64", "Age >= 65",
                 "Non-smoker", "Smoker/ex-smoker",
                 "Diabetic", "Non-diabetic",
                 "Hypertensive", "Non-hypertensive",
                 "APOB mutation"),
    N_Wales = c("UPDATE", "UPDATE", "UPDATE",
                "UPDATE", "UPDATE", "UPDATE",
                "UPDATE", "UPDATE",
                "UPDATE", "UPDATE",
                "UPDATE", "UPDATE",
                "UPDATE"),
    Delta_AUC = c("+0.0674", "+0.0506", "+0.0402",
                  "UPDATE", "UPDATE", "UPDATE",
                  "UPDATE", "UPDATE",
                  "UPDATE", "UPDATE",
                  "UPDATE", "UPDATE",
                  "UPDATE"),
    P_value = c("0.0008", "<0.0001", "UPDATE",
                "UPDATE", "UPDATE", "UPDATE",
                "UPDATE", "UPDATE",
                "UPDATE", "UPDATE",
                "UPDATE", "UPDATE",
                "UPDATE"),
    DRAGON3_Delta = c("RUN 09", "RUN 09", "RUN 09",
                      "RUN 09", "RUN 09", "RUN 09",
                      "RUN 09", "RUN 09",
                      "RUN 09", "RUN 09",
                      "RUN 09", "RUN 09",
                      "RUN 09"),
    stringsAsFactors = FALSE
)

write.csv(table5, paste0(TAB_DIR, "Table5_subgroups.csv"), row.names = FALSE)
cat(sprintf("  Saved: Table5_subgroups.csv\n"))

cat("\n")
cat(sprintf("  %-25s %8s %10s %10s %12s\n",
            "Subgroup", "N", "Delta", "P", "DRAGON3"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in 1:nrow(table5)) {
    cat(sprintf("  %-25s %8s %10s %10s %12s\n",
                table5$Subgroup[i], table5$N_Wales[i],
                table5$Delta_AUC[i], table5$P_value[i],
                table5$DRAGON3_Delta[i]))
}

# =============================================================================
# TABLE 6: FROZEN v2-A COEFFICIENTS (for clinical deployment)
# =============================================================================

cat("\n===================================================================\n")
cat("TABLE 6: FROZEN v2-A COEFFICIENTS\n")
cat("===================================================================\n\n")

table6 <- data.frame(
    Variable = c("(Intercept)", "Age (years)", "Sex (1=male)",
                 "LDL-C (mmol/L)", "HDL-C (mmol/L)",
                 "log(TG/HDL)", "Current/ex-smoker (1=yes)",
                 "Diabetes (1=yes)", "Hypertension (1=yes)"),
    Coefficient = c(-4.201192, 0.063350, 0.151549,
                    0.031590, -0.905700,
                    0.000000, 0.413440,
                    0.683894, 0.622197),
    Odds_Ratio = c(NA, exp(0.063350), exp(0.151549),
                   exp(0.031590), exp(-0.905700),
                   NA, exp(0.413440),
                   exp(0.683894), exp(0.622197)),
    Status = c("-", "Active", "Active",
               "Active", "Active",
               "Zeroed by EN", "Active",
               "Active", "Active"),
    Recal_Intercept = c(-4.201192 + 0.3655, NA, NA,
                        NA, NA,
                        NA, NA,
                        NA, NA),
    stringsAsFactors = FALSE
)

# Format nicely
table6$Coefficient <- sprintf("%.6f", table6$Coefficient)
table6$Odds_Ratio <- ifelse(is.na(table6$Odds_Ratio), "-",
                            sprintf("%.3f", as.numeric(gsub("NA", "0", table6$Odds_Ratio))))
# Fix the OR column
or_vals <- c(NA, exp(0.063350), exp(0.151549), exp(0.031590), exp(-0.905700),
             NA, exp(0.413440), exp(0.683894), exp(0.622197))
table6$Odds_Ratio <- ifelse(is.na(or_vals), "-", sprintf("%.3f", or_vals))
table6$Recal_Intercept <- ifelse(is.na(table6$Recal_Intercept), "-",
                                 sprintf("%.6f", as.numeric(table6$Recal_Intercept)))
# fix the recal col
recal_vals <- c(-4.201192 + 0.3655, rep(NA, 8))
table6$Recal_Intercept <- ifelse(is.na(recal_vals), "-", sprintf("%.6f", recal_vals))

write.csv(table6, paste0(TAB_DIR, "Table6_frozen_coefficients.csv"), row.names = FALSE)
cat(sprintf("  Saved: Table6_frozen_coefficients.csv\n"))

cat("\n")
cat(sprintf("  %-28s %12s %8s %14s %14s\n",
            "Variable", "Beta", "OR", "Status", "Recal.Int"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in 1:nrow(table6)) {
    cat(sprintf("  %-28s %12s %8s %14s %14s\n",
                table6$Variable[i], table6$Coefficient[i],
                table6$Odds_Ratio[i], table6$Status[i],
                table6$Recal_Intercept[i]))
}

# =============================================================================
# TABLE 7: DECISION CURVE ANALYSIS SUMMARY
# =============================================================================

cat("\n===================================================================\n")
cat("TABLE 7: DECISION CURVE ANALYSIS (Wales, key thresholds)\n")
cat("===================================================================\n\n")

table7 <- data.frame(
    Threshold = c("5%", "10%", "15%", "20%", "25%", "30%", "40%", "50%"),
    CALON2_v2A_Net_Benefit = c("UPDATE", "UPDATE", "UPDATE", "UPDATE",
                                "UPDATE", "UPDATE", "UPDATE", "UPDATE"),
    SAFEHEART_Net_Benefit = c("UPDATE", "UPDATE", "UPDATE", "UPDATE",
                              "UPDATE", "UPDATE", "UPDATE", "UPDATE"),
    Treat_All = c("UPDATE", "UPDATE", "UPDATE", "UPDATE",
                  "UPDATE", "UPDATE", "UPDATE", "UPDATE"),
    v2A_Advantage = c("UPDATE", "higher", "higher", "higher",
                      "higher", "higher", "higher", "UPDATE"),
    stringsAsFactors = FALSE
)

write.csv(table7, paste0(TAB_DIR, "Table7_DCA.csv"), row.names = FALSE)
cat(sprintf("  Saved: Table7_DCA.csv\n"))
cat("  NOTE: v2-A had higher net benefit at thresholds 10%%-48%% on Wales.\n")

# =============================================================================
# SUPPLEMENTARY TABLE S1: SAFEHEART-RE BENCHMARK COMPARISON
# =============================================================================

cat("\n===================================================================\n")
cat("SUPP TABLE S1: SAFEHEART-RE PUBLISHED vs CALON-2 v2-A\n")
cat("===================================================================\n\n")

table_s1 <- data.frame(
    Feature = c("Number of variables", "Active predictors",
                "Modelling method", "Development cohort",
                "Includes BMI", "Includes sex (OR weight)",
                "Uses LDL-C directly", "Uses HDL-C",
                "Collinearity (max VIF)",
                "UKB CV-AUC", "Wales AUC",
                "DRAGON3 AUC",
                "IDI vs comparator"),
    SAFEHEART_RE = c("6", "6",
                     "Logistic regression", "SAFEHEART (UK)",
                     "Yes (beta=0.025)", "Yes (OR=2.17)",
                     "Yes (beta=0.109)", "No",
                     "Not reported",
                     "0.7095", "0.7211",
                     "RUN 09",
                     "(Reference)"),
    CALON2_v2A = c("8 (7 active)", "7",
                   "Elastic Net (alpha=0.5)", "UK Biobank FH",
                   "No", "Yes (beta=0.152)",
                   "Yes (beta=0.032)", "Yes (beta=-0.906)",
                   "2.30",
                   "0.7247", "0.7551",
                   "RUN 09",
                   "+0.0333, P<0.0001"),
    stringsAsFactors = FALSE
)

write.csv(table_s1, paste0(TAB_DIR, "TableS1_safeheart_comparison.csv"), row.names = FALSE)
cat(sprintf("  Saved: TableS1_safeheart_comparison.csv\n"))

cat("\n")
cat(sprintf("  %-32s %-28s %-28s\n", "Feature", "SAFEHEART-RE", "CALON-2 v2-A"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in 1:nrow(table_s1)) {
    cat(sprintf("  %-32s %-28s %-28s\n",
                table_s1$Feature[i], table_s1$SAFEHEART_RE[i],
                table_s1$CALON2_v2A[i]))
}

# =============================================================================
# FIGURE GUIDANCE
# =============================================================================

cat("\n===================================================================\n")
cat("FIGURE GUIDANCE (for paper)\n")
cat("===================================================================\n\n")

cat("  Figure 1: Study flow diagram (TRIPOD)\n")
cat("    - UKB → development → internal 10-fold CV\n")
cat("    - Wales → external validation 1\n")
cat("    - DRAGON3 → external validation 2\n\n")

cat("  Figure 2: ROC curves (Wales external validation)\n")
cat("    - CALON-2 v2-A (AUC=0.7551) vs SAFEHEART-RE (AUC=0.7211)\n")
cat("    - Include 95%% CI shading, DeLong P-value annotation\n\n")

cat("  Figure 3: Calibration plots (Wales)\n")
cat("    - Panel A: v2-A original (intercept=-0.3655)\n")
cat("    - Panel B: v2-A recalibrated (intercept≈0)\n")
cat("    - Panel C: SAFEHEART-RE\n\n")

cat("  Figure 4: Decision Curve Analysis (Wales)\n")
cat("    - Net benefit: v2-A vs SAFEHEART vs treat-all vs treat-none\n")
cat("    - Threshold range: 1%%-50%%\n\n")

cat("  Figure 5: Forest plot — Subgroup analysis\n")
cat("    - Delta AUC by subgroup with 95%% CI\n\n")

cat("  Supplementary Figure S1: ROC curves (DRAGON3)\n")
cat("  Supplementary Figure S2: Calibration plots (DRAGON3)\n")
cat("  Supplementary Figure S3: Bootstrap distribution of delta AUC\n")

# =============================================================================
# DONE
# =============================================================================

cat("\n===================================================================\n")
cat("DONE — PUBLICATION TABLES\n")
cat("===================================================================\n\n")

cat("  OUTPUT FILES:\n")
cat(sprintf("    %sTable2_discrimination.csv\n", TAB_DIR))
cat(sprintf("    %sTable3_calibration.csv\n", TAB_DIR))
cat(sprintf("    %sTable4_headtohead.csv\n", TAB_DIR))
cat(sprintf("    %sTable5_subgroups.csv\n", TAB_DIR))
cat(sprintf("    %sTable6_frozen_coefficients.csv\n", TAB_DIR))
cat(sprintf("    %sTable7_DCA.csv\n", TAB_DIR))
cat(sprintf("    %sTableS1_safeheart_comparison.csv\n", TAB_DIR))

cat("\n  NOTE: Cells marked 'UPDATE' require values from the Posit Cloud run.\n")
cat("        Cells marked 'RUN 09' will be filled after running script 09.\n")
cat("        After both runs, update this script with final values and re-run.\n")
cat(sprintf("\n  Completed at: %s\n\n", Sys.time()))
