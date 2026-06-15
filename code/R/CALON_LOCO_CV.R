#!/usr/bin/env Rscript
# CALON-Lite LOCO-CV Pipeline
# LOCO-CV for ASCVD Prediction in FH (3 cohorts)
# Uses median imputation for missing values

cat("=== CALON-Lite LOCO-CV Pipeline ===
")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "

")

pkgs <- c("tidyverse", "pROC", "boot", "ggplot2", "gridExtra", "cowplot")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(tidyverse); library(pROC); library(boot)
  library(ggplot2); library(gridExtra); library(cowplot)
})

if (dir.exists("/cloud/project")) {
  data_dir <- "/cloud/project"
} else if (dir.exists("C:/Users/nader/Downloads/calon_ukb_pipeline")) {
  data_dir <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
} else {
  data_dir <- getwd()
}
cat("Data directory:", data_dir, "

")
