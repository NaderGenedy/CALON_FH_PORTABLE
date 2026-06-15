################################################################################
# AF3 MANUSCRIPT FIGURES - Nature/Lancet Calibre
# AlphaFold3 Structural Biology in Familial Hypercholesterolemia
# 15 Multi-Panel Figures
# Author: Dr Nader Genedy
# Date: March 2026
################################################################################

# ==============================================================================
# SECTION 0: PACKAGES, THEME, PALETTE, DATA LOADING, HELPERS
# ==============================================================================

required_pkgs <- c("ggplot2", "patchwork", "dplyr", "tidyr", "pROC", "ggrepel",
                   "scales", "RColorBrewer", "grid", "gridExtra", "data.table",
                   "broom", "forcats", "stringr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# --- Paths ---
BASE     <- "C:/Users/nader/Downloads/calon_ukb_pipeline"
AF_DIR   <- file.path(BASE, "alphafold/analysis")
FIG_DIR  <- file.path(BASE, "figures")
NMR_DIR  <- file.path(BASE, "alphafold")
MRI_DIR  <- file.path(BASE, "alphafold")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Nature theme ---
SINGLE_COL <- 88   # mm
DOUBLE_COL <- 180  # mm
mm2in <- function(mm) mm / 25.4

nature_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    text             = element_text(size = 7, family = "Arial"),
    axis.text        = element_text(size = 6, color = "black"),
    axis.title       = element_text(size = 7, face = "bold"),
    plot.title       = element_text(size = 8, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 6, color = "grey30"),
    legend.text      = element_text(size = 6),
    legend.title     = element_text(size = 6, face = "bold"),
    legend.key.size  = unit(3, "mm"),
    strip.text       = element_text(size = 6, face = "bold"),
    strip.background = element_blank(),
    panel.grid       = element_blank(),
    plot.margin      = margin(2, 2, 2, 2, "mm"),
    legend.margin    = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )
theme_set(nature_theme)

# --- Colour palettes ---
gene_pal   <- c("LDLR" = "#2166AC", "APOB" = "#B2182B", "PCSK9" = "#35978F")
domain_pal <- c(
  "Signal peptide"          = "#FDB863",
  "Ligand-binding"          = "#2166AC",
  "Ligand-binding R1"       = "#2166AC",
  "Ligand-binding R2"       = "#4393C3",
  "Ligand-binding R3"       = "#92C5DE",
  "Ligand-binding R4"       = "#4393C3",
  "Ligand-binding R5"       = "#2166AC",
  "Ligand-binding R6"       = "#92C5DE",
  "Ligand-binding R7"       = "#D1E5F0",
  "EGF-precursor-homology"  = "#B2182B",
  "EGF-A"                   = "#D6604D",
  "EGF-B"                   = "#F4A582",
  "Beta-propeller"          = "#FDDBC7",
  "EGF-C"                   = "#D6604D",
  "O-linked sugar"          = "#35978F",
  "Transmembrane"           = "#762A83",
  "Cytoplasmic"             = "#9970AB",
  "Linker (EGF-C/O-linked)" = "#66C2A5",
  "ApoB-RBD"               = "#B2182B",
  "Prodomain"               = "#FDB863",
  "Catalytic"               = "#2166AC",
  "V-domain"                = "#35978F",
  "Cys-His-rich"            = "#D6604D"
)

plddt_pal <- c("very_low" = "#FF7043", "low" = "#FFA726",
               "high" = "#66BB6A", "very_high" = "#1565C0")

ddg_class_pal <- c("Stabilising" = "#2166AC", "Neutral" = "#BDBDBD",
                   "Destabilising" = "#EF8A62", "Highly destabilising" = "#B2182B")

sss_tert_pal <- c("Low" = "#2166AC", "Medium" = "#FDB863", "High" = "#B2182B")

# --- Helper: safe column finder ---
find_col <- function(df, candidates) {
  cn <- colnames(df)
  for (c in candidates) {
    m <- grep(paste0("^", gsub("\\.", "\\\\.", c), "$"), cn, ignore.case = TRUE, value = TRUE)
    if (length(m) > 0) return(m[1])
    m2 <- grep(gsub("[_.]", "[_.]", c), cn, ignore.case = TRUE, value = TRUE)
    if (length(m2) > 0) return(m2[1])
  }
  return(NA_character_)
}

# --- Helper: save figure as PDF + PNG ---
save_figure <- function(p, name, w_mm = DOUBLE_COL, h_mm = 160) {
  w_in <- mm2in(w_mm)
  h_in <- mm2in(h_mm)
  pdf_path <- file.path(FIG_DIR, paste0(name, ".pdf"))
  png_path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(pdf_path, plot = p, width = w_in, height = h_in, dpi = 300, device = cairo_pdf)
  ggsave(png_path, plot = p, width = w_in, height = h_in, dpi = 600, device = "png",
         type = "cairo")
  message("  Saved: ", pdf_path)
  message("  Saved: ", png_path)
}

# --- Helper: Wilson confidence interval ---
wilson_ci <- function(x, n, z = 1.96) {
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  data.frame(est = p, lo = pmax(centre - margin, 0), hi = pmin(centre + margin, 1))
}

# --- Helper: add SSS tertiles ---
add_sss_tertile <- function(df, col = "sss") {
  vals <- df[[col]]
  brks <- quantile(vals, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if (brks[2] == brks[3]) brks[3] <- brks[3] + .Machine$double.eps
  df$sss_tertile <- cut(vals, breaks = brks, include.lowest = TRUE,
                        labels = c("Low", "Medium", "High"))
  df
}

# --- Helper: classify ddG ---
classify_ddg <- function(x) {
  ifelse(is.na(x), "Neutral",
         ifelse(x < -1, "Stabilising",
                ifelse(x <= 1, "Neutral",
                       ifelse(x <= 4, "Destabilising", "Highly destabilising"))))
}

# ==============================================================================
# SECTION 0b: LOAD CORE DATA
# ==============================================================================
message("Loading core datasets...")

comp_sss <- tryCatch(read.csv(file.path(AF_DIR, "comprehensive_sss_analysis.csv"),
                               stringsAsFactors = FALSE),
                     error = function(e) { warning("comp_sss: ", e$message); NULL })

nobel <- tryCatch(read.csv(file.path(AF_DIR, "nobel_multimodal_integration.csv"),
                           stringsAsFactors = FALSE),
                  error = function(e) { warning("nobel: ", e$message); NULL })

sss_scores <- tryCatch(read.csv(file.path(AF_DIR, "structural_severity_scores.csv"),
                                stringsAsFactors = FALSE),
                       error = function(e) { warning("sss_scores: ", e$message); NULL })

foldx <- tryCatch(read.csv(file.path(AF_DIR, "foldx_ddg_results.csv"),
                           stringsAsFactors = FALSE),
                  error = function(e) { warning("foldx: ", e$message); NULL })

llt_response <- tryCatch(read.csv(file.path(AF_DIR, "variant_llt_response.csv"),
                                  stringsAsFactors = FALSE),
                         error = function(e) { warning("llt_response: ", e$message); NULL })

loco_cv <- tryCatch(read.csv(file.path(AF_DIR, "loco_cv_pooled_sss_results.csv"),
                             stringsAsFactors = FALSE),
                    error = function(e) { warning("loco_cv: ", e$message); NULL })

plddt_ldlr <- tryCatch(read.csv(file.path(AF_DIR, "LDLR_wildtype_plddt_per_residue.csv"),
                                stringsAsFactors = FALSE),
                       error = function(e) { warning("plddt_ldlr: ", e$message); NULL })

plddt_pcsk9 <- tryCatch(read.csv(file.path(AF_DIR, "PCSK9_wildtype_plddt_per_residue.csv"),
                                 stringsAsFactors = FALSE),
                        error = function(e) { warning("plddt_pcsk9: ", e$message); NULL })

plddt_apob <- tryCatch(read.csv(file.path(AF_DIR, "ApoB_RBD_plddt_per_residue.csv"),
                                stringsAsFactors = FALSE),
                       error = function(e) { warning("plddt_apob: ", e$message); NULL })

subgroup <- tryCatch(read.csv(file.path(BASE, "Subgroup_analysis_results.csv"),
                              stringsAsFactors = FALSE),
                     error = function(e) { warning("subgroup: ", e$message); NULL })

# Harmonise domain column for pLDDT files without domain
if (!is.null(plddt_pcsk9) && !"domain" %in% colnames(plddt_pcsk9)) {
  plddt_pcsk9$domain <- ifelse(plddt_pcsk9$residue_num <= 152, "Prodomain",
                        ifelse(plddt_pcsk9$residue_num <= 449, "Catalytic",
                        ifelse(plddt_pcsk9$residue_num <= 452, "Hinge",
                               "V-domain")))
}

if (!is.null(plddt_apob) && !"domain" %in% colnames(plddt_apob)) {
  plddt_apob$domain <- "ApoB-RBD"
}

# Add gene to llt_response if missing
if (!is.null(llt_response) && (!"gene" %in% colnames(llt_response) || all(is.na(llt_response$gene)))) {
  llt_response$gene <- ifelse(grepl("^LDLR", llt_response$variant), "LDLR",
                       ifelse(grepl("^APOB", llt_response$variant), "APOB",
                       ifelse(grepl("^PCSK9", llt_response$variant), "PCSK9", NA)))
}

message("Data loading complete.")

# ==============================================================================
# FIGURE 1: Study Design & Cohort
# ==============================================================================
figure_1 <- function() {
  message("Generating Figure 1: Study Design & Cohort...")

  # --- Panel a: Study flow diagram ---
  boxes <- data.frame(
    label = c("FH Cohort\nn = 1,362", "Genetic Testing\nPathogenic variants", "Gene Distribution\nLDLR / APOB / PCSK9",
              "SSS Mapping\nn = 1,126", "Analytic Cohort\nn = 418"),
    x = c(1, 2, 3, 4, 5), y = rep(0, 5),
    stringsAsFactors = FALSE
  )
  arrows_df <- data.frame(x = 1:4, xend = 2:5, y = rep(0, 4), yend = rep(0, 4))

  p1a <- ggplot() +
    geom_rect(data = boxes, aes(xmin = x - 0.42, xmax = x + 0.42, ymin = -0.25, ymax = 0.25),
              fill = c("#E3F2FD", "#BBDEFB", "#90CAF9", "#64B5F6", "#42A5F5"),
              color = "grey30", linewidth = 0.3) +
    geom_segment(data = arrows_df, aes(x = x + 0.42, xend = xend - 0.42, y = y, yend = yend),
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                 linewidth = 0.4, color = "grey30") +
    geom_text(data = boxes, aes(x = x, y = y, label = label), size = 2, lineheight = 0.85) +
    coord_cartesian(xlim = c(0.3, 5.7), ylim = c(-0.5, 0.5)) +
    theme_void(base_size = 7, base_family = "Arial") +
    labs(title = "Study flow") +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0))

  # --- Panel b: Gene distribution bar chart ---
  gene_tab <- comp_sss %>%
    count(gene) %>%
    mutate(pct = round(100 * n / sum(n), 1),
           lbl = paste0(n, " (", pct, "%)"))

  p1b <- ggplot(gene_tab, aes(x = reorder(gene, -n), y = n, fill = gene)) +
    geom_col(width = 0.6, color = "grey30", linewidth = 0.2) +
    geom_text(aes(label = lbl), vjust = -0.4, size = 2) +
    scale_fill_manual(values = gene_pal, guide = "none") +
    labs(x = NULL, y = "Number of patients", title = "Gene distribution") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  # --- Panel c: LDLR domain lollipop plot ---
  ldlr_vars <- sss_scores %>% filter(gene == "LDLR", !is.na(protein_position), !is.na(sss))

  domain_rects <- data.frame(
    xmin = c(1, 22, 353, 693, 750, 771),
    xmax = c(21, 352, 692, 749, 770, 860),
    label = c("Signal\npeptide", "Ligand-\nbinding", "EGF-precursor\nhomology",
              "O-linked\nsugar", "TM", "Cyto"),
    fill = c("#FDB863", "#2166AC", "#B2182B", "#35978F", "#762A83", "#9970AB"),
    stringsAsFactors = FALSE
  )

  p1c <- ggplot() +
    geom_rect(data = domain_rects,
              aes(xmin = xmin, xmax = xmax, ymin = -0.05, ymax = -0.01),
              fill = domain_rects$fill, alpha = 0.5) +
    geom_text(data = domain_rects,
              aes(x = (xmin + xmax) / 2, y = -0.08, label = label),
              size = 1.5, lineheight = 0.8) +
    geom_segment(data = ldlr_vars, aes(x = protein_position, xend = protein_position,
                                        y = 0, yend = sss), color = "grey50", linewidth = 0.2) +
    geom_point(data = ldlr_vars, aes(x = protein_position, y = sss, color = domain),
               size = 0.6, alpha = 0.7) +
    scale_color_manual(values = domain_pal, name = "Domain") +
    labs(x = "Protein position (aa)", y = "SSS", title = "LDLR variant map") +
    scale_x_continuous(breaks = seq(0, 860, 100)) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 1.5)))

  # --- Panel d: SSS distribution histogram ---
  p1d <- ggplot(comp_sss, aes(x = sss, fill = domain)) +
    geom_histogram(bins = 30, color = "white", linewidth = 0.15, alpha = 0.85) +
    scale_fill_manual(values = domain_pal, name = "Domain") +
    labs(x = "Structural Severity Score", y = "Count", title = "SSS distribution") +
    guides(fill = guide_legend(ncol = 2, override.aes = list(alpha = 1)))

  combined <- (p1a | p1b) / (p1c | p1d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_01_cohort", h_mm = 140)
  message("  Figure 1 done.")
}

# ==============================================================================
# FIGURE 2: AlphaFold3 Structural Quality
# ==============================================================================
figure_2 <- function() {
  message("Generating Figure 2: AlphaFold3 Structural Quality...")

  # --- Panel a: LDLR pLDDT ribbon ---
  p2a <- ggplot(plddt_ldlr, aes(x = residue_num, y = plddt, color = domain)) +
    geom_line(linewidth = 0.3) +
    geom_hline(yintercept = c(50, 70, 90), linetype = "dashed", linewidth = 0.2, color = "grey50") +
    annotate("text", x = max(plddt_ldlr$residue_num) * 0.98,
             y = c(50, 70, 90) + 2, label = c("50", "70", "90"),
             size = 1.8, hjust = 1, color = "grey40") +
    scale_color_manual(values = domain_pal, name = "Domain") +
    labs(x = "Residue number", y = "pLDDT", title = "LDLR per-residue confidence") +
    scale_x_continuous(breaks = seq(0, 860, 100)) +
    guides(color = guide_legend(ncol = 2, override.aes = list(linewidth = 1)))

  # --- Panel b: PCSK9 pLDDT ribbon ---
  p2b <- ggplot(plddt_pcsk9, aes(x = residue_num, y = plddt, color = domain)) +
    geom_line(linewidth = 0.3) +
    geom_hline(yintercept = c(50, 70, 90), linetype = "dashed", linewidth = 0.2, color = "grey50") +
    annotate("text", x = max(plddt_pcsk9$residue_num) * 0.98,
             y = c(50, 70, 90) + 2, label = c("50", "70", "90"),
             size = 1.8, hjust = 1, color = "grey40") +
    scale_color_manual(values = domain_pal, name = "Domain") +
    labs(x = "Residue number", y = "pLDDT", title = "PCSK9 per-residue confidence") +
    guides(color = guide_legend(ncol = 1, override.aes = list(linewidth = 1)))

  # --- Panel c: Domain-averaged pLDDT barplot ---
  all_plddt <- bind_rows(
    plddt_ldlr %>% mutate(protein = "LDLR"),
    plddt_pcsk9 %>% mutate(protein = "PCSK9"),
    plddt_apob %>% mutate(protein = "ApoB-RBD")
  )
  dom_avg <- all_plddt %>%
    group_by(protein, domain) %>%
    summarise(mean_plddt = mean(plddt, na.rm = TRUE),
              sd_plddt = sd(plddt, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    mutate(se = sd_plddt / sqrt(n))

  p2c <- ggplot(dom_avg, aes(x = reorder(domain, mean_plddt), y = mean_plddt, fill = protein)) +
    geom_col(width = 0.65, color = "grey30", linewidth = 0.2) +
    geom_errorbar(aes(ymin = pmax(mean_plddt - se, 0), ymax = pmin(mean_plddt + se, 100)),
                  width = 0.2, linewidth = 0.25) +
    scale_fill_manual(values = c("LDLR" = "#2166AC", "PCSK9" = "#35978F", "ApoB-RBD" = "#B2182B"),
                      name = "Protein") +
    coord_flip() +
    labs(x = NULL, y = "Mean pLDDT", title = "Domain-averaged confidence")

  # --- Panel d: pLDDT vs LDL-C scatter ---
  merged_d <- sss_scores %>%
    filter(!is.na(plddt)) %>%
    inner_join(comp_sss %>% select(mutation, last_ldl, ldl1) %>% distinct(),
               by = c("variant_id" = "mutation")) %>%
    mutate(ldl_use = ifelse(!is.na(last_ldl), last_ldl, ldl1)) %>%
    filter(!is.na(ldl_use))

  if (nrow(merged_d) >= 5) {
    ct <- cor.test(merged_d$plddt, merged_d$ldl_use, method = "pearson")
    lbl <- sprintf("r = %.2f, p = %.3f", ct$estimate, ct$p.value)
    p2d <- ggplot(merged_d, aes(x = plddt, y = ldl_use)) +
      geom_point(alpha = 0.4, size = 0.8, color = "#2166AC") +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "#B2182B", fill = "#FDDBC7") +
      annotate("text", x = min(merged_d$plddt, na.rm = TRUE),
               y = max(merged_d$ldl_use, na.rm = TRUE), label = lbl,
               hjust = 0, vjust = 1, size = 2.2, fontface = "italic") +
      labs(x = "pLDDT", y = "LDL-C (mmol/L)", title = "Structural confidence vs LDL-C")
  } else {
    p2d <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 3) +
      theme_void() + labs(title = "pLDDT vs LDL-C")
  }

  combined <- (p2a | p2b) / (p2c | p2d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_02_structural_quality", h_mm = 150)
  message("  Figure 2 done.")
}

# ==============================================================================
# FIGURE 3: FoldX Thermodynamic Stability
# ==============================================================================
figure_3 <- function() {
  message("Generating Figure 3: FoldX Thermodynamic Stability...")

  ddg_col <- find_col(foldx, c("ddG_kcal_mol", "ddG", "ddG_kcal"))
  if (is.na(ddg_col)) { warning("No ddG column in FoldX"); return(invisible(NULL)) }

  foldx$ddG_val <- as.numeric(foldx[[ddg_col]])
  foldx$ddg_class <- factor(classify_ddg(foldx$ddG_val),
                            levels = c("Stabilising", "Neutral", "Destabilising", "Highly destabilising"))

  # --- Panel a: ddG by domain violin ---
  p3a <- ggplot(foldx, aes(x = domain, y = ddG_val, fill = domain)) +
    geom_violin(alpha = 0.5, linewidth = 0.2, scale = "width") +
    geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
    geom_hline(yintercept = c(-1, 1, 4), linetype = "dashed", linewidth = 0.2, color = "grey50") +
    scale_fill_manual(values = domain_pal, guide = "none") +
    labs(x = NULL, y = expression(Delta*Delta*"G (kcal/mol)"), title = "FoldX stability by domain") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- Panel b: ddG vs Last LDL ---
  comp_ddg <- comp_sss %>% filter(!is.na(ddG), !is.na(last_ldl))
  if (nrow(comp_ddg) < 5) comp_ddg <- comp_sss %>% filter(!is.na(ddG), !is.na(ldl1)) %>%
    mutate(last_ldl = ldl1)

  if (nrow(comp_ddg) >= 5) {
    ct_b <- cor.test(comp_ddg$ddG, comp_ddg$last_ldl, method = "pearson")
    lbl_b <- sprintf("r = %.2f, p = %.3f\nn = %d", ct_b$estimate, ct_b$p.value, nrow(comp_ddg))
  } else {
    lbl_b <- "Insufficient data"
  }

  p3b <- ggplot(comp_ddg, aes(x = ddG, y = last_ldl)) +
    geom_point(alpha = 0.4, size = 0.8, color = "#2166AC") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "#B2182B", fill = "#FDDBC7") +
    annotate("text", x = min(comp_ddg$ddG, na.rm = TRUE),
             y = max(comp_ddg$last_ldl, na.rm = TRUE), label = lbl_b,
             hjust = 0, vjust = 1, size = 2, fontface = "italic") +
    labs(x = expression(Delta*Delta*"G (kcal/mol)"), y = "Last LDL-C (mmol/L)",
         title = expression(Delta*Delta*"G vs LDL-C"))

  # --- Panel c: ddG vs ApoB ---
  comp_apob <- comp_sss %>% filter(!is.na(ddG), !is.na(apob))
  if (nrow(comp_apob) >= 5) {
    ct_c <- cor.test(comp_apob$ddG, comp_apob$apob, method = "pearson")
    lbl_c <- sprintf("r = %.2f, p = %.3f\nn = %d", ct_c$estimate, ct_c$p.value, nrow(comp_apob))
  } else {
    lbl_c <- "Insufficient data"
  }

  p3c <- ggplot(comp_apob, aes(x = ddG, y = apob)) +
    geom_point(alpha = 0.4, size = 0.8, color = "#35978F") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "#B2182B", fill = "#FDDBC7") +
    annotate("text", x = min(comp_apob$ddG, na.rm = TRUE),
             y = max(comp_apob$apob, na.rm = TRUE), label = lbl_c,
             hjust = 0, vjust = 1, size = 2, fontface = "italic") +
    labs(x = expression(Delta*Delta*"G (kcal/mol)"), y = "ApoB (g/L)",
         title = expression(Delta*Delta*"G vs ApoB"))

  # --- Panel d: ddG classification stacked bar ---
  class_tab <- foldx %>% count(ddg_class) %>% mutate(pct = n / sum(n))
  p3d <- ggplot(class_tab, aes(x = "FoldX variants", y = n, fill = ddg_class)) +
    geom_col(width = 0.5, color = "grey30", linewidth = 0.2) +
    geom_text(aes(label = paste0(n, " (", round(pct * 100), "%)")),
              position = position_stack(vjust = 0.5), size = 2) +
    scale_fill_manual(values = ddg_class_pal, name = "Classification") +
    labs(x = NULL, y = "Count", title = "Thermodynamic classification") +
    coord_flip()

  combined <- (p3a | p3b) / (p3c | p3d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_03_foldx_stability", h_mm = 140)
  message("  Figure 3 done.")
}

# ==============================================================================
# FIGURE 4: SSS Construction & Validation
# ==============================================================================
figure_4 <- function() {
  message("Generating Figure 4: SSS Construction & Validation...")

  # --- Panel a: SSS components conceptual diagram ---
  comps <- data.frame(
    component = c("Domain severity", "ddG (FoldX)", "pLDDT", "Interface distance", "Clinical preset"),
    weight = c(0.30, 0.25, 0.20, 0.15, 0.10),
    stringsAsFactors = FALSE
  )
  comps$component <- factor(comps$component, levels = rev(comps$component))

  p4a <- ggplot(comps, aes(x = component, y = weight, fill = component)) +
    geom_col(width = 0.6, color = "grey30", linewidth = 0.2, show.legend = FALSE) +
    geom_text(aes(label = paste0(weight * 100, "%")), hjust = -0.2, size = 2.2) +
    scale_fill_brewer(palette = "Set2") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Weight", title = "SSS component weights")

  # --- Panel b: SSS by gene boxplot ---
  p4b <- ggplot(comp_sss, aes(x = gene, y = sss, fill = gene)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
    scale_fill_manual(values = gene_pal, guide = "none") +
    labs(x = NULL, y = "Structural Severity Score", title = "SSS by gene") +
    stat_summary(fun = median, geom = "text",
                 aes(label = sprintf("%.2f", after_stat(y))),
                 vjust = -0.8, size = 2, color = "black")

  # --- Panel c: SSS tertile -> ASCVD OR forest plot ---
  df4c <- add_sss_tertile(comp_sss)
  df4c$sss_tertile <- relevel(df4c$sss_tertile, ref = "Low")

  fit_c <- tryCatch(glm(ascvd ~ sss_tertile + age + sex, data = df4c, family = "binomial"),
                    error = function(e) NULL)
  if (!is.null(fit_c)) {
    or_tab <- broom::tidy(fit_c, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(grepl("sss_tertile", term)) %>%
      mutate(term = gsub("sss_tertile", "", term))

    p4c <- ggplot(or_tab, aes(x = estimate, y = term)) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "grey50") +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth = 0.4) +
      geom_point(size = 2, color = "#B2182B") +
      geom_text(aes(label = sprintf("OR %.2f (%.2f-%.2f)", estimate, conf.low, conf.high)),
                hjust = -0.1, size = 2, vjust = -0.8) +
      labs(x = "Odds Ratio (95% CI)", y = "SSS tertile (vs Low)",
           title = "ASCVD risk by SSS tertile") +
      scale_x_continuous(trans = "log2")
  } else {
    p4c <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Model failed", size = 3) +
      theme_void()
  }

  # --- Panel d: ROC curves ---
  roc_df <- comp_sss %>% filter(!is.na(ascvd), !is.na(sss), !is.na(age), !is.na(sex))
  roc_base <- tryCatch(pROC::roc(ascvd ~ age + sex, data = roc_df, quiet = TRUE),
                       error = function(e) {
                         pred_base <- predict(glm(ascvd ~ age + sex, data = roc_df, family = "binomial"), type = "response")
                         pROC::roc(roc_df$ascvd, pred_base, quiet = TRUE)
                       })
  pred_sss <- predict(glm(ascvd ~ age + sex + sss, data = roc_df, family = "binomial"), type = "response")
  roc_sss <- pROC::roc(roc_df$ascvd, pred_sss, quiet = TRUE)

  roc_plot_df <- bind_rows(
    data.frame(sens = roc_base$sensitivities, spec = 1 - roc_base$specificities,
               model = paste0("Base (AUC=", round(pROC::auc(roc_base), 3), ")")),
    data.frame(sens = roc_sss$sensitivities, spec = 1 - roc_sss$specificities,
               model = paste0("+SSS (AUC=", round(pROC::auc(roc_sss), 3), ")"))
  )

  p4d <- ggplot(roc_plot_df, aes(x = spec, y = sens, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey60") +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("#BDBDBD", "#B2182B"), name = NULL) +
    labs(x = "1 - Specificity", y = "Sensitivity", title = "ROC: ASCVD discrimination") +
    coord_equal() +
    theme(legend.position = c(0.7, 0.25), legend.background = element_blank())

  # --- Panel e: 5 SSS strategies from LOCO-CV ---
  sss_cols <- grep("^sss_", colnames(loco_cv), value = TRUE)
  if (length(sss_cols) >= 2 && "ascvd" %in% colnames(loco_cv)) {
    auc_list <- lapply(sss_cols, function(sc) {
      dd <- loco_cv[!is.na(loco_cv[[sc]]) & !is.na(loco_cv$ascvd), ]
      if (nrow(dd) < 20 || length(unique(dd$ascvd)) < 2) return(NULL)
      r <- pROC::roc(dd$ascvd, dd[[sc]], quiet = TRUE)
      data.frame(strategy = gsub("^sss_", "", sc), auc = as.numeric(pROC::auc(r)))
    })
    auc_df <- bind_rows(auc_list)
    p4e <- ggplot(auc_df, aes(x = reorder(strategy, auc), y = auc)) +
      geom_col(fill = "#2166AC", width = 0.6, alpha = 0.8, color = "grey30", linewidth = 0.2) +
      geom_text(aes(label = sprintf("%.3f", auc)), hjust = -0.15, size = 2) +
      coord_flip(ylim = c(0.45, max(auc_df$auc, na.rm = TRUE) * 1.08)) +
      labs(x = NULL, y = "AUC", title = "LOCO-CV SSS strategies")
  } else {
    p4e <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "LOCO-CV data unavailable", size = 3) +
      theme_void()
  }

  combined <- (p4a | p4b | p4c) / (p4d | p4e) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_04_sss_validation", h_mm = 150)
  message("  Figure 4 done.")
}

# ==============================================================================
# FIGURE 5: Domain-Specific Clinical Phenotypes
# ==============================================================================
figure_5 <- function() {
  message("Generating Figure 5: Domain-Specific Clinical Phenotypes...")

  top_domains <- comp_sss %>% count(domain) %>% filter(n >= 5) %>% pull(domain)
  df5 <- comp_sss %>% filter(domain %in% top_domains)

  # --- Panel a: Domain phenotype radar ---
  dom_summary <- df5 %>%
    group_by(domain) %>%
    summarise(
      LDL = mean(ldl1, na.rm = TRUE),
      ApoB = mean(apob, na.rm = TRUE),
      Discordance = mean(apob / ldl1, na.rm = TRUE),
      Lp_a = mean(lpa, na.rm = TRUE),
      ASCVD = mean(ascvd, na.rm = TRUE) * 100,
      Xanthomata = mean(xanthomata, na.rm = TRUE) * 100,
      .groups = "drop"
    )

  # Normalize 0-1
  radar_long <- dom_summary %>%
    pivot_longer(-domain, names_to = "metric", values_to = "value") %>%
    group_by(metric) %>%
    mutate(norm = (value - min(value, na.rm = TRUE)) /
             (max(value, na.rm = TRUE) - min(value, na.rm = TRUE) + 1e-10)) %>%
    ungroup()

  p5a <- ggplot(radar_long, aes(x = metric, y = norm, group = domain, color = domain)) +
    geom_polygon(fill = NA, linewidth = 0.4, alpha = 0.3) +
    geom_point(size = 0.8) +
    coord_polar() +
    scale_color_manual(values = domain_pal, name = "Domain") +
    labs(title = "Domain phenotype profile", x = NULL, y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    guides(color = guide_legend(ncol = 2))

  # --- Panel b: Domain ASCVD rates with Wilson CI ---
  ascvd_rates <- df5 %>%
    group_by(domain) %>%
    summarise(n = n(), events = sum(ascvd, na.rm = TRUE), .groups = "drop") %>%
    filter(n >= 3) %>%
    rowwise() %>%
    mutate(wilson_ci(events, n)) %>%
    ungroup()

  p5b <- ggplot(ascvd_rates, aes(x = reorder(domain, est), y = est * 100)) +
    geom_col(aes(fill = domain), width = 0.6, alpha = 0.8, color = "grey30", linewidth = 0.2,
             show.legend = FALSE) +
    geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), width = 0.2, linewidth = 0.3) +
    geom_text(aes(label = paste0(events, "/", n)), vjust = -0.5, size = 1.8) +
    scale_fill_manual(values = domain_pal) +
    coord_flip() +
    labs(x = NULL, y = "ASCVD prevalence (%)", title = "ASCVD by domain")

  # --- Panel c: Domain treatment response ---
  llt_dom <- llt_response %>%
    filter(!is.na(domain), !is.na(ldl_pct_reduction)) %>%
    group_by(domain) %>%
    filter(n() >= 3) %>%
    summarise(
      mean_red = mean(ldl_pct_reduction, na.rm = TRUE),
      se = sd(ldl_pct_reduction, na.rm = TRUE) / sqrt(n()),
      n = n(), .groups = "drop"
    )

  p5c <- ggplot(llt_dom, aes(x = reorder(domain, mean_red), y = mean_red)) +
    geom_col(aes(fill = domain), width = 0.6, alpha = 0.8, color = "grey30", linewidth = 0.2,
             show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_red - se, ymax = mean_red + se), width = 0.2, linewidth = 0.3) +
    scale_fill_manual(values = domain_pal) +
    coord_flip() +
    labs(x = NULL, y = "Mean LDL % reduction", title = "Treatment response by domain")

  # --- Panel d: LDLR domain schematic ---
  doms <- data.frame(
    xmin = c(1, 22, 353, 693, 750, 771),
    xmax = c(21, 352, 692, 749, 770, 860),
    ymin = rep(0, 6), ymax = rep(1, 6),
    label = c("SP", "Ligand-binding\n(R1-R7)", "EGF-precursor\nhomology", "O-linked", "TM", "Cyto"),
    fill = c("#FDB863", "#2166AC", "#B2182B", "#35978F", "#762A83", "#9970AB"),
    stringsAsFactors = FALSE
  )

  p5d <- ggplot() +
    geom_rect(data = doms, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = doms$fill, alpha = 0.7, color = "grey30", linewidth = 0.3) +
    geom_text(data = doms, aes(x = (xmin + xmax) / 2, y = 0.5, label = label),
              size = 1.8, lineheight = 0.8, color = "white", fontface = "bold") +
    geom_text(data = doms, aes(x = (xmin + xmax) / 2, y = -0.15,
                                label = paste0(xmin, "-", xmax)),
              size = 1.5, color = "grey30") +
    coord_cartesian(ylim = c(-0.3, 1.3)) +
    labs(title = "LDLR domain architecture", x = "Amino acid position", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

  combined <- (p5a | p5b) / (p5c | p5d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_05_domain_phenotypes", h_mm = 155)
  message("  Figure 5 done.")
}

# ==============================================================================
# FIGURE 6: Treatment Response & Pharmacogenomics
# ==============================================================================
figure_6 <- function() {
  message("Generating Figure 6: Treatment Response & Pharmacogenomics...")

  llt <- llt_response
  llt$ldl_red <- llt$ldl_pct_reduction
  llt <- add_sss_tertile(llt, "sss")

  # --- Panel a: LDL reduction by SSS tertile ---
  llt_a <- llt %>% filter(!is.na(ldl_red), !is.na(sss_tertile))
  p6a <- ggplot(llt_a, aes(x = sss_tertile, y = ldl_red, fill = sss_tertile)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7, width = 0.6) +
    scale_fill_manual(values = sss_tert_pal, guide = "none") +
    labs(x = "SSS tertile", y = "LDL-C % reduction", title = "Treatment response by SSS")

  # --- Panel b: Statin-specific response ---
  statin_col <- find_col(llt, c("statin_type", "statin", "drug"))
  if (!is.na(statin_col)) {
    llt_b <- llt %>%
      filter(!is.na(ldl_red), !is.na(sss_tertile), !is.na(.data[[statin_col]]),
             .data[[statin_col]] != "", .data[[statin_col]] != "none") %>%
      group_by(.data[[statin_col]]) %>% filter(n() >= 10) %>% ungroup()

    if (nrow(llt_b) > 0) {
      p6b <- ggplot(llt_b, aes(x = sss_tertile, y = ldl_red, fill = sss_tertile)) +
        geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7, width = 0.6) +
        facet_wrap(as.formula(paste("~", statin_col)), scales = "free_y") +
        scale_fill_manual(values = sss_tert_pal, guide = "none") +
        labs(x = "SSS tertile", y = "LDL-C % reduction", title = "Statin-specific response")
    } else {
      p6b <- ggplot() + annotate("text", x = .5, y = .5, label = "No statin subgroups with n>=10", size = 2.5) + theme_void()
    }
  } else {
    p6b <- ggplot() + annotate("text", x = .5, y = .5, label = "No statin type column", size = 2.5) + theme_void()
  }

  # --- Panel c: 13 treatment-resistant variants forest plot ---
  var_response <- llt %>%
    filter(!is.na(ldl_red)) %>%
    group_by(variant) %>%
    filter(n() >= 3) %>%
    summarise(mean_red = mean(ldl_red, na.rm = TRUE),
              se = sd(ldl_red, na.rm = TRUE) / sqrt(n()),
              n = n(), .groups = "drop") %>%
    arrange(mean_red) %>%
    head(13)

  p6c <- ggplot(var_response, aes(x = mean_red, y = reorder(variant, -mean_red))) +
    geom_errorbarh(aes(xmin = mean_red - 1.96 * se, xmax = mean_red + 1.96 * se),
                   height = 0.2, linewidth = 0.3) +
    geom_point(size = 1.5, color = "#B2182B") +
    geom_text(aes(label = paste0("n=", n)), hjust = -0.3, size = 1.8) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2) +
    labs(x = "Mean LDL-C % reduction", y = NULL, title = "Treatment-resistant variants")

  # --- Panel d: Null vs missense vs splicing ---
  if (!is.null(sss_scores) && "variant_type" %in% colnames(sss_scores)) {
    type_merge <- llt %>%
      inner_join(sss_scores %>% select(variant_id, variant_type) %>% distinct(),
                 by = c("variant" = "variant_id")) %>%
      filter(!is.na(ldl_red), !is.na(variant_type))

    if (nrow(type_merge) >= 10) {
      p6d <- ggplot(type_merge, aes(x = variant_type, y = ldl_red, fill = variant_type)) +
        geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7, width = 0.6) +
        scale_fill_brewer(palette = "Set2", guide = "none") +
        labs(x = "Variant type", y = "LDL-C % reduction", title = "Response by variant type")
    } else {
      p6d <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient merged data", size = 2.5) + theme_void()
    }
  } else {
    p6d <- ggplot() + annotate("text", x = .5, y = .5, label = "No variant_type column", size = 2.5) + theme_void()
  }

  # --- Panel e: PCSK9i user SSS distribution ---
  pcsk9i_col <- find_col(llt, c("pcsk9i"))
  if (!is.na(pcsk9i_col)) {
    llt$pcsk9i_use <- as.character(llt[[pcsk9i_col]]) %in% c("True", "TRUE", "1", "Yes")
    pcsk9i_dat <- llt %>% filter(!is.na(sss))

    p6e <- ggplot(pcsk9i_dat, aes(x = sss, fill = pcsk9i_use)) +
      geom_histogram(bins = 25, alpha = 0.6, position = "identity", color = "white", linewidth = 0.15) +
      scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#B2182B"),
                        labels = c("No PCSK9i", "PCSK9i"), name = NULL) +
      labs(x = "SSS", y = "Count", title = "PCSK9i use by SSS") +
      theme(legend.position = c(0.8, 0.8), legend.background = element_blank())
  } else {
    p6e <- ggplot() + annotate("text", x = .5, y = .5, label = "No pcsk9i column", size = 2.5) + theme_void()
  }

  combined <- (p6a | p6b) / (p6c | p6d | p6e) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_06_treatment_response", h_mm = 165)
  message("  Figure 6 done.")
}

# ==============================================================================
# FIGURE 7: ApoB/LDL-C Discordance
# ==============================================================================
figure_7 <- function() {
  message("Generating Figure 7: ApoB/LDL-C Discordance...")

  df7 <- comp_sss %>% filter(!is.na(apob), !is.na(ldl1))
  med_apob <- median(df7$apob, na.rm = TRUE)
  med_ldl  <- median(df7$ldl1, na.rm = TRUE)
  df7$concordance <- ifelse(
    (df7$apob >= med_apob & df7$ldl1 >= med_ldl) |
    (df7$apob < med_apob & df7$ldl1 < med_ldl),
    "Concordant", "Discordant"
  )

  # --- Panel a: ApoB vs LDL-C scatter ---
  p7a <- ggplot(df7, aes(x = ldl1, y = apob, color = concordance)) +
    geom_hline(yintercept = med_apob, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_vline(xintercept = med_ldl, linetype = "dashed", linewidth = 0.2, color = "grey50") +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Concordant" = "#2166AC", "Discordant" = "#B2182B"), name = NULL) +
    labs(x = "LDL-C (mmol/L)", y = "ApoB (g/L)", title = "ApoB/LDL-C concordance") +
    annotate("text", x = max(df7$ldl1, na.rm = TRUE) * 0.95,
             y = max(df7$apob, na.rm = TRUE) * 0.95,
             label = sprintf("n = %d\nDiscordant: %d (%.0f%%)",
                             nrow(df7), sum(df7$concordance == "Discordant"),
                             100 * mean(df7$concordance == "Discordant")),
             hjust = 1, vjust = 1, size = 2) +
    theme(legend.position = c(0.15, 0.9), legend.background = element_blank())

  # --- Panel b: Discordance ratio by domain ---
  df7$disc_ratio <- df7$apob / df7$ldl1
  top_dom <- df7 %>% count(domain) %>% filter(n >= 5) %>% pull(domain)
  df7b <- df7 %>% filter(domain %in% top_dom)

  p7b <- ggplot(df7b, aes(x = reorder(domain, disc_ratio, FUN = median), y = disc_ratio, fill = domain)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7, width = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = domain_pal) +
    coord_flip() +
    labs(x = NULL, y = "ApoB/LDL-C ratio", title = "Discordance by domain")

  # --- Panel c: Discordance tertile -> ASCVD forest plot ---
  df7c <- df7 %>% filter(!is.na(ascvd))
  brks <- quantile(df7c$disc_ratio, c(0, 1/3, 2/3, 1), na.rm = TRUE)
  if (brks[2] == brks[3]) brks[3] <- brks[3] + 1e-6
  df7c$disc_tert <- cut(df7c$disc_ratio, breaks = brks, include.lowest = TRUE,
                        labels = c("Low", "Medium", "High"))
  df7c$disc_tert <- relevel(df7c$disc_tert, ref = "Low")

  fit7 <- tryCatch(glm(ascvd ~ disc_tert + age + sex, data = df7c, family = "binomial"),
                   error = function(e) NULL)
  if (!is.null(fit7)) {
    or7 <- broom::tidy(fit7, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(grepl("disc_tert", term)) %>%
      mutate(term = gsub("disc_tert", "", term))

    p7c <- ggplot(or7, aes(x = estimate, y = term)) +
      geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "grey50") +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth = 0.4) +
      geom_point(size = 2, color = "#B2182B") +
      geom_text(aes(label = sprintf("OR %.2f", estimate)), hjust = -0.2, vjust = -0.8, size = 2) +
      labs(x = "Odds Ratio (95% CI)", y = "Discordance tertile", title = "Discordance -> ASCVD") +
      scale_x_continuous(trans = "log2")
  } else {
    p7c <- ggplot() + annotate("text", x = .5, y = .5, label = "Model failed", size = 3) + theme_void()
  }

  # --- Panel d: Top 10 discordant variants from nobel ---
  disc_col <- find_col(nobel, c("mean_discordance"))
  if (!is.na(disc_col)) {
    top10 <- nobel %>% arrange(desc(.data[[disc_col]])) %>% head(10)
    p7d <- ggplot(top10, aes(x = .data[[disc_col]], y = reorder(variant, .data[[disc_col]]))) +
      geom_col(fill = "#B2182B", alpha = 0.8, width = 0.6, color = "grey30", linewidth = 0.2) +
      geom_text(aes(label = sprintf("%.2f", .data[[disc_col]])), hjust = -0.15, size = 2) +
      labs(x = "Mean discordance", y = NULL, title = "Top discordant variants") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.2)))
  } else {
    p7d <- ggplot() + annotate("text", x = .5, y = .5, label = "No discordance column", size = 3) + theme_void()
  }

  combined <- (p7a | p7b) / (p7c | p7d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_07_discordance", h_mm = 145)
  message("  Figure 7 done.")
}

# ==============================================================================
# FIGURE 8: Lp(a) x SSS Interaction
# ==============================================================================
figure_8 <- function() {
  message("Generating Figure 8: Lp(a) x SSS Interaction...")

  df8 <- comp_sss %>% filter(!is.na(lpa))

  # --- Panel a: Lp(a) distribution ---
  p8a <- ggplot(df8, aes(x = lpa)) +
    geom_histogram(bins = 40, fill = "#2166AC", color = "white", linewidth = 0.15, alpha = 0.8) +
    geom_vline(xintercept = 143, linetype = "dashed", linewidth = 0.4, color = "#B2182B") +
    annotate("text", x = 145, y = Inf, label = "143 nmol/L\nthreshold", hjust = -0.1,
             vjust = 1.5, size = 2, color = "#B2182B") +
    labs(x = "Lp(a) (nmol/L)", y = "Count", title = "Lp(a) distribution")

  # --- Panel b: Lp(a) by SSS tertile ---
  df8t <- add_sss_tertile(df8)
  p8b <- ggplot(df8t, aes(x = sss_tertile, y = lpa, fill = sss_tertile)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7, width = 0.6) +
    scale_fill_manual(values = sss_tert_pal, guide = "none") +
    labs(x = "SSS tertile", y = "Lp(a) (nmol/L)", title = "Lp(a) by SSS tertile")

  # --- Panel c: 2x2 heatmap high/low SSS x high/low Lp(a) -> ASCVD% ---
  df8c <- comp_sss %>% filter(!is.na(lpa), !is.na(sss), !is.na(ascvd))
  df8c$sss_cat <- ifelse(df8c$sss >= median(df8c$sss, na.rm = TRUE), "High SSS", "Low SSS")
  df8c$lpa_cat <- ifelse(df8c$lpa >= 143, "High Lp(a)", "Low Lp(a)")

  heat_tab <- df8c %>%
    group_by(sss_cat, lpa_cat) %>%
    summarise(ascvd_pct = round(100 * mean(ascvd, na.rm = TRUE), 1),
              n = n(), .groups = "drop")

  p8c <- ggplot(heat_tab, aes(x = lpa_cat, y = sss_cat, fill = ascvd_pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(ascvd_pct, "%\n(n=", n, ")")), size = 2.5, fontface = "bold") +
    scale_fill_gradient(low = "#FDDBC7", high = "#B2182B", name = "ASCVD %") +
    labs(x = NULL, y = NULL, title = "SSS x Lp(a) interaction") +
    theme(axis.text = element_text(size = 6, face = "bold"))

  # --- Panel d: Top Lp(a) variants from nobel ---
  lpa_col <- find_col(nobel, c("mean_lpa"))
  if (!is.na(lpa_col)) {
    top_lpa <- nobel %>% filter(!is.na(.data[[lpa_col]])) %>%
      arrange(desc(.data[[lpa_col]])) %>% head(10)
    p8d <- ggplot(top_lpa, aes(x = .data[[lpa_col]], y = reorder(variant, .data[[lpa_col]]))) +
      geom_col(fill = "#35978F", alpha = 0.8, width = 0.6, color = "grey30", linewidth = 0.2) +
      geom_text(aes(label = round(.data[[lpa_col]], 0)), hjust = -0.2, size = 2) +
      labs(x = "Mean Lp(a) (nmol/L)", y = NULL, title = "Top Lp(a) variants") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.2)))
  } else {
    p8d <- ggplot() + annotate("text", x = .5, y = .5, label = "No Lp(a) in nobel", size = 3) + theme_void()
  }

  combined <- (p8a | p8b) / (p8c | p8d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_08_lpa_interaction", h_mm = 140)
  message("  Figure 8 done.")
}

# ==============================================================================
# FIGURE 9: NMR Metabolomic Cascade (STREAMING)
# ==============================================================================
figure_9 <- function() {
  message("Generating Figure 9: NMR Metabolomic Cascade (streaming)...")

  # Build file list
  nmr_files <- c()
  for (i in 1:5) {
    for (j in letters[1:5]) {
      f <- file.path(NMR_DIR, paste0("calon_batch_nmr", i, j, ".csv"))
      if (file.exists(f)) nmr_files <- c(nmr_files, f)
    }
  }
  if (length(nmr_files) == 0) {
    warning("No NMR batch files found"); return(invisible(NULL))
  }
  message("  Found ", length(nmr_files), " NMR batch files")

  # Map field codes to biomarker names
  nmr_map <- c(
    "participant.p23400_i0" = "Total_cholesterol",
    "participant.p23405_i0" = "LDL_cholesterol",
    "participant.p23406_i0" = "LDL_ester_chol",
    "participant.p23407_i0" = "HDL_cholesterol",
    "participant.p23444_i0" = "Remnant_cholesterol",
    "participant.p23474_i0" = "GlycA",
    "participant.p23449_i0" = "ApoB",
    "participant.p23450_i0" = "ApoA1",
    "participant.p23454_i0" = "Total_FA",
    "participant.p23455_i0" = "SFA",
    "participant.p23457_i0" = "MUFA",
    "participant.p23458_i0" = "PUFA"
  )

  # Determine available columns from first file header
  hdr <- colnames(read.csv(nmr_files[1], nrows = 1, check.names = FALSE))
  avail_codes <- intersect(names(nmr_map), hdr)
  if (length(avail_codes) < 2) {
    # Try with make.names
    hdr_safe <- make.names(hdr)
    avail_codes_safe <- intersect(make.names(names(nmr_map)), hdr_safe)
    if (length(avail_codes_safe) >= 2) {
      names(nmr_map) <- make.names(names(nmr_map))
      avail_codes <- avail_codes_safe
      message("  Using make.names() column matching for NMR")
    } else {
      warning("  Fewer than 2 NMR columns matched. Trying partial match...")
      avail_codes <- hdr[grep("p234", hdr)]
      if (length(avail_codes) < 2) {
        warning("  Cannot find NMR columns"); return(invisible(NULL))
      }
      # Create dynamic mapping
      nmr_map <- setNames(paste0("NMR_", seq_along(avail_codes)), avail_codes)
    }
  }

  ldl_col_nmr <- find_col(data.frame(setNames(rep(1, length(hdr)), make.names(hdr))),
                           c("participant.p23405_i0", "participant.p23405.i0"))
  if (is.na(ldl_col_nmr)) {
    ldl_col_nmr <- grep("p23405", make.names(hdr), value = TRUE)[1]
  }
  if (is.na(ldl_col_nmr)) {
    warning("  Cannot find LDL NMR column"); return(invisible(NULL))
  }

  cols_to_read <- unique(c("participant.eid", ldl_col_nmr, avail_codes))
  cols_to_read_safe <- make.names(cols_to_read)

  # --- STREAMING: accumulate running sums per quintile ---
  message("  Streaming through NMR files...")
  all_ldl_vals <- c()

  # Pass 1: collect all LDL values for quintile breaks
  for (f in nmr_files) {
    chunk <- tryCatch({
      d <- data.table::fread(f, select = make.names(cols_to_read), showProgress = FALSE)
      as.numeric(d[[make.names(ldl_col_nmr)]])
    }, error = function(e) numeric(0))
    chunk <- chunk[!is.na(chunk)]
    all_ldl_vals <- c(all_ldl_vals, chunk)
    rm(chunk); gc(verbose = FALSE)
  }

  if (length(all_ldl_vals) < 100) {
    warning("  Too few LDL values from NMR"); return(invisible(NULL))
  }

  q_breaks <- quantile(all_ldl_vals, probs = seq(0, 1, 0.2), na.rm = TRUE)
  rm(all_ldl_vals); gc(verbose = FALSE)

  bio_cols_safe <- make.names(avail_codes)
  bio_cols_safe <- setdiff(bio_cols_safe, make.names(ldl_col_nmr))

  # Pass 2: accumulate sums and counts per tier
  tier_sums <- list()
  tier_counts <- list()
  for (bc in bio_cols_safe) {
    tier_sums[[bc]]   <- rep(0, 5)
    tier_counts[[bc]] <- rep(0, 5)
  }

  for (f in nmr_files) {
    tryCatch({
      d <- data.table::fread(f, select = unique(c(make.names(ldl_col_nmr), bio_cols_safe)),
                             showProgress = FALSE)
      ldl_v <- as.numeric(d[[make.names(ldl_col_nmr)]])
      tier <- cut(ldl_v, breaks = q_breaks, include.lowest = TRUE, labels = FALSE)

      for (bc in bio_cols_safe) {
        vals <- as.numeric(d[[bc]])
        for (t in 1:5) {
          idx <- which(tier == t & !is.na(vals))
          tier_sums[[bc]][t]   <- tier_sums[[bc]][t] + sum(vals[idx])
          tier_counts[[bc]][t] <- tier_counts[[bc]][t] + length(idx)
        }
      }
      rm(d); gc(verbose = FALSE)
    }, error = function(e) message("    Warning reading ", basename(f), ": ", e$message))
  }

  # Compute means
  tier_means <- data.frame(tier = 1:5)
  for (bc in bio_cols_safe) {
    tier_means[[bc]] <- ifelse(tier_counts[[bc]] > 0,
                               tier_sums[[bc]] / tier_counts[[bc]], NA)
  }

  # Compute fold-change vs tier 1
  tier_fc <- tier_means
  for (bc in bio_cols_safe) {
    ref <- tier_fc[[bc]][1]
    if (!is.na(ref) && ref != 0) {
      tier_fc[[bc]] <- tier_fc[[bc]] / ref
    }
  }

  # Map back to biomarker names
  rename_vec <- nmr_map[avail_codes]
  names(rename_vec) <- make.names(avail_codes)
  rename_vec <- rename_vec[names(rename_vec) %in% bio_cols_safe]

  # --- Panel a: Heatmap ---
  fc_long <- tier_fc %>%
    select(tier, all_of(names(rename_vec))) %>%
    pivot_longer(-tier, names_to = "code", values_to = "fc") %>%
    mutate(biomarker = rename_vec[code])

  p9a <- ggplot(fc_long, aes(x = factor(tier), y = biomarker, fill = fc)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", fc)), size = 1.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 1,
                         name = "Fold\nchange") +
    labs(x = "LDL quintile", y = NULL, title = "NMR biomarker cascade")

  # --- Panel b: LDL subclass stacked bar (use available LDL-related cols) ---
  ldl_related <- bio_cols_safe[grep("p23405|p23406|p23407|p23400", bio_cols_safe)]
  if (length(ldl_related) >= 2) {
    ldl_sub <- tier_means %>%
      select(tier, all_of(ldl_related)) %>%
      pivot_longer(-tier, names_to = "code", values_to = "value") %>%
      mutate(biomarker = ifelse(code %in% names(rename_vec), rename_vec[code], code))

    p9b <- ggplot(ldl_sub, aes(x = factor(tier), y = value, fill = biomarker)) +
      geom_col(position = "stack", color = "grey30", linewidth = 0.15) +
      scale_fill_brewer(palette = "Blues", name = NULL) +
      labs(x = "LDL quintile", y = "Concentration", title = "Lipid subclass composition")
  } else {
    p9b <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient LDL subclass data", size = 2.5) + theme_void()
  }

  # --- Panel c: GlycA gradient ---
  glyca_col <- bio_cols_safe[grep("p23474", bio_cols_safe)]
  if (length(glyca_col) >= 1) {
    glyca_col <- glyca_col[1]
    glyca_df <- data.frame(tier = 1:5, glyca = tier_means[[glyca_col]])

    p9c <- ggplot(glyca_df, aes(x = tier, y = glyca)) +
      geom_line(linewidth = 0.6, color = "#B2182B") +
      geom_point(size = 2, color = "#B2182B") +
      labs(x = "LDL quintile", y = "GlycA (mmol/L)", title = "Inflammation gradient")
  } else {
    p9c <- ggplot() + annotate("text", x = .5, y = .5, label = "No GlycA column", size = 2.5) + theme_void()
  }

  # --- Panel d: FA composition grouped bar ---
  fa_cols <- bio_cols_safe[grep("p23454|p23455|p23457|p23458", bio_cols_safe)]
  if (length(fa_cols) >= 2) {
    fa_df <- tier_means %>%
      select(tier, all_of(fa_cols)) %>%
      pivot_longer(-tier, names_to = "code", values_to = "value") %>%
      mutate(biomarker = ifelse(code %in% names(rename_vec), rename_vec[code], code))

    p9d <- ggplot(fa_df, aes(x = factor(tier), y = value, fill = biomarker)) +
      geom_col(position = "dodge", color = "grey30", linewidth = 0.15, width = 0.7) +
      scale_fill_brewer(palette = "Set2", name = NULL) +
      labs(x = "LDL quintile", y = "Concentration", title = "Fatty acid composition")
  } else {
    p9d <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient FA data", size = 2.5) + theme_void()
  }

  combined <- (p9a | p9b) / (p9c | p9d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_09_nmr_cascade", h_mm = 150)
  message("  Figure 9 done.")
}

# ==============================================================================
# FIGURE 10: Cardiac MRI Imaging (STREAMING)
# ==============================================================================
figure_10 <- function() {
  message("Generating Figure 10: Cardiac MRI Imaging (streaming)...")

  mri_files <- list(
    lv     = file.path(MRI_DIR, "calon_extra_mri_lv (1).csv"),
    la_rv  = file.path(MRI_DIR, "calon_extra_mri_la_rv.csv"),
    aorta  = file.path(MRI_DIR, "calon_extra_mri_aorta (1).csv"),
    strain = file.path(MRI_DIR, "calon_extra_mri_strain.csv"),
    cacs   = file.path(MRI_DIR, "calon_extra_mri_cat162.csv")
  )

  # Helper: stream a MRI file, compute quintile means for given cols
  stream_mri <- function(filepath, target_cols) {
    if (!file.exists(filepath)) { warning("  MRI file not found: ", filepath); return(NULL) }
    hdr <- make.names(colnames(read.csv(filepath, nrows = 1, check.names = FALSE)))
    avail <- intersect(make.names(target_cols), hdr)
    if (length(avail) == 0) {
      # Try partial match
      for (tc in target_cols) {
        code <- gsub("participant\\.p", "", tc)
        code <- gsub("_i.*", "", code)
        m <- grep(code, hdr, value = TRUE)
        if (length(m) > 0) avail <- c(avail, m[1])
      }
    }
    if (length(avail) == 0) { warning("  No target cols found in ", basename(filepath)); return(NULL) }

    d <- data.table::fread(filepath, select = avail, showProgress = FALSE)
    result <- list()
    for (col in avail) {
      vals <- as.numeric(d[[col]])
      vals <- vals[!is.na(vals)]
      if (length(vals) < 50) next
      brks <- quantile(vals, probs = seq(0, 1, 0.2), na.rm = TRUE)
      if (any(duplicated(brks))) {
        brks <- unique(brks)
        if (length(brks) < 3) next
      }
      tier <- cut(as.numeric(d[[col]]), breaks = brks, include.lowest = TRUE, labels = FALSE)
      tier_df <- data.frame(tier = tier, value = as.numeric(d[[col]])) %>%
        filter(!is.na(tier), !is.na(value))
      result[[col]] <- tier_df
    }
    rm(d); gc(verbose = FALSE)
    return(result)
  }

  # --- Panel a: Aortic distensibility ---
  aorta_targets <- c("participant.p24118_i2", "participant.p24119_i2")
  aorta_data <- stream_mri(mri_files$aorta, aorta_targets)
  if (!is.null(aorta_data) && length(aorta_data) > 0) {
    ad <- aorta_data[[1]]
    ad$tier <- factor(ad$tier)
    p10a <- ggplot(ad, aes(x = tier, y = value, fill = tier)) +
      geom_boxplot(outlier.size = 0.2, linewidth = 0.3, alpha = 0.7, width = 0.6, show.legend = FALSE) +
      scale_fill_brewer(palette = "Blues") +
      labs(x = "Quintile", y = "Aortic distensibility", title = "Aortic distensibility")
  } else {
    p10a <- ggplot() + annotate("text", x = .5, y = .5, label = "Aorta data unavailable", size = 2.5) + theme_void()
  }

  # --- Panel b: LV mass ---
  lv_targets <- c("participant.p24100_i2", "participant.p24101_i2", "participant.p24102_i2",
                   "participant.p24103_i2", "participant.p24104_i2", "participant.p24105_i2")
  lv_data <- stream_mri(mri_files$lv, lv_targets)
  if (!is.null(lv_data) && length(lv_data) > 0) {
    # Use last available LV metric
    lv_d <- lv_data[[length(lv_data)]]
    lv_d$tier <- factor(lv_d$tier)
    col_name <- names(lv_data)[length(lv_data)]
    p10b <- ggplot(lv_d, aes(x = tier, y = value, fill = tier)) +
      geom_boxplot(outlier.size = 0.2, linewidth = 0.3, alpha = 0.7, width = 0.6, show.legend = FALSE) +
      scale_fill_brewer(palette = "Reds") +
      labs(x = "Quintile", y = "LV parameter", title = paste0("LV: ", col_name))
  } else {
    p10b <- ggplot() + annotate("text", x = .5, y = .5, label = "LV data unavailable", size = 2.5) + theme_void()
  }

  # --- Panel c: GLS (strain) ---
  strain_targets <- c("participant.p24140_i2", "participant.p24157_i2", "participant.p24174_i2",
                       "participant.p24181_i2")
  strain_data <- stream_mri(mri_files$strain, strain_targets)
  if (!is.null(strain_data) && length(strain_data) > 0) {
    gls <- strain_data[[1]]
    gls$tier <- factor(gls$tier)
    p10c <- ggplot(gls, aes(x = tier, y = value, fill = tier)) +
      geom_boxplot(outlier.size = 0.2, linewidth = 0.3, alpha = 0.7, width = 0.6, show.legend = FALSE) +
      scale_fill_brewer(palette = "Greens") +
      labs(x = "Quintile", y = "GLS (%)", title = "Global longitudinal strain")
  } else {
    p10c <- ggplot() + annotate("text", x = .5, y = .5, label = "Strain data unavailable", size = 2.5) + theme_void()
  }

  # --- Panel d: CACS ---
  cacs_targets <- c("participant.p31063_i2", "participant.p31060_i2", "participant.p31075_i2",
                     "participant.p31085_i2")
  cacs_data <- stream_mri(mri_files$cacs, cacs_targets)
  if (!is.null(cacs_data) && length(cacs_data) > 0) {
    cacs_d <- cacs_data[[1]]
    cacs_d <- cacs_d %>% filter(value > 0)
    if (nrow(cacs_d) > 10) {
      p10d <- ggplot(cacs_d, aes(x = log10(value + 1))) +
        geom_histogram(bins = 40, fill = "#762A83", color = "white", linewidth = 0.15, alpha = 0.8) +
        labs(x = "log10(CACS + 1)", y = "Count", title = "Coronary artery calcification")
    } else {
      p10d <- ggplot() + annotate("text", x = .5, y = .5, label = "Too few CACS values", size = 2.5) + theme_void()
    }
  } else {
    p10d <- ggplot() + annotate("text", x = .5, y = .5, label = "CACS data unavailable", size = 2.5) + theme_void()
  }

  combined <- (p10a | p10b) / (p10c | p10d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_10_cardiac_mri", h_mm = 140)
  message("  Figure 10 done.")
}

# ==============================================================================
# FIGURE 11: Multi-Modal Integration Cascade (full page)
# ==============================================================================
figure_11 <- function() {
  message("Generating Figure 11: Multi-Modal Integration Cascade...")

  # --- Panel a: 9-level cascade diagram ---
  levels_df <- data.frame(
    level = 1:9,
    y = seq(9, 1, -1),
    label = c("1. Variant\nIdentification",
              "2. Structure\n(pLDDT)",
              "3. Stability\n(ddG)",
              "4. Function\n(SSS)",
              "5. Metabolism\n(NMR)",
              "6. Particles\n(ApoB/LDL)",
              "7. Vessel\n(MRI)",
              "8. Inflammation\n(GlycA)",
              "9. Events\n(ASCVD)"),
    fill = c("#E3F2FD", "#BBDEFB", "#90CAF9", "#64B5F6", "#42A5F5",
             "#1E88E5", "#1565C0", "#0D47A1", "#B71C1C"),
    text_col = c(rep("black", 5), rep("white", 4)),
    stringsAsFactors = FALSE
  )

  arrows_df <- data.frame(
    x = rep(3, 8), xend = rep(3, 8),
    y = levels_df$y[1:8] - 0.35, yend = levels_df$y[2:9] + 0.35
  )

  p11a <- ggplot() +
    geom_rect(data = levels_df,
              aes(xmin = 1.5, xmax = 4.5, ymin = y - 0.35, ymax = y + 0.35),
              fill = levels_df$fill, color = "grey30", linewidth = 0.3) +
    geom_text(data = levels_df, aes(x = 3, y = y, label = label),
              size = 2, lineheight = 0.8, color = levels_df$text_col, fontface = "bold") +
    geom_segment(data = arrows_df, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                 linewidth = 0.4, color = "grey30") +
    coord_cartesian(xlim = c(0.5, 5.5), ylim = c(0.3, 9.7)) +
    theme_void(base_size = 7, base_family = "Arial") +
    labs(title = "Multi-modal integration cascade") +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))

  # --- Panel b: Correlation matrix heatmap from nobel ---
  cor_vars <- c("sss", "ddG", "plddt", "mean_ldl", "mean_apob",
                "mean_discordance", "mean_lpa", "ascvd_pct")
  avail_vars <- intersect(cor_vars, colnames(nobel))
  if (length(avail_vars) >= 3) {
    cor_df <- nobel[, avail_vars, drop = FALSE]
    cor_df <- cor_df[complete.cases(cor_df), ]
    if (nrow(cor_df) >= 5) {
      cm <- cor(cor_df, use = "pairwise.complete.obs", method = "pearson")
      cm_long <- as.data.frame(as.table(cm))
      colnames(cm_long) <- c("Var1", "Var2", "r")

      p11b <- ggplot(cm_long, aes(x = Var1, y = Var2, fill = r)) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = sprintf("%.2f", r)), size = 1.8) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, limits = c(-1, 1), name = "r") +
        labs(x = NULL, y = NULL, title = "Correlation matrix") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
              axis.text.y = element_text(size = 5))
    } else {
      p11b <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient data", size = 3) + theme_void()
    }
  } else {
    p11b <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient variables", size = 3) + theme_void()
  }

  combined <- (p11a | p11b) +
    plot_layout(widths = c(1, 1.5)) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_11_multimodal_cascade", w_mm = DOUBLE_COL, h_mm = 200)
  message("  Figure 11 done.")
}

# ==============================================================================
# FIGURE 12: Compound Risk Stratification
# ==============================================================================
figure_12 <- function() {
  message("Generating Figure 12: Compound Risk Stratification...")

  df12 <- comp_sss %>% filter(!is.na(sss), !is.na(ascvd))

  # --- Panel a: 2x2 risk heatmap (SSS x Lp(a)) ---
  df12a <- df12 %>% filter(!is.na(lpa))
  df12a$sss_cat <- ifelse(df12a$sss >= median(df12a$sss, na.rm = TRUE), "High SSS", "Low SSS")
  df12a$lpa_cat <- ifelse(df12a$lpa >= 143, "High Lp(a)", "Low Lp(a)")

  # If apob available, add discordance dimension
  if ("apob" %in% colnames(df12a) && sum(!is.na(df12a$apob)) > 20) {
    df12a$disc <- df12a$apob / df12a$ldl1
    df12a$disc_cat <- ifelse(df12a$disc >= median(df12a$disc, na.rm = TRUE), "High Disc", "Low Disc")

    heat12 <- df12a %>%
      group_by(sss_cat, lpa_cat) %>%
      summarise(ascvd_pct = round(100 * mean(ascvd, na.rm = TRUE), 1),
                n = n(), .groups = "drop")
  } else {
    heat12 <- df12a %>%
      group_by(sss_cat, lpa_cat) %>%
      summarise(ascvd_pct = round(100 * mean(ascvd, na.rm = TRUE), 1),
                n = n(), .groups = "drop")
  }

  p12a <- ggplot(heat12, aes(x = lpa_cat, y = sss_cat, fill = ascvd_pct)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = paste0(ascvd_pct, "%\n(n=", n, ")")), size = 2.5, fontface = "bold") +
    scale_fill_gradient(low = "#FDDBC7", high = "#67001F", name = "ASCVD %") +
    labs(x = NULL, y = NULL, title = "Compound risk stratification")

  # --- Panel b: Incremental AUC ---
  roc_df <- df12 %>% filter(!is.na(age), !is.na(sex))

  models <- list()
  models[["Age + Sex"]] <- glm(ascvd ~ age + sex, data = roc_df, family = "binomial")
  models[["+ SSS"]] <- glm(ascvd ~ age + sex + sss, data = roc_df, family = "binomial")

  if (sum(!is.na(roc_df$lpa)) > 50) {
    roc_lpa <- roc_df %>% filter(!is.na(lpa))
    models[["+ Lp(a)"]] <- glm(ascvd ~ age + sex + sss + lpa, data = roc_lpa, family = "binomial")
  }
  if (sum(!is.na(roc_df$apob) & !is.na(roc_df$ldl1)) > 50) {
    roc_disc <- roc_df %>% filter(!is.na(apob), !is.na(ldl1))
    roc_disc$disc <- roc_disc$apob / roc_disc$ldl1
    models[["+ Discordance"]] <- glm(ascvd ~ age + sex + sss + disc, data = roc_disc, family = "binomial")
  }

  auc_list <- lapply(names(models), function(nm) {
    m <- models[[nm]]
    pred <- predict(m, type = "response")
    r <- pROC::roc(m$model[[1]], pred, quiet = TRUE)
    data.frame(model = nm, auc = as.numeric(pROC::auc(r)))
  })
  auc_df <- bind_rows(auc_list)
  auc_df$model <- factor(auc_df$model, levels = auc_df$model)

  p12b <- ggplot(auc_df, aes(x = model, y = auc, fill = model)) +
    geom_col(width = 0.6, color = "grey30", linewidth = 0.2, alpha = 0.8, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.3, size = 2) +
    scale_fill_brewer(palette = "Blues") +
    scale_y_continuous(limits = c(0.45, max(auc_df$auc) * 1.08), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "AUC", title = "Incremental discrimination") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  # --- Panel c: Risk group comparison ---
  if ("disc_cat" %in% colnames(df12a)) {
    risk_groups <- df12a %>%
      mutate(risk_group = case_when(
        sss_cat == "High SSS" & lpa_cat == "High Lp(a)" & disc_cat == "High Disc" ~ "Triple High",
        sss_cat == "Low SSS" & lpa_cat == "Low Lp(a)" & disc_cat == "Low Disc" ~ "Triple Low",
        TRUE ~ "Intermediate"
      )) %>%
      group_by(risk_group) %>%
      summarise(ascvd_pct = 100 * mean(ascvd, na.rm = TRUE), n = n(), .groups = "drop") %>%
      filter(risk_group != "Intermediate")
  } else {
    risk_groups <- df12a %>%
      mutate(risk_group = case_when(
        sss_cat == "High SSS" & lpa_cat == "High Lp(a)" ~ "Double High",
        sss_cat == "Low SSS" & lpa_cat == "Low Lp(a)" ~ "Double Low",
        TRUE ~ "Intermediate"
      )) %>%
      group_by(risk_group) %>%
      summarise(ascvd_pct = 100 * mean(ascvd, na.rm = TRUE), n = n(), .groups = "drop") %>%
      filter(risk_group != "Intermediate")
  }

  p12c <- ggplot(risk_groups, aes(x = risk_group, y = ascvd_pct, fill = risk_group)) +
    geom_col(width = 0.5, color = "grey30", linewidth = 0.2, alpha = 0.8) +
    geom_text(aes(label = paste0(round(ascvd_pct, 1), "%\n(n=", n, ")")),
              vjust = -0.3, size = 2) +
    scale_fill_manual(values = c("Triple High" = "#B2182B", "Triple Low" = "#2166AC",
                                  "Double High" = "#B2182B", "Double Low" = "#2166AC"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "ASCVD prevalence (%)", title = "Extreme risk groups")

  # --- Panel d: Clinical decision tree ---
  tree_nodes <- data.frame(
    x = c(3, 1.5, 4.5, 0.75, 2.25, 3.75, 5.25),
    y = c(3, 2, 2, 1, 1, 1, 1),
    label = c("FH Patient\n(SSS available)", "SSS < 0.4\nLow risk", "SSS >= 0.4\nHigh risk",
              "Standard\nstatin", "Intensify\ntherapy", "Add\nezetimibe", "Consider\nPCSK9i"),
    fill = c("#E3F2FD", "#C8E6C9", "#FFCDD2", "#A5D6A7", "#FFF9C4", "#FFCC80", "#EF9A9A"),
    stringsAsFactors = FALSE
  )
  tree_edges <- data.frame(
    x = c(3, 3, 1.5, 1.5, 4.5, 4.5),
    xend = c(1.5, 4.5, 0.75, 2.25, 3.75, 5.25),
    y = c(2.7, 2.7, 1.7, 1.7, 1.7, 1.7),
    yend = c(2.3, 2.3, 1.3, 1.3, 1.3, 1.3)
  )

  p12d <- ggplot() +
    geom_segment(data = tree_edges, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 0.3, color = "grey40") +
    geom_rect(data = tree_nodes,
              aes(xmin = x - 0.55, xmax = x + 0.55, ymin = y - 0.25, ymax = y + 0.25),
              fill = tree_nodes$fill, color = "grey30", linewidth = 0.3) +
    geom_text(data = tree_nodes, aes(x = x, y = y, label = label),
              size = 1.6, lineheight = 0.8) +
    coord_cartesian(xlim = c(-0.2, 6.2), ylim = c(0.5, 3.5)) +
    theme_void(base_size = 7) +
    labs(title = "Clinical decision framework") +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))

  combined <- (p12a | p12b) / (p12c | p12d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_12_compound_risk", h_mm = 155)
  message("  Figure 12 done.")
}

# ==============================================================================
# FIGURE 13: Variant Phenotype Clustering
# ==============================================================================
figure_13 <- function() {
  message("Generating Figure 13: Variant Phenotype Clustering...")

  # --- Panel a: Cluster scatter ---
  if (!"cluster" %in% colnames(nobel)) {
    warning("No cluster column in nobel"); return(invisible(NULL))
  }
  nob <- nobel %>% filter(!is.na(sss), !is.na(mean_ldl), !is.na(cluster))

  p13a <- ggplot(nob, aes(x = sss, y = mean_ldl, color = cluster)) +
    geom_point(size = 1.5, alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = variant), size = 1.4, max.overlaps = 15,
                              segment.size = 0.2, segment.color = "grey60") +
    scale_color_brewer(palette = "Set1", name = "Cluster") +
    labs(x = "Structural Severity Score", y = "Mean LDL-C (mmol/L)",
         title = "Variant phenotype clusters")

  # --- Panel b: Cluster characteristics ---
  char_vars <- c("sss", "ddG", "mean_ldl", "mean_apob", "ascvd_pct", "mean_lpa")
  avail_cv <- intersect(char_vars, colnames(nob))
  if (length(avail_cv) >= 2) {
    clust_means <- nob %>%
      group_by(cluster) %>%
      summarise(across(all_of(avail_cv), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      pivot_longer(-cluster, names_to = "metric", values_to = "value")

    p13b <- ggplot(clust_means, aes(x = metric, y = value, fill = cluster)) +
      geom_col(position = "dodge", width = 0.7, color = "grey30", linewidth = 0.15) +
      scale_fill_brewer(palette = "Set1", name = "Cluster") +
      labs(x = NULL, y = "Mean value", title = "Cluster characteristics") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      facet_wrap(~metric, scales = "free_y", nrow = 1)
  } else {
    p13b <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient variables", size = 3) + theme_void()
  }

  # --- Panel c: Cluster ASCVD rates ---
  if ("ascvd_pct" %in% colnames(nob) && "n" %in% colnames(nob)) {
    clust_ascvd <- nob %>%
      group_by(cluster) %>%
      summarise(mean_ascvd = mean(ascvd_pct, na.rm = TRUE),
                se = sd(ascvd_pct, na.rm = TRUE) / sqrt(n()),
                n = n(), .groups = "drop")

    p13c <- ggplot(clust_ascvd, aes(x = cluster, y = mean_ascvd, fill = cluster)) +
      geom_col(width = 0.6, color = "grey30", linewidth = 0.2, alpha = 0.8) +
      geom_errorbar(aes(ymin = pmax(mean_ascvd - 1.96 * se, 0),
                        ymax = mean_ascvd + 1.96 * se), width = 0.2, linewidth = 0.3) +
      geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 2) +
      scale_fill_brewer(palette = "Set1", guide = "none") +
      labs(x = "Cluster", y = "Mean ASCVD %", title = "ASCVD by cluster") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  } else {
    p13c <- ggplot() + annotate("text", x = .5, y = .5, label = "No ascvd_pct column", size = 3) + theme_void()
  }

  # --- Panel d: Exemplar variant table-plot ---
  top_per_cluster <- nob %>%
    group_by(cluster) %>%
    arrange(desc(sss)) %>%
    slice_head(n = 3) %>%
    ungroup() %>%
    select(cluster, variant, sss, mean_ldl) %>%
    mutate(row = row_number())

  p13d <- ggplot(top_per_cluster, aes(y = reorder(variant, -row))) +
    geom_text(aes(x = 0.1, label = cluster), size = 1.8, hjust = 0, color = "grey30") +
    geom_text(aes(x = 0.4, label = variant), size = 1.8, hjust = 0) +
    geom_text(aes(x = 0.7, label = sprintf("%.2f", sss)), size = 1.8, hjust = 0, color = "#2166AC") +
    geom_text(aes(x = 0.9, label = sprintf("%.1f", mean_ldl)), size = 1.8, hjust = 0, color = "#B2182B") +
    annotate("text", x = c(0.1, 0.4, 0.7, 0.9),
             y = rep(nrow(top_per_cluster) + 0.8, 4),
             label = c("Cluster", "Variant", "SSS", "LDL"),
             size = 2, fontface = "bold", hjust = 0) +
    xlim(0, 1.1) +
    theme_void(base_size = 7) +
    labs(title = "Exemplar variants") +
    theme(plot.title = element_text(size = 7, face = "bold"))

  combined <- (p13a | p13b) / (p13c | p13d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_13_clustering", h_mm = 160)
  message("  Figure 13 done.")
}

# ==============================================================================
# FIGURE 14: Treatment Paradox & Age Interaction
# ==============================================================================
figure_14 <- function() {
  message("Generating Figure 14: Treatment Paradox & Age Interaction...")

  df14 <- add_sss_tertile(comp_sss)

  # --- Panel a: SSS x Age interaction ---
  df14a <- df14 %>% filter(!is.na(age), !is.na(ascvd), !is.na(sss_tertile))
  df14a$age_decade <- cut(df14a$age, breaks = c(0, 30, 40, 50, 60, 70, 100),
                          labels = c("<30", "30-39", "40-49", "50-59", "60-69", "70+"),
                          include.lowest = TRUE)

  age_sss <- df14a %>%
    group_by(age_decade, sss_tertile) %>%
    summarise(ascvd_pct = 100 * mean(ascvd, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= 3)

  p14a <- ggplot(age_sss, aes(x = age_decade, y = ascvd_pct, color = sss_tertile, group = sss_tertile)) +
    geom_line(linewidth = 0.5) +
    geom_point(aes(size = n), alpha = 0.7) +
    scale_color_manual(values = sss_tert_pal, name = "SSS tertile") +
    scale_size_continuous(range = c(1, 3), guide = "none") +
    labs(x = "Age decade", y = "ASCVD prevalence (%)", title = "Age x SSS interaction")

  # --- Panel b: SSS x Statin interaction ---
  df14b <- df14 %>% filter(!is.na(on_statin), !is.na(ascvd), !is.na(sss_tertile))
  df14b$statin_grp <- ifelse(df14b$on_statin == 1, "On statin", "No statin")

  statin_sss <- df14b %>%
    group_by(sss_tertile, statin_grp) %>%
    summarise(ascvd_pct = 100 * mean(ascvd, na.rm = TRUE), n = n(), .groups = "drop")

  p14b <- ggplot(statin_sss, aes(x = sss_tertile, y = ascvd_pct, fill = statin_grp)) +
    geom_col(position = "dodge", width = 0.6, color = "grey30", linewidth = 0.2, alpha = 0.8) +
    geom_text(aes(label = paste0("n=", n)), position = position_dodge(0.6),
              vjust = -0.3, size = 1.8) +
    scale_fill_manual(values = c("On statin" = "#2166AC", "No statin" = "#BDBDBD"), name = NULL) +
    labs(x = "SSS tertile", y = "ASCVD prevalence (%)", title = "SSS x statin status") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  # --- Panel c: Cumulative LDL-year burden ---
  df14c <- df14 %>%
    filter(!is.na(age), !is.na(sss_tertile)) %>%
    mutate(ldl_use = ifelse(!is.na(matched_ldl), matched_ldl, ldl1),
           ldl_years = age * ldl_use)

  if (sum(!is.na(df14c$ldl_years)) >= 20) {
    p14c <- ggplot(df14c, aes(x = sss_tertile, y = ldl_years, fill = sss_tertile)) +
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.7, width = 0.6) +
      scale_fill_manual(values = sss_tert_pal, guide = "none") +
      labs(x = "SSS tertile", y = "Cumulative LDL-years\n(age x LDL-C)", title = "LDL burden by SSS")
  } else {
    p14c <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient LDL-year data", size = 3) + theme_void()
  }

  # --- Panel d: Age-stratified AUC ---
  df14d <- df14 %>% filter(!is.na(age), !is.na(ascvd), !is.na(sss))
  df14d$age_grp <- cut(df14d$age, breaks = c(0, 45, 60, 100),
                       labels = c("<45", "45-60", ">60"), include.lowest = TRUE)

  age_auc <- df14d %>%
    group_by(age_grp) %>%
    summarise(
      auc = tryCatch({
        pred <- predict(glm(ascvd ~ sss, family = "binomial"), type = "response")
        as.numeric(pROC::auc(pROC::roc(ascvd, pred, quiet = TRUE)))
      }, error = function(e) NA_real_),
      n = n(), .groups = "drop"
    )

  # Compute AUC properly for each age group
  auc_results <- lapply(levels(df14d$age_grp), function(g) {
    sub <- df14d %>% filter(age_grp == g)
    if (nrow(sub) < 20 || length(unique(sub$ascvd)) < 2) {
      return(data.frame(age_grp = g, auc = NA_real_, n = nrow(sub)))
    }
    pred <- predict(glm(ascvd ~ sss + age + sex, data = sub, family = "binomial"), type = "response")
    r <- pROC::roc(sub$ascvd, pred, quiet = TRUE)
    data.frame(age_grp = g, auc = as.numeric(pROC::auc(r)), n = nrow(sub))
  })
  auc_age <- bind_rows(auc_results) %>% filter(!is.na(auc))

  if (nrow(auc_age) >= 1) {
    p14d <- ggplot(auc_age, aes(x = age_grp, y = auc, fill = age_grp)) +
      geom_col(width = 0.5, color = "grey30", linewidth = 0.2, alpha = 0.8, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.3f\n(n=%d)", auc, n)), vjust = -0.2, size = 2) +
      scale_fill_brewer(palette = "Oranges") +
      scale_y_continuous(limits = c(0.4, max(auc_age$auc, na.rm = TRUE) * 1.1),
                         expand = expansion(mult = c(0, 0.1))) +
      labs(x = "Age group", y = "AUC (SSS + age + sex)", title = "Age-stratified discrimination")
  } else {
    p14d <- ggplot() + annotate("text", x = .5, y = .5, label = "AUC computation failed", size = 3) + theme_void()
  }

  combined <- (p14a | p14b) / (p14c | p14d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_14_age_interaction", h_mm = 145)
  message("  Figure 14 done.")
}

# ==============================================================================
# FIGURE 15: CALON-Lite vs SAFEHEART
# ==============================================================================
figure_15 <- function() {
  message("Generating Figure 15: CALON-Lite vs SAFEHEART...")

  if (is.null(subgroup)) {
    warning("Subgroup data not available"); return(invisible(NULL))
  }

  # --- Panel a: Master forest plot ---
  sg <- subgroup
  sg$label <- paste0(sg$direction, " | ", sg$subgroup, ": ", sg$value)

  p15a <- ggplot(sg, aes(y = reorder(label, delta))) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
    geom_point(aes(x = auc_lite), color = "#2166AC", size = 1.5, shape = 16) +
    geom_point(aes(x = auc_safe), color = "#B2182B", size = 1.5, shape = 17) +
    geom_segment(aes(x = auc_lite, xend = auc_safe, yend = reorder(label, delta)),
                 linewidth = 0.3, color = "grey50") +
    geom_text(aes(x = pmax(auc_lite, auc_safe) + 0.005,
                  label = sprintf("d=%.3f", delta)), size = 1.8, hjust = 0) +
    labs(x = "AUC", y = NULL, title = "CALON-Lite vs SAFEHEART") +
    annotate("text", x = c(0.55, 0.55), y = c(1.5, 0.8),
             label = c("CALON-Lite", "SAFEHEART"),
             color = c("#2166AC", "#B2182B"), size = 2, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

  # --- Panel b: Calibration plot ---
  cal_df <- comp_sss %>% filter(!is.na(ascvd), !is.na(sss), !is.na(age), !is.na(sex))
  if (nrow(cal_df) >= 50) {
    cal_model <- glm(ascvd ~ sss + age + sex, data = cal_df, family = "binomial")
    cal_df$pred <- predict(cal_model, type = "response")
    cal_df$decile <- ntile(cal_df$pred, 10)

    cal_summary <- cal_df %>%
      group_by(decile) %>%
      summarise(pred_mean = mean(pred), obs_mean = mean(ascvd), n = n(), .groups = "drop")

    p15b <- ggplot(cal_summary, aes(x = pred_mean, y = obs_mean)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
      geom_point(size = 2, color = "#2166AC") +
      geom_line(linewidth = 0.4, color = "#2166AC") +
      geom_text(aes(label = decile), vjust = -0.7, size = 1.8, color = "grey40") +
      labs(x = "Predicted probability", y = "Observed proportion",
           title = "Calibration plot") +
      coord_equal(xlim = c(0, max(cal_summary$pred_mean, cal_summary$obs_mean) * 1.1),
                  ylim = c(0, max(cal_summary$pred_mean, cal_summary$obs_mean) * 1.1))
  } else {
    p15b <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient data", size = 3) + theme_void()
  }

  # --- Panel c: Decision curve analysis ---
  if (nrow(cal_df) >= 50) {
    thresholds <- seq(0.01, 0.50, 0.01)
    prevalence <- mean(cal_df$ascvd)

    dca_results <- lapply(thresholds, function(pt) {
      # Treat all
      nb_all <- prevalence - (1 - prevalence) * pt / (1 - pt)
      # Model
      tp <- sum(cal_df$pred >= pt & cal_df$ascvd == 1) / nrow(cal_df)
      fp <- sum(cal_df$pred >= pt & cal_df$ascvd == 0) / nrow(cal_df)
      nb_model <- tp - fp * pt / (1 - pt)

      data.frame(threshold = pt, treat_all = nb_all, model = nb_model, treat_none = 0)
    })
    dca_df <- bind_rows(dca_results) %>%
      pivot_longer(-threshold, names_to = "strategy", values_to = "net_benefit")

    p15c <- ggplot(dca_df, aes(x = threshold, y = net_benefit, color = strategy)) +
      geom_line(linewidth = 0.5) +
      scale_color_manual(values = c("treat_all" = "grey60", "model" = "#B2182B", "treat_none" = "grey30"),
                         labels = c("SSS model", "Treat all", "Treat none"), name = NULL) +
      labs(x = "Threshold probability", y = "Net benefit", title = "Decision curve analysis") +
      coord_cartesian(ylim = c(-0.05, max(dca_df$net_benefit, na.rm = TRUE) * 1.1)) +
      theme(legend.position = c(0.7, 0.8), legend.background = element_blank())
  } else {
    p15c <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient data", size = 3) + theme_void()
  }

  # --- Panel d: NRI bar ---
  if (nrow(cal_df) >= 50) {
    base_pred <- predict(glm(ascvd ~ age + sex, data = cal_df, family = "binomial"), type = "response")
    sss_pred <- predict(glm(ascvd ~ age + sex + sss, data = cal_df, family = "binomial"), type = "response")

    # Compute continuous NRI
    events <- cal_df$ascvd == 1
    nonevents <- cal_df$ascvd == 0

    nri_events <- mean(sss_pred[events] > base_pred[events]) -
      mean(sss_pred[events] < base_pred[events])
    nri_nonevents <- mean(sss_pred[nonevents] < base_pred[nonevents]) -
      mean(sss_pred[nonevents] > base_pred[nonevents])
    nri_total <- nri_events + nri_nonevents

    nri_df <- data.frame(
      component = c("Events", "Non-events", "Total NRI"),
      nri = c(nri_events, nri_nonevents, nri_total)
    )
    nri_df$component <- factor(nri_df$component, levels = c("Events", "Non-events", "Total NRI"))

    p15d <- ggplot(nri_df, aes(x = component, y = nri, fill = component)) +
      geom_col(width = 0.5, color = "grey30", linewidth = 0.2, alpha = 0.8) +
      geom_text(aes(label = sprintf("%.3f", nri)), vjust = ifelse(nri_df$nri >= 0, -0.3, 1.3), size = 2.2) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      scale_fill_manual(values = c("Events" = "#2166AC", "Non-events" = "#35978F", "Total NRI" = "#B2182B"),
                        guide = "none") +
      labs(x = NULL, y = "Net Reclassification Index", title = "Continuous NRI (+SSS vs base)")
  } else {
    p15d <- ggplot() + annotate("text", x = .5, y = .5, label = "Insufficient data", size = 3) + theme_void()
  }

  combined <- (p15a | p15b) / (p15c | p15d) +
    plot_annotation(tag_levels = "a",
                    theme = theme(plot.tag = element_text(size = 8, face = "bold")))
  save_figure(combined, "Figure_15_calon_safeheart", h_mm = 155)
  message("  Figure 15 done.")
}

# ==============================================================================
# SECTION 16: EXECUTE ALL FIGURES
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("GENERATING ALL 15 FIGURES")
message(paste(rep("=", 70), collapse = ""), "\n")

fig_functions <- list(
  figure_1, figure_2, figure_3, figure_4, figure_5,
  figure_6, figure_7, figure_8, figure_9, figure_10,
  figure_11, figure_12, figure_13, figure_14, figure_15
)

for (i in seq_along(fig_functions)) {
  tryCatch({
    fig_functions[[i]]()
  }, error = function(e) {
    message(sprintf("  ERROR in Figure %d: %s", i, e$message))
  })
}

message("\n", paste(rep("=", 70), collapse = ""))
message("ALL FIGURES COMPLETE")
message("Output directory: ", FIG_DIR)
message(paste(rep("=", 70), collapse = ""))
