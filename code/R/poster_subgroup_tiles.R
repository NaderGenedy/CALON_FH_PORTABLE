# =====================================================================
# CALON-FH poster — 13-tile subgroup-wins infographic (R)
# =====================================================================
# Renders the Direction A subgroup grid with EXACT numbers baked in.
# Each tile shows: subgroup label, mini AUC bar comparison (CALON vs SAFEHEART),
# ΔAUC, and NRI% (or just ΔAUC if NRI is noisy due to <10 events).
#
# Style: Nature/Cell figure-audit compliant. Three-tone palette only.
# Output: ./poster_figures/p5_subgroup_tiles.{pdf,png}
# =====================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(forcats); library(patchwork); library(scales); library(stringr)
})

OUT     <- "C:/Users/nader/Downloads/calon_ukb_pipeline/output/v2"
FIG_DIR <- file.path(OUT, "poster_figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

COL_CALON   <- "#DA291C"
COL_SRE     <- "#5A6772"
COL_NEUTRAL <- "#1B1F23"
COL_NOISY   <- "#9CA3AF"
COL_LIGHT   <- "#E8EAED"

sg <- read_csv(file.path(OUT, "CALON_FINAL_v7_subgroup_NRI.csv"), show_col_types = FALSE)
sgA <- sg |> filter(direction == "A_Wales_to_UKB")

# Tag noisy subgroups (event count too small for stable NRI)
sgA <- sgA |>
  mutate(noisy = events < 10,
         delta_lbl  = sprintf("ΔAUC %+.3f", delta_auc),
         nri_lbl    = ifelse(noisy,
                              sprintf("NRI %+.0f%% †", nri_total*100),
                              sprintf("NRI %+.0f%%", nri_total*100)),
         tile_order = factor(level,
                              levels = c("Female","Male",
                                          "<50y","50-65y",">=65y",
                                          "LDL<4.14","LDL 4.14-6.5","LDL>=6.5",
                                          "No DM","T2DM",
                                          "Non-smoker","Smoker",
                                          "No HTN","HTN")))

# Per-tile data: two-row tibble for the mini bar
tiles_bars <- sgA |>
  select(tile_order, auc_calon, auc_sre, noisy) |>
  pivot_longer(c(auc_calon, auc_sre), names_to = "model", values_to = "auc") |>
  mutate(model = recode(model, "auc_calon" = "CALON-FH", "auc_sre" = "SAFEHEART-RE"))

# Single ggplot per tile (we'll wrap with facet_wrap)
p_tiles <- ggplot(tiles_bars, aes(x = model, y = auc, fill = model)) +
  geom_col(width = 0.65, colour = COL_NEUTRAL, linewidth = 0.3, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f", auc), y = auc),
            vjust = -0.3, size = 2.6, colour = COL_NEUTRAL, fontface = "bold") +
  scale_fill_manual(values = c("CALON-FH" = COL_CALON, "SAFEHEART-RE" = COL_SRE),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 0.95), expand = expansion(0)) +
  geom_text(data = sgA,
            aes(x = 1.5, y = 0.88, label = delta_lbl),
            inherit.aes = FALSE, size = 3, fontface = "bold",
            colour = COL_CALON) +
  geom_text(data = sgA,
            aes(x = 1.5, y = 0.79, label = nri_lbl, colour = noisy),
            inherit.aes = FALSE, size = 2.7, fontface = "bold",
            show.legend = FALSE) +
  geom_text(data = sgA,
            aes(x = 1.5, y = 0.05,
                label = sprintf("n=%d, ev=%d", n, events)),
            inherit.aes = FALSE, size = 2.2, colour = COL_NEUTRAL) +
  scale_colour_manual(values = c("TRUE" = COL_NOISY, "FALSE" = COL_NEUTRAL),
                      guide = "none") +
  facet_wrap(~ tile_order, ncol = 7, nrow = 2,
             strip.position = "top") +
  labs(title = "CALON-FH wins in every prespecified subgroup",
       subtitle = "Direction A external validation (Wales-clean trained → UKB tested) — ΔAUC and NRI total per stratum",
       caption  = paste("Each tile: CALON-FH (red) vs SAFEHEART-RE (grey) external AUC; ΔAUC and NRI total (3-tier, cuts 0.05/0.20).",
                         "\n† = sparse-event subgroup (events < 10) — NRI estimate noisy; ΔAUC remains the primary metric."),
       x = NULL, y = "External AUC") +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, colour = COL_NEUTRAL),
    plot.subtitle    = element_text(size = 10, colour = COL_NEUTRAL),
    plot.caption     = element_text(size = 8, colour = COL_NEUTRAL, hjust = 0,
                                     lineheight = 1.2),
    strip.text       = element_text(face = "bold", colour = COL_NEUTRAL, size = 10),
    strip.background = element_rect(fill = COL_LIGHT, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = COL_LIGHT, linewidth = 0.2),
    panel.border     = element_rect(colour = COL_NEUTRAL, fill = NA, linewidth = 0.4),
    axis.text.x      = element_text(angle = 25, hjust = 1, size = 7,
                                     colour = COL_NEUTRAL),
    axis.text.y      = element_text(size = 7, colour = COL_NEUTRAL),
    axis.title.y     = element_text(size = 9, colour = COL_NEUTRAL),
    axis.line        = element_line(colour = COL_NEUTRAL, linewidth = 0.3),
    legend.position  = "top",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 9, colour = COL_NEUTRAL),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.spacing.x  = unit(0.4, "lines"),
    panel.spacing.y  = unit(0.6, "lines")
  )

ggsave(file.path(FIG_DIR, "p5_subgroup_tiles.png"), p_tiles,
       width = 16, height = 7.5, dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, "p5_subgroup_tiles.pdf"), p_tiles,
       width = 16, height = 7.5, device = cairo_pdf, bg = "white")

cat("\n[DONE] Subgroup tile infographic written to:\n")
cat("  ", file.path(FIG_DIR, "p5_subgroup_tiles.png"), "\n")
cat("  ", file.path(FIG_DIR, "p5_subgroup_tiles.pdf"), "\n")
