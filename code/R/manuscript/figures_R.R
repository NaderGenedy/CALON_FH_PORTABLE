# ============================================================================
# CALON-FH Manuscript: Nature-Calibre R Figures
# Author: Dr Nader Genedy, 2026
# ============================================================================

library(ggplot2)
library(dplyr)
library(cowplot)
library(viridis)
library(scales)

# --- Nature theme ---
theme_nature <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Arial"),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      plot.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 8, face = "bold")
    )
}

# ============================================================================
# FIGURE 2a: Lp(a) Quintile Gradients (Bar chart, 6 outcomes)
# ============================================================================
fig2a_data <- data.frame(
  Outcome = rep(c("MI", "Aortic\nStenosis", "Carotid", "PVD", "Stroke", "TIA"), each = 5),
  Quintile = rep(c("Q1", "Q2", "Q3", "Q4", "Q5"), 6),
  Rate = c(
    3.11, 3.27, 3.36, 3.63, 4.51,  # MI
    1.59, 1.62, 1.67, 1.84, 2.21,  # AS
    0.54, 0.51, 0.57, 0.68, 0.80,  # Carotid
    3.46, 3.51, 3.52, 3.71, 4.11,  # PVD
    2.30, 2.27, 2.40, 2.40, 2.51,  # Stroke
    1.20, 1.15, 1.20, 1.17, 1.28   # TIA
  ),
  RR = rep(c(1.45, 1.39, 1.48, 1.19, 1.09, 1.07), each = 5)
)
fig2a_data$Quintile <- factor(fig2a_data$Quintile, levels = c("Q1","Q2","Q3","Q4","Q5"))
fig2a_data$Outcome <- factor(fig2a_data$Outcome,
  levels = c("MI", "Carotid", "Aortic\nStenosis", "PVD", "Stroke", "TIA"))

p2a <- ggplot(fig2a_data, aes(x = Quintile, y = Rate, fill = Quintile)) +
  geom_col(width = 0.7) +
  facet_wrap(~Outcome, scales = "free_y", nrow = 2) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  labs(x = "Lp(a) Quintile", y = "Event Rate (%)",
       title = "Lp(a) quintile gradients across vascular beds") +
  theme_nature() +
  theme(legend.position = "none")

ggsave("manuscript/fig2a_quintile_gradients.pdf", p2a, width = 183, height = 120, units = "mm", dpi = 300)

# ============================================================================
# FIGURE 3a: LVEF vs Lp(a) Quintile (Box plot)
# ============================================================================
fig3_data <- data.frame(
  Quintile = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  LVEF = c(59.6, 59.6, 59.7, 59.8, 59.5),
  LVmass = c(86.0, 86.2, 85.0, 84.2, 86.0),
  ASCVD = c(6.6, 7.0, 7.2, 7.6, 9.2)
)
fig3_data$Quintile <- factor(fig3_data$Quintile, levels = c("Q1","Q2","Q3","Q4","Q5"))

p3a <- ggplot(fig3_data, aes(x = Quintile, y = LVEF)) +
  geom_col(fill = "#3B528B", width = 0.6) +
  geom_hline(yintercept = 55, linetype = "dashed", colour = "red", linewidth = 0.4) +
  annotate("text", x = 5, y = 55.5, label = "Lower limit normal", size = 2.5, colour = "red") +
  coord_cartesian(ylim = c(50, 65)) +
  labs(x = "Lp(a) Quintile", y = "LVEF (%)",
       title = "LVEF preserved across Lp(a) range (n = 60,947)") +
  annotate("text", x = 3, y = 63, label = "rho = 0.003, P = 0.49", size = 2.5, fontface = "italic") +
  theme_nature()

p3b <- ggplot(fig3_data, aes(x = Quintile, y = ASCVD)) +
  geom_col(fill = "#440154", width = 0.6) +
  labs(x = "Lp(a) Quintile", y = "ASCVD Rate (%)",
       title = "ASCVD increases despite preserved LV function") +
  theme_nature()

p3 <- plot_grid(p3a, p3b, labels = c("a", "b"), label_size = 12, nrow = 1)
ggsave("manuscript/fig3_cardiac_mri.pdf", p3, width = 183, height = 80, units = "mm", dpi = 300)

# ============================================================================
# FIGURE 4a: Correlation Heatmap (FH vs Non-FH)
# ============================================================================
corr_fh <- data.frame(
  Variable = c("ApoB", "LDL-C", "TC", "HDL-C", "TG", "ApoB/LDL", "BMI", "CRP"),
  FH = c(0.029, 0.031, 0.007, -0.006, -0.083, 0.005, -0.001, -0.031),
  NonFH = c(0.080, 0.081, 0.069, 0.021, -0.048, -0.010, NA, NA)
)

library(tidyr)
corr_long <- corr_fh %>%
  pivot_longer(cols = c(FH, NonFH), names_to = "Group", values_to = "rho") %>%
  mutate(Variable = factor(Variable, levels = rev(corr_fh$Variable)),
         Group = factor(Group, levels = c("FH", "NonFH"), labels = c("FH", "Non-FH")))

p4a <- ggplot(corr_long, aes(x = Group, y = Variable, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(rho), "", sprintf("%.3f", rho))), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-0.1, 0.1),
                        name = "Spearman rho", na.value = "grey90") +
  labs(x = "", y = "", title = "Lp(a) correlation with metabolic markers") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 0))

ggsave("manuscript/fig4a_correlation_heatmap.pdf", p4a, width = 89, height = 100, units = "mm", dpi = 300)

# ============================================================================
# FIGURE 7: Menopause Effect (Line plot)
# ============================================================================
fig7_data <- data.frame(
  Age = rep(c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"), 2),
  Sex = rep(c("Women", "Men"), each = 6),
  Lpa = c(19.7, 20.1, 21.3, 22.6, 22.9, 23.3,  # Women
          19.7, 19.5, 19.4, 19.5, 19.9, 20.0)   # Men
)
fig7_data$Age <- factor(fig7_data$Age, levels = c("40-44","45-49","50-54","55-59","60-64","65-69"))

p7a <- ggplot(fig7_data, aes(x = Age, y = Lpa, colour = Sex, group = Sex)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 3, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  annotate("text", x = 3.2, y = 24, label = "Menopause", size = 2.5,
           colour = "grey40", hjust = 0) +
  scale_colour_manual(values = c("Women" = "#D7191C", "Men" = "#2C7BB6")) +
  labs(x = "Age group", y = "Median Lp(a) (nmol/L)",
       title = "Menopausal Lp(a) rise: +14% (P < 0.0001)") +
  theme_nature()

# FH women: no effect
fig7b_data <- data.frame(
  Age = c("40-49", "50-54", "55-64", "65-74"),
  FH = c(27.6, 25.7, 27.8, 27.1),
  NonFH = c(20.0, 21.2, 22.8, 23.5)
)
fig7b_data$Age <- factor(fig7b_data$Age, levels = fig7b_data$Age)

library(tidyr)
fig7b_long <- fig7b_data %>% pivot_longer(cols = c(FH, NonFH), names_to = "Group", values_to = "Lpa")

p7b <- ggplot(fig7b_long, aes(x = Age, y = Lpa, colour = Group, group = Group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("FH" = "#D7191C", "NonFH" = "#2C7BB6"),
                       labels = c("FH Women", "Non-FH Women")) +
  labs(x = "Age group", y = "Median Lp(a) (nmol/L)",
       title = "No menopausal effect in FH (P = 0.65)") +
  theme_nature()

p7 <- plot_grid(p7a, p7b, labels = c("a", "b"), label_size = 12, nrow = 1)
ggsave("manuscript/fig7_menopause.pdf", p7, width = 183, height = 80, units = "mm", dpi = 300)

# ============================================================================
# FIGURE 8a: Forest Plot (Pan-Vascular RR)
# ============================================================================
fig8_data <- data.frame(
  Outcome = c("MI", "Carotid", "Aortic Stenosis", "Composite ASCVD", "PVD", "Stroke", "TIA"),
  RR = c(1.45, 1.48, 1.39, 1.23, 1.19, 1.09, 1.07),
  P = c(1.4e-45, 9.3e-10, 2.2e-18, 4.9e-50, 5.0e-11, 1.1e-2, 0.16)
)
fig8_data$Outcome <- factor(fig8_data$Outcome, levels = rev(fig8_data$Outcome))
fig8_data$sig <- ifelse(fig8_data$P < 0.001, "***", ifelse(fig8_data$P < 0.01, "**",
                  ifelse(fig8_data$P < 0.05, "*", "NS")))

p8 <- ggplot(fig8_data, aes(x = RR, y = Outcome)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_point(aes(size = -log10(P), colour = RR), shape = 18) +
  geom_text(aes(label = sprintf("RR=%.2f %s", RR, sig)), hjust = -0.2, size = 2.5) +
  scale_colour_viridis_c(option = "C", begin = 0.2, end = 0.9, name = "RR") +
  scale_size_continuous(range = c(3, 8), name = "-log10(P)") +
  coord_cartesian(xlim = c(0.9, 1.7)) +
  labs(x = "Relative Risk (Q5 vs Q1)", y = "",
       title = "Lp(a) Q5/Q1 relative risk by vascular bed") +
  theme_nature() +
  theme(legend.position = "right")

ggsave("manuscript/fig8_forest_panvascular.pdf", p8, width = 183, height = 100, units = "mm", dpi = 300)

cat("All figures saved to manuscript/ directory\n")
