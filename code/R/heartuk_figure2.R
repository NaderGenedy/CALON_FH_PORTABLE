# =============================================================================
# HEART UK Abstract - Figure 2
# ApoB/LDL-C Discordance: FH vs Non-FH Prognostic Amplification
# =============================================================================

library(ggplot2)
library(grid)
library(gridExtra)

# ----- Panel A: OR comparison (paired bar with CI whiskers) -----
or_df <- data.frame(
  Population = factor(c("FH\n(n=1,788)", "Non-FH\n(n=423,002)"),
                      levels=c("FH\n(n=1,788)", "Non-FH\n(n=423,002)")),
  OR = c(3.57, 2.88),
  CI_low = c(2.76, 2.83),
  CI_high = c(4.62, 2.93)
)

panelA <- ggplot(or_df, aes(x=Population, y=OR, fill=Population)) +
  geom_col(width=0.55, colour="grey20", linewidth=0.3) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.15, linewidth=0.5) +
  geom_text(aes(label=sprintf("OR %.2f", OR)), vjust=-0.8, size=2.5, fontface="bold") +
  scale_fill_manual(values=c("#1B4F72", "#85C1E9")) +
  scale_y_continuous(limits=c(0, 5.5), breaks=seq(0, 5, 1),
                     name="Odds Ratio for ASCVD\n(ApoB/LDL-C \u22650.31 threshold)") +
  labs(title="A. Prognostic Power Amplification",
       subtitle="Same ratio, stronger signal in FH") +
  theme_minimal(base_size=9) +
  theme(
    text = element_text(family="sans"),
    plot.title = element_text(face="bold", size=8),
    plot.subtitle = element_text(size=6.5, colour="grey40"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7),
    axis.text = element_text(size=7),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 8, 5, 5, "pt")
  )

# ----- Panel B: ASCVD Event Rate by Discordance Threshold -----
event_df <- data.frame(
  Population = factor(rep(c("FH", "Non-FH"), each=2),
                      levels=c("FH", "Non-FH")),
  Threshold = factor(rep(c("Ratio <0.31\n(Concordant)", "Ratio \u22650.31\n(Discordant)"), 2),
                     levels=c("Ratio <0.31\n(Concordant)", "Ratio \u22650.31\n(Discordant)")),
  Rate = c(18.2, 44.3, 22.0, 38.4),
  label = c("18.2%", "44.3%", "22.0%", "38.4%")
)

panelB <- ggplot(event_df, aes(x=Threshold, y=Rate, fill=Population)) +
  geom_col(position=position_dodge(width=0.65), width=0.55,
           colour="grey20", linewidth=0.3) +
  geom_text(aes(label=label, group=Population),
            position=position_dodge(width=0.65), vjust=-0.5,
            size=2.3, fontface="bold") +
  scale_fill_manual(values=c("FH"="#1B4F72", "Non-FH"="#85C1E9")) +
  scale_y_continuous(limits=c(0, 55), breaks=seq(0, 50, 10),
                     name="ASCVD Event Rate (%)") +
  labs(title="B. ASCVD Event Rate by Discordance",
       subtitle="2.4-fold enrichment in discordant FH") +
  theme_minimal(base_size=9) +
  theme(
    text = element_text(family="sans"),
    plot.title = element_text(face="bold", size=8),
    plot.subtitle = element_text(size=6.5, colour="grey40"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7),
    axis.text = element_text(size=7),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=6.5),
    legend.key.size = unit(0.3, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 8, "pt")
  )

# Combine panels
combined <- grid.arrange(panelA, panelB, ncol=2, widths=c(1, 1.15))

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_fig2.jpg",
       combined, width=10.5, height=6, units="cm", dpi=300, bg="white")

cat("Figure 2 saved successfully\n")
