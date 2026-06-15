# =============================================================================
# HEART UK Abstract - TUDOR Figure 2
# Gene-Specific Discrimination + ApoB Augmentation
# Forest plot: AUC by gene (LDLR, APOB) in both cohorts + ApoB models
# =============================================================================

library(ggplot2)

df <- data.frame(
  Label = c(
    "Wales: LDLR (n=724)",
    "Wales: APOB (n=96)",
    "UKB: LDLR (n=515)",
    "UKB: APOB (n=213)",
    "",
    "UKB: TUDOR base",
    "UKB: + ApoB/LDL-C",
    "UKB: + ApoB + ratio",
    "UKB: + ApoB + ratio + Lp(a)"
  ),
  AUC = c(0.839, 0.841, 0.717, 0.830,
          NA,
          0.753, 0.760, 0.771, 0.776),
  CI_lo = c(0.817, 0.790, 0.693, 0.802,
            NA,
            0.733, 0.740, 0.751, 0.755),
  CI_hi = c(0.861, 0.892, 0.741, 0.858,
            NA,
            0.772, 0.780, 0.791, 0.797),
  Group = c(rep("Gene-specific", 4), NA, rep("ApoB augmentation", 4)),
  stringsAsFactors = FALSE
)

# Remove spacer
df_plot <- df[!is.na(df$AUC),]
df_plot$Label <- factor(df_plot$Label, levels=rev(df_plot$Label))
df_plot$Group <- factor(df_plot$Group, levels=c("Gene-specific", "ApoB augmentation"))

cols <- c("Gene-specific"="#1B4F72", "ApoB augmentation"="#E67E22")

p <- ggplot(df_plot, aes(x=AUC, y=Label, colour=Group)) +
  geom_vline(xintercept=c(0.7, 0.8), linetype="dotted", colour="grey70", linewidth=0.3) +
  geom_pointrange(aes(xmin=CI_lo, xmax=CI_hi), linewidth=0.5, size=0.4) +
  geom_text(aes(x=CI_hi + 0.008, label=sprintf("%.3f", AUC)),
            hjust=0, size=1.9, show.legend=FALSE, fontface="bold") +
  scale_colour_manual(values=cols) +
  scale_x_continuous(limits=c(0.65, 0.92), breaks=seq(0.65, 0.90, 0.05),
                     name="AUC (95% CI)") +
  labs(title="Gene-Specific Validation & ApoB Augmentation",
       colour=NULL) +
  theme_minimal(base_size=7) +
  theme(
    plot.title = element_text(face="bold", size=7, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=5.5),
    axis.title.x = element_text(size=6.5),
    axis.text.x = element_text(size=5.5),
    legend.position = "bottom",
    legend.text = element_text(size=5.5),
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(-3, 0, 0, 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 10, 2, 3, "pt")
  )

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_fig2.png",
       p, width=10.5, height=5, units="cm", dpi=300, bg="white")

cat("TUDOR Figure 2 saved\n")
