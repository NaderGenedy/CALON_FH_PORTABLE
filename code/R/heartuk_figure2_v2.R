# =============================================================================
# HEART UK Abstract - Figure 2 (v2 - refined dual panel)
# ApoB/LDL-C Discordance: FH vs Non-FH
# =============================================================================

library(ggplot2)
library(gridExtra)

# ----- Panel A: OR comparison -----
or_df <- data.frame(
  Pop = factor(c("FH\n(n=1,788)", "Non-FH\n(n=423,002)"),
               levels=c("FH\n(n=1,788)", "Non-FH\n(n=423,002)")),
  OR = c(3.57, 2.88),
  CI_low = c(2.76, 2.83),
  CI_high = c(4.62, 2.93)
)

pA <- ggplot(or_df, aes(x=Pop, y=OR, fill=Pop)) +
  geom_col(width=0.5, colour="grey30", linewidth=0.25) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0.12, linewidth=0.4) +
  geom_text(aes(label=sprintf("%.2f", OR)), vjust=-0.6, size=2.2, fontface="bold") +
  scale_fill_manual(values=c("#1B4F72", "#85C1E9")) +
  scale_y_continuous(limits=c(0, 5.2), breaks=seq(0, 5, 1),
                     name="OR for ASCVD") +
  labs(title="A. OR at \u22650.31 threshold") +
  theme_minimal(base_size=7) +
  theme(
    plot.title = element_text(face="bold", size=7, hjust=0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=6),
    axis.text = element_text(size=5.5),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 4, 2, 3, "pt")
  )

# ----- Panel B: ASCVD Event Rates -----
ev_df <- data.frame(
  Pop = factor(rep(c("FH", "Non-FH"), each=2), levels=c("FH", "Non-FH")),
  Thr = factor(rep(c("<0.31", "\u22650.31"), 2), levels=c("<0.31", "\u22650.31")),
  Rate = c(18.2, 44.3, 22.0, 38.4),
  lab = c("18.2%", "44.3%", "22.0%", "38.4%")
)

pB <- ggplot(ev_df, aes(x=Thr, y=Rate, fill=Pop)) +
  geom_col(position=position_dodge(width=0.6), width=0.5,
           colour="grey30", linewidth=0.25) +
  geom_text(aes(label=lab, group=Pop),
            position=position_dodge(width=0.6), vjust=-0.4,
            size=1.9, fontface="bold") +
  scale_fill_manual(values=c("FH"="#1B4F72", "Non-FH"="#85C1E9")) +
  scale_y_continuous(limits=c(0, 52), breaks=seq(0, 50, 10),
                     name="ASCVD rate (%)") +
  labs(title="B. ASCVD rate by discordance") +
  theme_minimal(base_size=7) +
  theme(
    plot.title = element_text(face="bold", size=7, hjust=0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=6),
    axis.text = element_text(size=5.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=5.5),
    legend.key.size = unit(0.2, "cm"),
    legend.margin = margin(-3, 0, 0, 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 3, 2, 4, "pt")
  )

combined <- grid.arrange(pA, pB, ncol=2, widths=c(1, 1.1))

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_fig2.jpg",
       combined, width=10.5, height=5, units="cm", dpi=300, bg="white")

cat("Figure 2 v2 saved\n")
