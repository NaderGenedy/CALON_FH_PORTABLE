# =============================================================================
# HEART UK Abstract - TUDOR Figure 1
# Dual-Validation AUC: TUDOR vs eDLCN vs LDL-C vs Trig_Filter
# Grouped bar chart with CI whiskers - two panels (Wales | UKB)
# =============================================================================

library(ggplot2)
library(gridExtra)

# --- Wales External Validation (TRIPOD 2b) ---
wales <- data.frame(
  Model = factor(c("TUDOR", "DLCN", "FAMCAT\napprox.", "LDL-C\nalone"),
                 levels=c("TUDOR", "DLCN", "FAMCAT\napprox.", "LDL-C\nalone")),
  AUC = c(0.842, 0.791, 0.775, 0.710),
  CI_lo = c(0.822, 0.770, 0.752, 0.688),
  CI_hi = c(0.863, 0.812, 0.798, 0.732)
)

cols4 <- c("TUDOR"="#1B4F72", "DLCN"="#E74C3C", "FAMCAT\napprox."="#F39C12", "LDL-C\nalone"="#95A5A6")

pW <- ggplot(wales, aes(x=Model, y=AUC, fill=Model)) +
  geom_col(width=0.6, colour="grey30", linewidth=0.25) +
  geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi), width=0.15, linewidth=0.4) +
  geom_text(aes(label=sprintf("%.3f", AUC)), vjust=-0.5, size=2, fontface="bold") +
  scale_fill_manual(values=cols4) +
  scale_y_continuous(limits=c(0, 0.95), breaks=seq(0, 0.9, 0.1), name="AUC (95% CI)") +
  labs(title="A. Wales Registry (n=7,253)") +
  theme_minimal(base_size=7) +
  theme(
    plot.title = element_text(face="bold", size=7, hjust=0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=5.5, lineheight=0.85),
    axis.title.y = element_text(size=6),
    axis.text.y = element_text(size=5.5),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 4, 2, 3, "pt")
  )

# --- UK Biobank External Validation (TRIPOD 4) ---
ukb <- data.frame(
  Model = factor(c("TUDOR", "eDLCN", "Trig_Filter\nalone", "LDL-C\nalone"),
                 levels=c("TUDOR", "eDLCN", "Trig_Filter\nalone", "LDL-C\nalone")),
  AUC = c(0.750, 0.636, 0.700, 0.689),
  CI_lo = c(0.731, 0.614, 0.680, 0.668),
  CI_hi = c(0.770, 0.658, 0.720, 0.710)
)

cols4b <- c("TUDOR"="#1B4F72", "eDLCN"="#E74C3C", "Trig_Filter\nalone"="#27AE60", "LDL-C\nalone"="#95A5A6")

pU <- ggplot(ukb, aes(x=Model, y=AUC, fill=Model)) +
  geom_col(width=0.6, colour="grey30", linewidth=0.25) +
  geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi), width=0.15, linewidth=0.4) +
  geom_text(aes(label=sprintf("%.3f", AUC)), vjust=-0.5, size=2, fontface="bold") +
  geom_text(data=data.frame(x=1.5, y=0.82, lab="DeLong p<10\u207B\u00B2\u2074"),
            aes(x=x, y=y, label=lab), inherit.aes=FALSE, size=1.8,
            fontface="italic", colour="grey30") +
  scale_fill_manual(values=cols4b) +
  scale_y_continuous(limits=c(0, 0.95), breaks=seq(0, 0.9, 0.1), name="AUC (95% CI)") +
  labs(title="B. UK Biobank Lipid Clinic (n=58,021)") +
  theme_minimal(base_size=7) +
  theme(
    plot.title = element_text(face="bold", size=7, hjust=0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=5.5, lineheight=0.85),
    axis.title.y = element_text(size=6),
    axis.text.y = element_text(size=5.5),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 3, 2, 4, "pt")
  )

combined <- grid.arrange(pW, pU, ncol=2, widths=c(1, 1))

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/tudor_fig1.png",
       combined, width=10.5, height=5.5, units="cm", dpi=300, bg="white")

cat("TUDOR Figure 1 saved\n")
