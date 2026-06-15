# =============================================================================
# HEART UK Abstract - Figure 1 (v2 - refined for 10.5cm constraint)
# LOCO-CV Discrimination: CALON-Lite vs SAFEHEART-RE
# =============================================================================

library(ggplot2)

df <- data.frame(
  Fold = rep(c("SW hold-out\n(n=418)",
               "UKB hold-out\n(n=1,625)",
               "Wales hold-out\n(n=2,405)"), each=2),
  Model = rep(c("CALON-Lite", "SRE refitted"), 3),
  AUC = c(0.790, 0.768,
          0.722, 0.700,
          0.855, 0.814),
  CI_low = c(0.741, 0.716,
             0.694, 0.671,
             0.837, 0.791),
  CI_high = c(0.839, 0.820,
              0.748, 0.729,
              0.872, 0.837),
  stringsAsFactors = FALSE
)

df$DeLong <- c("p=0.64", "", "p<0.001", "", "p=0.04", "")
df$Fold <- factor(df$Fold, levels=rev(unique(df$Fold)))
df$Model <- factor(df$Model, levels=c("SRE refitted", "CALON-Lite"))

pd <- position_dodge(width=0.55)
cols <- c("CALON-Lite" = "#1B4F72", "SRE refitted" = "#C0392B")

p1 <- ggplot(df, aes(x=AUC, y=Fold, colour=Model, shape=Model)) +
  geom_vline(xintercept=c(0.7, 0.8), linetype="dotted", colour="grey70", linewidth=0.3) +
  geom_pointrange(aes(xmin=CI_low, xmax=CI_high),
                  position=pd, linewidth=0.6, size=0.6) +
  geom_text(aes(x=CI_high + 0.008, label=DeLong),
            hjust=0, size=2, colour="grey30", fontface="italic",
            position=pd, show.legend=FALSE) +
  scale_colour_manual(values=cols) +
  scale_shape_manual(values=c("CALON-Lite"=16, "SRE refitted"=17)) +
  scale_x_continuous(limits=c(0.64, 0.92), breaks=seq(0.65, 0.90, 0.05),
                     name="C-statistic (95% CI)") +
  labs(title="CALON-Lite vs SRE: LOCO-CV (N=4,448; 835 events)",
       colour=NULL, shape=NULL) +
  theme_minimal(base_size=8) +
  theme(
    plot.title = element_text(face="bold", size=7.5, hjust=0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=6.5, lineheight=0.85),
    axis.title.x = element_text(size=7),
    axis.text.x = element_text(size=6),
    legend.position = "bottom",
    legend.text = element_text(size=6),
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(-5, 0, 0, 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(3, 8, 2, 3, "pt")
  )

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_fig1.jpg",
       p1, width=10.5, height=4.5, units="cm", dpi=300, bg="white")

cat("Figure 1 v2 saved\n")
