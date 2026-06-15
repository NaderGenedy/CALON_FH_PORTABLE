# =============================================================================
# HEART UK Abstract - Figure 1
# LOCO-CV Discrimination: CALON-Lite vs SAFEHEART-RE (Forest Plot)
# =============================================================================

library(ggplot2)
library(grid)

# -- Data from manuscript LOCO-CV results --
df <- data.frame(
  Fold = rep(c("Fold 1\n(Hold-out: South Wales, n=418)",
               "Fold 2\n(Hold-out: UK Biobank, n=1,625)",
               "Fold 3\n(Hold-out: Wales, n=2,405)"), each=2),
  Model = rep(c("CALON-Lite", "SRE (refitted)"), 3),
  AUC = c(0.790, 0.768,   # SW holdout
          0.722, 0.700,   # UKB holdout
          0.855, 0.814),  # Wales holdout
  CI_low = c(0.741, 0.716,
             0.694, 0.671,
             0.837, 0.791),
  CI_high = c(0.839, 0.820,
              0.748, 0.729,
              0.872, 0.837),
  DeLong_p = c("p=0.644", "p=0.040", "p<0.001"),
  stringsAsFactors = FALSE
)

# Reverse factor order so Fold 1 is at top
df$Fold <- factor(df$Fold, levels=rev(unique(df$Fold)))
df$Model <- factor(df$Model, levels=c("SRE (refitted)", "CALON-Lite"))

# Dodge position
pd <- position_dodge(width=0.6)

# Colour palette
cols <- c("CALON-Lite" = "#1B4F72", "SRE (refitted)" = "#E74C3C")

p1 <- ggplot(df, aes(x=AUC, y=Fold, colour=Model, shape=Model)) +
  geom_vline(xintercept=0.7, linetype="dashed", colour="grey60", linewidth=0.4) +
  geom_vline(xintercept=0.8, linetype="dashed", colour="grey60", linewidth=0.4) +
  geom_pointrange(aes(xmin=CI_low, xmax=CI_high),
                  position=pd, linewidth=0.7, size=0.8, fatten=4) +
  geom_text(data=df[df$Model=="CALON-Lite",],
            aes(x=CI_high + 0.012, label=DeLong_p),
            hjust=0, size=2.5, colour="grey30", fontface="italic",
            position=pd, show.legend=FALSE) +
  scale_colour_manual(values=cols) +
  scale_shape_manual(values=c("CALON-Lite"=16, "SRE (refitted)"=17)) +
  scale_x_continuous(limits=c(0.62, 0.92), breaks=seq(0.65, 0.90, 0.05),
                     name="C-statistic (95% CI)") +
  labs(title="CALON-Lite vs SAFEHEART-RE: LOCO-CV Discrimination",
       subtitle="Leave-one-cohort-out cross-validation across 3 genetically confirmed FH cohorts (N=4,448)",
       colour=NULL, shape=NULL) +
  theme_minimal(base_size=9) +
  theme(
    text = element_text(family="sans"),
    plot.title = element_text(face="bold", size=9, hjust=0),
    plot.subtitle = element_text(size=7, colour="grey40", hjust=0),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=7.5, lineheight=0.9),
    axis.title.x = element_text(size=8),
    legend.position = "bottom",
    legend.text = element_text(size=7),
    legend.key.size = unit(0.35, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 12, 5, 5, "pt")
  )

ggsave("C:/Users/nader/Downloads/calon_ukb_pipeline/heartuk_fig1.jpg",
       p1, width=10.5, height=5.5, units="cm", dpi=300, bg="white")

cat("Figure 1 saved successfully\n")
