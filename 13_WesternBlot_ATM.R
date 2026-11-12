#-----------------------------------------
# binding Data
#-----------------------------------------
celllines <- c("HCC78", "HCC827", "HCC44", "NCI-H23", "A549")
# ATM density, measured on Fiji
WB_intens_ATM_raw <- data.frame(
  Cellline = celllines,
  Exp1 = c(0.9118, 1.5197, 0.9773, 0.0221, 0.4787),
  Exp2 = c(0.7158, 1.1931, 0.7672, 0.0174, 0.3758),
  Exp3 = c(1.3632, 2.2720, 1.4611, 0.0330, 0.7156),
  Exp4 = c(0.6907, 1.0261, 0.6569, 0.0139, 0.2467)
)

WB_intens_ATM_long <- WB_intens_ATM_raw %>%
  pivot_longer(
    cols = starts_with("Exp"),
    names_to = "Experiment",
    values_to = "Intensity"
  ) %>%
  mutate(
    Cellline = factor(Cellline, levels = celllines),
    Experiment = factor(Experiment)
  )

# Mean, SD und SEM
WB_intens_ATM_summary <- WB_intens_ATM_long %>%
  group_by(Cellline) %>%
  summarise(
    mean_ATM = mean(Intensity),
    sd_ATM = sd(Intensity),
    sem_ATM = sd_ATM / sqrt(n()),
    .groups = "drop"
  )
#-----------------------------------------
# import SF4Gy and merging to density data
#-----------------------------------------
setwd(file.path(dirs$results, "04_ColonyFormingAssay_8Gy"))

CFA_8Gy_plot_data <- read.csv2("04_CFA_5Celllines_8Gy_means.csv",sep=",")

WB_SF_all_4Gy<-CFA_8Gy_plot_data%>%
  filter(
    Dose==4
  )%>%
  mutate(
    SF_mean=ifelse(is.na(SF_mean),"",SF_mean)
  )

WB_ATM_SF_Cor <- data.frame(
  Cellline = WB_SF_all_4Gy$Cellline,
  SF = as.numeric(WB_SF_all_4Gy$SF_mean),
  SEM_SF = as.numeric(WB_SF_all_4Gy$SF_SEM)
) %>%
  left_join(
    WB_intens_ATM_summary,
    by = "Cellline"
  ) %>%
  mutate(
    Cellline = factor(Cellline, levels = celllines)
  )
#-----------------------------------------
# statistical analysis
#-----------------------------------------
WB_intens_ATM_long$Cellline <- relevel(
  WB_intens_ATM_long$Cellline,
  ref = "NCI-H23"
)

WB_mod_ATM_vs_SF <- nlme::lme(
  Intensity ~ Cellline,
  random = ~ 1 | Experiment,
  data = WB_intens_ATM_long,
  method = "REML"
)

WB_ATM_vs_SF_anova_all<-anova(WB_mod_ATM_vs_SF)
summary(WB_mod_ATM_vs_SF)
WB_emm_ATM_vs_SF <- emmeans(
  WB_mod_ATM_vs_SF,
  ~ Cellline
)

WB_ATM_vs_SF_pairs <- pairs(
  WB_emm_ATM_vs_SF,
  adjust = "tukey"
)

WB_ATM_vs_SF_post_hoc_all<-summary(WB_ATM_vs_SF_pairs)
#-----------------------------------------
# statistical analysis vs NCI-H23
#-----------------------------------------
WB_ATM_vs_SF_Dunnett <- contrast(
  WB_emm_ATM_vs_SF,
  method = "trt.vs.ctrl",
  ref = "NCI-H23",
  adjust = "dunnett"
)
#-----------------------------------------
# correlation ATM ~ SF_4Gy
#-----------------------------------------
WB_ATM_vs_SF_spearman <- cor.test(
  WB_ATM_SF_Cor$SF,
  WB_ATM_SF_Cor$mean_ATM,
  method = "spearman",
  exact = FALSE
)

WB_ATM_vs_SF_lm <- lm(
  mean_ATM ~ SF,
  data = WB_ATM_SF_Cor
)

WB_ATM_vs_SF_stats <- data.frame(
  rho = unname(WB_ATM_vs_SF_spearman$estimate),
  p = WB_ATM_vs_SF_spearman$p.value,
  R2 = summary(WB_ATM_vs_SF_lm)$r.squared,
  n = nrow(WB_ATM_SF_Cor)
)
#-----------------------------------------
# preparing for plotting
#-----------------------------------------
WB_ATM_colours <- c(
  "HCC44"   = "#EA4025",
  "A549"    = "#5F8F29",
  "HCC78"   = "#4294F7",
  "HCC827"  = "#F09837",
  "NCI-H23" = "#872111"
)

WB_ATM_significance <- as.data.frame(
  summary(WB_ATM_vs_SF_Dunnett)
) %>%
  mutate(
    Cellline = sub(" - \\(NCI-H23\\)", "", contrast),
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  ) %>%
  dplyr::select(Cellline, significance) %>%
  bind_rows(
    data.frame(
      Cellline = "NCI-H23",
      significance = "ref."
    )
  ) %>%
  mutate(
    Cellline = factor(
      Cellline,
      levels = c("HCC78", "HCC827", "HCC44", "NCI-H23", "A549")
    )
  )
#-----------------------------------------
# bar plot ATM intensity ~ Cell line
#----------------------------------------
WB_p_ATM_vs_SF_intensity <- ggplot(
  WB_intens_ATM_summary,
  aes(
    x = Cellline,
    y = mean_ATM,
    fill = Cellline
  )
) +
  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.75
  ) +
  geom_errorbar(
    aes(
      ymin = mean_ATM - sem_ATM,
      ymax = mean_ATM + sem_ATM
    ),
    width = 0.2,
    linewidth = 0.75,
    colour = "black"
  ) +
  geom_point(
    data = WB_intens_ATM_long,
    aes(
      x = Cellline,
      y = Intensity
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.5,
    stroke = 0.4,
    fill = "black",
    colour = "black",
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 123
    )
  ) +
  scale_fill_manual(
    values = WB_ATM_colours
  ) +
  scale_y_continuous(
    breaks = c(0, 1, 2),
    minor_breaks = c(0.5, 1.5),
    limits = c(0, 2.4),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Normalized ATM intensity (a.u.)"
  ) +
   theme_prism_xcell() +
  theme(
    axis.title = element_text(
      size = 8,
      face = "bold"
    ),
    axis.text = element_text(
      size = 7,
      colour = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "none"
  )

WB_p_ATM_vs_SF_intensity <- WB_p_ATM_vs_SF_intensity +
  geom_text(
    data = WB_ATM_significance,
    aes(
      x = Cellline,
      y = 2.3,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.5
  )
#----------------------------------------
# correlation plot ATM content ~ SF_4Gy
#----------------------------------------
WB_p_ATM_vs_SF_correlation <- ggplot(
  WB_ATM_SF_Cor,
  aes(
    x = SF,
    y = mean_ATM
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "black",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = SF - SEM_SF,
      xmax = SF + SEM_SF
    ),
    width = 0.1,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_ATM - sem_ATM,
      ymax = mean_ATM + sem_ATM
    ),
    width = 0.1,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_point(
    aes(fill = Cellline),
    shape = 21,
    size = 3,
    stroke = 0.6,
    colour = "black"
  ) +
  geom_text(
    aes(label = Cellline),
    size = 2.5,
    hjust = -0.15,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 0.01,
    y = 1.7,
    label = paste0(
      "\u03c1 = ", sprintf("%.2f", WB_ATM_vs_SF_stats$rho),
      "\np = ", format.pval(WB_ATM_vs_SF_stats$p, digits = 2)
    ),
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  scale_fill_manual(
    values = WB_ATM_colours
  ) +
  scale_x_continuous(
    breaks = c(0, 0.1,0.2,0.3),
    minor_breaks = c(seq(0,0.3,0.05)),
    limits = c(0, 0.30),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1, 1.5),
    minor_breaks = seq(0.25, 1.75, 0.5),
    limits = c(0, 1.8),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Surviving fraction (4 Gy)",
    y = "Normalized ATM intensity (a.u.)"
  ) +
  theme_prism_xcell() +
  theme(
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      colour = "black"
    ),
    legend.position = "none"
  )
#----------------------------------------
# saving results
#----------------------------------------
dir.create(
  file.path(dirs$results, "13_WesternBlot_ATM"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "13_WesternBlot_ATM"))

#1. saving graph barplot intensity
ggsave(
  filename = "13_WesternBlot_ATM_intensity.jpeg",
  plot = WB_p_ATM_vs_SF_intensity,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "13_WesternBlot_ATM_intensity.pdf",
  plot = WB_p_ATM_vs_SF_intensity,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. saving graphs correlation
ggsave(
  filename = "13_WesternBlot_ATM_4Gy_correlation.jpeg",
  plot = WB_p_ATM_vs_SF_correlation,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "13_WesternBlot_ATM_4Gy_correlation.pdf",
  plot = WB_p_ATM_vs_SF_correlation,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#3. saving graph combined intensity and correlation plot
WB_p_ATM_vs_SF_final <- (
  WB_p_ATM_vs_SF_correlation +
    WB_p_ATM_vs_SF_intensity
) +
  plot_layout(
    widths = c(2, 1)
  )

ggsave(
  filename = "13_WesternBlot_ATM_4Gy_intensity_final.jpeg",
  plot = WB_p_ATM_vs_SF_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "13_WesternBlot_ATM_4Gy_intensity_final.pdf",
  plot = WB_p_ATM_vs_SF_final,
  width = 21*2/3,
  height = 7,
  units = "cm"
)  

#4. saving aov
write.csv(
  WB_ATM_vs_SF_anova_all,
  file = file.path(dirs$results, "13_WesternBlot_ATM", "13_WesternBlot_ATM_vs_SF_anova_all.csv"),
  row.names = FALSE)
  
#5. saving Duunnet pairs with NCI-H23 as reference
write.csv(
  WB_ATM_vs_SF_Dunnett,
  file = file.path(dirs$results, "13_WesternBlot_ATM", "13_WesternBlot_ATM_vs_SF_Dunnett.csv"),
  row.names = FALSE
)

#5. saving post hoc analysis
write.csv(
  WB_ATM_vs_SF_post_hoc_all,
  file = file.path(dirs$results, "13_WesternBlot_ATM", "13_WesternBlot_ATM_vs_SF_post_hoc_all.csv"),
  row.names = FALSE
)

#5. saving post hoc analysis
write.csv(
  WB_intens_ATM_summary,
  file = file.path(dirs$results, "13_WesternBlot_ATM", "13_WesternBlot_ATM_intensity_summary.csv"),
  row.names = FALSE
)
#-----------------------------------------------
setwd(dirs$code)
