#------------------------------------------
# Set directory
#------------------------------------------
setwd(file.path(base_dir,"functional_assays","CFA","Exp1_DoseResponse8Gy","raw")) #replace with path on your own PC
Dose <-c(0,2,4,6,8)
#------------------------------------------
# preparing NCI-H23
#------------------------------------------
H23_CFA_8Gy<- read.csv("CFA_1_8Gy_H23.csv",sep=",")

H23_rep_CFA_8Gy <- list(
  analyse_CFA_rep(H23_CFA_8Gy[2:31,1:3])
)
#------------------------------------------
# preparing HCC44
#------------------------------------------
HCC44_CFA_8Gy<- read.csv("CFA_1_8Gy_HCC44.csv")

HCC44_rep_CFA_8Gy <- list(
  analyse_CFA_rep(HCC44_CFA_8Gy[2:31,1:3]),
  analyse_CFA_rep(HCC44_CFA_8Gy[2:31,7:9]),
  analyse_CFA_rep(HCC44_CFA_8Gy[2:31,13:15]),
  analyse_CFA_rep(HCC44_CFA_8Gy[2:31,19:21])
)
#------------------------------------------
# preparing A594 
#------------------------------------------
A549_CFA_8Gy<- read.csv("CFA_1_8Gy_A549.csv")

A549_rep_CFA_8Gy <- list(
  analyse_CFA_rep(A549_CFA_8Gy[2:21,1:3]),
  analyse_CFA_rep(A549_CFA_8Gy[2:21,7:9]),
  analyse_CFA_rep(A549_CFA_8Gy[2:21,13:15]),
  analyse_CFA_rep(A549_CFA_8Gy[2:21,13:15])
)
#------------------------------------------
# preparing HCC827 
#------------------------------------------
HCC827_CFA_8Gy<- read.csv("CFA_1_8Gy_HCC827.csv")

HCC827_rep_CFA_8Gy <- list(
  analyse_CFA_rep(HCC827_CFA_8Gy[2:31,1:3]),
  analyse_CFA_rep(HCC827_CFA_8Gy[2:31,7:9]),
  analyse_CFA_rep(HCC827_CFA_8Gy[2:31,13:15])
)
#------------------------------------------
# preparing HCC78
#------------------------------------------
HCC78_CFA_8Gy<- read.csv("CFA_1_8Gy_HCC78.csv")

HCC78_rep_CFA_8Gy <- list(
  analyse_CFA_rep(HCC78_CFA_8Gy[2:12,1:3]),
  analyse_CFA_rep(HCC78_CFA_8Gy[2:31,7:9]),
  analyse_CFA_rep(HCC78_CFA_8Gy[2:31,13:15]),
  analyse_CFA_rep(HCC78_CFA_8Gy[2:31,19:21])
)
#------------------------------------------
# binding data
#------------------------------------------
stat_data_CFA_8Gy <- bind_rows(
  data.frame(
    Cellline = "NCI-H23",
    Dose = c(0, 2, 4, 6, 8),
    SF = H23_rep_CFA_8Gy[[1]][["summary"]][["SF"]],
    experiment = factor(1)
  ),
  make_stat_CFA(HCC44_rep_CFA_8Gy,  "HCC44"),
  make_stat_CFA(A549_rep_CFA_8Gy,   "A549"),
  make_stat_CFA(HCC78_rep_CFA_8Gy,  "HCC78"),
  make_stat_CFA(HCC827_rep_CFA_8Gy, "HCC827")
) %>%
  mutate(
    Cellline = factor(Cellline),
    sf_adj = pmax(SF, 1e-8),
    log_sf_adj = log(sf_adj),
    dose_f = as.numeric(Dose)
  )
#------------------------------------------
# ANOVA
#------------------------------------------
anova_data_CFA_8Gy <- stat_data_CFA_8Gy %>%
  filter(Dose != 6) %>%
  mutate(
    Cellline = relevel(Cellline, ref = "HCC44"),
    exp_id = interaction(Cellline, experiment, drop = TRUE)
  )

mod_simple_CFA_8Gy <- lme(
  log_sf_adj ~ dose_f * Cellline,
  random = ~ 1 | exp_id,
  data = anova_data_CFA_8Gy,
  method = "REML"
)

summary(mod_simple_CFA_8Gy)
anova(mod_simple_CFA_8Gy)

significance_CFA_8Gy<-emmeans(mod_simple_CFA_8Gy, pairwise ~ Cellline | dose_f, adjust = "tukey")
#------------------------------------------
# linear-quadratic model
#------------------------------------------
mean_data_CFA_8Gy <- stat_data_CFA_8Gy %>%
  group_by(Cellline, Dose) %>%
  summarise(
    mean_log_sf = mean(log_sf_adj, na.rm = TRUE),
    .groups = "drop"
  )

lq_fits_CFA_8Gy <- mean_data_CFA_8Gy %>%
  group_split(Cellline) %>%
  setNames(unique(mean_data_CFA_8Gy$Cellline)) %>%
  lapply(fit_lq_model)
#------------------------------------------
# preparing graph
#------------------------------------------
CFA_8Gy_colours <- c(
  "HCC44"  = "#4294F7",
  "NCI-H23"    = "#5F8F29",
  "HCC827" = "#F09837",
  "A549"   = "#EA4025",
  "HCC78"  = "#872111"
)

CFA_8Gy_shapes <- c(
  "HCC44"  = 22,
  "NCI-H23"    = 23,
  "HCC827" = 25,
  "A549"   = 25,
  "HCC78"  = 25
)

CFA_8Gy_plot_data <- stat_data_CFA_8Gy %>%
group_by(Cellline, Dose) %>%
summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    .groups = "drop"
  )

dseq <- seq(0, 8, length.out = 200)
#------------------------------------------
# plotting graph
#------------------------------------------
CFA_8Gy_plot <- ggplot(
  CFA_8Gy_plot_data,
  aes(x = Dose, y = SF_mean, group = Cellline)
) +
  geom_point(
    aes(fill = Cellline),
    shape = 21,
    size = 2,
    colour = "black"
  ) +
  # error bars
  geom_errorbar(
    data = CFA_8Gy_plot_data,
    aes(
      x = Dose,
      ymin = pmax(SF_mean - SF_SEM, 1e-4),
      ymax = SF_mean + SF_SEM
    ),
    width = 0.15,
    linewidth = 0.5,
    colour = "black"
  ) +
  # surving fraction
  geom_point(
    data = CFA_8Gy_plot_data,
    aes(
      x = Dose,
      y = SF_mean,
      fill = Cellline,
      shape = Cellline
    ),
    size = 2.5,
    stroke = 0.6,
    colour = "black"
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ 0 + x + I(x^2),
    se = FALSE,
    colour = "black",
    linewidth = 0.6,
    linetype = "dashed",
    n = 200
  )  +
  # colors and shapes
  scale_fill_manual(
    values = CFA_8Gy_colours,
    breaks = c("HCC44", "NCI-H23", "HCC827", "A549", "HCC78"),
    labels = c(
      "HCC44 - ref",
      "NCI-H23 - **",
      "HCC827 - n.s.",
      "A549 - ***",
      "HCC78 - ***"
    )
  ) +
  scale_shape_manual(
    values = CFA_8Gy_shapes,
    breaks = c("HCC44", "NCI-H23", "HCC827", "A549", "HCC78"),
    labels = c(
      "HCC44 - ref",
      "NCI-H23 - **",
      "HCC827 - n.s.",
      "A549 - ***",
      "HCC78 - ***"
    )
  ) +
  # logarithmic y axes
  scale_y_log10(
    limits = c(0.001, 2),
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c(expression(10^-3), "0.01", "0.1", "1")
  ) +
  scale_x_continuous(
    limits = c(0, 8.5),
    breaks = c(0, 2, 4, 6, 8),
    minor_breaks = seq(0, 8, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Dose [Gy]",
    y = "Surviving Fraction"
  ) +
ggprism::theme_prism(base_size = 6) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.spacing.y = unit(0.1, "mm"),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )+
  guides(
  fill = guide_legend(
    ncol = 3,
    byrow = TRUE,
    override.aes = list(
      size = 2,
      stroke = 0.25
    )
  ),
  shape = guide_legend(
    ncol = 3,
    byrow = TRUE,
    override.aes = list(
      size = 2,
      stroke = 0.5
    )
  )
)
#-------------------------------------------
# saving results
#-------------------------------------------
dir.create(
  file.path(dirs$results, "04_CFA_8Gy"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "04_CFA_8Gy"))

# 1. saving mean values used for plotting
write.csv(
  mean_data_CFA_8Gy,
  file = file.path(dirs$results, "04_CFA_8Gy", "CFA_8Gy_means.csv"),
  row.names = FALSE
)

# 2. saving lq modell results
write.csv(
  lq_fits_CFA_8Gy,
  file = file.path(dirs$results, "04_CFA_8Gy", "CFA_8Gy_LQ_coefficients.csv"),
  row.names = FALSE
)

# 3. saving anova results
ANOVA_results_CFA_8Gy <- as.data.frame(
  anova(mod_simple_CFA_8Gy)
)

ANOVA_results_CFA_8Gy$Effect <- rownames(ANOVA_results_CFA_8Gy)
rownames(ANOVA_results_CFA_8Gy) <- NULL

ANOVA_results_CFA_8Gy <- ANOVA_results_CFA_8Gy %>%
  dplyr::select(Effect, everything())

write.csv(
  ANOVA_results_CFA_8Gy,
  file = file.path(dirs$results, "04_CFA_8Gy", "CFA_8Gy_ANOVA.csv"),
  row.names = FALSE
)

# 5. Post-hoc-comparison
emm_CFA_8Gy <- emmeans(
  mod_simple_CFA_8Gy,
    ~ Cellline | dose_f
  )

posthoc_CFA_8Gy <- pairs(
  emm_CFA_8Gy,
  adjust = "tukey"
)

posthoc_results_CFA_8Gy <- as.data.frame(
  summary(posthoc_CFA_8Gy)
)

write.csv(
  posthoc_results_CFA_8Gy,
  file = file.path(dirs$results, "04_CFA_8Gy", "CFA_8Gy_posthoc.csv"),
  row.names = FALSE
)

# 6. saving figure
ggsave(
  filename = file.path(dirs$results, "04_CFA_8Gy", "CFA_8Gy_plot.jpeg"),
  plot = CFA_8Gy_plot,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(dirs$results, "04_CFA_8Gy","CFA_8Gy_plot.pdf"),
  plot = CFA_8Gy_plot,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)
#--------------------------------------------
setwd(dirs$code)
