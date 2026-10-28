#------------------------------------------
#Set directory
#------------------------------------------
setwd(file.path(base_dir,"functional_assays","CFA","Exp2_DoseResponse2Gy","raw")) #replace with path on your own PC
#------------------------------------------
#loading raw data
#------------------------------------------
HCC44_2Gy_files <- c(
  "CFA_1_2Gy_HCC44_rep1.csv",
  "CFA_1_2Gy_HCC44_rep2.csv",
  "CFA_1_2Gy_HCC44_rep3.csv"
)

HCC44_2Gy_rep <- lapply(
  HCC44_2Gy_files,
  function(x) analyse_CFA_rep(read.csv(x))
)

H23_2Gy_files <- c(
  "CFA_1_2Gy_H23_rep1.csv",
  "CFA_1_2Gy_H23_rep2.csv",
  "CFA_1_2Gy_H23_rep3.csv",
  "CFA_1_2Gy_H23_rep4.csv"
)

H23_2Gy_rep <- lapply(
  H23_2Gy_files,
  function(x) analyse_CFA_rep(read.csv(x))
)
#------------------------------------------
# preparing for modelling
#------------------------------------------
stat_data_CFA_2Gy <- bind_rows(
    make_CFA_data(H23_2Gy_rep, "NCI-H23") %>%
      mutate(Cellline = "NCI-H23"),
    
    make_CFA_data(HCC44_2Gy_rep, "HCC44") %>%
      mutate(Cellline = "HCC44")
  ) %>%
  filter(!is.na(SF), !is.nan(SF)) %>%
  mutate(
    Cellline = factor(Cellline),
    experiment = factor(experiment),  
    exp_id = interaction(Cellline, experiment, drop = TRUE),
    sf_adj = pmax(SF, 1e-8),
    log_sf_adj = log(sf_adj),
    dose_f = factor(Dose)
  )

stat_data_CFA_2Gy$Cellline <- relevel(
  stat_data_CFA_2Gy$Cellline,
  ref = "HCC44"
)
#------------------------------------------
# ANOVA
#------------------------------------------
mod_CFA_2Gy <- lme(
  fixed = log_sf_adj ~ dose_f + Cellline,
  random = ~ 1 | exp_id,
  data = stat_data_CFA_2Gy,
  method = "REML"
)

summary(mod_CFA_2Gy)

anova_CFA_2Gy <- anova(mod_CFA_2Gy)
#------------------------------------------
# Post-hoc comparisons 
#------------------------------------------
emm_CFA_2Gy <- emmeans(
  mod_CFA_2Gy,
  ~ Cellline | dose_f
)

posthoc_CFA_2Gy <- pairs(
  emm_CFA_2Gy,
  adjust = "tukey"
)
#------------------------------------------
# linear-quadratic fit
#------------------------------------------
mean_data_CFA_2Gy <- stat_data_CFA_2Gy %>%
  group_by(Cellline, Dose) %>%
  summarise(
    mean_log_sf = mean(log_sf_adj, na.rm = TRUE),
    .groups = "drop"
  )

coefficients_CFA_2Gy <- mean_data_CFA_2Gy %>%
  group_split(Cellline) %>%
  setNames(unique(mean_data_CFA_2Gy$Cellline)) %>%
  lapply(fit_lq_model)
#------------------------------------------
# binding for plotting
#------------------------------------------
CFA_plot_data_2Gy <- bind_rows(
  make_stat_CFA(H23_2Gy_rep, "NCI-H23"),
  make_stat_CFA(HCC44_2Gy_rep, "HCC44")
) %>%
  group_by(Cellline, Dose) %>%
  summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    .groups = "drop"
  )

dseq <- seq(0, 2, length.out = 200)

LQ_data_CFA_2Gy <- bind_rows(
  data.frame(
    Cellline = "HCC44",
    Dose = dseq,
    SF = exp(
      -(coefficients_CFA_2Gy[["HCC44"]]["alpha"] * dseq +
          coefficients_CFA_2Gy[["HCC44"]]["beta"] * dseq^2)
    )
  ),
  data.frame(
    Cellline = "H23",
    Dose = dseq,
    SF = exp(
      -(coefficients_CFA_2Gy[["NCI-H23"]]["alpha"] * dseq +
          coefficients_CFA_2Gy[["NCI-H23"]]["beta"] * dseq^2)
    )
  )
)

CFA_colours_2Gy <- c(
  "HCC44"  = "#4294F7",
  "NCI-H23"    = "#5F8F29"
)

CFA_shapes_2Gy <- c(
  "HCC44" = 22,
  "NCI-H23"   = 21
)
#------------------------------------------
# plotting 
#------------------------------------------
pd <- position_dodge(width = 0.10)

CFA_plot_2Gy <- ggplot() +
  # Linear-quadratic fits
  geom_line(
    data = LQ_data_CFA_2Gy,
    aes(
      x = Dose,
      y = SF,
      group = Cellline
    ),
    colour = "black",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  # SEM
  geom_errorbar(
    data = CFA_plot_data_2Gy,
    aes(
      x = Dose,
      ymin = pmax(SF_mean - SF_SEM, 0.001),
      ymax = SF_mean + SF_SEM,
      group = Cellline
    ),
    position = pd,
    width = 0.04,
    linewidth = 0.5,
    colour = "black"
  ) +
  # mean value
  geom_point(
    data = CFA_plot_data_2Gy,
    aes(
      x = Dose,
      y = SF_mean,
      fill = Cellline,
      shape = Cellline,
      group = Cellline
    ),
    position = pd,
    size = 2.5,
    stroke = 0.6,
    colour = "black"
  ) +
  scale_fill_manual(
    values = CFA_colours_2Gy,
    breaks = c("HCC44", "NCI-H23"),
    labels = c("HCC44", "NCI-H23")
  ) +
  scale_shape_manual(
    values = CFA_shapes_2Gy,
    breaks = c("HCC44", "NCI-H23"),
    labels = c("HCC44", "NCI-H23")
  ) +
  # logarithmic y axis
  scale_y_log10(
    limits = c(0.01, 2),
    breaks = c(0.01, 0.1, 1),
    labels = scales::label_math(10^.x)
  ) +
  scale_x_continuous(
    limits = c(-0.05, 2.1),
    breaks = c(0, 0.5, 1, 1.5, 2),
    minor_breaks = seq(0, 2, 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = "Dose [Gy]",
    y = "Surviving Fraction"
  ) +
  ggprism::theme_prism(base_size = 8) +
  theme(
    axis.title = element_text(
      size = 8,
      face = "bold"
    ),
    axis.text = element_text(
      size = 7,
      colour = "black"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 2)
    ),
    shape = guide_legend(
      override.aes = list(size = 2)
    )
  )

#------------------------------------------
# saving results
#------------------------------------------
dir.create(
  file.path(dirs$results, "05_CFA_2Gy"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "05_CFA_2Gy"))

# 1. saving figure
ggsave(
  filename = file.path(dirs$results, "05_CFA_2Gy", "CFA_2Gy_gplot.jpeg"),
  plot = CFA_plot_2Gy,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(dirs$results, "05_CFA_2Gy", "CFA_2Gy_gplot.pdf"),
  plot = CFA_plot_2Gy,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 2. saving ANOVA 
ANOVA_results_CFA_2Gy <- as.data.frame(
  anova(mod_CFA_2Gy)
)

ANOVA_results_CFA_2Gy$Effect <- rownames(ANOVA_results_CFA_2Gy)
rownames(ANOVA_results_CFA_2Gy) <- NULL

ANOVA_results_CFA_2Gy <- ANOVA_results_CFA_2Gy %>%
  dplyr::select(Effect, dplyr::everything())

write.csv(
  ANOVA_results_CFA_2Gy,
  file = file.path(dirs$results, "05_CFA_2Gy", "CFA_2Gy_ANOVA.csv"),
  row.names = FALSE
)

# 3. saving Post-hoc-Analysis 
posthoc_results_CFA_2Gy <- as.data.frame(
  summary(posthoc_CFA_2Gy)
)

write.csv(
  posthoc_results_CFA_2Gy,
  file = file.path(dirs$results, "05_CFA_2Gy", "CFA_2Gy_posthoc.csv"),
  row.names = FALSE
)

# 4. saving LQ-coeffizient 
LQ_results_CFA_2Gy <- do.call(
  rbind,
  lapply(names(coefficients_CFA_2Gy), function(cell) {
    
    data.frame(
      Cellline = cell,
      alpha = unname(coefficients_CFA_2Gy[[cell]]["alpha"]),
      beta  = unname(coefficients_CFA_2Gy[[cell]]["beta"])
    )
  })
)

rownames(LQ_results_CFA_2Gy) <- NULL

write.csv(
  LQ_results_CFA_2Gy,
  file = file.path(dirs$results, "05_CFA_2Gy", "CFA_2Gy_LQ_coefficients.csv"),
  row.names = FALSE
)
#------------------------------------------
setwd(dirs$code)
