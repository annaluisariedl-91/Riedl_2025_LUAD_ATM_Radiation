#------------------------------------------
#Set directory
#------------------------------------------
setwd(file.path(dirs$CFA,"Exp2_DoseResponse2Gy","raw")) #replace with path on your own PC
#------------------------------------------
#loading raw data
#------------------------------------------
CFA_2Gy_HCC44_files <- c(
  "CFA_1_2Gy_HCC44_rep1.csv",
  "CFA_1_2Gy_HCC44_rep2.csv",
  "CFA_1_2Gy_HCC44_rep3.csv"
)

CFA_2Gy_HCC44_rep <- lapply(
  CFA_2Gy_HCC44_files,
  function(x) analyse_CFA_rep(read.csv(x))
)

CFA_2Gy_H23_files <- c(
  "CFA_1_2Gy_H23_rep1.csv",
  "CFA_1_2Gy_H23_rep2.csv",
  "CFA_1_2Gy_H23_rep3.csv",
  "CFA_1_2Gy_H23_rep4.csv"
)

CFA_2Gy_H23_rep <- lapply(
  CFA_2Gy_H23_files,
  function(x) analyse_CFA_rep(read.csv(x))
)
#------------------------------------------
# preparing for modelling
#------------------------------------------
CFA_2Gy_stat_data <- bind_rows(
    make_CFA_data(CFA_2Gy_H23_rep, "NCI-H23") %>%
      mutate(Cellline = "NCI-H23"),
    
    make_CFA_data(CFA_2Gy_HCC44_rep, "HCC44") %>%
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

CFA_2Gy_stat_data$Cellline <- relevel(
  CFA_2Gy_stat_data$Cellline,
  ref = "HCC44"
)
#------------------------------------------
# linear-quadratic fit
#------------------------------------------
CFA_2Gy_mean_data<- CFA_2Gy_stat_data %>%
  group_by(Cellline, Dose) %>%
  summarise(
    mean_log_sf = mean(log_sf_adj, na.rm = TRUE),
    .groups = "drop"
  )

CFA_2Gy_coefficients <- CFA_2Gy_mean_data %>%
  group_split(Cellline) %>%
  setNames(unique(CFA_2Gy_mean_data$Cellline)) %>%
  lapply(fit_lq_model)
#------------------------------------------
# Welsh t Test
#------------------------------------------
CFA_2Gy_t.test <- CFA_2Gy_stat_data %>%
  filter(Dose == 2, sf_adj > 0) %>%
  mutate(log_SF = log(sf_adj)) %>%
  {
    t.test(log_SF ~ Cellline, data = .)
  }

CFA_2Gy_res_p <- CFA_2Gy_t.test$p.value

CFA_2Gy_sig_label <- paste0(
  "SF(2Gy): HCC44 vs. NCI-H23,\n ",
  dplyr::case_when(
    CFA_2Gy_res_p < 0.001 ~ "***",
    CFA_2Gy_res_p < 0.01  ~ "**",
    CFA_2Gy_res_p < 0.05  ~ "*",
    TRUE      ~ "n.s."
  ),
  " (p = ", format.pval(CFA_2Gy_res_p, digits = 2, eps = 0.001), ")"
)
#------------------------------------------
# binding for plotting
#------------------------------------------
CFA_2Gy_plot_data <- bind_rows(
  make_stat_CFA(CFA_2Gy_H23_rep, "NCI-H23"),
  make_stat_CFA(CFA_2Gy_HCC44_rep, "HCC44")
) %>%
  group_by(Cellline, Dose) %>%
  summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    .groups = "drop"
  )

CFA_2Gy_dseq <- seq(0, 2, length.out = 200)
CFA_2Gy_LQ_data <- bind_rows(
  data.frame(
    Cellline = "HCC44",
    Dose = CFA_2Gy_dseq,
    SF = exp(
      -(CFA_2Gy_coefficients[["HCC44"]]["alpha"] * CFA_2Gy_dseq +
          CFA_2Gy_coefficients[["HCC44"]]["beta"] * CFA_2Gy_dseq^2)
    )
  ),
  data.frame(
    Cellline = "H23",
    Dose = CFA_2Gy_dseq,
    SF = exp(
      -(CFA_2Gy_coefficients[["NCI-H23"]]["alpha"] * CFA_2Gy_dseq +
          CFA_2Gy_coefficients[["NCI-H23"]]["beta"] * CFA_2Gy_dseq^2)
    )
  )
)

CFA_2Gy_colours <- c(
  "HCC44"  = "#0196FF",
  "NCI-H23"    = "#4C8F00")

CFA_2Gy_shapes <- c(
  "HCC44" = 22,
  "NCI-H23"   = 21
)
#------------------------------------------
# plotting 
#------------------------------------------
pd <- position_dodge(width = 0.10)

CFA_2Gy_plot <- ggplot() +
  # mean value
  geom_point(
    data = CFA_2Gy_plot_data,
    aes(
      x = Dose,
      y = SF_mean,
      fill = Cellline,
      shape = Cellline,
      group = Cellline
    ),
    position = pd,
    size = 2,
    stroke = 0.6,
    colour = "black"
  ) +
  # Linear-quadratic fits
  geom_line(
    data = CFA_2Gy_LQ_data,
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
    data = CFA_2Gy_plot_data,
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
  scale_fill_manual(
    values = CFA_2Gy_colours,
    breaks = c("HCC44", "NCI-H23"),
    labels = c("HCC44", "NCI-H23")
  ) +
  scale_shape_manual(
    values = CFA_2Gy_shapes,
    breaks = c("HCC44", "NCI-H23"),
    labels = c("HCC44", "NCI-H23")
  ) +
  # logarithmic y axis
  scale_y_log10(
    limits = c(0.01, 2),
    breaks = c(1e-0,1e-1, 1e-2, 1e-3),
    labels = expression(10^0,10^-1, 10^-2, 10^-3)
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
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      colour = "black"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(
      face = "plain",
      size = 10,
      colour = "black"
    ),
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

CFA_2Gy_plot <- CFA_2Gy_plot +
  annotate(
    "text",
    x = 0.1,
    y = 0.02,
    label = CFA_2Gy_sig_label,
    size = 2.5,
    hjust=0
  )
#------------------------------------------
# saving results
#------------------------------------------
dir.create(
  file.path(dirs$results, "05_ColonyFormingAssay_2Gy"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "05_ColonyFormingAssay_2Gy"))

# 1. saving figure
ggsave(
  filename = "05_CFA_2Gy_HCC44_H23_plot.jpeg",
  plot = CFA_2Gy_plot,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename ="05_CFA_2Gy_HCC44_H23_plot.pdf",
  plot = CFA_2Gy_plot,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 2. Welsh t-test
CFA_2Gy_t_test_results <- data.frame(
  Comparison = "HCC44 vs NCI-H23",
  t_value = unname(CFA_2Gy_t.test$statistic),
  df = unname(CFA_2Gy_t.test$parameter),
  p_value = CFA_2Gy_t.test$p.value,
  CI_lower = CFA_2Gy_t.test$conf.int[1],
  CI_upper = CFA_2Gy_t.test$conf.int[2],
  Mean_HCC44 = unname(CFA_2Gy_t.test$estimate[1]),
  Mean_H23 = unname(CFA_2Gy_t.test$estimate[2])
)
write.csv(
  CFA_2Gy_t_test_results,
  file = file.path(dirs$results, "05_ColonyFormingAssay_2Gy", "05_CFA_2Gy_HCC44_H23_ANOVA.csv"),
  row.names = FALSE
)

# 3. saving LQ-coeffizient 
CFA_2Gy_LQ_results <- do.call(
  rbind,
  lapply(names(CFA_2Gy_coefficients), function(cell) {
    data.frame(
      Cellline = cell,
      alpha = unname(CFA_2Gy_coefficients[[cell]]["alpha"]),
      beta  = unname(CFA_2Gy_coefficients[[cell]]["beta"])
    )
  })
)

rownames(CFA_2Gy_LQ_results) <- NULL

write.csv(
  CFA_2Gy_LQ_results,
  file = file.path(dirs$results, "05_ColonyFormingAssay_2Gy", "05_CFA_2Gy_HCC44_H23_LQ_coefficients.csv"),
  row.names = FALSE
)
#------------------------------------------
setwd(dirs$code)
