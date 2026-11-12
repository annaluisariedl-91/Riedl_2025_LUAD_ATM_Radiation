#------------------------------------------
# Set directory
#------------------------------------------
setwd(file.path(dirs$CFA,"Exp1_DoseResponse8Gy","raw")) #replace with path on your own PC
#------------------------------------------
# preparing NCI-H23
#------------------------------------------
CFA_8Gy_H23<- read.csv("CFA_1_8Gy_H23.csv",sep=",")

CFA_8Gy_H23_rep <- list(
  analyse_CFA_rep(CFA_8Gy_H23[2:31,1:3])
)
#------------------------------------------
# preparing HCC44
#------------------------------------------
CFA_8Gy_HCC44<- read.csv("CFA_1_8Gy_HCC44.csv")

CFA_8Gy_HCC44_rep <- list(
  analyse_CFA_rep(CFA_8Gy_HCC44[2:31,1:3]),
  analyse_CFA_rep(CFA_8Gy_HCC44[2:31,7:9]),
  analyse_CFA_rep(CFA_8Gy_HCC44[2:31,13:15]),
  analyse_CFA_rep(CFA_8Gy_HCC44[2:31,19:21])
)
#------------------------------------------
# preparing A594 
#------------------------------------------
CFA_8Gy_A549<- read.csv("CFA_1_8Gy_A549.csv")

CFA_A549_rep_8Gy <- list(
  analyse_CFA_rep(CFA_8Gy_A549[2:21,1:3]),
  analyse_CFA_rep(CFA_8Gy_A549[2:21,7:9]),
  analyse_CFA_rep(CFA_8Gy_A549[2:21,13:15]),
  analyse_CFA_rep(CFA_8Gy_A549[2:21,13:15])
)
#------------------------------------------
# preparing HCC827 
#------------------------------------------
CFA_8Gy_HCC827<- read.csv("CFA_1_8Gy_HCC827.csv")

CFA_8Gy_HCC827_rep <- list(
  analyse_CFA_rep(CFA_8Gy_HCC827[2:31,1:3]),
  analyse_CFA_rep(CFA_8Gy_HCC827[2:31,7:9]),
  analyse_CFA_rep(CFA_8Gy_HCC827[2:31,13:15])
)
#------------------------------------------
# preparing HCC78
#------------------------------------------
CFA_8Gy_HCC78<- read.csv("CFA_1_8Gy_HCC78.csv")

CFA_8Gy_HCC78_rep <- list(
  analyse_CFA_rep(CFA_8Gy_HCC78[2:12,1:3]),
  analyse_CFA_rep(CFA_8Gy_HCC78[2:31,7:9]),
  analyse_CFA_rep(CFA_8Gy_HCC78[2:31,13:15]),
  analyse_CFA_rep(CFA_8Gy_HCC78[2:31,19:21])
)
#------------------------------------------
# binding data
#------------------------------------------
CFA_8Gy_stat_data <- bind_rows(
  data.frame(
    Cellline = "NCI-H23",
    Dose = c(0, 2, 4, 6, 8),
    SF = CFA_8Gy_H23_rep[[1]][["summary"]][["SF"]],
    experiment = factor(1)
  ),
  make_stat_CFA(CFA_8Gy_HCC44_rep,  "HCC44"),
  make_stat_CFA(CFA_A549_rep_8Gy,   "A549"),
  make_stat_CFA(CFA_8Gy_HCC78_rep,  "HCC78"),
  make_stat_CFA(CFA_8Gy_HCC827_rep, "HCC827")
) %>%
  mutate(
    Cellline = factor(Cellline),
    sf_adj = pmax(SF, 1e-8),
    log_sf_adj = log(sf_adj),
    dose_f = as.numeric(Dose)
  )
#------------------------------------------
# linear-quadratic model
#------------------------------------------
CFA_8Gy_mean_data <- CFA_8Gy_stat_data %>%
  group_by(Cellline, Dose) %>%
  summarise(
    mean_log_sf = mean(log_sf_adj, na.rm = TRUE),
    .groups = "drop"
  )

CFA_8Gy_lq_fits <- CFA_8Gy_mean_data %>%
  group_split(Cellline) %>%
  setNames(
    CFA_8Gy_mean_data %>%
      distinct(Cellline) %>%
      pull(Cellline) %>%
      as.character()
  ) %>%
  lapply(fit_lq_model)

CFA_8Gy_LQ_curves <- purrr::imap_dfr(
  CFA_8Gy_lq_fits,
  function(coef, Cellline) {
    Dose <- seq(0, 8, length.out = 200)
    data.frame(
      Cellline = Cellline,
      Dose = Dose,
      SF_pred = exp(
        -coef["alpha"] * Dose -
          coef["beta"] * Dose^2
      )
    )
  }
)
#------------------------------------------
# ANOVA
#------------------------------------------
CFA_8Gy_anova_data <- CFA_8Gy_stat_data %>%
  filter(Dose != 6) %>%
  mutate(
    Cellline = relevel(Cellline, ref = "NCI-H23"),
    exp_id = interaction(Cellline, experiment, drop = TRUE)
  )

CFA_8Gy_mod_simple <- lme(
  log_sf_adj ~ dose_f * Cellline,
  random = ~ 1 | exp_id,
  data = CFA_8Gy_anova_data,
  method = "REML"
)

CFA_8Gy_significance<-emmeans(CFA_8Gy_mod_simple, pairwise ~ Cellline | dose_f, adjust = "tukey")
#------------------------------------------
# preparing graph
#------------------------------------------
CFA_8Gy_colours <- c(
  "HCC44"  = "#FE2500",
  "NCI-H23"    = "#921000",
  "HCC827" = "#FE9301",
  "A549"   = "#4C8F00",
  "HCC78"  = "#0196FF"
)

CFA_8Gy_shapes <- c(
  "HCC44"  = 22,
  "NCI-H23"    = 23,
  "HCC827" = 25,
  "A549"   = 25,
  "HCC78"  = 25
)

CFA_8Gy_order <- c(
  "HCC44", "NCI-H23", "HCC827", "A549", "HCC78")
  
CFA_8Gy_plot_data <- CFA_8Gy_stat_data %>%
group_by(Cellline, Dose) %>%
summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    .groups = "drop"
  )%>%
  dplyr::mutate(
    SF_mean = if_else(
      Cellline == "NCI-H23" & Dose == 2,
      0.001,
      SF_mean
    )
  )

CFA_8Gy_sig <- as.data.frame(
  CFA_8Gy_significance$contrasts
) %>%
  dplyr::filter(grepl("NCI-H23", contrast)) %>%
  dplyr::mutate(
    Cellline = sub("^\\(NCI-H23\\) - ", "", contrast),
    Cellline = sub("^\\s+", "", Cellline),
    Significance = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "n.s."
    ))
  
CFA_8Gy_sig_labels <- setNames(
  paste0(
    CFA_8Gy_sig$Cellline,
    " - ",
    CFA_8Gy_sig$Significance
  ),
  CFA_8Gy_sig$Cellline
)

CFA_8Gy_sig_labels <- c(
  "NCI-H23" = "NCI-H23 - ref",
  CFA_8Gy_sig_labels
)

CFA_8Gy_sig_labels <- CFA_8Gy_sig_labels[
  match(CFA_8Gy_order, names(CFA_8Gy_sig_labels))
]
#------------------------------------------
# plotting graph
#------------------------------------------
CFA_8Gy_plot <- ggplot(
  CFA_8Gy_plot_data,
  aes(
    x = Dose,
    y = SF_mean,
    group = Cellline
  )
) +
  # Error bars
  geom_errorbar(
    aes(
      ymin = pmax(SF_mean - SF_SEM, 1e-4),
      ymax = SF_mean + SF_SEM
    ),
    width = 0.15,
    linewidth = 0.5,
    colour = "black"
  ) +
  # Surviving fraction
  geom_point(
    aes(
      fill = Cellline,
      shape = Cellline
    ),
    size = 2.5,
    stroke = 0.6,
    colour = "black"
  ) +
  # LQ curves
  geom_line(
    data = CFA_8Gy_LQ_curves,
    aes(
      x = Dose,
      y = SF_pred,
      group = Cellline
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.6,
    linetype = "dashed"
  )+
  # Colours with significance labels
  scale_fill_manual(
    values = CFA_8Gy_colours,
    breaks = CFA_8Gy_order,
    labels =  unname(CFA_8Gy_sig_labels),
    name = NULL
  ) +
  # Shapes with significance labels
  scale_shape_manual(
    values = CFA_8Gy_shapes,
    breaks = CFA_8Gy_order,
    labels =  unname(CFA_8Gy_sig_labels),
    name = NULL
  ) +
  # Logarithmic y-axis
  scale_y_log10(
    limits = c(0.001, 3),
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c(
      expression(10^-3),
      "0.01",
      "0.1",
      "1"
    )
  ) +
  # X-axis
  scale_x_continuous(
    limits = c(0, 8.5),
    breaks = c(0, 2, 4, 6, 8),
    minor_breaks = seq(0, 8, 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Dose (Gy)",
    y = "Surviving fraction"
  )  +
  ggprism::theme_prism(base_size = 6) +
  theme(
    axis.title = element_text(
    size = 12,
    face = "bold"
  ),
  axis.text = element_text(
    size = 10,
    colour = "black"
  ),
    legend.position ="bottom",
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
        stroke = 0.5
      )
    ),
    shape = guide_legend(
      ncol = 3,
      byrow = TRUE,
      override.aes = list(
        size = 1.5,
        stroke = 0.4
        )
    )
  )
#-------------------------------------------
# saving results
#-------------------------------------------
dir.create(
  file.path(dirs$results, "04_ColonyFormingAssay_8Gy"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "04_ColonyFormingAssay_8Gy"))

# 1. saving mean values used for plotting
write.csv(
  CFA_8Gy_plot_data,
  file = file.path(dirs$results, "04_ColonyFormingAssay_8Gy", "04_CFA_5Celllines_8Gy_means.csv"),
  row.names = FALSE
)

# 2. saving lq modell results
write.csv(
  CFA_8Gy_lq_fits,
  file = file.path(dirs$results, "04_ColonyFormingAssay_8Gy", "04_CFA_5Celllines_8Gy_LQ_coefficients.csv"),
  row.names = FALSE
)

# 3. saving anova results
openxlsx::write.xlsx(
  list(
    Pairwise = as.data.frame(CFA_8Gy_significance$contrasts),
    EMMs = as.data.frame(CFA_8Gy_significance$emmeans)
  ),
  file = "04_CFA_5Celllines_8Gy_significance.xlsx"
)

# 4. saving figure
ggsave(
  filename = file.path(dirs$results, "04_ColonyFormingAssay_8Gy", "04_CFA_5Celllines_8Gy_plot.jpeg"),
  plot = CFA_8Gy_plot,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(dirs$results, "04_ColonyFormingAssay_8Gy","04_CFA_5Celllines_8Gy_plot.pdf"),
  plot = CFA_8Gy_plot,
  width = 21/3,
  height = 7,
  units = "cm",
  dpi = 300
)
#--------------------------------------------
setwd(dirs$code)
