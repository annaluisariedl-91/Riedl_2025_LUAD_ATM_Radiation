#-----------------------------------------
# setting path
#-----------------------------------------
setwd(file.path(dirs$MTS,"Exp2_IC50_AZD1390","processed")) #replace with path on your own PC

MTS_raw_IC50_HCC44 <- read_IC50("HCC44")
MTS_raw_IC50_H23   <- read_IC50("NCI-H23")
#-----------------------------------------
# binding and normalization 
#-----------------------------------------
MTS_IC50_raw <- bind_rows(
  MTS_raw_IC50_HCC44,
  MTS_raw_IC50_H23
)

MTS_IC50_norm_control <- MTS_IC50_raw %>%
  group_by(Cellline, BioReplicate, Block) %>%
  mutate(
    Control_plate = mean(
      Intensity[Concentration == 0],
      na.rm = TRUE
    ),
    rel_Intensity = Intensity / Control_plate,
    Dose_log = Concentration
  ) %>%
  ungroup()

MTS_IC50_plate <- MTS_IC50_norm_control %>%
  group_by(Cellline, BioReplicate, Block, Concentration) %>%
  summarise(
    n_tech = sum(!is.na(rel_Intensity)),
    mean_rel_Intensity = mean(rel_Intensity, na.rm = TRUE),
    .groups = "drop"
  )

MTS_IC50_biorep <- MTS_IC50_plate %>%
  group_by(Cellline, BioReplicate, Concentration) %>%
  summarise(
    n_plate = sum(!is.na(mean_rel_Intensity)),
    mean_rel_Intensity = mean(mean_rel_Intensity, na.rm = TRUE),
    .groups = "drop"
  )

MTS_IC50_summary <- MTS_IC50_biorep %>%
  group_by(Cellline, Concentration) %>%
  summarise(
    n = sum(!is.na(mean_rel_Intensity)),
    mean_rel_Intensity2 = mean(mean_rel_Intensity, na.rm = TRUE),
    sd_rel_Intensity = sd(mean_rel_Intensity, na.rm = TRUE),
    sem_rel_Intensity = sd_rel_Intensity / sqrt(n),
    .groups = "drop"
  )
#-----------------------------------------
# IC50 modelling
#-----------------------------------------
MTS_IC50_global_model <- MTS_IC50_raw %>%
  group_by(Cellline) %>%
  mutate(
    Dose_log = Concentration, 
    viability_norm = Intensity / mean(Intensity[Dose_log == 0]) * 100,
    Dose_log = ifelse(Dose_log == 0, 1e-6, Dose_log)
  ) %>%
  ungroup()

MTS_IC50_global_model_IC50 <- drc::drm(
  viability_norm ~ Dose_log,
  data = MTS_IC50_global_model %>%
    dplyr::filter(Cellline == "HCC44"),
  fct = drc::LL.4()
)

MTS_IC50_global_model_values <- drc::ED(
  MTS_IC50_global_model_IC50,
  50,
  type = "absolute",
  interval = "delta"
) %>%
  as.data.frame() %>%
  dplyr::mutate(
    Cellline = "HCC44",
    .before = 1
  )
#-----------------------------------------
# Prediction for curve plotting
#-----------------------------------------
MTS_IC50_predictions <- MTS_IC50_summary %>%
  filter(Concentration > 0) %>%
  group_by(Cellline) %>%
  group_modify(~ {
    model <- loess(
      mean_rel_Intensity2 ~ log10(Concentration),
      data = .x,
      span = 1
    )
    newdata <- data.frame(
      Concentration = 10^seq(
        log10(min(.x$Concentration)),
        log10(max(.x$Concentration)),
        length.out = 200
      )
    )
    newdata$prediction <- predict(
      model,
      newdata = newdata
    )
    newdata
  }) %>%
  ungroup()

MTS_IC50_text <- MTS_IC50_global_model_values %>%
  mutate(
    label = ifelse(
      is.na(Estimate),
      "Not reached",
      sprintf("%.2f µM", Estimate)
    )
  ) %>%
  dplyr::select(Cellline, label) %>%
  bind_rows(
    data.frame(
      Cellline = "NCI-H23",
      label = "Not reached"
    )
  ) %>%
  summarise(
    text = paste(
      paste0(Cellline, ": ", label),
      collapse = "\n"
    )
  ) %>%
  pull(text)

MTS_IC50_text <- paste0(
  "IC50 (unirradiated)\n",
  MTS_IC50_text,
  "\n---- 50% viability"
)
#-----------------------------------------
# statistics
#-----------------------------------------
MTS_IC50_stats <- MTS_IC50_biorep %>%
  group_by(
    Cellline,
    BioReplicate,
    #Block,
    Concentration
  ) %>%
  summarise(
    viability = mean(mean_rel_Intensity, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Cellline = factor(Cellline),
    Concentration = factor(
      Concentration,
      levels = c(0, 0.01, 0.1, 1, 10, 100)
    ),
    BioReplicate = factor(BioReplicate),
    #Block = factor(Block),
    SampleID = interaction(
      Cellline,
      BioReplicate,
      #Block,
      drop = TRUE
    )
  )

MTS_IC50_model <- nlme::lme(
  viability ~ Cellline * Concentration,
  random = ~ 1 | BioReplicate,#/Block,
  data = MTS_IC50_stats,
  method = "REML"
)

MTS_IC50_anova<-anova(MTS_IC50_model)

MTS_IC50_emm <- emmeans(
  MTS_IC50_model,
  ~ Concentration | Cellline
)

MTS_IC50_Dunnett <- contrast(
  MTS_IC50_emm,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "dunnett"
)

MTS_IC50_equivalence <- test(
  MTS_IC50_Dunnett,
  delta = 10,
  side = "equivalence",
  adjust = "none"
)
#-----------------------------------------
# prepare plotting
#-----------------------------------------
MTS_IC50_plotting <- MTS_IC50_norm_control %>%
  group_by(
    Cellline,
    BioReplicate,
    Concentration
  ) %>%
  summarise(
    mean_rel_Intensity = mean(rel_Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(
    Cellline,
    Concentration
  ) %>%
  summarise(
    mean_rel_Intensity_m = mean(mean_rel_Intensity) * 100,
    sem_rel = sd(mean_rel_Intensity) /
      sqrt(sum(!is.na(mean_rel_Intensity))) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Cellline = factor(
      Cellline,
      levels = c("NCI-H23", "HCC44")
    )
  )

MTS_IC50_colours <- c(
  "NCI-H23" = "#5F8F29",
  "HCC44"   = "#4294F7"
)
#-----------------------------------------
# Plotting IC50
#-----------------------------------------
MTS_IC50_p_unirradiated <- ggplot(
  MTS_IC50_summary %>%
    filter(Concentration > 0),
  aes(
    x = Concentration,
    y = mean_rel_Intensity2*100,
    colour = Cellline,
    group = Cellline
  )
) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    span = 1,
    se = FALSE,
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = (mean_rel_Intensity2*100)  - (sem_rel_Intensity*100),
      ymax = (mean_rel_Intensity2*100)  + (sem_rel_Intensity*100)
    ),
    width =0.5,
    linewidth = 0.6,
    colour = "black"
  ) +
  geom_point(
    aes(fill = Cellline),
    shape = 21,
    size = 2,
    stroke = 0.6,
    colour = "black",
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 100,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey50"
  ) +
  geom_hline(
    yintercept = 50,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey50"
  ) +
  scale_colour_manual(
    values = MTS_IC50_colours,
    breaks = c("HCC44", "NCI-H23"),
    name = NULL
  ) +
  scale_fill_manual(
    values = MTS_IC50_colours,
    guide = "none"
  ) +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
    labels = c("0.001", "0.01", "0.1", "1", "10", "100"),
    limits = c(0.001, 150),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(0, 50, 100),
    minor_breaks = c(25, 75),
    limits = c(0, 120),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "AZD1390 [µM]",
    y = "Relative viability [%]"
  ) +
  theme_prism_xcell() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7.5),
    legend.key.width = grid::unit(6, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 5,
      l = 10
    ),
      axis.title = element_text(
        size = 12,
        face = "bold"
      ),
      axis.text = element_text(
        size = 10,
        colour = "black"
      ),
      axis.text.x = element_text(
        face = "plain",
        size = 10,
        colour = "black"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
       ) +
  guides(
    fill = "none",
    colour = guide_legend(
      title = NULL,
      nrow = 1,
      override.aes = list(
        linewidth = 0.5
      )
    )
  )

MTS_IC50_p_unirradiated <-MTS_IC50_p_unirradiated +
  geom_vline(
    data = MTS_IC50_global_model_values %>%
      filter(!is.na(Estimate)),
    aes(
      xintercept = Estimate,
      colour = Cellline
    ),
    linetype = "dotted",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 0.002,
    y = 45,
    label = paste0(
      "Cell line (IC50)\n",
      MTS_IC50_text
    ),
    hjust = 0,
    vjust = 1,
    size = 2.5,
    lineheight = 1.1
  )
#----------------------------------------
# saving results
#-----------------------------------------
dir.create(
  file.path(dirs$results, "12_MTS_IC50_AZD1390"),
  showWarnings = FALSE,
  recursive = TRUE
)

setwd(file.path(dirs$results, "12_MTS_IC50_AZD1390"))

#1. saving dose response curve
ggsave(
  filename = "12_MTS_IC50_AZD1390.jpeg",
  plot = MTS_IC50_p_unirradiated,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "12_MTS_IC50_AZD1390.pdf",
  plot = MTS_IC50_p_unirradiated,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#2. equivalence test
write.csv(
  MTS_IC50_equivalence,
  file = file.path(dirs$results, "12_MTS_IC50_AZD1390", "12_MTS_IC50_equivalence.csv"),
  row.names = FALSE
)

#3. ANOVA
write.csv(
  MTS_IC50_anova,
  file = file.path(dirs$results, "12_MTS_IC50_AZD1390", "12_MTS_IC50_anova.csv"),
  row.names = FALSE
)

#-----------------------------------------------
setwd(dirs$code)
