#-----------------------------------------
# setting path
#-----------------------------------------
setwd(file.path(base_dir,"functional_assays","MTS assays","Exp2_IC50_AZD1390","processed")) #replace with path on your own PC

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
  group_by(Cellline, BioReplicate) %>%
  mutate(
    rel_Intensity = Intensity /
      mean(Intensity[Concentration == 0], na.rm = TRUE),
    Dose_log = Concentration 
  ) %>%
  ungroup()

MTS_IC50_summary <- MTS_IC50_norm_control %>%
  group_by(Cellline, Concentration) %>%
  summarise(
    n = sum(!is.na(rel_Intensity)),
    mean_rel_Intensity = mean(rel_Intensity, na.rm = TRUE),
    sd_rel_Intensity = sd(rel_Intensity, na.rm = TRUE),
    sem_rel_Intensity = sd_rel_Intensity / sqrt(n),
    .groups = "drop"
  )
#-----------------------------------------
# IC50 modelling
#-----------------------------------------
MTS_global_IC50_model <- MTS_IC50_raw %>%
  group_by(Cellline) %>%
  mutate(
    Dose_log = Concentration, 
    viability_norm = Intensity / mean(Intensity[Dose_log == 0]) * 100,
    Dose_log = ifelse(Dose_log == 0, 1e-6, Dose_log)
  ) %>%
  ungroup()

MTS_model_IC50 <- drc::drm(
  viability_norm ~ Dose_log,
  data = MTS_global_IC50_model %>%
    filter(Cellline == "HCC44"),
  fct = drc::LL.4()
)

MTS_IC50_values <- as.data.frame(
  drc::ED(
    MTS_model_IC50,
    50,
    type = "absolute",
    interval = "delta"
  )
)

MTS_IC50_values$Cellline <- "HCC44"

MTS_IC50_predictions <- MTS_IC50_summary %>%
  filter(Concentration > 0) %>%
  group_by(Cellline) %>%
  group_modify(~ {
    model <- loess(
      mean_rel_Intensity ~ log10(Concentration),
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

MTS_IC50_text <- MTS_IC50_values %>%
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
MTS_IC50_stats <- MTS_IC50_norm_control %>%
  group_by(
    Cellline,
    BioReplicate,
    Block,
    Concentration
  ) %>%
  summarise(
    viability = mean(rel_Intensity, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    Cellline = factor(Cellline),
    Concentration = factor(
      Concentration,
      levels = c(0, 0.01, 0.1, 1, 10, 100)
    ),
    BioReplicate = factor(BioReplicate),
    Block = factor(Block),
    SampleID = interaction(
      Cellline,
      BioReplicate,
      Block,
      drop = TRUE
    )
  )

MTS_model_IC50 <- nlme::lme(
  viability ~ Cellline * Concentration,
  random = ~ 1 | BioReplicate/Block,
  data = MTS_IC50_stats,
  method = "REML"
)

MTS_anova_IC50<-anova(MTS_model_IC50)

MTS_emm_IC50 <- emmeans(
  MTS_model_IC50,
  ~ Concentration | Cellline
)

MTS_Dunnett_IC50 <- contrast(
  MTS_emm_IC50,
  method = "trt.vs.ctrl",
  ref = 1,
  adjust = "dunnett"
)

summary(MTS_Dunnett_IC50)

MTS_equivalence_IC50 <- test(
  MTS_Dunnett_IC50,
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
MTS_p_IC50_unirradiated <- ggplot(
  MTS_IC50_plotting %>%
    filter(Concentration > 0),
  aes(
    x = Concentration,
    y = mean_rel_Intensity_m,
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
      ymin = mean_rel_Intensity_m - sem_rel,
      ymax = mean_rel_Intensity_m + sem_rel
    ),
    width = 0.5,
    linewidth = 0.5,
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
    colour = "grey60"
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
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 5,
      l = 10
    )
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

MTS_p_IC50_unirradiated <-MTS_p_IC50_unirradiated +
  geom_vline(
    data = MTS_IC50_values %>%
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
    y = 25,
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
  filename = "MTS_IC50_AZD1390.jpeg",
  plot = MTS_p_IC50_unirradiated,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "MTS_IC50_AZD1390.pdf",
  plot = MTS_p_IC50_unirradiated,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#2. equivalence test
write.csv(
  MTS_equivalence_IC50,
  file = file.path(dirs$results, "12_MTS_IC50_AZD1390", "equivalence_IC50_testresults.csv"),
  row.names = FALSE
)

#3. Dunnett pairs
write.csv(
  MTS_Dunnett_IC50,
  file = file.path(dirs$results, "12_MTS_IC50_AZD1390", "Dunnett_IC50_testresults.csv"),
  row.names = FALSE
)

#4. post hoc analysis
write.csv(
  MTS_emm_IC50,
  file = file.path(dirs$results, "12_MTS_IC50_AZD1390", "emm_IC50_testresults.csv"),
  row.names = FALSE
)

#5. ANOVA
write.csv(
  MTS_anova_IC50,
  file = file.path(dirs$results, "12_MTS_IC50_AZD1390", "anova_IC50_testresults.csv"),
  row.names = FALSE
)

#-----------------------------------------------
setwd(dirs$code)
