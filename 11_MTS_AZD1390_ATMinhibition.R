#------------------------------------------
# setting path
#------------------------------------------
setwd(file.path(dirs$MTS,"Exp3_ATMinhibition","processed")) #replace with path on your own PC
#time <- c(0, 24, 48, 72, 96, 120)
#------------------------------------------
# importing data
#------------------------------------------
MTS_AZD1390_raw_HCC44 <- read_replicates_MTS_AZD1390("HCC44")
MTS_AZD1390_raw_H23   <- read_replicates_MTS_AZD1390("H23")

MTS_AZD1390_raw <- bind_rows(
  prepare_MTS_ATMi(MTS_AZD1390_raw_HCC44, "HCC44",   "DMSO",    "0 Gy"),
  prepare_MTS_ATMi(MTS_AZD1390_raw_HCC44, "HCC44",   "DMSO",    "2 Gy"),
  prepare_MTS_ATMi(MTS_AZD1390_raw_HCC44, "HCC44",   "AZD1390", "0 Gy"),
  prepare_MTS_ATMi(MTS_AZD1390_raw_HCC44, "HCC44",   "AZD1390", "2 Gy"),
  
  prepare_MTS_ATMi(MTS_AZD1390_raw_H23,   "NCI-H23", "DMSO",    "0 Gy"),
  prepare_MTS_ATMi(MTS_AZD1390_raw_H23,   "NCI-H23", "DMSO",    "2 Gy"),
  prepare_MTS_ATMi(MTS_AZD1390_raw_H23,   "NCI-H23", "AZD1390", "0 Gy"),
  prepare_MTS_ATMi(MTS_AZD1390_raw_H23,   "NCI-H23", "AZD1390", "2 Gy")
) %>%
  mutate(
    Cellline = factor(Cellline),
    Treatment = factor(
      Treatment,
      levels = c("DMSO", "AZD1390")
    ),
    Dose = factor(
      Dose,
      levels = c("0 Gy", "2 Gy")
    ),
    BioReplicate = factor(BioReplicate),
    log_rel =log10(Intensity)
  )
#------------------------------------------
# raw data
#------------------------------------------
MTS_AZD1390_raw_mean <- MTS_AZD1390_raw %>%
  summarise(
    mean_Intensity = mean(Intensity, na.rm = TRUE),
    .by = c(Cellline, Dose, Treatment, BioReplicate, Time)
  ) 

MTS_AZD1390_raw_mean_plot<- MTS_AZD1390_raw_mean%>%
  summarise(
    mean_Intensity_raw = mean(mean_Intensity, na.rm = TRUE),
    SD_abs = sd(mean_Intensity, na.rm = TRUE),
    n = sum(!is.na(mean_Intensity)),
    .by = c(Cellline, Dose, Treatment, Time)
  ) %>%
  mutate(
    SEM_abs= SD_abs / sqrt(n),
    Cellline = relevel(factor(Cellline), ref = "HCC44"),
    Dose = relevel(factor(Dose), ref = "0 Gy"),
    Treatment = relevel(factor(Treatment), ref = "DMSO"),
    Group = factor(
      paste(Dose, Treatment),
      levels = c(
        "0 Gy DMSO",
        "0 Gy AZD1390",
        "2 Gy DMSO",
        "2 Gy AZD1390"
      )
    )
  )
#------------------------------------------
# normalization
#------------------------------------------
MTS_AZD1390_norm <- MTS_AZD1390_raw %>%
  group_by(Cellline, Dose, Treatment, BioReplicate) %>%
  mutate(
    ref_0h = mean(Intensity[Time == 0], na.rm = TRUE),
    rel_Intensity = Intensity / ref_0h
  ) %>%
  ungroup()

MTS_AZD1390_norm_mean <- MTS_AZD1390_norm %>%
  summarise(
    mean_Intensity = mean(rel_Intensity, na.rm = TRUE),
    .by = c(Cellline, Dose, Treatment, BioReplicate, Time)
  ) %>%
  mutate(
    log_rel = if_else(
      mean_Intensity > 0,
      log(mean_Intensity),
      NA_real_
    )
  )

MTS_AZD1390_norm_mean_plot<- MTS_AZD1390_norm_mean%>%
  summarise(
    mean_Intensity_norm = mean(mean_Intensity, na.rm = TRUE),
    SD_rel = sd(mean_Intensity, na.rm = TRUE),
    n = sum(!is.na(mean_Intensity)),
    .by = c(Cellline, Dose, Treatment, Time)
  ) %>%
  mutate(
    SEM_rel = SD_rel / sqrt(n),
    
    log_rel = if_else(
      mean_Intensity_norm > 0,
      log(mean_Intensity_norm),
      NA_real_
    ),
    
    Cellline = relevel(factor(Cellline), ref = "HCC44"),
    Dose = relevel(factor(Dose), ref = "0 Gy"),
    Treatment = relevel(factor(Treatment), ref = "DMSO"),
    
    Group = factor(
      paste(Dose, Treatment),
      levels = c(
        "0 Gy DMSO",
        "0 Gy AZD1390",
        "2 Gy DMSO",
        "2 Gy AZD1390"
      )
    )
  )
#------------------------------------------
# Growth rates
#------------------------------------------
MTS_AZD1390_growthrate <- MTS_AZD1390_norm_mean %>%
  filter(
    is.finite(mean_Intensity),
    is.finite(Time)
  ) %>%
  group_by(
    BioReplicate,
    Cellline,
    Dose,
    Treatment
  ) %>%
  filter(n_distinct(Time) >= 2) %>%
  group_modify(
    ~ tibble(
      GrowthRate = as.numeric(coef(
        lm(log_rel ~ Time, data = .x)
      )[["Time"]]
    ))
  ) %>%
  ungroup()

MTS_AZD1390_growthrate_summary <- MTS_AZD1390_growthrate %>%
  group_by(Cellline, Dose, Treatment) %>%
  summarise(
    n = sum(is.finite(GrowthRate)),
    GrowthRate_mean = mean(GrowthRate, na.rm = TRUE),
    SD = sd(GrowthRate, na.rm = TRUE),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )%>%
  mutate(
    label = sprintf(
      "k = %.4f ± %.4f h⁻¹",
      GrowthRate_mean,
      SEM
    )
  )

MTS_AZD1390_growthrate_label <- MTS_AZD1390_growthrate_summary %>%
    filter(
    Treatment=="AZD1390"
  )%>%
  mutate(
    label = sprintf(
      "%s %s: %.4f ± %.4f",
      Treatment,
      Dose,
      GrowthRate_mean,
      SEM
    )
  ) %>%
  group_by(Cellline) %>%
  summarise(
    label = paste(
      c("Growth rate k [h⁻¹]:", label),
      collapse = "\n"
    ),
    .groups = "drop"
  )
#------------------------------------------
# Anova viability
#------------------------------------------
MTS_AZD1390_mod_norm <- nlme::lme(
  log_rel ~ Time * Dose * Treatment * Cellline,
  random = ~ 1 | BioReplicate,
  data = MTS_AZD1390_raw
)

MTS_AZD1390_res_aov_norm <- anova(MTS_AZD1390_mod_norm)
#
MTS_AZD1390_mod_aov_growthrate <- nlme::lme(
  GrowthRate ~ Dose * Cellline * Treatment,
  random = ~ 1 | BioReplicate,
  data = MTS_AZD1390_growthrate
)

MTS_AZD1390_res_aov_growthrate <- anova(MTS_AZD1390_mod_aov_growthrate)
#------------------------------------------
# Mixed linear models
#------------------------------------------
MTS_AZD1390_mod_norm <- nlme::lme(
  log_rel ~ Time * Dose * Treatment * Cellline,
  random = ~ 1 | BioReplicate,
  data = MTS_AZD1390_norm_mean
)

# Treatment: DMSO vs. AZD1390
MTS_AZD1390_mod_norm$Treatment <- factor(
  MTS_AZD1390_mod_norm$Treatment,
  levels = c("DMSO", "AZD1390")
)

MTS_AZD1390_res_mod_DMSO_AZD1390 <- emmeans::emmeans(
  MTS_AZD1390_mod_norm,
  ~ Treatment | Cellline * Dose* Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

MTS_AZD1390_res_sig_DMSO_AZD1390 <- pairs(
  MTS_AZD1390_res_mod_DMSO_AZD1390,
  adjust = "BH"
)
# Dose: 0 Gy vs. 2 Gy
MTS_AZD1390_res_mod_0Gy_2Gy <- emmeans::emmeans(
  MTS_AZD1390_mod_norm,
  ~ Dose | Cellline * Treatment * Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

MTS_AZD1390_res_sig_0Gy_2Gy <- pairs(
  MTS_AZD1390_res_mod_0Gy_2Gy,
  adjust = "BH"
)
# Cellline: H23 vs. HCC44 
MTS_AZD1390_res_mod_H23_HCC44 <- emmeans::emmeans(
  MTS_AZD1390_mod_norm,
  ~ Cellline | Dose * Treatment * Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

MTS_AZD1390_res_sig_H23_HCC44 <- pairs(
  MTS_AZD1390_res_mod_H23_HCC44,
  adjust = "BH"
)
#------------------------------------------
# Delta growth rates
#------------------------------------------
MTS_AZD1390_growthrate_delta_dose <- MTS_AZD1390_growthrate %>%
  dplyr::summarise(
    mean_GrowthRate = mean(GrowthRate, na.rm = TRUE),
    .by = c(Cellline, Treatment, Dose)
  ) %>%
  tidyr::pivot_wider(
    names_from = Dose,
    values_from = mean_GrowthRate
  ) %>%
  dplyr::mutate(
    Comparison = "0 Gy - 2 Gy"
  )

MTS_AZD1390_growthrate_delta_dose_label <- MTS_AZD1390_growthrate_delta_dose %>%
  dplyr::transmute(
    Cellline,
    Treatment,
    Comparison,
    mean_GrowthRate =  `2 Gy`-`0 Gy`
  )
#
MTS_AZD1390_growthrate_delta_Cellline <- MTS_AZD1390_growthrate %>%
  dplyr::summarise(
    mean_GrowthRate = mean(GrowthRate, na.rm = TRUE),
    .by = c(Cellline, Treatment, Dose)
  ) %>%
  tidyr::pivot_wider(
    names_from = Cellline,
    values_from = mean_GrowthRate
  ) %>%
  dplyr::mutate(
    Comparison = "NCI-H23 - HCC44"
  )

MTS_AZD1390_growthrate_delta_Cellline_label <- MTS_AZD1390_growthrate_delta_Cellline %>%
  dplyr::transmute(
    Dose,
    Treatment,
    Comparison,
    mean_GrowthRate =  `NCI-H23` - HCC44
  )
#
MTS_AZD1390_growthrate_delta_Treatment <- MTS_AZD1390_growthrate %>%
  dplyr::summarise(
    mean_GrowthRate = mean(GrowthRate, na.rm = TRUE),
    .by = c(Cellline, Treatment, Dose)
  ) %>%
  tidyr::pivot_wider(
    names_from = Treatment,
    values_from = mean_GrowthRate
  )%>%
  dplyr::mutate(
    Comparison = "DMSO - AZD1390"
  )

MTS_AZD1390_growthrate_delta_Treatment_label <- MTS_AZD1390_growthrate_delta_Treatment %>%
  dplyr::transmute(
    Cellline,
    Dose,
    Comparison,
    mean_GrowthRate =  AZD1390 - DMSO
  )
#
MTS_AZD1390_growthrate_delta_all <- dplyr::bind_rows(
  MTS_AZD1390_growthrate %>%
    tidyr::pivot_wider(
      id_cols = c(Cellline, Treatment),
      names_from = Dose,
      values_from = GrowthRate,
      values_fn = mean
    ) %>%
    dplyr::transmute(
      Cellline,
      Treatment,
      Dose = NA_character_,
      Comparison = "0 Gy - 2 Gy",
      delta_k = `0 Gy` - `2 Gy`
    ),
  MTS_AZD1390_growthrate %>%
    tidyr::pivot_wider(
      id_cols = c(Cellline, Dose),
      names_from = Treatment,
      values_from = GrowthRate,
      values_fn = mean
    ) %>%
    dplyr::transmute(
      Cellline,
      Treatment = NA_character_,
      Dose,
      Comparison = "DMSO - AZD1390",
      delta_k = DMSO - AZD1390
    ),
  MTS_AZD1390_growthrate %>%
    tidyr::pivot_wider(
      id_cols = c(Treatment, Dose),
      names_from = Cellline,
      values_from = GrowthRate,
      values_fn = mean
    ) %>%
    dplyr::transmute(
      Cellline = NA_character_,
      Treatment,
      Dose,
      Comparison = "NCI-H23 - HCC44",
      delta_k = `NCI-H23` - HCC44
    )
) %>%
  tidyr::unite(
    Group,
    Cellline, Treatment, Dose,
    sep = " ",
    remove = FALSE,
    na.rm = TRUE
  )
#------------------------------------------
# prepare for plotting
#------------------------------------------
MTS_ATMi_colours <- c(
  "0 Gy DMSO"    = "#5F8F29",
  "0 Gy AZD1390" = "#4294F7",
  "2 Gy DMSO"    = "#F09837",
  "2 Gy AZD1390" = "#EA4025"
)

MTS_ATMi_slope_colours <- c(
  "HCC44 0 Gy DMSO"    = "#5F8F29",
  "HCC44 0 Gy AZD1390" = "#4294F7",
  "HCC44 2 Gy DMSO"    = "#F09837",
  "HCC44 2 Gy AZD1390" = "#EA4025",
  "NCI-H23 0 Gy DMSO"    = "#5F8F29",
  "NCI-H23 0 Gy AZD1390" = "#4294F7",
  "NCI-H23 2 Gy DMSO"    = "#F09837",
  "NCI-H23 2 Gy AZD1390" = "#EA4025"
)

MTS_AZD1390_sig_norm <- pairs(
  MTS_AZD1390_res_mod_DMSO_AZD1390,
  adjust = "BH"
) %>%
  as.data.frame() %>%
  mutate(
    significance = case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001  ~ "***",
      p.value < 0.01   ~ "**",
      p.value < 0.05   ~ "*",
      TRUE             ~ "ns"
    )
  )

MTS_AZD1390_sig_norm <- MTS_AZD1390_sig_norm %>%
  filter(Time >= 72,
         Dose =="2 Gy")
#------------------------------------------
# plotting normalized data
#------------------------------------------
MTS_AZD1390_plot_norm <- ggplot(
  MTS_AZD1390_norm_mean_plot,
  aes(
    x = Time,
    y = mean_Intensity_norm ,
    colour = Group,
    group = Group
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean_Intensity_norm - SEM_rel,
      ymax = mean_Intensity_norm + SEM_rel
    ),
    width = 2,
    linewidth = 0.45,
    colour = "black"
  ) +
  facet_wrap(
    ~ Cellline,
    nrow = 1
  ) +
  scale_colour_manual(
    values = MTS_ATMi_colours
  ) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(
    xlim = c(0, 130),
    ylim= c(0,6.25)
  )+ 
  labs(
    x = "Time [h]",
    y = "Viability [a.u.]"
  ) +
  theme_MTS_growth()+
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
      r = 5,
      b = 5,
      l = 5
    ),
    strip.text = element_text(
      size = 10,
      face = "bold",
      colour = "black"
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

MTS_AZD1390_plot_norm <- MTS_AZD1390_plot_norm +
  geom_text(
    data = MTS_AZD1390_sig_norm,
    aes(
      x = Time,
      y = 6,
      label = significance
    ),
    inherit.aes = FALSE,
    size =2.5
  ) +
  geom_text(
    data = MTS_AZD1390_sig_norm,
    aes(
      x =38,
      y = 6,
      label = contrast
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_text(
    data = MTS_AZD1390_growthrate_label,
    aes(
      x = 55,
      y = c(1.8),
      label =label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.5,
    lineheight = 0.9
  )+
 geom_text(
    data=MTS_AZD1390_growthrate_delta_dose_label,
    aes(
      x = 55,
      y = 0.3,
      label = paste("(ΔΔk: ):",round(mean_GrowthRate,4) )
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.5,
    colour = "black"
  )
#------------------------------------------
# Additional Graph - plotting raw data
#------------------------------------------
MTS_AZD1390_plot_raw <- ggplot(
  MTS_AZD1390_raw_mean_plot,
  aes(
    x = Time,
    y = mean_Intensity_raw ,
    colour = Group,
    group = Group
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean_Intensity_raw  - SEM_abs,
      ymax = mean_Intensity_raw  + SEM_abs
    ),
    width = 2,
    linewidth = 0.45,
    colour = "black"
  ) +
  facet_wrap(
    ~ Cellline,
    nrow = 1
  ) +
  scale_colour_manual(
    values = MTS_ATMi_colours
  ) +
  coord_cartesian(
    xlim = c(0, 130)
  )+ 
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Time [h]",
    y = "Viability [a.u.]"
  ) +
  theme_MTS_growth()+
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
      r = 5,
      b = 5,
      l = 5
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
# Additional Graph - plotting slope
#------------------------------------------
MTS_AZD1390_p_growthrate_delta<- ggplot(
  MTS_AZD1390_growthrate_delta_all,
  aes(
    x = Group,
    y = delta_k,
    fill = Comparison
  )
) +
  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ Comparison,
    scales = "free_x",
    ncol = 1
  ) +
  labs(
    x = NULL,
    y = expression(Delta * " Growth rate (" * k * ")")
  ) +
  ggprism::theme_prism(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    strip.text = element_text(
      face = "bold"
    )
  )

#-----------------------------------------
# saving results
#-----------------------------------------
dir.create(
  file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R"))
       
#1. saving raw graph
ggsave(
  filename = "11_MTS_AZD1390__raw.jpeg",
  plot = MTS_AZD1390_plot_raw,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)
  
ggsave(
  filename = "11_MTS_AZD1390__raw.pdf",
  plot = MTS_AZD1390_plot_raw,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. saving normalized graph
ggsave(
  filename = "11_MTS_AZD1390__norm.jpeg",
  plot = MTS_AZD1390_plot_norm,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "11_MTS_AZD1390_norm.pdf",
  plot = MTS_AZD1390_plot_norm,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#3. saving slope graph
ggsave(
  filename = "11_MTS_AZD1390_slope.jpeg",
  plot = MTS_AZD1390_p_growthrate_delta,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "11_MTS_AZD1390_slope.pdf",
  plot = MTS_AZD1390_p_growthrate_delta,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#4. saving aov with norm data
write.csv(
  MTS_AZD1390_res_aov_norm,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_res_aov_norm.csv"),
  row.names = FALSE
)

#5. saving aov of growth rate
write.csv(
  MTS_AZD1390_res_aov_growthrate,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_res_aov_growthrate.csv"),
  row.names = FALSE
)

#6. saving growth rate
write.csv(
  MTS_AZD1390_growthrate_summary,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_growthrate_summary.csv"),
  row.names = FALSE
)

#7. saving significance treatment
write.csv(
  MTS_AZD1390_res_sig_DMSO_AZD1390,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_res_sig_DMSO_AZD1390.csv"),
  row.names = FALSE
)

#8. saving significance Dose
write.csv(
  MTS_AZD1390_res_sig_H23_HCC44,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_res_sig_H23_HCC44.csv"),
  row.names = FALSE
)

#9. saving significance cell line
write.csv(
  MTS_AZD1390_res_sig_0Gy_2Gy,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_res_sig_0Gy_2Gy.csv"),
  row.names = FALSE
)

#10. saving growth rate difference
write.csv(
  MTS_AZD1390_growthrate_delta_Cellline_label,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_growthrate_delta_Cellline.csv"),
  row.names = FALSE
)

#11. saving Growth rate treatment
write.csv(
  MTS_AZD1390_growthrate_delta_Treatment_label,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_growthrate_delta_Treatment.csv"),
  row.names = FALSE
)

#12. saving Growth rate dose
write.csv(
  MTS_AZD1390_growthrate_delta_dose_label,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_growthrate_delta_dose.csv"),
  row.names = FALSE
)

#13. saving Growth rate Cell line
write.csv(
  MTS_AZD1390_growthrate_delta_all,
  file = file.path(dirs$results, "11_MTS_AZD1390_ATMinhibition.R", "11_MTS_AZD1390_significant_slope.csv"),
  row.names = FALSE
)
#-----------------------------------------------
setwd(dirs$code)
