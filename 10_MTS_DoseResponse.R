#------------------------------------------
# Setting path
#-----------------------------------------
setwd(file.path(dirs$MTS,"Exp1_DoseResponse","processed")) #replace with path on your own PC

mts_doseResponse_raw_HCC44 <- lapply(
  paste0("DoseResponse_HCC44_Rep", 1:3, ".csv"),
  read.csv
)

mts_doseResponse_raw_H23 <- lapply(
  paste0("DoseResponse_H23_Rep", 1:3, ".csv"),
  read.csv
)
#-----------------------------------------
# Importing and Binding Data
#-----------------------------------------
mts_doseresponse_raw <- bind_rows(
  prepare_growth_data(
    mts_doseResponse_raw_HCC44,
    "HCC44",
    "0 Gy",
    2:7
  ),
  prepare_growth_data(
    mts_doseResponse_raw_HCC44,
    "HCC44",
    "2 Gy",
    10:15
  ),
  prepare_growth_data(
    mts_doseResponse_raw_H23,
    "NCI-H23",
    "0 Gy",
    2:7
  ),
  prepare_growth_data(
    mts_doseResponse_raw_H23,
    "NCI-H23",
    "2 Gy",
    10:15
  )
)

mts_mts_doseresponse_raw_mean <- mts_doseresponse_raw %>%
  dplyr::group_by(
    Cellline,
    Dose,
    Time
  ) %>%
  dplyr::summarise(
    mean_Intensity = mean(Intensity, na.rm = TRUE),
    SD = sd(Intensity, na.rm = TRUE),
    n = sum(!is.na(Intensity)),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )  %>%
  mutate(
    Group = paste(Cellline, Dose)
  )
#-----------------------------------------
# normalization
#-----------------------------------------
# Normalisation to 0 h
mts_doseresponse_norm <- mts_doseresponse_raw %>%
  dplyr::group_by(
    Cellline,
    Dose,
    BioReplicate
  ) %>%
  dplyr::mutate(
    ref_0h = mean(Intensity[Time == 0], na.rm = TRUE),
    rel_Intensity = dplyr::if_else(
      is.finite(ref_0h) & ref_0h > 0,
      Intensity / ref_0h,
      NA_real_
    )
  ) %>%
  dplyr::ungroup()


# Mean of technical replicates
mts_doseresponse_biorep_norm <- mts_doseresponse_norm %>%
  dplyr::group_by(
    Cellline,
    Dose,
    BioReplicate,
    Time
  ) %>%
  dplyr::summarise(
    mean_Intensity = mean(Intensity, na.rm = TRUE),
    rel_Intensity = mean(rel_Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    log_rel = log(rel_Intensity)
  )

mts_doseresponse_biorep_norm_mean <- mts_doseresponse_biorep_norm %>%
  dplyr::group_by(
    Cellline,
    Dose,
    Time
  ) %>%
  dplyr::summarise(
    mean_rel_Intensity = mean(rel_Intensity, na.rm = TRUE),
    mean_Intensity = mean(mean_Intensity, na.rm = TRUE),
    SD = sd(rel_Intensity, na.rm = TRUE),
    n = sum(!is.na(rel_Intensity)),
    SEM = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Group = paste(Cellline, Dose)
  )
#-----------------------------------------
# mixed linear models
#-----------------------------------------
# ANOVA
mts_doseresponse_norm_aov_mod <- mts_doseresponse_norm %>%
  dplyr::filter(
    !is.na(rel_Intensity),
    !is.na(Time),
    !is.na(Dose),
    !is.na(Cellline),
    !is.na(BioReplicate),
    !is.na(Plate)
  )

mts_doseresponse_norm_aov <- nlme::lme(
  Intensity ~ Time * Dose * Cellline,
    random = ~ 1 | BioReplicate/Plate,
    data = mts_doseresponse_norm_aov_mod
  )

mts_doseresponse_res_norm_aov<- anova(mts_doseresponse_norm_aov)

# HCC44 vs NCI-H2
mts_doseresponse_norm_mod <- nlme::lme(
  Intensity ~ Time * Dose * Cellline,
  random = ~ 1 | BioReplicate/Plate,
  data = mts_doseresponse_norm_aov_mod
)

mts_doseresponse_emm_72h_cell <- emmeans::emmeans(
  mts_doseresponse_norm_mod,
    ~ Cellline | Dose,
    at = list(Time = 72)
  )

mts_doseresponse_sig_72h_cell<-pairs(
  mts_doseresponse_emm_72h_cell,
  adjust = "tukey"
)

# HCC44 vs NCI-H23
mts_doseresponse_emm_dose <- emmeans::emmeans(
  mts_doseresponse_norm_mod,
    ~ Cellline |  Dose* Time,
    at = list(
      Time = c(24, 48, 72, 96, 120)
    )
  )

mts_doseresponse_sig_cell<-pairs(
  mts_doseresponse_emm_dose,
  adjust = "BH"
)

# 0 Gy vs. 2 Gy - time series
mts_doseresponse_norm_emm_time <- emmeans::emmeans(
  mts_doseresponse_emm_dose,
  ~ Dose | Cellline * Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

mts_growth_significance_dose<-pairs(
  mts_doseresponse_norm_emm_time,
  adjust = "tukey"
)
#-----------------------------------------
# slope (Growth rate)
#-----------------------------------------
# growth rate with raw data
mts_doseresponse_raw_biorep_slopes <- mts_doseresponse_raw %>%
  group_by(
    BioReplicate,
    Plate,
    Cellline,
    Dose
  ) %>%
  group_modify(
    ~ tibble(
      slope = coef(
        lm(Intensity   ~ Time, data = .x)
      )[["Time"]]
    )
  ) %>%
  ungroup()

mts_doseresponse_raw_slopes <- mts_doseresponse_raw_biorep_slopes %>%
  group_by(
    BioReplicate,
    Cellline,
    Dose
  ) %>%
  summarise(
    slope = mean(slope, na.rm = TRUE),
    .groups = "drop"
  )

mts_doseresponse_raw_slopes_mean <- mts_doseresponse_raw_slopes %>%
  group_by(
    Cellline,
    Dose
  ) %>%
  summarise(
    mean_Slope = mean(slope, na.rm = TRUE),
    SD = sd(slope, na.rm = TRUE),
    n = sum(!is.na(slope)),
    SEM = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(Cellline, Dose)
  )

mts_doseresponse_mod_raw_slopes <- nlme::lme(
  slope ~ Dose * Cellline,
  random = ~ 1 | BioReplicate,
  data = mts_doseresponse_raw_slopes
)

mts_doseresponse_aov_res_raw_slopes <- anova(
  mts_doseresponse_mod_raw_slopes
)

# growth rate with normaliased Data
mts_doseresponse_biorep_norm_slopes <- mts_doseresponse_biorep_norm %>%
  group_by(
    BioReplicate,
    Cellline,
    Dose
  ) %>%
  group_modify(
    ~ tibble(
      slope = coef(
        lm(mean_Intensity  ~ Time, data = .x)
      )[["Time"]]
    )
  ) %>%
  ungroup()

mts_doseresponse_norm_slopes <- mts_doseresponse_biorep_norm_slopes %>%
  group_by(
    BioReplicate,
    Cellline,
    Dose
  ) %>%
  summarise(
    slope_mean = mean(slope, na.rm = TRUE),
    SD = sd(slope, na.rm = TRUE),
    n = sum(!is.na(slope)),
    SEM = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(Cellline, Dose)
  )

mts_doseresponse_norm_slopes_data <- mts_doseresponse_biorep_norm_slopes %>%
  group_by(
    Cellline,
    Dose
  ) %>%
  summarise(
    slope_mean = mean(slope, na.rm = TRUE),
    SD = sd(slope, na.rm = TRUE),
    n = sum(!is.na(slope)),
    SEM = SD / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(Cellline, Dose)
  )

# modelling
mts_doseresponse_mod_norm_slopes <- nlme::lme(
  slope_mean ~ Dose * Cellline,
  random = ~ 1 | BioReplicate,
  data = mts_doseresponse_norm_slopes
)

mts_doseresponse_aov_res_norm_slopes <- anova(
  mts_doseresponse_mod_norm_slopes
)

# Delta Growth rate
mts_doseresponse_norm_delta_slopes <- mts_doseresponse_norm_slopes %>%
  tidyr::pivot_wider(
    names_from = Dose,
    values_from = slope_mean
  ) %>%
  mutate(
    delta_growth_rate = `0 Gy` - `2 Gy`
  )

mts_doseresponse_norm_delta_slopes_res <- mts_doseresponse_norm_delta_slopes %>%
  group_by(Cellline) %>%
  summarise(
    mean_delta = mean(delta_growth_rate, na.rm = TRUE),
    SD = sd(delta_growth_rate, na.rm = TRUE),
    n = sum(!is.na(delta_growth_rate)),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

mts_doseresponse_label_slope_label <- mts_doseresponse_norm_delta_slopes_res %>%
  mutate(
    label = sprintf(
      "%s: %.4f ± %.4f h⁻¹",
      Cellline,
      mean_delta,
      SEM
    )
  ) %>%
  summarise(
    label = paste(
      c("ΔGrowth rate (0 vs. 2 Gy):", label),
      collapse = "\n"
    )
  )
#-----------------------------------------
# prepare plotting
#-----------------------------------------
mts_growth_colours <- c(
  "HCC44 0 Gy"    = "#4294F7",
  "HCC44 2 Gy"    = "#EA4025",
  "NCI-H23 0 Gy"  = "#5F8F29",
  "NCI-H23 2 Gy"  = "#F09837"
)

mts_doseresponse_raw <- mts_doseresponse_raw %>%
  dplyr::filter(
    !is.na(Intensity),
    !is.na(Time),
    !is.na(Dose),
    !is.na(Cellline),
    !is.na(BioReplicate),
    !is.na(Plate)
  )

mts_doseresponse_raw_sig_plot <- nlme::lme(
  Intensity ~ Time * Dose * Cellline,
  random = ~ 1 | BioReplicate/Plate,
  data = mts_doseresponse_raw
)

mts_doseresponse_raw_emm_plot <- emmeans::emmeans(
  mts_doseresponse_raw_sig_plot,
  ~ Dose | Cellline* Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

mts_doseresponse_raw_contrast_interaction <- emmeans::contrast(
  mts_doseresponse_raw_emm_plot,
  interaction = c("pairwise", "pairwise"),
  by = "Time"
)

mts_doseresponse_sig_interaction <- mts_doseresponse_raw_contrast_interaction %>%
  as.data.frame() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    significance = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    )
  )

mts_doseresponse_sig_interaction <- mts_doseresponse_sig_interaction%>%
  filter(Time>=48)
#-----------------------------------------
# plotting normalized data
#-----------------------------------------
mts_doseresponse_p_norm <- ggplot(
  mts_doseresponse_biorep_norm_mean,
  aes(
    x = Time,
    y = mean_rel_Intensity,
    colour = Group,
    group = Group
  )
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean_rel_Intensity - SEM,
      ymax = mean_rel_Intensity + SEM
    ),
    width = 2,
    linewidth = 0.45,
    colour = "black"
  ) +
  scale_colour_manual(
    values = mts_growth_colours
  ) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Time [h]",
    y = "Relative Cell Index [a.u.]"
  )+
  coord_cartesian(
    xlim = c(0, 130),
    ylim= c(0,6.5)
  ) +
  theme_MTS_growth() +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text = element_text(
      size = 8,
      colour = "black"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(
      face = "plain",
      size = 8,
      colour = "black"
    ),
    plot.margin = margin(
      t = 10,
      r = 5,
      b = 5,
      l = 5
    )
    )+
  guides(
    fill = guide_legend(
      override.aes = list(size = 2)
    ),
    shape = guide_legend(
      override.aes = list(size = 2)
    )
  )

mts_doseresponse_p_norm <- mts_doseresponse_p_norm +
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  geom_text(
    data = mts_doseresponse_sig_interaction,
    aes(
      x = Time,
      y =6,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 2.5,
    hjust=0
  )  +
  annotate(
    "text",
    x = 48,
    y = 6.45,
    label = "Cell line × Dose:",
    size = 2.5,
    fontface = "bold",
    hjust=0
  ) +
  geom_text(
    data = mts_doseresponse_label_slope_label,
    aes(
      x = 60,
      y = 1.4,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2,
    colour = "black"
  )
#-----------------------------------------
# Additional Graphes - plotting raw data
#-----------------------------------------
mts_doseresponse_p_raw <- ggplot(
  mts_mts_doseresponse_raw_mean,
  aes(
    x = Time,
    y = mean_Intensity,
    colour = Group,
    group = Group
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_errorbar(
    aes(
      ymin = mean_Intensity - SEM,
      ymax = mean_Intensity + SEM
    ),
    width = 2,
    linewidth = 0.45,
    colour = "black"
  ) +
  scale_colour_manual(
    values = mts_growth_colours
  ) +
  scale_x_continuous(
    breaks = c(0, 24, 48, 72, 96, 120),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Time [h]",
    y = "Viability [a.u.]"
  ) +
  theme_MTS_growth()
#-----------------------------------------
# Additional Graphes - plotting Growth rate
#-----------------------------------------
mts_doseresponse_p_slope <- ggplot(
  mts_doseresponse_norm_slopes_data,
  aes(
    x = Group,
    y = slope_mean,
    fill = Group
  )
) +
  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = slope_mean - SEM,
      ymax = slope_mean + SEM
    ),
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_point(
    data = mts_doseresponse_norm_slopes_data,
    aes(
      x = Group,
      y = slope_mean
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2,
    fill = "white",
    colour = "black",
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  scale_fill_manual(
    values = mts_growth_colours
  ) +
  labs(
    x = NULL,
    y = expression(
      "Growth rate " * italic(k) * " [h"^{-1} * "]"
    )
  ) +
  guides(fill = "none") +
  theme_MTS_growth() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )
#-----------------------------------------
# saving results
#-----------------------------------------
dir.create(
  file.path(dirs$results, "10_MTS_DoseResponse"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "10_MTS_DoseResponse"))
        
#1. saving raw graph
ggsave(
  filename = "10_MTS_DoseRespons_raw.jpeg",
  plot = mts_doseresponse_p_raw,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "10_MTS_DoseRespons_raw.pdf",
  plot = mts_doseresponse_p_raw,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#2. saving normalized graph
ggsave(
  filename = "10_MTS_DoseRespons_norm.jpeg",
  plot = mts_doseresponse_p_norm,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "10_MTS_DoseRespons_norm.pdf",
  plot = mts_doseresponse_p_norm,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3. saving slope graph
ggsave(
  filename = "10_MTS_DoseRespons_slope.jpeg",
  plot = mts_doseresponse_p_slope,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "10_MTS_DoseRespons_slope.pdf",
  plot = mts_doseresponse_p_slope,
  width = 21*2/3,
  height = 7,
  units = "cm"
)
#4. saving aov of growth -raw
write.csv(
  mts_doseresponse_aov_res_raw_slopes,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_aov_growthrate_raw.csv"),
  row.names = FALSE
)

#5. saving aov of growth - normalized
write.csv(
  mts_doseresponse_res_norm_aov,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_aov_growthrate_norem.csv"),
  row.names = FALSE
) 

#6. summary of normalised data
write.csv(
  mts_doseresponse_biorep_norm_mean,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_summary_norm.csv"),
  row.names = FALSE
)

#7. mts_doseresponse_norm_slopes_data
write.csv(
  mts_doseresponse_aov_res_norm_slopes,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_summary_growthrate.csv"),
  row.names = FALSE
)

#8. saving significance for plot
write.csv(
  mts_doseresponse_raw_contrast_interaction,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_raw_interaction_plot.csv"),
  row.names = FALSE
) 

#9. saving significance for difference growth between cell lines
write.csv(
  mts_doseresponse_sig_72h_cell,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_growth_72h.csv"),
  row.names = FALSE
)

#9. saving significance
write.csv(
  mts_doseresponse_sig_cell,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_growth_Time_Cellline.csv"),
  row.names = FALSE
)            

#10. saving ANOVA slopes - raw
write.csv(
  mts_doseresponse_aov_res_raw_slopes,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_slopes_raw_AOV.csv"),
  row.names = FALSE
)  

#11. saving ANOVA slopes - normalised
write.csv(
  mts_doseresponse_aov_res_norm_slopes,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_slopes_norm_AOV.csv"),
  row.names = FALSE
)  
     
#12. saving summary Delta slopes
 write.csv(
   mts_doseresponse_norm_delta_slopes_res,
   file = file.path(dirs$results, "10_MTS_DoseResponse", "10_MTS_DoseRespons_slopes_delta.csv"),
   row.names = FALSE
 )  
#-----------------------------------------------
setwd(dirs$code)
