#------------------------------------------
# setting path
#------------------------------------------
setwd(file.path(base_dir,"functional_assays","MTS assays","Exp3_ATMinhibition","processed")) #replace with path on your own PC
time <- c(0, 24, 48, 72, 96, 120)
#------------------------------------------
# importing data
#------------------------------------------
raw_HCC44 <- read_replicates_MTS_AZD1390("HCC44")
raw_H23   <- read_replicates_MTS_AZD1390("H23")

MTS_growth_raw_ATMi <- bind_rows(
  prepare_MTS_ATMi(raw_HCC44, "HCC44",   "DMSO",    "0 Gy"),
  prepare_MTS_ATMi(raw_HCC44, "HCC44",   "DMSO",    "2 Gy"),
  prepare_MTS_ATMi(raw_HCC44, "HCC44",   "AZD1390", "0 Gy"),
  prepare_MTS_ATMi(raw_HCC44, "HCC44",   "AZD1390", "2 Gy"),
  
  prepare_MTS_ATMi(raw_H23,   "NCI-H23", "DMSO",    "0 Gy"),
  prepare_MTS_ATMi(raw_H23,   "NCI-H23", "DMSO",    "2 Gy"),
  prepare_MTS_ATMi(raw_H23,   "NCI-H23", "AZD1390", "0 Gy"),
  prepare_MTS_ATMi(raw_H23,   "NCI-H23", "AZD1390", "2 Gy")
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
    BioReplicate = factor(BioReplicate)
  )
#------------------------------------------
# normalization
#------------------------------------------
MTS_data_stat_mean_ATMi <- MTS_growth_raw_ATMi %>%
  group_by(
    Cellline,
    Dose,
    Treatment,
    BioReplicate
  ) %>%
  mutate(
    ref_0h = mean(Intensity[Time == 0], na.rm = TRUE),
    rel_Intensity = Intensity / ref_0h
  ) %>%
  group_by(
    Cellline,
    Dose,
    Treatment,
    BioReplicate,
    Time
  ) %>%
  summarise(
    mean_Intensity = mean(Intensity, na.rm = TRUE),
    mean_rel_Intensity = mean(rel_Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Cellline = relevel(factor(Cellline), ref = "HCC44"),
    Dose = relevel(factor(Dose), ref = "0 Gy"),
    Treatment = relevel(factor(Treatment), ref = "DMSO"),
    BioReplicate = factor(BioReplicate),
    log_rel = log(mean_rel_Intensity)
  )

plot_MTS_summary_ATMi <- MTS_data_stat_mean_ATMi %>%
  group_by(
    Cellline,
    Dose,
    Treatment,
    Time
  ) %>%
  summarise(
    mean2_Intensity = mean(mean_Intensity, na.rm = TRUE),
    SD_raw = sd(mean_Intensity, na.rm = TRUE),
    SEM_raw = SD_raw / sqrt(sum(!is.na(mean_Intensity))),
    
    mean2_rel_Intensity = mean(mean_rel_Intensity, na.rm = TRUE),
    SD_rel = sd(mean_rel_Intensity, na.rm = TRUE),
    SEM_rel = SD_rel / sqrt(sum(!is.na(mean_rel_Intensity))),
    
    n = sum(!is.na(mean_rel_Intensity)),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(Dose, Treatment),
    Group = factor(
      Group,
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
slopes_MTS_ATMi <- MTS_data_stat_mean_ATMi %>%
  group_by(
    BioReplicate,
    Cellline,
    Dose,
    Treatment
  ) %>%
  group_modify(
    ~ tibble(
      slope = coef(
        lm(log_rel ~ Time, data = .x)
      )[["Time"]]
    )
  ) %>%
  ungroup()
#------------------------------------------
# Mixed linear models
#------------------------------------------
mod_raw_MTS_ATMi <- nlme::lme(
  log_rel ~ Time * Dose * Treatment * Cellline,
  random = ~ 1 | BioReplicate,
  data = MTS_data_stat_mean_ATMi
)

summary(mod_raw_MTS_ATMi)
MTS_raw_ATMi_aov <- anova(mod_raw_MTS_ATMi)

mod_mts_slope_ATMi <- nlme::lme(
  slope ~ Dose * Cellline * Treatment,
  random = ~ 1 | BioReplicate,
  data = slopes_MTS_ATMi
)

summary(mod_mts_slope_ATMi)
slope_mts_ATMi_aov <- anova(mod_mts_slope_ATMi)

emm_mts_ATMi_dose <- emmeans::emmeans(
  mod_mts_slope_ATMi,
  ~ Treatment | Cellline * Dose
)

significance_mts_ATMi <- pairs(
  emm_mts_ATMi_dose,
  adjust = "BH"
)

summary(significance_mts_ATMi)
#------------------------------------------
# slope
#------------------------------------------
slopes_MTS_mean_ATMi <- slopes_MTS_ATMi %>%
  group_by(
    Cellline,
    Dose,
    Treatment
  ) %>%
  summarise(
    mean_Slope = mean(slope, na.rm = TRUE),
    SD = sd(slope, na.rm = TRUE),
    n = sum(!is.na(slope)),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

slopes_MTS_mean_ATMi <- slopes_MTS_mean_ATMi %>%
  mutate(
    Group = paste(Cellline, Dose, Treatment),
    Group = factor(
      Group,
      levels = c(
        "HCC44 0 Gy DMSO",
        "HCC44 0 Gy AZD1390",
        "HCC44 2 Gy DMSO",
        "HCC44 2 Gy AZD1390",
        "NCI-H23 0 Gy DMSO",
        "NCI-H23 0 Gy AZD1390",
        "NCI-H23 2 Gy DMSO",
        "NCI-H23 2 Gy AZD1390"
      )
    )
  )
#------------------------------------------
# Anova viability
#------------------------------------------
# Treatment: DMSO vs. AZD1390
emm_mts_treatment_ATMi <- emmeans::emmeans(
  mod_mts_slope_ATMi,
  ~ Treatment | Cellline * Dose,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

sig_mts_treatment_ATMi <- pairs(
  emm_mts_treatment_ATMi,
  adjust = "BH"
)
# Dose: 0 Gy vs. 2 Gy
emm_mts_dose_ATMi <- emmeans::emmeans(
  mod_raw_MTS_ATMi,
  ~ Dose | Cellline * Treatment * Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

sig_mts_dose_ATMi <- pairs(
  emm_mts_dose_ATMi,
  adjust = "BH"
)
# Cellline: H23 vs. HCC44 
emm_mts_cell_ATMi <- emmeans::emmeans(
  mod_raw_MTS_ATMi,
  ~ Cellline | Dose * Treatment * Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

sig_mts_cell_ATMi <- pairs(
  emm_mts_cell_ATMi,
  adjust = "BH"
)
#------------------------------------------
# Growth Rates:
#------------------------------------------  
# Growth rate: Dose
emm_mts_slope_dose_ATMi <- emmeans::emmeans(
  mod_mts_slope_ATMi,
  ~ Dose | Cellline * Treatment
)

sig_mts_slope_dose_ATMi <- pairs(
  emm_mts_slope_dose_ATMi,
  adjust = "tukey"
)
# Growth rate: Treatment 
emm_mts_slope_treatment_ATMi <- emmeans::emmeans(
  mod_mts_slope_ATMi,
  ~ Treatment | Cellline * Dose
)

sig_mts_slope_treatment_ATMi <- pairs(
  emm_mts_slope_treatment_ATMi,
  adjust = "tukey"
)
#------------------------------------------
# Delta growth rates
#------------------------------------------
delta_MTS_rad_ATMi <- slopes_MTS_mean_ATMi %>%
  dplyr::select(
    Cellline,
    Treatment,
    Dose,
    mean_Slope
  ) %>%
  tidyr::pivot_wider(
    names_from = Dose,
    values_from = mean_Slope
  ) %>%
  mutate(
    delta_growth_rate = `0 Gy` - `2 Gy`
  )

delta_MTS_treatment_ATMi <- slopes_MTS_mean_ATMi %>%
  dplyr::select(
    Cellline,
    Dose,
    Treatment,
    mean_Slope
  ) %>%
  tidyr::pivot_wider(
    names_from = Treatment,
    values_from = mean_Slope
  ) %>%
  mutate(
    delta_growth_rate = DMSO - AZD1390
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


sig_norm_mts_ATMi <- pairs(
  emm_mts_ATMi_dose,
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

sig_norm_mts_ATMi <- sig_norm_mts_ATMi %>%
  left_join(
    plot_MTS_summary_ATMi %>%
      group_by(Cellline, Dose) %>%
      summarise(
        y_norm = max(mean2_rel_Intensity + SEM_rel, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("Cellline", "Dose")
  ) %>%
  left_join(
    plot_MTS_summary_ATMi %>%
      group_by(Cellline, Dose) %>%
      summarise(
        y_raw = max(mean2_Intensity + SEM_raw, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("Cellline", "Dose")
  ) %>%
  mutate(
    y_norm = y_norm + 0.2,
    y_raw  = y_raw + 0.2
  )

sig_norm_plot_mts_ATMi <- sig_norm_mts_ATMi %>%
  filter(significance != "ns")
#------------------------------------------
# plotting raw data
#------------------------------------------
p_mts_raw_ATMi <- ggplot(
  plot_MTS_summary_ATMi,
  aes(
    x = Time,
    y = mean2_Intensity,
    colour = Group,
    group = Group
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 1.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean2_Intensity - SEM_raw,
      ymax = mean2_Intensity + SEM_raw
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
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Time [h]",
    y = "Viability [a.u.]"
  ) +
  theme_MTS_growth()

p_mts_raw_ATMi <- p_mts_raw_ATMi +
  geom_text(
    data = sig_norm_plot_mts_ATMi,
    aes(
      x = 95,
      y = y_raw - 0.475,
      label = paste0(contrast, "\n", significance),
      vjust=0.25
    ),
    inherit.aes = FALSE,
    size = 2.5
  )
#------------------------------------------
# plotting normalized data
#------------------------------------------
p_mts_norm_ATMi <- ggplot(
  plot_MTS_summary_ATMi,
  aes(
    x = Time,
    y = mean2_rel_Intensity,
    colour = Group,
    group = Group
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 1.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean2_rel_Intensity - SEM_raw,
      ymax = mean2_rel_Intensity + SEM_raw
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
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Time [h]",
    y = "Viability [a.u.]"
  ) +
  theme_MTS_growth()

p_mts_norm_ATMi <- p_mts_norm_ATMi +
  geom_text(
    data = sig_norm_plot_mts_ATMi,
    aes(
      x = 100,
      y = y_norm - 1.45,
      label = paste0(contrast, "\n", significance)
    ),
    inherit.aes = FALSE,
    size = 2.5
  )
#------------------------------------------
# plotting slope
#------------------------------------------
p_MTS_slope_ATMi <- ggplot(
  slopes_MTS_mean_ATMi,
  aes(
    x = Group,
    y = mean_Slope,
    fill = Group
  )
) +
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.8,
    colour = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_Slope - SEM,
      ymax = mean_Slope + SEM
    ),
    width = 0.15,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_point(
    data = slopes_MTS_mean_ATMi,
    aes(
      x = Group,
      y = mean_Slope,
      fill = Group
    ),
    shape = 21,
    size = 1.8,
    stroke = 0.5,
    colour = "black",
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    inherit.aes = FALSE
  ) +
    facet_wrap(
    ~ Cellline,
    nrow = 1
  ) +
    scale_fill_manual(
    values = MTS_ATMi_slope_colours
  ) +
    labs(
    x = NULL,
    y = expression(
      "Growth rate " * italic(k) * " [h"^{-1} * "]"
    )
  ) +
    guides(
    fill = guide_legend(
      title = NULL,
      override.aes = list(
        size = 2
      )
    )
  ) +
    theme_MTS_growth() +
    theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
#-----------------------------------------
# saving results
#-----------------------------------------
dir.create(
  file.path(dirs$results, "11_MTS_ATMInhibition"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "11_MTS_ATMInhibition"))

#1. saving raw graph
ggsave(
  filename = "mts_ATMi_raw.jpeg",
  plot = p_mts_raw_ATMi,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "mts_ATMi_raw.pdf",
  plot = p_mts_raw_ATMi,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. saving normalized graph
ggsave(
  filename = "mts_ATMi_norm.jpeg",
  plot = p_mts_norm_ATMi,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "mts_ATMi_norm.pdf",
  plot = p_mts_norm_ATMi,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#3. saving slope graph
ggsave(
  filename = "mts_ATMi_slope.jpeg",
  plot = p_MTS_slope_ATMi,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "mts_ATMi_slope.pdf",
  plot = p_MTS_slope_ATMi,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#4. saving aov of ATMi
write.csv(
  delta_MTS_rad_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "mts_ATMi_aov_delta_growthrate_rad.csv"),
  row.names = FALSE
)

#5. saving aov of ATMi
write.csv(
  delta_MTS_treatment_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "mts_ATMi_aov_delta_growthrate_AZD1390.csv"),
  row.names = FALSE
)

#6. saving aov of slope
write.csv(
  sig_mts_slope_treatment_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "mts_ATMi_aov_slope.csv"),
  row.names = FALSE
)

#7. saving linear ATMi model
MTS_ATMi_results <- as.data.frame(
  summary(mod_mts_slope_ATMi)$tTable
)

MTS_ATMi_results$term <- rownames(MTS_ATMi_results)

write.csv(
  MTS_ATMi_results,
  file.path(
    dirs$results,
    "11_MTS_ATMInhibition",
    "mts_ATMi_linearModel.csv"
  ),
  row.names = FALSE
)
  
#8. saving significance treatment
write.csv(
  sig_mts_treatment_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "slope_mts_ATMi_significant_treatment.csv"),
  row.names = FALSE
)

#9. saving significance Dose
write.csv(
  sig_mts_dose_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "slope_mts_ATMi_significant_dose.csv"),
  row.names = FALSE
)

#10. saving significance cell line
write.csv(
  sig_mts_cell_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "slope_mts_ATMi_significant_cellline.csv"),
  row.names = FALSE
)

#11. saving significance slope
write.csv(
  sig_mts_slope_dose_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "slope_mts_ATMi_significant_slope.csv"),
  row.names = FALSE
)
   
#12. savinglope
write.csv(
  slopes_MTS_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "slopes_MTS_ATMi.csv"),
  row.names = FALSE
)

#13. saving AOV
write.csv(
  MTS_raw_ATMi_aov,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "MTS_raw_ATMi_aov.csv"),
  row.names = FALSE
)

#11. saving significance
write.csv(
  significance_mts_ATMi,
  file = file.path(dirs$results, "11_MTS_ATMInhibition", "significance_mts_ATMi.csv"),
  row.names = FALSE
)
#-----------------------------------------------
setwd(dirs$code)
