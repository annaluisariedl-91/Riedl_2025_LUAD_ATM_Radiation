#------------------------------------------
# Setting path
#-----------------------------------------
setwd(file.path(base_dir,"functional_assays","MTS assays","Exp1_DoseResponse","processed")) #replace with path on your own PC

raw_HCC44 <- lapply(
  paste0("DoseResponse_HCC44_Rep", 1:3, ".csv"),
  read.csv
)

raw_H23 <- lapply(
  paste0("DoseResponse_H23_Rep", 1:3, ".csv"),
  read.csv
)
#-----------------------------------------
# Importing and Binding Data
#-----------------------------------------
growth_raw_mts_growth <- bind_rows(
  prepare_growth_data(
    raw_HCC44,
    "HCC44",
    "0 Gy",
    2:7
  ),
  prepare_growth_data(
    raw_HCC44,
    "HCC44",
    "2 Gy",
    10:15
  ),
  prepare_growth_data(
    raw_H23,
    "NCI-H23",
    "0 Gy",
    2:7
  ),
  prepare_growth_data(
    raw_H23,
    "NCI-H23",
    "2 Gy",
    10:15
  )
)
#-----------------------------------------
# normalization
#-----------------------------------------
growth_plate_mts_growth <- growth_raw_mts_growth %>%
  dplyr::group_by(
    Cellline,
    Dose,
    BioReplicate,
    Plate,
    Time
  )  %>%
  dplyr::summarise(
    mean_Intensity = mean(Intensity, na.rm = TRUE),
    SD = sd(Intensity, na.rm = TRUE),
    n = sum(!is.na(Intensity)),
    SEM = SD / sqrt(n),
    .groups = "drop" )



norm_plate_mts_growth <- growth_plate_mts_growth %>%
  dplyr::group_by(
    Cellline,
    BioReplicate,
    Plate
  ) %>%
  dplyr::mutate(
    max_Intensity = max(mean_Intensity, na.rm = TRUE),
    rel_Intensity = dplyr::if_else(
      is.finite(max_Intensity) & max_Intensity > 0,
      mean_Intensity / max_Intensity,
      NA_real_
    )
  ) %>%
  dplyr::ungroup()

norm_biorep_mts_growth <- norm_plate_mts_growth %>%
  dplyr::group_by(
    Cellline,
    Dose,
    BioReplicate,
    Time
  ) %>%
  dplyr::summarise(
    rel_Intensity = mean(rel_Intensity, na.rm = TRUE),
    mean_Intensity =mean(mean_Intensity, na.rm = TRUE),
    .groups = "drop"
  )

norm_sum_biorep_mts_growth <- norm_biorep_mts_growth %>%
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
mod_data_mts_growth <- growth_raw_mts_growth %>%
  dplyr::filter(
    !is.na(Intensity),
    !is.na(Time),
    !is.na(Dose),
    !is.na(Cellline),
    !is.na(BioReplicate),
    !is.na(Plate)
  )

mod_mts_raw_growth <- nlme::lme(
  Intensity ~ Time * Dose * Cellline,
    random = ~ 1 | BioReplicate/Plate,
    data = mod_data_mts_growth
  )

mod_data_mts_growth_summary<-summary(mod_mts_raw_growth)
mod_data_mts_growth_aov<- anova(mod_mts_raw_growth)

emm_mts_growth_72h_cell <- emmeans::emmeans(
  mod_mts_raw_growth,
    ~ Cellline | Dose,
    at = list(Time = 72)
  )

significance_mts_growth_72h<-pairs(
  emm_mts_growth_72h_cell,
  adjust = "tukey"
)

emm_mts_growth_dose <- emmeans::emmeans(
  mod_mts_raw_growth,
    ~ Dose | Cellline * Time,
    at = list(
      Time = c(24, 48, 72, 96, 120)
    )
  )

significance_mts_growth_dose<-pairs(
  emm_mts_growth_dose,
  adjust = "BH"
)
#-----------------------------------------
# slope
#-----------------------------------------
slopes_mts_growth <- growth_raw_mts_growth %>%
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

slopes_biorep_mts_growth <- slopes_mts_growth %>%
  group_by(
    BioReplicate,
    Cellline,
    Dose
  ) %>%
  summarise(
    slope = mean(slope, na.rm = TRUE),
    .groups = "drop"
  )

slopes_mean_mts_growth <- slopes_biorep_mts_growth %>%
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

mod_slope_mts_growth <- nlme::lme(
  slope ~ Dose * Cellline,
  random = ~ 1 | BioReplicate,
  data = slopes_mts_growth
)

summary(mod_slope_mts_growth)

slope_mts_growth_aov <- anova(
  mod_slope_mts_growth
)
#-----------------------------------------
# Post-hoc
#-----------------------------------------
emm_slope_mts_growth <- emmeans::emmeans(
  mod_slope_mts_growth,
  ~ Dose | Cellline
)  

slope_mts_growth_posthoc<-pairs(
  emm_slope_mts_growth,
  adjust = "tukey"
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

emm_norm_mts_growth <- emmeans::emmeans(
  mod_mts_raw_growth,
  ~ Dose | Cellline * Time,
  at = list(
    Time = c(24, 48, 72, 96, 120)
  )
)

sig_norm_mts_growth <- pairs(
  emm_norm_mts_growth,
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

sig_norm_mts_growth <- sig_norm_mts_growth %>%
  left_join(
    norm_sum_biorep_mts_growth %>%
      filter(Dose == "2 Gy") %>%
      group_by(Cellline, Time) %>%
      summarise(
        y_norm = max(
          if_else(
            Cellline == "HCC44",
            mean_rel_Intensity + SEM+0.12,
            mean_rel_Intensity - SEM-0.05
          ),
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = c("Cellline", "Time")
  ) %>%
  left_join(
    norm_sum_biorep_mts_growth %>%
      filter(Dose == "2 Gy") %>%
      group_by(Cellline, Time) %>%
      summarise(
        y_raw = max(
          if_else(
            Cellline == "HCC44",
            mean_Intensity + SEM-0.1,
            mean_Intensity - SEM-0.025
          ),
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = c("Cellline", "Time")
  )

sig_norm_plot_mts_growth <- sig_norm_mts_growth %>%
  filter(significance != "ns")
#-----------------------------------------
# plotting raw data
#-----------------------------------------
p_raw_mts_growth <- ggplot(
  norm_sum_biorep_mts_growth,
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

p_raw_mts_growth <- p_raw_mts_growth +
  geom_text(
    data = sig_norm_plot_mts_growth,
    aes(
      x = Time,
      y = y_raw-0.085,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3
  )+
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
#-----------------------------------------
# plotting normalized data
#-----------------------------------------
p_norm_mts_growth <- ggplot(
  norm_sum_biorep_mts_growth,
  aes(
    x = Time,
    y = mean_rel_Intensity,
    colour = Group,
    group = Group
  )
) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.5) +
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
  ) +
  theme_MTS_growth()

p_norm_mts_growth <- p_norm_mts_growth +
  geom_text(
    data = sig_norm_plot_mts_growth,
    aes(
      x = Time,
      y = y_norm,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3
  )+
  guides(
    colour = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
#-----------------------------------------
# plotting slope
#-----------------------------------------
p_slope_mts_growth <- ggplot(
  slopes_mean_mts_growth,
  aes(
    x = Group,
    y = mean_Slope,
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
      ymin = mean_Slope - SEM,
      ymax = mean_Slope + SEM
    ),
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_point(
    data = slopes_mean_mts_growth,
    aes(
      x = Group,
      y = mean_Slope
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
  filename = "mts_growth_raw.jpeg",
  plot = p_raw_mts_growth,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "mts_growth_raw.pdf",
  plot = p_raw_mts_growth,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#2. saving normalized graph
ggsave(
  filename = "mts_growth_norm.jpeg",
  plot = p_norm_mts_growth,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "mts_growth_norm.pdf",
  plot = p_norm_mts_growth,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3. saving slope graph
ggsave(
  filename = "mts_growth_slope.jpeg",
  plot = p_slope_mts_growth,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "mts_growth_slope.pdf",
  plot = p_slope_mts_growth,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#4. saving aov of growth
write.csv(
  mod_data_mts_growth_aov,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "mts_growth_aov_growthrate.csv"),
  row.names = FALSE
)

#5. saving aov of slope
write.csv(
  slope_mts_growth_aov,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "mts_growth_aov_slope.csv"),
  row.names = FALSE
)

#6. saving significance
write.csv(
  sig_norm_mts_growth,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "significance_mts_norm_growth.csv"),
  row.names = FALSE
)

#7. saving significance
write.csv(
  slope_mts_growth_posthoc,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "slope_mts_growth_posthoc.csv"),
  row.names = FALSE
) 

#8. saving significance
write.csv(
  significance_mts_growth_72h,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "emm_mts_growth_72h.csv"),
  row.names = FALSE
)

#9. saving significance
write.csv(
  significance_mts_growth_dose,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "emm_mts_growth_dose.csv"),
  row.names = FALSE
)            

#10. saving significance
write.csv(
  sig_norm_mts_growth,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "emm_sig_norm_mts_growth.csv"),
  row.names = FALSE
)  

#11. saving slopes
write.csv(
  slopes_mean_mts_growth,
  file = file.path(dirs$results, "10_MTS_DoseResponse", "emm_sig_norm_mts_growth_slopes.csv"),
  row.names = FALSE
)  
#-----------------------------------------------
setwd(dirs$code)
