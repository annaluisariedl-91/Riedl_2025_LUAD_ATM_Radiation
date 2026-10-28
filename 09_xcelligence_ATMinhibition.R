#------------------------------------------
# Setting path
#------------------------------------------
setwd(file.path(base_dir,"functional_assays","xCelligence","Exp2_ATMinhibition","processed")) #replace with path on your own PC
#-----------------------------------------
# importing and binding of data
#-----------------------------------------
xCell_H23_AZD <- bind_rows(
  prepare_xcell_AZD(
    "extracted_mean_H23_rep1.csv",
    "NCI-H23", 1,
    cols = c(2, 9:3),
    doses = c("2Gy", "05Gy", "0Gy")
  ),
  prepare_xcell_AZD(
    "extracted_mean_H23_rep2.csv",
    "NCI-H23", 2,
    cols = 2:9,
    doses = c("2Gy", "05Gy", "0Gy")
  ),
  prepare_xcell_AZD(
    "extracted_mean_H23_rep3.csv",
    "NCI-H23", 3,
    cols = c(2, 4:9, 3),
    doses = c("2Gy", "05Gy", "0Gy")
  ),
  prepare_xcell_AZD(
    "extracted_mean_H23_rep4.csv",
    "NCI-H23", 4,
    cols = 2:9,
    doses = c("2Gy", "05Gy", "0Gy")
  )
)

xCell_HCC44_AZD <- bind_rows(
  prepare_xcell_AZD(
    "extracted_mean_HCC44_rep1.csv",
    "HCC44", 1,
    cols = c(2, 9:3),
    doses = c("4Gy", "2Gy", "0Gy")
  ),
  prepare_xcell_AZD(
    "extracted_mean_HCC44_rep2.csv",
    "HCC44", 2,
    cols = 2:9,
    doses = c("4Gy", "2Gy", "0Gy")
  ),
  prepare_xcell_AZD(
    "extracted_mean_HCC44_rep3.csv",
    "HCC44", 3,
    cols = 2:9,
    doses = c("4Gy", "2Gy", "0Gy")
  )
)

xcell_AZD_raw <- bind_rows(
  xCell_H23_AZD,
  xCell_HCC44_AZD
)

xcell_AZD_raw <- xcell_AZD_raw%>%
  filter(
    CellIndex >= -10,
    CellIndex <= 10,
  ) 

xcell_AZD_mean<-xcell_AZD_raw  %>%
  group_by(Cellline, Dose, Treatment, Timepoint) %>%
  summarise(
    n = sum(!is.na(CellIndex)),
    mean_CellIndex = mean(CellIndex, na.rm = TRUE),
    sd_CellIndex = sd(CellIndex, na.rm = TRUE),
    sem_CellIndex = sd_CellIndex / sqrt(n),
    Time=mean(Time, na.rm = TRUE),
    .groups = "drop"
  )
#------------------------------------------
# Benchmark
#------------------------------------------
# Maximal Cell Index
xCell_max_ci_AZD1390 <- xcell_AZD_mean %>%
  group_by(Cellline, Dose, Treatment) %>%
  summarise(
    maxCI = max(mean_CellIndex, na.rm = TRUE),
    .groups = "drop"
  )

# time point od maximal Cell Index
xCell_tmax_AZD1390  <- xcell_AZD_mean %>%
  group_by(Cellline, Dose, Treatment) %>%
  slice_max(
    mean_CellIndex,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()
#------------------------------------------
# growth rate
#------------------------------------------
xcell_growth_data_AZD1390 <-xcell_AZD_mean %>%
  mutate(
    Dose_num = as.numeric(gsub(" Gy", "", as.character(Dose)))
  ) %>%
  filter(
    Time >= 24,
    Time <= 72,
    mean_CellIndex > 0
  ) %>%
  group_by(Cellline, Dose, Dose_num, Treatment) %>%
  group_modify(~ {
    fit <- lm(log(mean_CellIndex ) ~ Time, data = .x)
    tibble(
      k = coef(fit)[["Time"]],
      R2 = summary(fit)$r.squared
    )
  }) %>%
  ungroup()

xcell_ATMi_growth_rate <- lm(
  k ~ Dose * Cellline,
  data = xcell_growth_data_AZD1390
)
#------------------------------------------
# normalization
#------------------------------------------
# Normalisation to 0 Gy DMSO
xcell_ATMi_norm_data_AZD1390 <- xcell_AZD_mean %>%
  group_by(Cellline, Time) %>%
  mutate(
    CI_norm = mean_CellIndex /
      mean_CellIndex[Dose == "0 Gy" & Treatment == "DMSO"][1]
  ) %>%
  ungroup()

# split cell line and dose
xcell_ATMi_split_data_AZD1390 <- xcell_ATMi_norm_data_AZD1390 %>%
  group_by(Cellline, Dose, Treatment) %>%
  group_split() %>%
  setNames(
    xcell_ATMi_norm_data_AZD1390 %>%
      distinct(Cellline, Dose, Treatment) %>%
      transmute(
        name = paste(Cellline, Dose, sep = "_")
      ) %>%
      pull(name)
  )

xcell_ATMi_lin_model_AZD1390 <- lm(
  mean_CellIndex ~ Time * Dose * Cellline*Treatment,
  data = xcell_AZD_mean
)

xcelligance_AOV_AZD1390<-anova(xcell_ATMi_lin_model_AZD1390)
xcelligance_sammary_AZD1390<-summary(xcell_ATMi_lin_model_AZD1390)
#------------------------------------------
# Area under the curve (AUC)
#-----------------------------------------
xcell_ATMi_control_data <- xcell_AZD_raw %>%
  dplyr::filter(Dose == "0 Gy") %>%
  dplyr::select(
    Cellline, Replicate, Treatment, Time,
    CI_control = CellIndex
  )

# Normalisation + AUC
xcell_ATMi_auc_data_AZD1390_all <- xcell_AZD_raw %>%
  dplyr::left_join(
    xcell_ATMi_control_data,
    by = c("Cellline", "Replicate", "Treatment", "Time")
  ) %>%
  dplyr::mutate(
    CI_norm = CellIndex / CI_control
  ) %>%
  dplyr::filter(
    !is.na(Time),
    !is.na(CI_norm)
  ) %>%
  dplyr::arrange(
    Cellline, Replicate, Treatment, Dose, Time
  ) %>%
  dplyr::group_by(
    Cellline, Replicate, Dose, Treatment
  ) %>%
  dplyr::summarise(
    AUC = pracma::trapz(Time, CI_norm),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    Cellline, Replicate, Treatment
  ) %>%
  dplyr::mutate(
    AUC_rel = AUC / AUC[Dose == "0 Gy"][1]
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Dose_num = as.numeric(
      gsub(" Gy", "", as.character(Dose))
    )
  )

# AUC model
xcell_ATMi_mod_auc_AZD1390_all <- lm(
  AUC ~ Dose * Cellline * Treatment,
  data = xcell_ATMi_auc_data_AZD1390_all
)

xcelligance_AUC_AOV_AZD1390 <- anova(
  xcell_ATMi_mod_auc_AZD1390_all
)

summary(xcell_ATMi_mod_auc_AZD1390_all)

#-----------------------------------------
# statistic
#-----------------------------------------
auc_xcell_ATMi_raw <- xcell_AZD_raw %>%
  dplyr::filter(
    Dose %in% c("0 Gy", "2 Gy"),
    !is.na(Time),
    !is.na(CellIndex)
  ) %>%
  dplyr::arrange(
    Cellline, Replicate, Dose, Treatment, Time
  ) %>%
  dplyr::group_by(
    Cellline, Replicate, Dose, Treatment
  ) %>%
  dplyr::summarise(
    AUC = pracma::trapz(Time, CellIndex),
    .groups = "drop"
  )

mod_auc_xcell_ATMi <- lm(
  AUC ~ Dose * Cellline * Treatment,
  data = auc_xcell_ATMi_raw
)

mod_AOV_xcell_ATMi<-anova(mod_auc_xcell_ATMi)

emm_xcell_ATMi <- emmeans(
  mod_auc_xcell_ATMi,
  ~ Treatment | Cellline * Dose
)

posthoc_xcell_ATMi<- pairs(
  emm_xcell_ATMi,
  adjust = "holm"
) %>%
  as.data.frame()
#-----------------------------------------
# preparing for plot
#-----------------------------------------
xcell_AZD_colours <- c(
  "HCC44 DMSO"           = "#4294F7",
  "HCC44 AZD1390"        = "#EA4025",      
  "NCI-H23 DMSO"         = "#5F8F29",
  "NCI-H23 AZD1390"      = "#F09837"
)

growth_colours <- c(
  "HCC44"     = "#3288BD",
  "NCI-H23"   = "#ABDDA4"
)
#-----------------------------------------
# plotting raw data
#-----------------------------------------
plot_data_xcell_ATMi <- xcell_AZD_mean %>%
  dplyr::filter(Dose %in% c("0 Gy", "2 Gy")) %>%
  dplyr::mutate(
    CellIndex = mean_CellIndex,
    Group = paste(Cellline, Treatment)
  )

plot_data_xcell_ATMi_0Gy <- plot_data_xcell_ATMi %>%
  dplyr::filter(Dose == "0 Gy") %>%
  mutate(
    CellIndex=mean_CellIndex,
    Group = paste(Cellline, Treatment)
  )
  

p_xcell_ATMi_0Gy_raw <- plot_xcell(
  data = plot_data_xcell_ATMi_0Gy,
  title = "AZD1390 and DMSO, 0 Gy",
  colours = xcell_AZD_colours,
  breaks = c(
    "NCI-H23 DMSO",
    "NCI-H23 AZD1390",
    "HCC44 DMSO",
    "HCC44 AZD1390"
  ),
  xmax = 125,
  ymax = 6.5,
  ymin = 0,
  x_breaks = c(0, 60, 120),
  x_minor = c(30,90),
  y_breaks = c(0, 3,6),
  y_minor = c(1.5,4.5)
)

plot_data_xcell_ATMi_2Gy <- xcell_AZD_mean %>%
  dplyr::filter(Dose == "2 Gy")%>%
  mutate(
    CellIndex=mean_CellIndex,
    Group = paste(Cellline, Treatment)
  )

p_xcell_ATMi_2Gy_raw <- plot_xcell(
  data = plot_data_xcell_ATMi_2Gy,
  title = "AZD1390 and DMSO, 2 Gy",
  colours = xcell_AZD_colours,
  breaks = c(
    "NCI-H23 DMSO",
    "NCI-H23 AZD1390",
    "HCC44 DMSO",
    "HCC44 AZD1390"
  ),
  xmax = 125,
  ymax = 6.5,
  ymin = 0,
  x_breaks = c(0, 60, 120),
  x_minor = c(30,90),
  y_breaks = c(0, 3,6),
  y_minor = c(1.5,4.5)
)

xCell_ATMi_raw_final <- (
  p_xcell_ATMi_0Gy_raw +
    p_xcell_ATMi_2Gy_raw +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    )
) &
  theme(
    legend.position = "bottom"
  )
#-----------------------------------------  
# plotting normalized data
#-----------------------------------------
plot_data_norm_xcell_ATMi_0Gy <- xcell_ATMi_norm_data_AZD1390 %>%
  dplyr::filter(Dose == "0 Gy")%>%
  mutate(
    CellIndex = CI_norm,
    Group = paste(Cellline, Treatment)
  )

p_xcell_norm_ATMi_0Gy <- plot_xcell(
  data = plot_data_norm_xcell_ATMi_0Gy,
  title = "AZD1390, 0 Gy, normalized Data",
  colours = xcell_AZD_colours,
  breaks = c(
    "NCI-H23 DMSO",
    "NCI-H23 AZD1390",
    "HCC44 DMSO",
    "HCC44 AZD1390"
  ),
  xmax = 125,
  ymax = 1.25,
  ymin = 0.75,
  x_breaks = c(0, 60, 120),
  x_minor = c(30,90),
  y_breaks = c(0.75,0,1.25),
  y_minor = c(0.8,0.9,1.1,1.2)
)

plot_data_norm_xcell_ATMi_2Gy <- xcell_ATMi_norm_data_AZD1390 %>%
  dplyr::filter(Dose == "2 Gy")%>%
  mutate(
    CellIndex = CI_norm,
    Group = paste(Cellline, Treatment)
  )

p_xcell_norm_ATMi_2Gy <- plot_xcell(
  data = plot_data_norm_xcell_ATMi_2Gy,
  title = "AZD1390, 2 Gy, normalized Data",
  colours = xcell_AZD_colours,
  breaks = c(
    "NCI-H23 DMSO",
    "NCI-H23 AZD1390",
    "HCC44 DMSO",
    "HCC44 AZD1390"
  ),
  xmax = 125,
  ymax = 1.25,
  ymin = 0,
  x_breaks = c(0, 60, 120),
  x_minor = c(30,90),
  y_breaks = c(0.5,0,1.25),
  y_minor = c(0.8,0.9,1.1,1.2)
)

p_xcell_norm_ATMi_0Gy <- p_xcell_norm_ATMi_0Gy +
  theme(legend.position = "right")

p_xcell_norm_ATMi_2Gy <- p_xcell_norm_ATMi_2Gy +
  theme(legend.position = "right")

xCell_norm_ATMi_final <- (
  p_xcell_norm_ATMi_0Gy +
    p_xcell_norm_ATMi_2Gy +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    )
) &
  theme(
    legend.position = "bottom"
  )
#-----------------------------------------
# growth rate
#-----------------------------------------
growth_data_xcell_ATMi <- xcell_growth_data_AZD1390 %>%
  mutate(
    Group = paste(Cellline, Treatment)
  )

p_xcell_growth_xcell_ATMi <-
  ggplot(
    growth_data_xcell_ATMi,
    aes(
      x = Dose_num,
      y = k,
      fill = Group
    )
  ) +
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.8,
    colour = "black"
  ) +
  scale_fill_manual(
    values = xcell_AZD_colours,
    breaks = c(
      "NCI-H23 DMSO",
      "NCI-H23 AZD1390",
      "HCC44 DMSO",
      "HCC44 AZD1390"
    )) +
  scale_x_continuous(
    breaks = c(0, 1, 2),
    minor_breaks = c(0.5, 1.5),
    limits = c(-0.1, 4.25),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(-0.01,0, 0.03),
    minor_breaks = c(-0.005, 0.005,0.01,0.015,0.02,0.025),
    limits = c(-0.01, 0.03),
    expand = expansion(mult = c(0, 0)),
    labels = function(x) {
      ifelse(x == 0, "0", sprintf("%.2f", x))
    }
  ) +
  labs(
    title = "Growth rate",
    x = "Dose [Gy]",
    y = "Growth rate [a.u.]",
    fill = NULL
  ) +
  theme_prism_xcell() +
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 3,
        shape = 21,
        colour = "black",
        stroke = 0.8
      )
    )
  )+
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(
      t = 5,
      r = 10,
      b = 5,
      l = 5
    )
  )
#-----------------------------------------
# plotting AUC
#-----------------------------------------
xcell_ATMi_auc_data_AZD1390_mean <- xcell_ATMi_auc_data_AZD1390_all %>%
  group_by(Cellline, Dose, Treatment) %>%
  summarise(
    mean_AUC_rel = mean(AUC_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Dose_num = as.numeric(gsub(" Gy", "", Dose)),
    Group = paste(Cellline, Treatment)
  )

p_xcell_ATMi_AUC <-
  ggplot(
    xcell_ATMi_auc_data_AZD1390_mean,
    aes(
      x = Dose_num,
      y = mean_AUC_rel,
      fill = Group
    )
  ) +
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.8,
    colour = "black"
  ) +
  scale_fill_manual(
    values = xcell_AZD_colours,
    breaks = c(
      "NCI-H23 DMSO",
      "NCI-H23 AZD1390",
      "HCC44 DMSO",
      "HCC44 AZD1390"
    ))+
  scale_x_continuous(
    breaks = c(0, 1, 2),
    minor_breaks = c(0.5, 1.5),
    limits = c(0, 2.25),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(0,1,2),
    minor_breaks = c(0.5,1.5),
    limits = c(0, 2),
    expand = expansion(mult = c(0, 0)),
    labels = function(x) {
      ifelse(x == 0, "0", sprintf("%.2f", x))
    }
  ) +
  labs(
    title = "Growth rate",
    x = "Dose [Gy]",
    y = "relative AUC [a.u.]",
    fill = NULL
  ) +
  theme_prism_xcell() +
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 3,
        shape = 21,
        colour = "black",
        stroke = 0.8
      )
    )
  )+
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = margin(
      t = 5,
      r = 10,
      b = 5,
      l = 5
    )
  )
#-----------------------------------------
# Saving results
#-----------------------------------------
dir.create(
  file.path(dirs$results, "09_xcelligence_ATMinhibition"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "09_xcelligence_ATMinhibition"))
        
#1. saving combined raw graph
ggsave(
  filename = "xCelligence_doseResponse_AZD1390_raw_all.jpeg",
  plot = xCell_ATMi_raw_final,
  width = 18,
  height = 8,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AZD1390_raw_all.pdf",
  plot = xCell_ATMi_raw_final,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. saving 0Gy raw graph
ggsave(
  filename = "xCelligence_doseResponse_AZD1390_0Gy_raw.jpeg",
  plot = p_xcell_ATMi_0Gy_raw,
  width = 18,
  height = 8,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AZD1390_0Gy_raw.pdf",
  plot = p_xcell_ATMi_0Gy_raw,
  width = 21*2/3,
  height = 7,
  units = "cm"
)
#3. saving 2 Gy raw graph
ggsave(
  filename = "xCelligence_doseResponse_AZD1390_2Gy_raw.jpeg",
  plot = p_xcell_ATMi_2Gy_raw,
  width = 18,
  height = 8,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AZD1390_2Gy_raw.pdf",
  plot = p_xcell_ATMi_2Gy_raw,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#4. saving graphes for growthrate
ggsave(
  filename = "xCelligence_doseResponse_AZD1390_growthrate.jpeg",
  plot = p_xcell_growth_xcell_ATMi,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AZD1390_growthrate.pdf",
  plot = p_xcell_growth_xcell_ATMi,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 
 
#5. saving graph combined with normalised data
ggsave(
  filename = "xCelligence_doseResponse_normData_AZD1390_all.jpeg",
  plot = xCell_norm_ATMi_final,
  width = 18,
  height = 8,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_normData_AZD1390_all.pdf",
  plot = xCell_norm_ATMi_final,
  width = 21*2/3,
  height = 7,
  units = "cm"
)  
#6. saving 2 Gy graph with normalized data
ggsave(
  filename = "xCelligence_doseResponse_normData_AZD1390_2Gy.jpeg",
  plot = p_xcell_norm_ATMi_2Gy,
  width = 18,
  height = 8,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_normData_AZD1390_2Gy.pdf",
  plot = p_xcell_norm_ATMi_2Gy,
  width = 21*2/3,
  height = 7,
  units = "cm"
) 
#7. saving 0Gy graph with normalised data
ggsave(
  filename = "xCelligence_doseResponse_normData_AZD1390_0Gy.jpeg",
  plot = p_xcell_norm_ATMi_0Gy,
  width = 18,
  height = 8,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_normData_AZD1390_0Gy.pdf",
  plot = p_xcell_norm_ATMi_0Gy,
  width = 21*2/3,
  height = 7,
  units = "cm"
) 
#8. saving graphes for AUC
ggsave(
  filename = "xCelligence_doseResponse_AZD1390_AUC.jpeg",
  plot = p_xcell_ATMi_AUC,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AZD1390_AUC.pdf",
  plot = p_xcell_ATMi_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#9. saving tmax
write.csv(
  xCell_tmax_AZD1390,
  file = file.path(dirs$results, "09_xcelligence_ATMinhibition", "xcelligence_growth_tmax_AZD1390.csv"),
  row.names = FALSE
)

#10. saving linear growth model
write.csv(
  broom::tidy(xcell_growth_data_AZD1390),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_growth_linearModel_AZD1390.csv"
  ),
  row.names = FALSE
)

#11. saving linear growth model
write.csv(
  broom::tidy(xcell_ATMi_growth_rate),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_growth_linearModel_AZD1390.csv"
  ),
  row.names = FALSE
)

#12. saving post hoc analysis
write.csv(posthoc_xcell_ATMi,
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_AZD1390_posthocanalysis.csv"
  ),
  row.names = FALSE
)

#13. saving ANOVA
write.csv(
  broom::tidy(xcelligance_AOV_AZD1390),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_AZD1390_ANOVA.csv"
  ),
  row.names = FALSE
)

#14. saving ANOVA
write.csv(
  broom::tidy(mod_AOV_xcell_ATMi),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_AZD1390_ANOVA.csv"
  ),
  row.names = FALSE
)

#15. saving analysis
write.csv(
  broom::tidy(xcelligance_sammary_AZD1390),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_AZD1390_analysis.csv"
  ),
  row.names = FALSE
)

#16. saving analysis AUC
write.csv(
  broom::tidy(xcell_ATMi_mod_auc_AZD1390_all),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_AZD1390_analysis_AUC.csv"
  ),
  row.names = FALSE
)  

#17. saving analysis AUC
write.csv(
  broom::tidy(xcelligance_AUC_AOV_AZD1390),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_AZD1390_analysis_AUC.csv"
  ),
  row.names = FALSE
) 
 
#18. saving ANOVA AUC
write.csv(
  broom::tidy(xcelligance_AUC_AOV_AZD1390),
  file = file.path(
    dirs$results,
    "09_xcelligence_ATMinhibition",
    "xcelligence_AZD1390_AOVs_AUC.csv"
  ),
  row.names = FALSE
)    

#------------------------------------------
setwd(dirs$code)
