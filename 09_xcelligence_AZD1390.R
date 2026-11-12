#------------------------------------------
# Setting path
#------------------------------------------
setwd(file.path(dirs$xCelligence,"Exp2_ATMinhibition","processed")) #replace with path on your own PC
#-----------------------------------------
# importing and binding of data
#-----------------------------------------
xCell_AZD1390_H23 <- bind_rows(
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

xCell_AZD1390_HCC44 <- bind_rows(
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

xcell_AZD1390_raw <- bind_rows(
  xCell_AZD1390_H23,
  xCell_AZD1390_HCC44
)

xcell_AZD1390_raw <- xcell_AZD1390_raw %>%
  group_by(Cellline, Treatment, Dose, Time) %>%
  mutate(
    median_CI = median(CellIndex, na.rm = TRUE),
    mad_CI = mad(CellIndex, na.rm = TRUE),
    CellIndex = ifelse(
      abs(CellIndex - median_CI) > 3 * mad_CI,
      NA_real_,
      CellIndex
    ),
    CellIndex = ifelse(
      CellIndex < (-0.01)|CellIndex >7,
      NA_real_,
      CellIndex
    )
  ) %>%
  select(-median_CI, -mad_CI) %>%
  ungroup()

xcell_AZD1390_mean<-xcell_AZD1390_raw  %>%
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
xCell_AZD1390_max_ci <- xcell_AZD1390_mean %>%
  group_by(Cellline, Dose, Treatment) %>%
  summarise(
    maxCI = max(mean_CellIndex, na.rm = TRUE),
    .groups = "drop"
  )

# time point od maximal Cell Index
xCell_AZD1390_tmax  <- xcell_AZD1390_mean %>%
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
xcell_AZD1390_growth_rate_data <- xcell_AZD1390_raw %>%
  mutate(
    Dose_num = as.numeric(gsub(" Gy", "", as.character(Dose)))
  ) %>%
  filter(
    Time >= 24,
    Time <= 72
  ) %>%
  group_by(Cellline, Dose, Dose_num, Treatment, Replicate) %>%
  group_modify(~ {
    fit <- lm(log(CellIndex) ~ Time, data = .x)
    tibble(
      k = coef(fit)[["Time"]],
      R2 = summary(fit)$r.squared
    )
  }) %>%
  ungroup()

xcell_AZD1390_res_aov_growth_rate <- aov( k ~ factor(Dose) * Treatment * Cellline, data = xcell_AZD1390_growth_rate_data )

xcell_AZD1390_res_aov_growth_rate_df <- as.data.frame(
  anova(xcell_AZD1390_res_aov_growth_rate)
)

xcell_AZD1390_res_aov_growth_rate$Effect <- rownames(xcell_AZD1390_res_aov_growth_rate)

xcell_AZD1390_res_aov_growth_rate_mean <- xcell_AZD1390_growth_rate_data %>%
  summarise(
    mean_k = mean(k, na.rm = TRUE),
    SD_k   = sd(k, na.rm = TRUE),
    n      = sum(!is.na(k)),
    SEM_k  = SD_k / sqrt(n),
    .by = c(Cellline, Dose, Treatment)
  )
#------------------------------------------
# normalization
#------------------------------------------
# Normalisation to 0 Gy DMSO
xcell_AZD1390_norm_data <- xcell_AZD1390_mean %>%
  group_by(Cellline, Time) %>%
  mutate(
    CI_norm = mean_CellIndex /
      mean_CellIndex[Dose == "0 Gy" & Treatment == "DMSO"][1]
  ) %>%
  ungroup()

xcell_AZD1390_lin_model <- lm(
  mean_CellIndex ~ Time * Dose * Cellline*Treatment,
  data = xcell_AZD1390_norm_data
)

xcell_AZD1390_res_lin_model<-anova(xcell_AZD1390_lin_model)
#------------------------------------------
# Area under the curve (AUC)
#-----------------------------------------
xcell_AZD1390_AUC_control_data <- xcell_AZD1390_raw %>%
  dplyr::filter(Dose == "0 Gy") %>%
  dplyr::select(
    Cellline, Replicate, Treatment, Time,Dose,
    CI_control = CellIndex
  )

# Normalisation + AUC
xcell_AZD1390_AUC_norm_all <- xcell_AZD1390_raw %>%
  dplyr::left_join(
    xcell_AZD1390_AUC_control_data %>%
      dplyr::select(-Dose),
    by = c("Cellline", "Replicate", "Treatment", "Time")
  ) %>%
  dplyr::mutate(
    CI_norm = CellIndex / CI_control
  ) %>%
  dplyr::filter(
    !is.na(Time),
    is.finite(CI_norm)
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
xcell_AZD1390_AUC_norm_all_lm <- lm(
  AUC ~ Dose * Cellline * Treatment,
  data = xcell_AZD1390_AUC_norm_all
)

xcell_AZD1390_AUC_norm_all_AOV_res<- anova(
  xcell_AZD1390_AUC_norm_all_lm
)
#-----------------------------------------
# statistic
#-----------------------------------------
xcell_AZD1390_mean_data <- xcell_AZD1390_raw %>%
  dplyr::filter(
    Dose %in% c("0 Gy", "2 Gy"),
    !is.na(Time)
  ) %>%
  dplyr::group_by(
    Cellline, Dose, Treatment, Time
  ) 

xcell_AZD1390_mod_data <- lm(
  CellIndex ~ Dose * Cellline * Treatment,
  data = xcell_AZD1390_mean_data
)

xcell_AZD1390_res_AOV_data<-anova(xcell_AZD1390_mod_data)

xcell_AZD1390_res_emm <- emmeans(
  xcell_AZD1390_mod_data,
  ~ Treatment |Dose*Cellline
)

xcell_AZD1390_res_emm_data<- pairs(
  xcell_AZD1390_res_emm,
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

xcell_AZD1390_significance <- as.data.frame(
  xcell_AZD1390_res_emm_data
) %>%
  mutate(
    label = case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001  ~ "***",
      p.value < 0.01   ~ "**",
      p.value < 0.05   ~ "*",
      TRUE             ~ "ns"
    ),
    x = 10,
    y = c(rep(5.2,2),rep(5.8,2))
  )
#-----------------------------------------
# plotting raw data
#-----------------------------------------
xcell_AZD1390_plot_data_AZD1390 <- xcell_AZD1390_mean %>%
  dplyr::filter(Dose %in% c("0 Gy", "2 Gy")) %>%
  dplyr::mutate(
    Group = paste(Cellline, Treatment)
  )%>%
  rename(
    Timepoint_2=Timepoint,
    CellIndex=mean_CellIndex,
    Timepoint=Time
  )

xcell_AZD1390_plot_data_0Gy <- xcell_AZD1390_plot_data_AZD1390 %>%
  dplyr::filter(Dose == "0 Gy") 

xcell_AZD1390_p_0Gy_raw <- plot_xcell(
  data = xcell_AZD1390_plot_data_0Gy,
  title = "AZD1390 and DMSO, 0 Gy",
  colours = xcell_AZD_colours,
  breaks = c(
    "NCI-H23 DMSO",
    "NCI-H23 AZD1390",
    "HCC44 DMSO",
    "HCC44 AZD1390"
  ),
  xmax = 135,
  ymax = 6.5,
  ymin = 0,
  x_breaks = c(0, 60, 120),
  x_minor = c(30,90),
  y_breaks = c(0, 3,6),
  y_minor = c(1.5,4.5)
)

xcell_AZD1390_p_0Gy_raw <- xcell_AZD1390_p_0Gy_raw +
  geom_text(
    data = xcell_AZD1390_significance %>%
      filter(Dose == "0 Gy"),
    aes(
      x = x-0.6,
      y = y,
      label =paste(Cellline,"|",contrast,": ", label),
      colour = "black",
      hjust=0
    ),
    inherit.aes = FALSE,
    size = 3
  )+ theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    colour = ggplot2::guide_legend(nrow = 2,
                                   override.aes = list(
                                     size = 3,
                                     stroke = 0.8
                                   )))

xcell_AZD1390_plot_data_2Gy <- xcell_AZD1390_plot_data_AZD1390 %>%
  dplyr::filter(Dose == "2 Gy") 

xcell_AZD1390_p_2Gy_raw <- plot_xcell(
  data = xcell_AZD1390_plot_data_2Gy,
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

xcell_AZD1390_p_2Gy_raw <- xcell_AZD1390_p_2Gy_raw +
  geom_text(
    data = xcell_AZD1390_significance %>%
      filter(Dose == "2 Gy"),
    aes(
      x = x-0.6,
      y = y,
      label =paste(Cellline,"|",contrast,": ", label),
      colour = "black",
      hjust=0
    ),
    inherit.aes = FALSE,
    size = 3
  )+ theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    colour = ggplot2::guide_legend(nrow = 2,
                                   override.aes = list(
                                     size = 3,
                                     stroke = 0.8
                                   )))

xcell_AZD1390_p_final_raw <- (
  xcell_AZD1390_p_0Gy_raw +
    xcell_AZD1390_p_2Gy_raw +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    )
) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(
      size = 8
    ),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm")
  )

#-----------------------------------------  
# Additional Graphes - plotting normalized data
#-----------------------------------------
xcell_AZD1390_plot_data_norm_0Gy <- xcell_AZD1390_norm_data %>%
  dplyr::filter(Dose == "0 Gy") %>%
  dplyr::mutate(
    CellIndex = CI_norm,
    Group = paste(Cellline, Treatment),
    CellIndex = ifelse(
      is.nan(CellIndex) | CellIndex > 2 | CellIndex < 0,
      0,
      CellIndex
    )
  )

xcell_AZD1390_p_norm_0Gy <- plot_xcell(
  data = xcell_AZD1390_plot_data_norm_0Gy,
  title = "AZD1390, 0 Gy, normalized Data",
  colours = xcell_AZD_colours,
  breaks = c(
    "NCI-H23 DMSO",
    "NCI-H23 AZD1390",
    "HCC44 DMSO",
    "HCC44 AZD1390"
  ),
  xmax = 150,
  ymax = 1.25,
  ymin = 0.75,
  x_breaks = c(0, 60, 120),
  x_minor = c(30,90),
  y_breaks = c(0.75,0,1.25),
  y_minor = c(0.8,0.9,1.1,1.2)
)

xcell_AZD1390_plot_data_norm_2Gy <- xcell_AZD1390_norm_data %>%
  dplyr::filter(Dose == "2 Gy")%>%
  dplyr::mutate(
    CellIndex = CI_norm,
    Group = paste(Cellline, Treatment),
    CellIndex = ifelse(
      is.nan(CellIndex) | CellIndex > 2 | CellIndex < 0,
      0,
      CellIndex
    )
  )

xcell_AZD1390_p_norm_2Gy <- plot_xcell(
  data = xcell_AZD1390_plot_data_norm_2Gy,
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
  ymin = 0.75,
  x_breaks = c(0, 60, 120),
  x_minor = c(30,90),
  y_breaks = c(0.5,0,1.25),
  y_minor = c(0.8,0.9,1.1,1.2)
)

xcell_AZD1390_p_norm_0Gy <- xcell_AZD1390_p_norm_0Gy +
  theme(legend.position = "bottom")

xcell_AZD1390_p_norm_2Gy <- xcell_AZD1390_p_norm_2Gy +
  theme(legend.position = "bottom")

xCell_AZD1390_p_norm_final <- (
  xcell_AZD1390_p_norm_0Gy +
    xcell_AZD1390_p_norm_2Gy +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    )
) &
  theme(
    legend.position = "bottom"
  )

#-----------------------------------------
# Additional Graphes - growth rate
#-----------------------------------------
xcell_AZD1390_plot_growth_data <- xcell_AZD1390_res_aov_growth_rate_mean %>%
  dplyr::mutate(
    Dose_num = as.numeric(gsub(" Gy", "", as.character(Dose))),
    Group = paste(Cellline, Treatment)
  )

xcell_AZD1390_p_growth_data <-
  ggplot(
    xcell_AZD1390_plot_growth_data,
    aes(
      x = Dose_num,
      y =  mean_k,
      fill = Group
    )
  ) +
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.8,
    colour = "black"
  )  +
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
# Additional Graphes - plotting AUC
#-----------------------------------------
xcell_AZD1390_plot_norm_AUC <- xcell_AZD1390_AUC_norm_all %>%
  group_by(Cellline, Dose, Treatment) %>%
  summarise(
    mean_AUC_rel = mean(AUC_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Dose_num = as.numeric(gsub(" Gy", "", Dose)),
    Group = paste(Cellline, Treatment)
  )

xcell_AZD1390_p_norm_AUC <-
  ggplot(
    xcell_AZD1390_plot_norm_AUC,
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
    expand = expansion(mult = c(0.05, 0))
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
  file.path(dirs$results, "09_xcelligence_AZD1390"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "09_xcelligence_AZD1390"))

#1. saving combined raw graph
ggsave(
  filename = "09_xCelligence_AZD1390_raw_all.jpeg",
  plot = xcell_AZD1390_p_final_raw,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_AZD1390_raw_all.pdf",
  plot = xcell_AZD1390_p_final_raw,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. saving 0Gy raw graph
ggsave(
  filename = "09_xCelligence_AZD1390_0Gy_raw.jpeg",
  plot = xcell_AZD1390_p_0Gy_raw,
  width = 21*3/6,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_AZD1390_0Gy_raw.pdf",
  plot = xcell_AZD1390_p_0Gy_raw,
  width = 21*3/6,
  height = 7,
  units = "cm"
)

#3. saving 2 Gy raw graph
ggsave(
  filename = "09_xCelligence_AZD1390_2Gy_raw.jpeg",
  plot = xcell_AZD1390_p_2Gy_raw,
  width = 21*3/6,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_AZD1390_2Gy_raw.pdf",
  plot = xcell_AZD1390_p_2Gy_raw,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#4. saving graphes for growthrate
ggsave(
  filename = "09_xCelligence_AZD1390_growthrate.jpeg",
  plot = xcell_AZD1390_p_growth_data,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_AZD1390_growthrate.pdf",
  plot = xcell_AZD1390_p_growth_data,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 
 
#5. saving graph combined with normalised data
ggsave(
  filename = "09_xCelligence_normData_AZD1390_all.jpeg",
  plot = xCell_AZD1390_p_norm_final,
  width = 21*3/6,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_normData_AZD1390_all.pdf",
  plot = xCell_AZD1390_p_norm_final,
  width = 21*2/3,
  height = 7,
  units = "cm"
)  

#6. saving 2 Gy graph with normalized data
ggsave(
  filename = "09_xCelligence_normData_AZD1390_2Gy.jpeg",
  plot = xcell_AZD1390_p_norm_2Gy,
  width = 21*3/6,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_normData_AZD1390_2Gy.pdf",
  plot = xcell_AZD1390_p_norm_2Gy,
  width = 21*2/3,
  height = 7,
  units = "cm"
) 
          
#7. saving 0Gy graph with normalised data
ggsave(
  filename = "09_xCelligence_normData_AZD1390_0Gy.jpeg",
  plot = xcell_AZD1390_p_norm_0Gy,
  width = 21*3/6,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_normData_AZD1390_0Gy.pdf",
  plot = xcell_AZD1390_p_norm_0Gy,
  width = 21*2/3,
  height = 7,
  units = "cm"
) 

#8. saving graphes for AUC
ggsave(
  filename = "09_xCelligence_AZD1390_AUC.jpeg",
  plot = xcell_AZD1390_p_norm_AUC,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "09_xCelligence_AZD1390_AUC.pdf",
  plot = xcell_AZD1390_p_norm_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#9. saving tmax
write.csv(
  xCell_AZD1390_tmax,
  file = file.path(dirs$results, "09_xcelligence_AZD1390", "09_xcelligence_growth_tmax_AZD1390.csv"),
  row.names = FALSE
)

#10. saving AOV with normaliased data
write.csv(
  broom::tidy(xcell_AZD1390_AUC_norm_all_AOV_res),
  file = file.path(
    dirs$results,
    "09_xcelligence_AZD1390",
    "09_xcelligence_AZD1390_AOV_norm.csv"
  ),
  row.names = FALSE
)

#11. saving AOV growth model
write.csv(
  (xcell_AZD1390_res_aov_growth_rate_df),
  file = file.path(
    dirs$results,
    "09_xcelligence_AZD1390",
    "09_xcelligence_AZD1390_AOV_growthrate.csv"
  ),
  row.names = FALSE
)

#12. saving growth rates
write.csv(xcell_AZD1390_res_aov_growth_rate_mean,
  file = file.path(
    dirs$results,
    "09_xcelligence_AZD1390",
    "09_xcelligence_AZD1390_growthrate.csv"
  ),
  row.names = FALSE
)

#13. saving growth rates statistics
write.csv(xcell_AZD1390_res_lin_model,
  file = file.path(
    dirs$results,
    "09_xcelligence_AZD1390",
    "09_xcelligence_AZD1390_growthrate_AOV.csv"
  ),
  row.names = FALSE
)
#------------------------------------------
setwd(dirs$code)
