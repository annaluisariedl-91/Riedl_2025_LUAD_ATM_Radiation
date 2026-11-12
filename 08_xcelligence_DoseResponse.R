#------------------------------------------
# Setting path
#------------------------------------------
setwd(file.path(dirs$xCelligence,"Exp1_Dose response","processed")) #replace with path on your own PC
#original xlxs format was transformed to .csv for easier importing, original is available at /raw
Doses <- c("2 Gy", "1.5 Gy", "1 Gy", "0.5 Gy", "0 Gy","Medium")
#------------------------------------------
# Binding and preparing data
#------------------------------------------
xCell_DoseResponse_Data_H23 <- prepare_xCell(
  "extracted_H23_rep1.csv"
)

xCell_DoseResponse_Data_HCC44 <- prepare_xCell(
  "extracted_HCC44_rep1.csv"
)

xCell_DoseResponse_stat_data_all <- bind_rows(
  xCell_DoseResponse_Data_H23   %>% mutate(Cellline = "NCI-H23"),
  xCell_DoseResponse_Data_HCC44 %>% mutate(Cellline = "HCC44")
) %>%
  pivot_longer(
    cols = all_of(Doses),
    names_to = "Dose",
    values_to = "CellIndex",
    values_transform = list(
      CellIndex = ~ as.numeric(gsub(",", ".", as.character(.x)))
    )
  ) %>%
  mutate(
    Cellline = factor(Cellline),
    Dose = factor(Dose, levels = Doses)
  )
#------------------------------------------
# Benchmark
#------------------------------------------
# Maximal Cell Index
xCell_DoseResponse_max_ci <- xCell_DoseResponse_stat_data_all %>%
  group_by(Cellline, Dose) %>%
  summarise(
    maxCI = max(CellIndex, na.rm = TRUE),
    .groups = "drop"
  )

# time point od maximal Cell Index
xCell_DoseResponse_tmax <- xCell_DoseResponse_stat_data_all %>%
  group_by(Cellline, Dose) %>%
  slice_max(
    CellIndex,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()
#------------------------------------------
# Growth rate
#------------------------------------------
xCell_DoseResponse_growth_rate_data <-xCell_DoseResponse_stat_data_all %>%
  mutate(
    Dose_num = as.numeric(gsub(" Gy", "", as.character(Dose)))
  ) %>%
  filter(
    Time >= 24,
    Time <= 72,
    CellIndex > 0
  ) %>%
  group_by(Cellline, Dose, Dose_num) %>%
  group_modify(~ {
    fit <- lm(log(CellIndex) ~ Time, data = .x)
    tibble(
      k = coef(fit)[["Time"]],
      R2 = summary(fit)$r.squared
    )
  }) %>%
  ungroup()

xCell_DoseResponse_Grwothrate_rate_lm <- lm(
  k ~ Dose_num * Cellline,
  data = xCell_DoseResponse_growth_rate_data
)

xCell_DoseResponse_significance_Growth_rate_lm<-emmeans(xCell_DoseResponse_Grwothrate_rate_lm, pairwise ~ Cellline | Dose_num, adjust = "tukey")
#------------------------------------------
# normalization
#------------------------------------------
# Normalisation to 0 Gy
xCell_DoseResponse_norm_data <- xCell_DoseResponse_stat_data_all %>%
  group_by(Cellline, Time) %>%
  mutate(
    CI_norm = CellIndex / CellIndex[Dose == "0 Gy"][1]
  ) %>%
  ungroup()

# split cell line and dose
xCell_DoseResponse_split_data <- xCell_DoseResponse_norm_data %>%
  group_by(Cellline, Dose) %>%
  group_split() %>%
  setNames(
    xCell_DoseResponse_norm_data %>%
      distinct(Cellline, Dose) %>%
      transmute(
        name = paste(Cellline, Dose, sep = "_")
      ) %>%
      pull(name)
  )

xCell_DoseResponse_norm_stat_lm <- lm(
  CellIndex ~ Time * Dose * Cellline,
  data = xCell_DoseResponse_stat_data_all
)

xCell_DoseResponse_significance_norm_lm<-emmeans(xCell_DoseResponse_norm_stat_lm, pairwise ~ Cellline | Dose, adjust = "tukey")
#------------------------------------------
# Area under the curve (AUC)
#------------------------------------------
xCell_DoseResponse_auc_data_all <- xCell_DoseResponse_norm_data %>%
  group_by(Cellline, Dose) %>%
  summarise(
    AUC = pracma::trapz(Time, CI_norm),
    .groups = "drop"
  )%>%
  filter(
    Dose!="Medium"
  )

# Model for AUC
xCell_DoseResponse_mod_auc_all <- lm(
  AUC ~ Dose * Cellline,
  data = xCell_DoseResponse_auc_data_all
)

xCell_DoseResponse_significance_auc_all<-emmeans(xCell_DoseResponse_mod_auc_all, pairwise ~ Cellline | Dose, adjust = "tukey")

# AUC relative to 0 Gy control
xCell_DoseResponse_auc_data_all <- xCell_DoseResponse_auc_data_all %>%
  group_by(Cellline) %>%
  mutate(
    AUC_rel = AUC / AUC[Dose == "0 Gy"][1]
  ) %>%
  ungroup()

xCell_DoseResponse_auc_data_all <- xCell_DoseResponse_auc_data_all%>%
  mutate(
    Dose_num = as.numeric(gsub(" Gy", "", as.character(Dose)))
  )
#-----------------------------------------
# preparing for plot
#-----------------------------------------
xcell_colours <- c(
  "0 Gy" = "#4294F7",
  "0.5 Gy"   = "#5F8F29",
  "1 Gy"     = "#F09837",
  "1.5 Gy"   = "#EA4025",
  "2 Gy"     = "#872111",
  "Medium"  = "grey70"
)

DoseResponse_colours <- c(
  "HCC44" = "#4294F7",
  "NCI-H23"   = "#5F8F29"
)
#-----------------------------------------  
# plotting raw data
#-----------------------------------------
xcell_DoseResponse_plot_data_H23 <- xCell_DoseResponse_stat_data_all %>%
  filter(Cellline == "NCI-H23") %>%
  mutate(
    Group = Dose,
    Timepoint=Time
  )

xcell_DoseResponse_p_raw_H23 <- plot_xcell(
  data = xcell_DoseResponse_plot_data_H23,
  title = "NCI-H23",
  colours = xcell_colours,
  breaks = c(
    "0 Gy",
    "0.5 Gy",
    "1 Gy",
    "1.5 Gy",
    "2 Gy"
  ),
  xmax = 175,
  ymax = 4.5,
  ymin = 0,
  x_breaks = c(0, 100, 175),
  x_minor = c(50, 150),
  y_breaks = c(0, 2, 4),
  y_minor = c(1, 3)
)

xcell_DoseResponse_plot_data_HCC44 <- xCell_DoseResponse_stat_data_all %>%
  filter(Cellline == "HCC44") %>%
  mutate(
    Group = Dose,
    Timepoint=Time
    )

xcell_DoseResponse_p_raw_HCC44 <- plot_xcell(
  data = xcell_DoseResponse_plot_data_HCC44,
  title = "HCC44",
  colours = xcell_colours,
  breaks = c(
    "0 Gy",
    "0.5 Gy",
    "1 Gy",
    "1.5 Gy",
    "2 Gy"
  ),
  xmax = 175,
  ymax = 4.5,
  ymin = 0,
  x_breaks = c(0, 100, 175),
  x_minor = c(50, 150),
  y_breaks = c(0, 2, 4),
  y_minor = c(1, 3)
)

xcell_DoseResponse_p_raw_H23 <- xcell_DoseResponse_p_raw_H23 +
  theme(legend.position = "right")

xcell_DoseResponse_p_raw_HCC44 <- xcell_DoseResponse_p_raw_HCC44 +
  theme(legend.position = "right")

xCell_DoseResponse_p_raw_final <- (
  xcell_DoseResponse_p_raw_H23 +
    xcell_DoseResponse_p_raw_HCC44 +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    )
) &
  theme(
    legend.position = "bottom"
  )
#-----------------------------------------  
# Additional Graphs - plotting normalized data
#-----------------------------------------
xcell_DoseResponse_plot_norm_data_H23 <- xCell_DoseResponse_norm_data %>%
  filter(Cellline == "NCI-H23") %>%
  mutate(
    Group = Dose,
    Timepoint=Time
  )
xcell_DoseResponse_plot_norm_data_H23<-xcell_DoseResponse_plot_norm_data_H23[,-5]
colnames(xcell_DoseResponse_plot_norm_data_H23)[5]<-"CellIndex"

xcell_DoseResponse_p_norm_H23 <- plot_xcell(
  data = xcell_DoseResponse_plot_norm_data_H23,
  title = "NCI-H23 - normalized Data",
  colours = xcell_colours,
  breaks = c(
    "0 Gy",
    "0.5 Gy",
    "1 Gy",
    "1.5 Gy",
    "2 Gy"
  ),
  xmax = 175,
  ymax = 1.5,
  ymin = 0,
  x_breaks = c(0, 100, 175),
  x_minor = c(50, 150),
  y_breaks = c(0, 2, 4),
  y_minor = c(1, 3)
)

xcell_DoseResponse_plot_norm_data_HCC44 <- xCell_DoseResponse_norm_data %>%
  filter(Cellline == "HCC44") %>%
  mutate(
    Group = Dose,
    Timepoint=Time
  )

xcell_DoseResponse_plot_norm_data_HCC44<-xcell_DoseResponse_plot_norm_data_HCC44[,-5]
colnames(xcell_DoseResponse_plot_norm_data_HCC44)[5]<-"CellIndex"

xcell_DoseResponse_p_norm_HCC44 <- plot_xcell(
  data = xcell_DoseResponse_plot_norm_data_HCC44,
  title = "HCC44 - normalized Data",
  colours = xcell_colours,
  breaks = c(
    "0 Gy",
    "0.5 Gy",
    "1 Gy",
    "1.5 Gy",
    "2 Gy"
  ),
  xmax = 175,
  ymax =1.5,
  ymin = 0,
  x_breaks = c(0, 100, 175),
  x_minor = c(50, 150),
  y_breaks = c(0, 2, 4),
  y_minor = c(1, 3)
)

xcell_DoseResponse_p_norm_H23 <- xcell_DoseResponse_p_norm_H23 +
  theme(legend.position = "right")

xcell_DoseResponse_p_norm_HCC44 <-xcell_DoseResponse_p_norm_HCC44 +
  theme(legend.position = "right")

xCell_DoseResponse_p_norm_final <- (
  xcell_DoseResponse_p_norm_H23 +
    xcell_DoseResponse_p_norm_HCC44 +
    patchwork::plot_layout(
      ncol = 2,
      guides = "collect"
    )
) &
  theme(
    legend.position = "bottom"
  )
#-----------------------------------------
# Additional Graphes - Growth rate
#-----------------------------------------
xCell_DoseResponse_p_Growth <-
  ggplot(
    xCell_DoseResponse_growth_rate_data,
    aes(
      x = Dose_num,
      y = k,
      fill = Cellline
    )
  ) +
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.8,
    colour = "black"
  ) +
  scale_fill_manual(
    values = DoseResponse_colours,
    breaks = c("HCC44", "NCI-H23")
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2),
    minor_breaks = c(0.5, 1.5),
    limits = c(0, 2.25),
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
xcell_DoseResponse_p_AUC <-
  ggplot(
    xCell_DoseResponse_auc_data_all,
    aes(
      x = Dose_num,
      y = AUC_rel,
      fill = Cellline
    )
  ) +
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.8,
    colour = "black"
  ) +
  scale_fill_manual(
    values = DoseResponse_colours,
    breaks = c("HCC44", "NCI-H23")
  ) +
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
    title = "DoseResponse rate",
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
  file.path(dirs$results, "08_xcelligance_DoseResponse"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "08_xcelligance_DoseResponse"))

#1. saving combined raw graph
ggsave(
  filename = "08_xCelligence_doseResponse_raw_all.jpeg",
  plot = xCell_DoseResponse_p_raw_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_raw_all.pdf",
  plot = xCell_DoseResponse_p_raw_final,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. saving NCI-H23 raw graph
ggsave(
  filename = "08_xCelligence_doseResponse_raw_H23.jpeg",
  plot = xcell_DoseResponse_p_raw_H23,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_raw_H23.pdf",
  plot = xcell_DoseResponse_p_raw_H23,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#3. saving HCC44 raw graph
ggsave(
  filename = "08_xCelligence_doseResponse_raw_HCC44.jpeg",
  plot = xcell_DoseResponse_p_raw_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_raw_HCC44.pdf",
  plot = xcell_DoseResponse_p_raw_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm"
)
            
#4. saving combined normalised graphes 
ggsave(
  filename = "08_xCelligence_doseResponse_norm_all.jpeg",
  plot = xCell_DoseResponse_p_norm_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_norm_all.pdf",
  plot = xCell_DoseResponse_p_norm_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#5. saving combined normalised H23 graphes 
ggsave(
  filename = "08_xCelligence_doseResponse_norm_H23.jpeg",
  plot = xcell_DoseResponse_p_norm_H23,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_norm_H23.pdf",
  plot = xcell_DoseResponse_p_norm_H23,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#6. saving HCC44 combined normalied graphes 
ggsave(
  filename = "08_xCelligence_doseResponse_norm_HCC44.jpeg",
  plot = xcell_DoseResponse_p_norm_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_norm_HCC44.pdf",
  plot = xcell_DoseResponse_p_norm_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#7. saving growth rate graph for data
ggsave(
  filename = "08_xCelligence_doseResponse_GrowthRate.jpeg",
  plot = xCell_DoseResponse_p_Growth,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_GrowthRate.pdf",
  plot = xCell_DoseResponse_p_Growth,
  width = 21*1/3,
  height = 7,
  units = "cm"
)  

#8. saving graphes for AUC
ggsave(
  filename = "08_xCelligence_doseResponse_AUC.jpeg",
  plot = xcell_DoseResponse_p_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "08_xCelligence_doseResponse_AUC.pdf",
  plot = xcell_DoseResponse_p_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#9. saving tmax
write.csv(
  xCell_DoseResponse_tmax,
  file = file.path(dirs$results, "08_xcelligance_DoseResponse", "xcelligence_DoseResponse_tmax.csv"),
  row.names = FALSE
)
xCell_DoseResponse_stat_data_all

#10. saving linear Growth model
write.csv(xCell_DoseResponse_significance_Growth_rate_lm,
  file = file.path(
    dirs$results,
    "08_xcelligance_DoseResponse",
    "xcelligence_DoseResponse_GrowthRate.csv"
  ),
  row.names = FALSE
)

#11. saving statistics 
write.csv(xCell_DoseResponse_stat_data_all,
          file = file.path(
            dirs$results,
            "08_xcelligance_DoseResponse",
            "xcelligence_DoseResponse_stats.csv"
          ),
          row.names = FALSE
)

#13. saving significance norm
write.csv(
  as.data.frame(xCell_DoseResponse_significance_norm_lm$contrasts),
  file = file.path(
    dirs$results,
    "08_xcelligance_DoseResponse",
    "xcelligence_DoseResponse_AOV_norm.csv"
  ),
  row.names = FALSE
)
#------------------------------------------
setwd(dirs$code)
