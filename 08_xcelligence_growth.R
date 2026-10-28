#------------------------------------------
# Setting path
#------------------------------------------
setwd(file.path(base_dir,"functional_assays","xCelligence","Exp1_Dose response","processed")) #replace with path on your own PC
#original xlxs format was transformed to .csv for easier importing, original is available at /raw
Doses <- c("2 Gy", "1.5 Gy", "1 Gy", "0.5 Gy", "0 Gy")
#------------------------------------------
# Binding and preparing data
#------------------------------------------
xCell_growth_Data_H23 <- prepare_xCell(
  "extracted_H23_rep1.csv"
)

xCell_growth_Data_HCC44 <- prepare_xCell(
  "extracted_HCC44_rep1.csv"
)

stat_growth_data_all <- bind_rows(
  xCell_growth_Data_H23   %>% mutate(Cellline = "NCI-H23"),
  xCell_growth_Data_HCC44 %>% mutate(Cellline = "HCC44")
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
max_ci_growth <- stat_growth_data_all %>%
  group_by(Cellline, Dose) %>%
  summarise(
    maxCI = max(CellIndex, na.rm = TRUE),
    .groups = "drop"
  )


# time point od maximal Cell Index
tmax_growth <- stat_growth_data_all %>%
  group_by(Cellline, Dose) %>%
  slice_max(
    CellIndex,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()
#------------------------------------------
# growth rate
#------------------------------------------
xCell_growth_data <-stat_growth_data_all %>%
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

xCell_growth_rate_lm <- lm(
  k ~ Dose * Cellline,
  data = xCell_growth_data
)
#------------------------------------------
# normalization
#------------------------------------------
# Normalisation to 0 Gy
norm_xCell_growth_data <- stat_growth_data_all %>%
  group_by(Cellline, Time) %>%
  mutate(
    CI_norm = CellIndex / CellIndex[Dose == "0 Gy"][1]
  ) %>%
  ungroup()

# split cell line and dose
split_xCell_growth_data <- norm_xCell_growth_data %>%
  group_by(Cellline, Dose) %>%
  group_split() %>%
  setNames(
    norm_xCell_growth_data %>%
      distinct(Cellline, Dose) %>%
      transmute(
        name = paste(Cellline, Dose, sep = "_")
      ) %>%
      pull(name)
  )

stat_xCell_growth_lin_model <- lm(
  CellIndex ~ Time * Dose * Cellline,
  data = stat_growth_data_all
)

xCell_growth_aov_lin_model<-anova(stat_xCell_growth_lin_model)
summary(stat_xCell_growth_lin_model)
#------------------------------------------
# Area under the curve (AUC)
#-----------------------------------------
auc_xCell_growth_data_all <- norm_xCell_growth_data %>%
  group_by(Cellline, Dose) %>%
  summarise(
    AUC = pracma::trapz(Time, CI_norm),
    .groups = "drop"
  )


# Model for AUC
mod_auc_xCell_growth_all <- lm(
  AUC ~ Dose * Cellline,
  data = auc_xCell_growth_data_all
)

xCell_growth_auc_lin_model<-anova(mod_auc_xCell_growth_all)
summary(mod_auc_xCell_growth_all)

# AUC relative to 0 Gy control
auc_xCell_growth_data_all <- auc_xCell_growth_data_all %>%
  group_by(Cellline) %>%
  mutate(
    AUC_rel = AUC / AUC[Dose == "0 Gy"][1]
  ) %>%
  ungroup()

auc_xCell_growth_data_all <- auc_xCell_growth_data_all%>%
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

growth_colours <- c(
  "HCC44" = "#3288BD",
  "NCI-H23"   = "#ABDDA4"
)
#-----------------------------------------  
# plotting raw data
#-----------------------------------------
plot_xcell_growth_data_H23 <- stat_growth_data_all %>%
  filter(Cellline == "NCI-H23") %>%
  mutate(
    Group = Dose
  )

p_xcell_growth_raw_H23 <- plot_xcell(
  data = plot_xcell_growth_data_H23,
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

plot_xcell_growth_data_HCC44 <- stat_growth_data_all %>%
  filter(Cellline == "HCC44") %>%
  mutate(
    Group = Dose
  )

p_xcell_growth_raw_HCC44 <- plot_xcell(
  data = plot_xcell_growth_data_HCC44,
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

p_xcell_growth_raw_H23 <- p_xcell_growth_raw_H23 +
  theme(legend.position = "right")

p_xcell_growth_raw_HCC44 <- p_xcell_growth_raw_HCC44 +
  theme(legend.position = "right")

xCell_growth_raw_final <- (
  p_xcell_growth_raw_H23 +
    p_xcell_growth_raw_HCC44 +
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
plot_xcell_growth_norm_data_H23 <- norm_xCell_growth_data %>%
  filter(Cellline == "NCI-H23") %>%
  mutate(
    Group = Dose
  )
plot_xcell_growth_norm_data_H23<-plot_xcell_growth_norm_data_H23[,-6]
colnames(plot_xcell_growth_norm_data_H23)[6]<-"CellIndex"

p_xcell_growth_norm_H23 <- plot_xcell(
  data = plot_xcell_growth_norm_data_H23,
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

plot_norm_data_HCC44 <- norm_xCell_growth_data %>%
  filter(Cellline == "HCC44") %>%
  mutate(
    Group = Dose
  )

plot_norm_data_HCC44<-plot_norm_data_HCC44[,-6]
colnames(plot_norm_data_HCC44)[6]<-"CellIndex"

p_xcell_growth_norm_HCC44 <- plot_xcell(
  data = plot_norm_data_HCC44,
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

p_xcell_growth_norm_H23 <- p_xcell_growth_norm_H23 +
  theme(legend.position = "right")

p_xcell_growth_norm_HCC44 <-p_xcell_growth_norm_HCC44 +
  theme(legend.position = "right")

xCell_growth_norm_final <- (
  p_xcell_growth_norm_H23 +
    p_xcell_growth_norm_HCC44 +
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
p_xcell_grwoth_GrowthRate <-
  ggplot(
    xCell_growth_data,
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
    values = growth_colours,
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
# plotting AUC
#-----------------------------------------
p_xcell_growth_AUC <-
  ggplot(
    auc_xCell_growth_data_all,
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
    values = growth_colours,
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
  file.path(dirs$results, "08_xcelligance_growth"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "08_xcelligance_growth"))

#1. saving combined raw graph
ggsave(
  filename = "xCelligence_doseResponse_all.jpeg",
  plot = xCell_growth_raw_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_all.pdf",
  plot = xCell_growth_raw_final,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. saving NCI-H23 raw graph
ggsave(
  filename = "xCelligence_doseResponse_H23.jpeg",
  plot = p_xcell_growth_raw_H23,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_H23.pdf",
  plot = p_xcell_growth_raw_H23,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#3. saving HCC44 raw graph
ggsave(
  filename = "xCelligence_doseResponse_HCC44.jpeg",
  plot = p_xcell_growth_raw_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_HCC44.pdf",
  plot = p_xcell_growth_raw_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#4. saving combined normalised graphes 
ggsave(
  filename = "xCelligence_doseResponse_norm.jpeg",
  plot = xCell_growth_norm_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_norm.pdf",
  plot = xCell_growth_norm_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#5. saving combined normalised H23 graphes 
ggsave(
  filename = "xCelligence_doseResponse_norm_H23.jpeg",
  plot = p_xcell_growth_norm_H23,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_norm_H23.pdf",
  plot = p_xcell_growth_norm_H23,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#6. saving HCC44 combined normalied graphes 
ggsave(
  filename = "xCelligence_doseResponse_norm_HCC44.jpeg",
  plot = p_xcell_growth_norm_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_norm_HCC44.pdf",
  plot = p_xcell_growth_norm_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#7. saving graph for data
ggsave(
  filename = "xCelligence_doseResponse_GrowthRate.jpeg",
  plot = p_xcell_grwoth_GrowthRate,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_GrowthRate.pdf",
  plot = p_xcell_grwoth_GrowthRate,
  width = 21*1/3,
  height = 7,
  units = "cm"
)  

#8. saving graphes for AUC
ggsave(
  filename = "xCelligence_doseResponse_AUC.jpeg",
  plot = p_xcell_growth_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AUC.pdf",
  plot = p_xcell_growth_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 


#9. saving graphes for AUC
ggsave(
  filename = "xCelligence_doseResponse_AUC.jpeg",
  plot = p_xcell_growth_AUC,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AUC.pdf",
  plot = p_xcell_growth_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#10. saving graphes for AUC
ggsave(
  filename = "xCelligence_doseResponse_AUC.jpeg",
  plot = p_xcell_growth_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
)

ggsave(
  filename = "xCelligence_doseResponse_AUC.pdf",
  plot = p_xcell_growth_AUC,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 600
) 

#11. saving tmax
write.csv(
  tmax_growth,
  file = file.path(dirs$results, "08_xcelligance_growth", "xcelligence_growth_tmax.csv"),
  row.names = FALSE
)

#12. saving linear growth model
write.csv(
  broom::tidy(xCell_growth_rate_lm),
  file = file.path(
    dirs$results,
    "08_xcelligance_growth",
    "xcelligence_growth_linearModel.csv"
  ),
  row.names = FALSE
)

#13. saving linear growth time model
write.csv(
  broom::tidy(stat_xCell_growth_lin_model),
  file = file.path(
    dirs$results,
    "08_xcelligance_growth",
    "xcelligence_growth_linearModel_time.csv"
  ),
  row.names = FALSE
)

#14. saving AOV
write.csv(
  broom::tidy(xCell_growth_aov_lin_model),
  file = file.path(
    dirs$results,
    "08_xcelligance_growth",
    "xcelligence_growth_AOV_time.csv"
  ),
  row.names = FALSE
)

#15. saving AUC
write.csv(
  broom::tidy(xCell_growth_auc_lin_model),
  file = file.path(
    dirs$results,
    "08_xcelligance_growth",
    "xcelligence_growth_AUC_time.csv"
  ),
  row.names = FALSE
)
#------------------------------------------
setwd(dirs$code)
