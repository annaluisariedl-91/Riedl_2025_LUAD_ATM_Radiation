#-------------------------------------------
# loading or renaming data, preparing colors
#------------------------------------------
load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)

mean_pHCC44_24h <- read.csv(
  file.path(dirs$results, "15_Proteome_DifferentialExpression", "15_Proteome_mean_pHCC44_24h.csv")
)

mean_pH23_24h <- read.csv(
  file.path(dirs$results, "15_Proteome_DifferentialExpression", "15_Proteome_mean_pH23_24h.csv")
)

for (cellline in c("H23", "HCC44")) {
  for (time in c("0h", "2h", "6h", "24h", "48h")) {
    for (treatment in c("c", "2Gy")) {
      name <- paste0(
        "raw_", cellline, "_", time, "_", treatment
      )
      assign(
        name,
        load_prot_raw(cellline, time, treatment)
      )
    }
  }
}

ann_colors <- list(
  Radiation = c(
    "ctrl" = "#4294F7",
    "2Gy"  = "#F09837"
  ),
  Cellline = c(
    "HCC44"   = "#EA4025",
    "NCI-H23" = "#872111"
  )
)

H23_colours <- c(
  "NCI-H23 ctrl" = "#5F8F29",
  "NCI-H23 2 Gy" = "#F09837"
)

HCC44_colours <- c(
  "HCC44 ctrl" ="#4294F7",
  "HCC44 2 Gy" = "#EA4025"
)

PCA_colours <- c(
  "HCC44 ctrl"   = "#4294F7",
  "HCC44 2 Gy"   =  "#EA4025",
  "NCI-H23 ctrl" = "#5F8F29",
  "NCI-H23 2 Gy" = "#F09837"
)

time_colours <- c(
  "0 h"  ="#4294F7",
  "2 h"  = "#5F8F29",
  "6 h"  = "#F09837",
  "24 h" = "#EA4025",
  "48 h" =  "#872111"
)

protein_colours <- c(
  "HCC44 0 h ctrl"      = "#4294F7",
  "HCC44 0 h 2 Gy"      = "#EA4025",
  "HCC44 24 h ctrl"     = "#4294F7",
  "HCC44 24 h 2 Gy"     = "#EA4025",
  "NCI-H23 0 h ctrl"    = "#5F8F29",
  "NCI-H23 0 h 2 Gy"    = "#F09837",
  "NCI-H23 24 h ctrl"   = "#5F8F29",
  "NCI-H23 24 h 2 Gy"   = "#F09837"
)

venn_colours <- c(
  "#4294F7",
  "#EA4025",
  "#5F8F29",
  "#F09837"
)

my_colors <-rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))
#--------------------------------------------
# Proteome - intensity binding data
#--------------------------------------------
# HCC44
Protein_intensity_HCC44 <- make_box_all("HCC44")

Protein_intensity_HCC44_values <- Protein_intensity_HCC44$box
Protein_intensity_HCC44_labels   <- Protein_intensity_HCC44$labels

# NCI-H23
Protein_intensity_result_H23 <- make_box_all("H23")

Protein_intensity_H23_values <- Protein_intensity_result_H23$box
Protein_intensity_H23_labels   <- Protein_intensity_result_H23$labels
#--------------------------------------------
# Intensity Bar Plot
#--------------------------------------------
times <- c("0h", "2h", "6h", "24h", "48h")
treatments <- c("c", "2Gy")
for (cellline in c("H23", "HCC44")) {
  names <- unlist(lapply(times, function(time) {
    paste0("raw_", cellline, "_", time, "_", treatments)
  }))
  combined <- do.call(
    cbind,
    lapply(names, function(name) {
      get(name)[, 5:8]
    })
  )
  colnames(combined) <- unlist(lapply(times, function(time) {
    unlist(lapply(treatments, function(treatment) {
      paste0(time, "_", treatment, "_", 1:4)
    }))
  }))
  assign(paste0("Protein_intensity_box_", cellline, "_all"), combined)
}

Protein_intensity_box_H23_long <-Protein_intensity_box_H23_all %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    Sample_ID = rep(
      seq_len(ncol(Protein_intensity_box_H23_all)),
      each = nrow(Protein_intensity_box_H23_all)
    ),
    Time = case_when(
      Sample_ID <= 8  ~ "0 h",
      Sample_ID <= 16 ~ "2 h",
      Sample_ID <= 24 ~ "6 h",
      Sample_ID <= 32 ~ "24 h",
      TRUE            ~ "48 h"
    ),
    Treatment = if_else(
      ((Sample_ID - 1) %% 8) < 4,
      "ctrl",
      "2 Gy"
    ),
    Group = paste("NCI-H23", Treatment)
  )

Protein_intensity_box_HCC44_long <- Protein_intensity_box_HCC44_all %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    Sample_ID = rep(
      seq_len(ncol(Protein_intensity_box_HCC44_all)),
      each = nrow(Protein_intensity_box_HCC44_all)
    ),
    Time = case_when(
      Sample_ID <= 8  ~ "0 h",
      Sample_ID <= 16 ~ "2 h",
      Sample_ID <= 24 ~ "6 h",
      Sample_ID <= 32 ~ "24 h",
      TRUE            ~ "48 h"
    ),
    Treatment = if_else(
      ((Sample_ID - 1) %% 8) < 4,
      "ctrl",
      "2 Gy"
    ),
    Group = paste("HCC44", Treatment)
  )

Protein_intensity_box_H23_long$Intensity<- log10(2^(Protein_intensity_box_H23_long$Intensity)+1)
Protein_intensity_box_HCC44_long$Intensity<- log10(2^(Protein_intensity_box_HCC44_long$Intensity)+1)

Protein_p_H23_intensity<- ggplot(
  Protein_intensity_box_H23_long,
    aes(
      x = factor(Sample_ID),
      y = Intensity,
      fill = Group
    )
  ) +
  geom_violin(
    width = 0.8,
    linewidth = 0.5,
    colour = "black",
    trim = TRUE,
    na.rm = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.5,
    stroke = 0.5,
    fill = "white",
    colour = "black",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = H23_colours
  ) +
  scale_y_continuous(
    breaks = c(0,5,10,15),
    minor_breaks = c(15,25,35,45),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_discrete(
    breaks = as.character(c(1, 9, 17, 25, 33)),
    labels = c("0 h", "2 h", "6 h", "24 h", "48 h")
  ) +
  coord_cartesian(
    ylim = c(-1, 10),
    clip = "off"
  ) +
  labs(
    title = "Proteome - NCI-H23",
    x = NULL,
    y = expression(log[10] * "(Intensity + 1)")
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "left",
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),
    plot.title = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )  +
  guides(
    fill = guide_legend(
      title = NULL,
      nrow = 1
    )
  )

Protein_p_HCC44_intensity <- ggplot(
  Protein_intensity_box_HCC44_long,
  aes(
    x = factor(Sample_ID),
    y = Intensity,
    fill = Group
  )
) +
  geom_violin(
    width = 0.8,
    linewidth = 0.5,
    colour = "black",
    trim = TRUE,
    na.rm = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.5,
    stroke = 0.5,
    fill = "white",
    colour = "black",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = HCC44_colours
  ) +
  scale_y_continuous(
    breaks = c(0, 5,10),
    minor_breaks = c(15,25,35,45),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_discrete(
    breaks = as.character(c(1, 9, 17, 25, 33)),
    labels = c("0 h", "2 h", "6 h", "24 h", "48 h")
  ) +
  coord_cartesian(
    ylim = c(-1, 10),
    clip = "off"
  ) +
  labs(
    title = "Proteome - HCC44",
    x = NULL,
    y = expression(log[10] * "(Intensity + 1)")
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "left",
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),
    plot.title = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )  +
  guides(
    fill = guide_legend(
      title = NULL,
      nrow = 1
    )
  )  

Protein_p_H23_HCC44_intensity_final <- Protein_p_H23_intensity / Protein_p_HCC44_intensity +
  plot_layout(
    guides = "collect"
  )  &
  theme(
    legend.position = "bottom"
  )
#--------------------------------------------
# Principal Components Analysis
#-------------------------------------------- 
Protein_PCA_data <- data.frame(
  log2(raw_H23_24h_c[, 5:8]),
  log2(raw_H23_24h_2Gy[, 5:8]),
  log2(raw_HCC44_24h_c[, 5:8]),
  log2(raw_HCC44_24h_2Gy[, 5:8])
)

Protein_PCA_data <- Protein_PCA_data %>%
  mutate(
    across(
      everything(),
      ~ replace(.x, !is.finite(.x), NA)
    )
  )

Protein_PCA_data <- Protein_PCA_data[
  complete.cases(Protein_PCA_data),
]

Protein_PCA_data <- t(Protein_PCA_data)

Protein_PCA_res <- prcomp(
  Protein_PCA_data,
  center = TRUE,
  scale. = TRUE
)

Protein_PCA_variance <- Protein_PCA_res$sdev^2 /
  sum(Protein_PCA_res$sdev^2) * 100

Protein_PC1_var <- Protein_PCA_variance[1]
Protein_PC2_var <- Protein_PCA_variance[2]

Protein_PCA_res_plot_data <- as.data.frame(Protein_PCA_res$x[, 1:2])

Protein_PCA_res_plot_data$Group <- factor(
  c(
    rep("NCI-H23 ctrl", 4),
    rep("NCI-H23 2 Gy", 4),
    rep("HCC44 ctrl", 4),
    rep("HCC44 2 Gy", 4)
  ),
  levels = c(
    "HCC44 ctrl",
    "HCC44 2 Gy",
    "NCI-H23 ctrl",
    "NCI-H23 2 Gy"
  )
)

Protein_p_PCA <- ggplot(
  Protein_PCA_res_plot_data,
  aes(
    x = PC1,
    y = PC2,
    fill = Group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey40"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey40"
  ) +
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.7,
    colour = "black"
  ) +
  scale_fill_manual(
    values = PCA_colours,
    breaks = c(
      "HCC44 ctrl",
      "HCC44 2 Gy",
      "NCI-H23 ctrl",
      "NCI-H23 2 Gy"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(-100,-50, 0, 50, 100),
    minor_breaks = c(-50, 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(-100,-50, 0, 50, 100),
    minor_breaks = c(-75,-25, 25,75),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(
    xlim = c(-100, 100),
    ylim = c(-100, 100),
    clip = "off"
  ) +
  labs(
    title = "Proteome, 0 h",
    x = paste0(
      "PCA1 (",
      round(Protein_PC1_var, 1),
      "%)"
    ),
    y = paste0(
      "PCA2 (",
      round(Protein_PC2_var, 1),
      "%)"
    )
  ) +
  theme_prism_xcell() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "left",
    
    legend.text = element_text(size = 8.5),
    
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),

      plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      colour = "black"
    ),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 5,
      l = 10
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = 21,
        size = 3,
        colour = "black"
      )
    )
  )
#--------------------------------------------
# Coefficient of Variation
#--------------------------------------------
Protein_cov_HCC44 <- bind_rows(
  #control
  calculate_cov(
    raw_HCC44_0h_c,
    "HCC44", "ctrl", "0 h"
  ),
  calculate_cov(
    raw_HCC44_2h_c,
    "HCC44", "ctrl", "2 h"
  ),
  calculate_cov(
    raw_HCC44_6h_c,
    "HCC44", "ctrl", "6 h"
  ),
  calculate_cov(
    raw_HCC44_24h_c,
    "HCC44", "ctrl", "24 h"
  ),
  calculate_cov(
    raw_HCC44_48h_c,
    "HCC44", "ctrl", "48 h"
  ),
  # 2 Gy
  calculate_cov(
    raw_HCC44_0h_2Gy,
    "HCC44", "2 Gy", "0 h"
  ),
  calculate_cov(
    raw_HCC44_2h_2Gy,
    "HCC44", "2 Gy", "2 h"
  ),
  calculate_cov(
    raw_HCC44_6h_2Gy,
    "HCC44", "2 Gy", "6 h"
  ),
  calculate_cov(
    raw_HCC44_24h_2Gy,
    "HCC44", "2 Gy", "24 h"
  ),
  calculate_cov(
    raw_HCC44_48h_2Gy,
    "HCC44", "2 Gy", "48 h"
  )
)

Protein_cov_H23 <- bind_rows(
  # control
  calculate_cov(
    raw_H23_0h_c,
    "NCI-H23", "ctrl", "0 h"
  ),
  calculate_cov(
    raw_H23_2h_c,
    "NCI-H23", "ctrl", "2 h"
  ),
  calculate_cov(
    raw_H23_6h_c,
    "NCI-H23", "ctrl", "6 h"
  ),
  calculate_cov(
    raw_H23_24h_c,
    "NCI-H23", "ctrl", "24 h"
  ),
  calculate_cov(
    raw_H23_48h_c,
    "NCI-H23", "ctrl", "48 h"
  ),
  # 2 Gy
  calculate_cov(
    raw_H23_0h_2Gy,
    "NCI-H23", "2 Gy", "0 h"
  ),
  calculate_cov(
    raw_H23_2h_2Gy,
    "NCI-H23", "2 Gy", "2 h"
  ),
  calculate_cov(
    raw_H23_6h_2Gy,
    "NCI-H23", "2 Gy", "6 h"
  ),
  calculate_cov(
    raw_H23_24h_2Gy,
    "NCI-H23", "2 Gy", "24 h"
  ),
  calculate_cov(
    raw_H23_48h_2Gy,
    "NCI-H23", "2 Gy", "48 h"
  )
)

Protein_cov_data <- bind_rows(
  Protein_cov_HCC44,
  Protein_cov_H23
) %>%
  mutate(
    Time = factor(
      Time,
      levels = c(
        "0 h",
        "2 h",
        "6 h",
        "24 h",
        "48 h"
      )
    ),
    Treatment = factor(
      Treatment,
      levels = c("ctrl", "2 Gy")
    )
  )

Protein_cov_density <- Protein_cov_data %>%
  group_by(Cellline, Treatment, Time) %>%
  group_modify(~ make_density(.x)) %>%
  ungroup()

Protein_p_cov_HCC44_ctrl <- plot_cov(
  Protein_cov_density,
  cellline = "HCC44",
  treatment = "ctrl",
  title = "HCC44 ctrl"
)

Protein_p_cov_HCC44_2Gy <- plot_cov(
  Protein_cov_density,
  cellline = "HCC44",
  treatment = "2 Gy",
  title = "HCC44 2 Gy"
)

Protein_p_cov_H23_ctrl <- plot_cov(
  Protein_cov_density,
  cellline = "NCI-H23",
  treatment = "ctrl",
  title = "NCI-H23 ctrl"
)

Protein_p_cov_H23_2Gy <- plot_cov(
  Protein_cov_density,
  cellline = "NCI-H23",
  treatment = "2 Gy",
  title = "NCI-H23 2 Gy"
)

Protein_p_cov_final <-
  (
    Protein_p_cov_HCC44_ctrl |
      Protein_p_cov_HCC44_2Gy
  ) /
  (
    Protein_p_cov_H23_ctrl |
      Protein_p_cov_H23_2Gy
  ) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    axis.text = element_text(
      size = 10,
      colour = "black"
    )
  )

Protein_p_cov_HCC44_ctrl <- Protein_p_cov_HCC44_ctrl +
  theme(
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 6)
  )
Protein_p_cov_HCC44_2Gy <- Protein_p_cov_HCC44_2Gy +
  theme(
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 6)
  )

Protein_p_cov_H23_ctrl <- Protein_p_cov_H23_ctrl +
  theme(
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 6)
  )

Protein_p_cov_H23_2Gy <- Protein_p_cov_H23_2Gy +
  theme(
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 6)
  )
#--------------------------------------------
# Bar plot identified Proteins
#--------------------------------------------
Protein_identified_proteins <- data.frame(
  Group = c(
    "HCC44 0 h ctrl",
    "HCC44 0 h 2 Gy",
    "HCC44 24 h ctrl",
    "HCC44 24 h 2 Gy",
    "NCI-H23 0 h ctrl",
    "NCI-H23 0 h 2 Gy",
    "NCI-H23 24 h ctrl",
    "NCI-H23 24 h 2 Gy"
  ),
  n = c(
    sum(!is.na(raw_HCC44_0h_c$mean)),
    sum(!is.na(raw_HCC44_0h_2Gy$mean)),
    sum(!is.na(raw_HCC44_24h_c$mean)),
    sum(!is.na(raw_HCC44_24h_2Gy$mean)),
    sum(!is.na(raw_H23_0h_c$mean)),
    sum(!is.na(raw_H23_0h_2Gy$mean)),
    sum(!is.na(raw_H23_24h_c$mean)),
    sum(!is.na(raw_H23_24h_2Gy$mean))
  )
)

Protein_identified_proteins$Group <- factor(
  Protein_identified_proteins$Group,
  levels = Protein_identified_proteins$Group
)

Protein_p_identified_proteins <- ggplot(
  Protein_identified_proteins,
  aes(
    x = Group,
    y = n,
    fill = Group
  )
) +
  geom_col(
    width = 0.7,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = n),
    vjust = -0.6,
    hjust = -0.05,
    angle = 35,
    size = 3.5
  ) +
  geom_hline(
    yintercept = 7557,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "black"
  ) +
  scale_fill_manual(
    values = protein_colours
  ) +
  scale_y_continuous(
    limits = c(0, 12000),
    breaks = c(0, 5000, 10000),
    labels = c("0", "5.000", "10.000"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(
      add = c(0.6, 2.95)
    )
  )+
  labs(
    x = NULL,
    y = "Quantified proteins (n)"
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    axis.ticks.length =
      grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 10,
      angle = 30,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "none",
    plot.margin = margin(
      t = 15,
      r = 10,
      b = 5,
      l = 10
    )
  )

Protein_p_identified_proteins <- Protein_p_identified_proteins +
 annotate(
    "text",
    x = 8.5,
    y = 8000,
    label = "Shared proteins:\n8870",
    hjust = 0,
    vjust = 1.6,
    size = 3
  )
#--------------------------------------------
# Heat map
#--------------------------------------------
Protein_hm_H23_HCC44_all <- data.frame(
  log2(raw_H23_24h_c[, 5:8]),
  log2(raw_H23_24h_2Gy[, 5:8]),
  log2(raw_HCC44_24h_c[, 5:8]),
  log2(raw_HCC44_24h_2Gy[, 5:8])
)

Protein_hm_H23_HCC44_all <- Protein_hm_H23_HCC44_all %>%
  mutate(across(everything(), ~ replace(.x, !is.finite(.x), NA)))

Protein_hm_H23_HCC44_all <- cbind(
  raw_H23_24h_c[, 1:3],
  Protein_hm_H23_HCC44_all
)

Protein_hm_H23_HCC44_all <-
  Protein_hm_H23_HCC44_all[
    complete.cases(Protein_hm_H23_HCC44_all[, 4:19]),
  ]

colnames(Protein_hm_H23_HCC44_all) <- c(
  "ENTREZ",
  "SYMBOL",
  "ID",

  "NCI-H23 c 24h.1",
  "NCI-H23 c 24h.2",
  "NCI-H23 c 24h.3",
  "NCI-H23 c 24h.4",
  
  "NCI-H23 2Gy 24h.1",
  "NCI-H23 2Gy 24h.2",
  "NCI-H23 2Gy 24h.3",
  "NCI-H23 2Gy 24h.4",
  
  "HCC44 c 24h.1",
  "HCC44 c 24h.2",
  "HCC44 c 24h.3",
  "HCC44 c 24h.4",
  
  "HCC44 2Gy 24h.1",
  "HCC44 2Gy 24h.2",
  "HCC44 2Gy 24h.3",
  "HCC44 2Gy 24h.4"
)

Protein_hm_matrix <- as.matrix(
  Protein_hm_H23_HCC44_all[, 4:19]
)

Protein_hm_matrix_z <- t(
  scale(
    t(Protein_hm_matrix)
  )
)

annotation_col <- data.frame(
  Radiation = factor(
    rep(c("ctrl", "2Gy"), each = 4, times = 2),
    levels = c("ctrl", "2Gy")
  ),
  `Cell line` = factor(
    rep(c("NCI-H23", "HCC44"), each = 8),
    levels = c("NCI-H23", "HCC44")
  ),
  check.names = FALSE
)

rownames(annotation_col) <- colnames(Protein_hm_matrix)

identical(
  colnames(Protein_hm_matrix_z),
  rownames(annotation_col)
)

Protein_p_hm_intensity_all <- pheatmap(
  Protein_hm_matrix_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = my_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  fontsize = 12
)
#--------------------------------------------
# Proteome Venn Diagramm
#--------------------------------------------
Protein_venn_proteome <- venn.diagram(
  x = list(
    "HCC44 ctrl"   = mean_pHCC44_24h[complete.cases(mean_pHCC44_24h[,5]),1],
    "HCC44 2 Gy"   = mean_pHCC44_24h[complete.cases(mean_pHCC44_24h[,7]),1],
    "NCI-H23 ctrl" = mean_pH23_24h[complete.cases(mean_pH23_24h[,5]),1],
    "NCI-H23 2 Gy" = mean_pH23_24h[complete.cases(mean_pH23_24h[,7]),1]
  ),
  filename = NULL,
  fill = venn_colours,
  alpha = rep(0.65, 4),
  col = rep("black", 4),
  lwd = 1.2,
  lty = "solid",
  print.mode = "raw",
  cex = 0.8,
  fontfamily = "sans",
  fontface = "plain",
  cat.cex = 0.9,
  cat.fontfamily = "sans",
  cat.fontface = "plain",
  cat.col = "black",
  cat.pos = c(335, 25, 340, 20),
  cat.dist = c(0.25, 0.25, 0.15, 0.15),
  scaled = FALSE,
  main = "",
  main.cex = 1,
  main.fontfamily = "sans",
  main.fontface = "bold",
  main.just = c(0,0.2),
  margin = 0.25
)


#--------------------------------------------
# saving results
#--------------------------------------------
dir.create(
  file.path(dirs$results, "19_Proteome_QualityControl"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "19_Proteome_QualityControl"))

#1. Saving Intensity violine plot
ggsave(
  "19_Proteome_Intensity.pdf",
  Protein_p_H23_HCC44_intensity_final,
  width = 28*1/2,
  height = 14,
  units = "cm"
)

ggsave(
  "19_Proteome_Intensity.jpg",
  Protein_p_H23_HCC44_intensity_final,
  width = 28*1/2,
  height = 14,
  units = "cm",
  dpi = 600
)

#1.a HCC44
ggsave(
  "19_Proteome_Intensity_HCC44.pdf",
  Protein_p_HCC44_intensity,
  width = 28*1/2,
  height = 7,
  units = "cm"
)

ggsave(
  "19_Proteome_Intensity_HCC44.jpg",
  Protein_p_HCC44_intensity,
  width = 28*1/2,
  height = 7,
  units = "cm",
  dpi = 600
)

#1.b NCI-H23
ggsave(
  "19_Proteome_Intensity_H23.pdf",
  Protein_p_H23_intensity,
  width = 28*1/2,
  height = 7,
  units = "cm"
)

ggsave(
  "19_Proteome_Intensity_H23.jpg",
  Protein_p_H23_intensity,
  width = 28*1/2,
  height = 7,
  units = "cm",
  dpi = 600
)

#2. Pricipal component Analysis
ggsave(
  filename = "19_Proteome_PCA_24h.pdf",
  plot = Protein_p_PCA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "19_Proteome_PCA_24h.jpg",
  plot = Protein_p_PCA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
  
#3. Coefficient of Variation
ggsave(
  filename = "19_Proteome_COV_all.pdf",
  plot = Protein_p_cov_final,
  width = 21*2/3,
  height = 14,
  units = "cm"
)

ggsave(
  filename = "19_Proteome_COV_all.jpg",
  plot = Protein_p_cov_final,
  width = 21*2/3,
  height = 14,
  units = "cm"
)

#3.a HCC44 control
ggsave(
  filename = "19_Proteome_COV_HCC44_ctrl.pdf",
  plot = Protein_p_cov_HCC44_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "19_Proteome_COV_HCC44_ctrl.jpg",
  plot = Protein_p_cov_HCC44_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3.b HCC44 2Gy
ggsave(
  filename = "19_Proteome_COV_HCC44_2Gy.pdf",
  plot = Protein_p_cov_HCC44_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "19_Proteome_COV_HCC44_2Gy.jpg",
  plot = Protein_p_cov_HCC44_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3.c NCI-H23 control
ggsave(
  filename = "19_Proteome_COV_H23_ctrl.pdf",
  plot = Protein_p_cov_H23_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "19_Proteome_COV_H23_ctrl.jpg",
  plot = Protein_p_cov_H23_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3.d NCI-H23 2 Gy
ggsave(
  filename = "19_Proteome_COV_H23_2Gy.pdf",
  plot = Protein_p_cov_H23_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "19_Proteome_COV_H23_2Gy.jpg",
  plot = Protein_p_cov_H23_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#5. Identified Proteins
ggsave(
  "19_Proteome_identifiedProteins.pdf",
  Protein_p_identified_proteins,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

ggsave(
  "19_Proteome_identifiedProteins.jpg",
  Protein_p_identified_proteins,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

#6. Heatmaps
pdf(
  "19_Proteome_Heatmap_all_unfiltered.pdf",
  width = 21*1/3,
  height = 3.5
)

grid::grid.newpage()
grid::grid.draw(Protein_p_hm_intensity_all$gtable)
dev.off()

#7. Venn Diagramm
pdf(
  "19_Proteome_VennDiagram",
  width = 9 / 2.54,
  height = 7 / 2.54,
  useDingbats = FALSE
)
grid.newpage()
grid.draw(Protein_venn_proteome)
grid.text(
  "Proteome",
  x = 0.5,
  y = 0.75,
  gp = gpar(
    fontsize = 14,
    fontface = "bold",
    fontfamily = "sans"
  )
)
dev.off()
#--------------------------------------------
setwd(dirs$code)
