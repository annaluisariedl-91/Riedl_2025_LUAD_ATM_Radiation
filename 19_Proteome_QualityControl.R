#-------------------------------------------
# loading or renaming data, preparing colors
#------------------------------------------
load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
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
  "ctrl" = "#5F8F29",
  "2 Gy" = "#F09837"
)

HCC44_colours <- c(
  "ctrl" ="#4294F7",
  "2 Gy" = "#EA4025"
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

my_colors <-rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))
#--------------------------------------------
# Proteome - intensity binding data
#--------------------------------------------
# HCC44
result_HCC44 <- make_box_all("HCC44")

box_HCC44_values <- result_HCC44$box
dat_HCC44_labels   <- result_HCC44$labels

# NCI-H23
result_H23 <- make_box_all("H23")

box_H23_values <- result_H23$box
dat_H23_labels   <- result_H23$labels
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
  assign(paste0("Proteome_box_", cellline, "_all"), combined)
}

Proteome_box_H23_long <- Proteome_box_H23_all %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    Sample_ID = rep(
      seq_len(ncol(Proteome_box_H23_all)),
      each = nrow(Proteome_box_H23_all)
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
    Group = paste(Time, Treatment)
  )

Proteome_box_HCC44_long <- Proteome_box_HCC44_all %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    Sample_ID = rep(
      seq_len(ncol(Proteome_box_HCC44_all)),
      each = nrow(Proteome_box_HCC44_all)
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
    Group = paste(Time, Treatment)
  )

Proteome_box_H23_long$Intensity<- log2(exp(Proteome_box_H23_long$Intensity)+1)
Proteome_box_HCC44_long$Intensity<- log2(exp(Proteome_box_HCC44_long$Intensity)+1)

p_H23_intensity_proteome <- ggplot(
  Proteome_box_HCC44_long,
    aes(
      x = factor(Sample_ID),
      y = Intensity,
      fill = Treatment
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
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "black"
  ) +
  scale_fill_manual(
    values = H23_colours
  ) +
  scale_y_continuous(
    breaks = c(0, 10,20,30,40,50),
    minor_breaks = c(15,25,35,45),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_x_discrete(
    breaks = as.character(c(1, 9, 17, 25, 33)),
    labels = c("0 h", "2 h", "6 h", "24 h", "48 h")
  ) +
  coord_cartesian(
    ylim = c(-2, 50),
    clip = "off"
  ) +
  labs(
    title = "Proteome - NCI-H23",
    x = NULL,
    y = expression(log[10] * "(Intensity + 1)")
  ) +
  theme_prism_xcell() +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "left",
    
    legend.text = element_text(size = 6.5),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    
    plot.margin = margin(10, 10, 5, 10),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    )
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      nrow = 1
    )
  )

p_HCC44_intensity_proteome <- ggplot(
  Proteome_box_HCC44_long,
  aes(
    x = factor(Sample_ID),
    y = Intensity,
    fill = Treatment
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
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "black"
  ) +
  scale_fill_manual(
    values = HCC44_colours
  ) +
  scale_y_continuous(
    breaks = c(0, 10,20,30,40,50),
    minor_breaks = c(15,25,35,45),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_x_discrete(
    breaks = as.character(c(1, 9, 17, 25, 33)),
    labels = c("0 h", "2 h", "6 h", "24 h", "48 h")
  ) +
  coord_cartesian(
    ylim = c(-2, 50),
    clip = "off"
  ) +
  labs(
    title = "Proteome - HCC44",
    x = NULL,
    y = expression(log[10] * "(Intensity + 1)")
  ) +
  theme_prism_xcell() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "left",
    
    legend.text = element_text(size = 6.5),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    
    plot.margin = margin(10, 10, 5, 10),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    )
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      nrow = 1
    )
  )  

p_proteome_intensity_qc <- p_H23_intensity_proteome / p_HCC44_intensity_proteome +
  plot_layout(
    ncol = 1,
    guides = "collect"
  )  
#--------------------------------------------
# Principal Components Analysis
#-------------------------------------------- 
pbase_pca_filt <- data.frame(
  log2(raw_H23_24h_c[, 5:8]),
  log2(raw_H23_24h_2Gy[, 5:8]),
  log2(raw_HCC44_24h_c[, 5:8]),
  log2(raw_HCC44_24h_2Gy[, 5:8])
)

pbase_pca_filt <- pbase_pca_filt %>%
  mutate(
    across(
      everything(),
      ~ replace(.x, !is.finite(.x), NA)
    )
  )

pbase_pca_filt <- pbase_pca_filt[
  complete.cases(pbase_pca_filt),
]

pbase_pca_filt <- t(pbase_pca_filt)

pc_pbase_pca_filt <- prcomp(
  pbase_pca_filt,
  center = TRUE,
  scale. = TRUE
)

proteome_pca_variance <- pc_pbase_pca_filt$sdev^2 /
  sum(pc_pbase_pca_filt$sdev^2) * 100

proteome_PC1_var <- proteome_pca_variance[1]
proteome_PC2_var <- proteome_pca_variance[2]

proteome_pca_data <- as.data.frame(pc_pbase_pca_filt$x[, 1:2])

proteome_pca_data$Group <- factor(
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

p_proteome_PCA <- ggplot(
  proteome_pca_data,
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
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = c(-100, 100),
    ylim = c(-100, 100),
    clip = "off"
  ) +
  labs(
    title = "Proteome, 24 h",
    x = paste0(
      "PCA1 (",
      round(proteome_PC1_var, 1),
      "%)"
    ),
    y = paste0(
      "PCA2 (",
      round(proteome_PC2_var, 1),
      "%)"
    )
  ) +
  theme_prism_xcell() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "left",
    
    legend.text = element_text(size = 6.5),
    
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),

      plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text = element_text(
      size = 9,
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
cov_proteome_HCC44 <- bind_rows(
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

cov_proteome_H23 <- bind_rows(
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

cov_proteome_data <- bind_rows(
  cov_proteome_HCC44,
  cov_proteome_H23
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

cov_proteome_density <- cov_proteome_data %>%
  group_by(
    Cellline,
    Treatment,
    Time
  ) %>%
  group_modify(
    ~ make_density(.x)
  ) %>%
  ungroup() %>%
  mutate(
    log_CV = log10(CV + 1)
  )

p_cov_HCC44_ctrl <- plot_cov(
  cov_proteome_density,
  cellline = "HCC44",
  treatment = "ctrl",
  title = "HCC44 ctrl"
)

p_cov_HCC44_2Gy <- plot_cov(
  cov_proteome_density,
  cellline = "HCC44",
  treatment = "2 Gy",
  title = "HCC44 2 Gy"
)

p_cov_H23_ctrl <- plot_cov(
  cov_proteome_density,
  cellline = "NCI-H23",
  treatment = "ctrl",
  title = "NCI-H23 ctrl"
)

p_cov_H23_2Gy <- plot_cov(
  cov_proteome_density,
  cellline = "NCI-H23",
  treatment = "2 Gy",
  title = "NCI-H23 2 Gy"
)

p_proteome_cov_final <-
  (
    p_cov_HCC44_ctrl |
    p_cov_HCC44_2Gy
  ) /
  (
    p_cov_H23_ctrl |
    p_cov_H23_2Gy
  ) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )
#--------------------------------------------
# Bar plot identified Proteins
#--------------------------------------------
identified_proteins <- data.frame(
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

identified_proteins$Group <- factor(
  identified_proteins$Group,
  levels = identified_proteins$Group
)

p_identified_proteins <- ggplot(
  identified_proteins,
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
    vjust = -0.4,
    size = 3
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
    limits = c(0, 11000),
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
      size = 8,
      angle = 30,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "none",
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

p_identified_proteins <- p_identified_proteins +
 annotate(
    "text",
    x = 8.5,
    y = 8400,
    label = "Shared proteins:\n8870",
    hjust = 0,
    size = 3
  )
#--------------------------------------------
# Heat map
#--------------------------------------------
pH23_pHCC44_24h_all_filt <- data.frame(
  log2(raw_H23_24h_c[, 5:8]),
  log2(raw_H23_24h_2Gy[, 5:8]),
  log2(raw_HCC44_24h_c[, 5:8]),
  log2(raw_HCC44_24h_2Gy[, 5:8])
)

pH23_pHCC44_24h_all_filt <- pH23_pHCC44_24h_all_filt %>%
  mutate(across(everything(), ~ replace(.x, !is.finite(.x), NA)))

pH23_pHCC44_24h_all_filt <- cbind(
  raw_H23_24h_c[, 1:3],
  pH23_pHCC44_24h_all_filt
)

pH23_pHCC44_24h_all_filt <-
  pH23_pHCC44_24h_all_filt[
    complete.cases(pH23_pHCC44_24h_all_filt[, 4:19]),
  ]

colnames(pH23_pHCC44_24h_all_filt) <- c(
  "ENTREZ",
  "SYMBOL",
  "ID",

  "NCI-H23 c 0h.1",
  "NCI-H23 c 0h.2",
  "NCI-H23 c 0h.3",
  "NCI-H23 c 0h.4",
  
  "NCI-H23 2Gy 0h.1",
  "NCI-H23 2Gy 0h.2",
  "NCI-H23 2Gy 0h.3",
  "NCI-H23 2Gy 0h.4",
  
  "HCC44 c 0h.1",
  "HCC44 c 0h.2",
  "HCC44 c 0h.3",
  "HCC44 c 0h.4",
  
  "HCC44 2Gy 0h.1",
  "HCC44 2Gy 0h.2",
  "HCC44 2Gy 0h.3",
  "HCC44 2Gy 0h.4"
)

proteome_matrix <- as.matrix(
  pH23_pHCC44_24h_all_filt[, 4:19]
)

proteome_z <- t(
  scale(
    t(proteome_matrix)
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

rownames(annotation_col) <- colnames(proteome_matrix)

identical(
  colnames(proteome_z),
  rownames(annotation_col)
)

Graph_QC_Proteom_2 <- pheatmap(
  proteome_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = my_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  border_color = NA,
  fontsize = 9
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
  "Proteome_qc_Intensity.pdf",
  p_proteome_intensity_qc,
  width = 28*1/2,
  height = 14,
  units = "cm"
)

ggsave(
  "Proteome_qc_Intensity.jpg",
  p_proteome_intensity_qc,
  width = 28*1/2,
  height = 14,
  units = "cm",
  dpi = 600
)

#1.a HCC44
ggsave(
  "Proteome_qc_Intensity_HCC44.pdf",
  p_HCC44_intensity_proteome,
  width = 28*1/2,
  height = 14,
  units = "cm"
)

ggsave(
  "Proteome_qc_Intensity_HCC44.jpg",
  p_HCC44_intensity_proteome,
  width = 28*1/2,
  height = 14,
  units = "cm",
  dpi = 600
)

#1.b NCI-H23
ggsave(
  "Proteome_qc_Intensity_H23.pdf",
  p_H23_intensity_proteome,
  width = 28*1/2,
  height = 14,
  units = "cm"
)

ggsave(
  "Proteome_qc_Intensity_H23.jpg",
  p_H23_intensity_proteome,
  width = 28*1/2,
  height = 14,
  units = "cm",
  dpi = 600
)

#2. Pricipal component Analysis
ggsave(
  filename = "Proteome_qc_PCA.pdf",
  plot = p_proteome_PCA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_qc_PCA.jpg",
  plot = p_proteome_PCA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3. Coefficient of Variation
ggsave(
  filename = "Proteome_qc_COV_all.pdf",
  plot = p_proteome_cov_final,
  width = 21*2/3,
  height = 14,
  units = "cm"
)

ggsave(
  filename = "Proteome_qc_COV_all.jpg",
  plot = p_proteome_cov_final,
  width = 21*2/3,
  height = 14,
  units = "cm"
)

#3.a HCC44 control
ggsave(
  filename = "Proteome_qc_COV_HCC44_ctrl.pdf",
  plot = p_cov_HCC44_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_qc_COV_HCC44_ctrl.jpg",
  plot = p_cov_HCC44_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3.b HCC44 2Gy
ggsave(
  filename = "Proteome_qc_COV_HCC44_2Gy.pdf",
  plot = p_cov_HCC44_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_qc_COV_HCC44_2Gy.jpg",
  plot = p_cov_HCC44_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3.c NCI-H23 control
ggsave(
  filename = "Proteome_qc_COV_H23_ctrl.pdf",
  plot = p_cov_H23_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_qc_COV_H23_ctrl.jpg",
  plot = p_cov_H23_ctrl,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3.d NCI-H23 2 Gy
ggsave(
  filename = "Proteome_qc_COV_H23_2Gy.pdf",
  plot = p_cov_H23_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_qc_COV_H23_2Gy.jpg",
  plot = p_cov_H23_2Gy,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#5. Identified Proteins
ggsave(
  "Proteome_qc_identifiedProteins.pdf",
  p_identified_proteins,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

ggsave(
  "Proteome_qc_identifiedProteins.jpg",
  p_identified_proteins,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 600
)

#6. Heatmaps
pdf(
  "Proteome_Heatmap.pdf",
  width = 21*1/3,
  height = 3.5
)

grid::grid.newpage()
grid::grid.draw(Graph_QC_Proteom_2$gtable)

dev.off()
#-----------------------------------------------
setwd(dirs$code)
