#-------------------------------------------
# loading or renaming data
#------------------------------------------
load_stat_data_R(
  file.path(dirs$results, "14_Transcriptome_DifferentialExpression"),
  mode = "ready"
)

RNA_RH23_24h   <- load_RNA("H23", "24h")
RNA_RHCC44_24h <- load_RNA("HCC44", "24h")

setwd(file.path(dirs$RNA))
RNA_RPKM<- data.frame(read.csv("RPKM_Bohnenberger.csv"))
RNA_TPM<- data.frame(read.csv("TPM_Bohnenberger.csv"))
#-------------------------------------------
# defining colors
#------------------------------------------
RPKM_colours <- c("#872111","#EA4025","#F09837","#5F8F29","#3F8D27","#4294F7")

RNA_colours <- c(
  "HCC44 ctrl"   = "#4294F7",
  "HCC44 2 Gy"    = "#EA4025",
  "NCI-H23 ctrl" = "#5F8F29",
  "NCI-H23 2 Gy"  ="#F09837")

venn_colours <- c(
  "#4294F7",
  "#EA4025",
  "#5F8F29",
  "#F09837"
)

ann_colors <- list(
  Radiation = c(
    "ctrl" = "#4294F7",
    "2Gy"  = "#F09837"
  ),
  "Cell line" = c(
    "HCC44"   = "#EA4025",
    "NCI-H23" = "#872111"
  )
)

heatmap_colors <- rev(
  colorRampPalette(
    brewer.pal(11, "RdYlBu")
  )(100)
)
#--------------------------------------------
# Violin plot of Intensities
#--------------------------------------------
RNA_intensity_HCC44_H23_24h <- cbind(
  RNA_RH23_24h[,c(1:2,5:11)],
  RNA_RHCC44_24h[, 5:12]
) |>
  setNames(c(
    "Genename", "ID",   "
    H23_ctrl_1",
    "H23_ctrl_2",
    "H23_ctrl_3",
    "H23_2Gy_1",
    "H23_2Gy_2",
    "H23_2Gy_3",
    "H23_2Gy_4",
    "HCC44_ctrl_1",
    "HCC44_ctrl_2",
    "HCC44_ctrl_3",
    "HCC44_ctrl_4",
    "HCC44_2Gy_1",
    "HCC44_2Gy_2",
    "HCC44_2Gy_3",
    "HCC44_2Gy_4"
  ))

RNA_intensity_data <- as.data.frame(
  log10((RNA_intensity_HCC44_H23_24h[, c(9:17, 3:8)]) + 1)
)

colnames(RNA_intensity_data) <- c(
  "HCC44_ctrl_1",
  "HCC44_ctrl_2",
  "HCC44_ctrl_3",
  "HCC44_ctrl_4",
  "HCC44_2Gy_1",
  "HCC44_2Gy_2",
  "HCC44_2Gy_3",
  "HCC44_2Gy_4",
  "H23_ctrl_1",
  "H23_ctrl_2",
  "H23_ctrl_3",
  "H23_2Gy_1",
  "H23_2Gy_2",
  "H23_2Gy_3",
  "H23_2Gy_4"
)

RNA_intensity_long <- RNA_intensity_data %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  mutate(
    Sample = factor(
      Sample,
      levels = colnames(RNA_intensity_data)
    ),
    Group = case_when(
      grepl("^HCC44_ctrl", Sample) ~ "HCC44 ctrl",
      grepl("^HCC44_2Gy", Sample)  ~ "HCC44 2 Gy",
      grepl("^H23_ctrl", Sample)   ~ "NCI-H23 ctrl",
      grepl("^H23_2Gy", Sample)    ~ "NCI-H23 2 Gy"
    ),
    Group = factor(
      Group,
      levels = c(
        "HCC44 ctrl",
        "HCC44 2 Gy",
        "NCI-H23 ctrl",
        "NCI-H23 2 Gy"
      )
    )
  )

RNA_p_transcriptome_Intensity <- ggplot(
  RNA_intensity_long,
  aes(
    x = Sample,
    y = Count,
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
    size = 2,
    stroke = 0.5,
    fill = "white",
    colour = "black",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = RNA_colours,
    breaks = c(
      "HCC44 ctrl",
      "HCC44 2 Gy",
      "NCI-H23 ctrl",
      "NCI-H23 2 Gy"
    )
  ) +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6),
    minor_breaks = c(1, 3, 5),
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_cartesian(
    ylim = c(0, 7),
    clip = "off"
  ) +
  labs(
    title = "Transcriptome - HCC44 & NCI-H23",
    x = NULL,
    y = expression(log[10] * "(Count + 1)")
  ) +
  theme_prism_xcell() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    axis.ticks.y = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.just = "left",
    legend.text = element_text(size = 9),
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
      b = 40,
      l = 10
    )
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      nrow = 1
    )
  )+
  annotate(
    "text",
    x = 1,
    y = -0.5,
    label = "HCC44 ctrl",
    angle = 30,
    hjust = 1,
    size = 3.5
  ) +
  annotate(
    "text",
    x =5,
    y = -0.5,
    label = "HCC44 2Gy",
    angle = 30,
    hjust = 1,
    size = 3.5
  ) +
  annotate(
    "text",
    x = 9,
    y = -0.5,
    label = "NCI-H23 ctrl",
    angle = 30,
    hjust = 1,
    size = 3.5
  ) +
  annotate(
    "text",
    x = 12,
    y = -0.5,
    label = "NCI-H23 2 Gy",
    angle = 30,
    hjust = 1,
    size =3.5
  )

RNA_p_transcriptome_Intensity <- RNA_p_transcriptome_Intensity +
  theme(
    legend.position = "none"
  )
#--------------------------------------------
# Principal Component Analysis
#--------------------------------------------
RNA_PCA_data <- data.frame(
log2(RNA_RH23_24h[, c(5:11)]),
log2(RNA_RHCC44_24h[, c(5:12)]))

RNA_PCA_data_filt <- RNA_PCA_data %>%
  mutate(
    across(
      everything(),
      ~ replace(.x, !is.finite(.x), NA)
    )
  )

RNA_PCA_data_filt <- RNA_PCA_data_filt[
  complete.cases(RNA_PCA_data_filt),
]

RNA_PCA_data_filt <- t(RNA_PCA_data_filt)

RNA_PCA_res <- prcomp(
  RNA_PCA_data_filt,
  center = TRUE,
  scale. = TRUE
)

RNA_PCA_variance <- RNA_PCA_res$sdev^2 /
  sum(RNA_PCA_res$sdev^2) * 100

RNA_PC1_var <- RNA_PCA_variance[1]
RNA_PC2_var <- RNA_PCA_variance[2]

RNA_PCA_plot_data <- as.data.frame(RNA_PCA_res$x[, 1:2])

RNA_PCA_plot_data$Group <- factor(
  c(
    rep("NCI-H23 ctrl", 3),
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

RNA_p_transcriptome_PCA <- ggplot(
  RNA_PCA_plot_data,
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
    values = RNA_colours,
    breaks = c(
      "HCC44 ctrl",
      "HCC44 2 Gy",
      "NCI-H23 ctrl",
      "NCI-H23 2 Gy"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(-100,0, 50, 100),
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
    title = "Transcriptome, 24 h",
    x = paste0(
      "PCA1 (",
      round(RNA_PC1_var, 1),
      "%)"
    ),
    y = paste0(
      "PCA2 (",
      round(RNA_PC2_var, 1),
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
# Venn Diagram Transcriptome
#--------------------------------------------
RNA_venn_transcriptome <- venn.diagram(
    x = list(
      "HCC44 ctrl"   = RNA_RHCC44_24h[complete.cases(RNA_RHCC44_24h[,13]),2],
      "HCC44 2 Gy"   = RNA_RHCC44_24h[complete.cases(RNA_RHCC44_24h[,15]),2],
      "NCI-H23 ctrl" = RNA_RH23_24h[complete.cases(RNA_RH23_24h[,12]),2],
      "NCI-H23 2 Gy" = RNA_RH23_24h[complete.cases(RNA_RH23_24h[,14]),2]
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
# Heatmap Transcriptome
#--------------------------------------------
RNA_hm_intensity_all <- data.frame(
  RNA_RH23_24h[, c(1:2,5:11)],
  RNA_RHCC44_24h[, 5:12]
)

colnames(RNA_hm_intensity_all) <- c(
  "Symbol",
  "ID",
  "H23 c1",
  "H23 c2",
  "H23 c3",
  "H23 2Gy1",
  "H23 2Gy2",
  "H23 2Gy3",
  "H23 2Gy4",
  
  "HCC44 c1",
  "HCC44 c2",
  "HCC44 c3",
  "HCC44 c4",
  "HCC44 2Gy1",
  "HCC44 2Gy2",
  "HCC44 2Gy3",
  "HCC44 2Gy4"
)

RNA_hm_intensity_all[, 3:17] <- RNA_hm_intensity_all[, 3:17] %>%
  mutate(
    across(everything(), as.numeric)
  )

RNA_hm_intensity_all[, 3:17] <- RNA_hm_intensity_all[, 3:17] %>%
  mutate(
    across(
      everything(),
      ~ replace(.x, !is.finite(.x), NA)
    )
  )

RNA_hm_intensity_all <- RNA_hm_intensity_all[
  complete.cases(RNA_hm_intensity_all[, 3:17]),
]

RNA_hm_intensity_matrix <- as.matrix(
  RNA_hm_intensity_all[, 3:17]
)

rownames(RNA_hm_intensity_matrix) <- RNA_hm_intensity_all$Symbol

RNA_hm_intensity_z <- t(
  scale(
    t(RNA_hm_intensity_matrix)
  )
)

annotation_col <- data.frame(
    Radiation = factor(
    c(
      rep("ctrl", 3),
      rep("2Gy",  4),
      rep("ctrl", 4),
      rep("2Gy",  4)
    ),
    levels = c("ctrl", "2Gy")
  ),
  `Cell line` = factor(
    c(
      rep("NCI-H23", 7),
      rep("HCC44", 8)
    ),
    levels = c("NCI-H23", "HCC44")
  ),
  check.names = FALSE
)

rownames(annotation_col) <- colnames(RNA_hm_intensity_z)

identical(
  colnames(RNA_hm_intensity_z),
  rownames(annotation_col)
)

RNA_p_hm_intensity_all <- pheatmap(
  RNA_hm_intensity_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 12
)
#--------------------------------------------
# Barplots of counts
#--------------------------------------------
RNA_n_Rcount <- data.frame(
  H23_24h_c = count(
    RNA_RH23_24h[complete.cases(RNA_RH23_24h$mean_C), ]
  ),
  H23_24h_2Gy = count(
    RNA_RH23_24h[complete.cases(RNA_RH23_24h$mean_2Gy), ]
  ),
  HCC44_24h_c = count(
    RNA_RHCC44_24h[complete.cases(RNA_RHCC44_24h$mean_C), ]
  ),
  HCC44_24h_2Gy = count(
    RNA_RHCC44_24h[complete.cases(RNA_RHCC44_24h$mean_2Gy), ]
  )
)

RNA_count_data <- data.frame(
  Group = factor(
    c(
      "HCC44 ctrl",
      "HCC44 2 Gy",
      "NCI-H23 ctrl",
      "NCI-H23 2 Gy"
    ),
    levels = c(
      "HCC44 ctrl",
      "HCC44 2 Gy",
      "NCI-H23 ctrl",
      "NCI-H23 2 Gy"
    )
  ),
  Count = as.numeric(
    RNA_n_Rcount[, c(3, 4, 1, 2)]
  )
)

RNA_p_identifiedRNAs <- ggplot(
  RNA_count_data,
  aes(
    x = Group,
    y = Count,
    fill = Group
  )
) +
  geom_hline(
    yintercept = 13751,
    linetype = "dotted",
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_col(
    width = 0.4,
    colour = "black",
    linewidth = 0.75
  ) +
  geom_text(
    aes(
      label = format(
        Count,
        big.mark = ".",
        decimal.mark = ",",
        scientific = FALSE
      )
    ),
    vjust = -0.6,
    hjust = -0.05,
    angle = 35,
    size = 3.5
  ) +
  annotate(
    "text",
    x = 4.3,
    y = 10251,
    label = "Shared RNAs:\n 13,751",
    hjust = 0,
    vjust = 0.1,
    size = 3.5
  ) +
  scale_fill_manual(
    values = RNA_colours
  ) +
  scale_y_continuous(
    limits = c(0, 17500),
    breaks = seq(0, 15000, by = 5000),
    minor_breaks = seq(2500, 20000, by = 2500),
    labels = function(x) {
      format(
        x,
        big.mark = ".",
        decimal.mark = ",",
        scientific = FALSE
      )
    },
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(
      add = c(0.6, 2.2)
    ))+
  labs(
    x = NULL,
    y = "Detected genes (n)"
  ) +
  theme_prism_xcell() +
  theme(
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(
      size =10,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      colour = "black",
      angle = 35,
      hjust = 1.1,
      vjust = 1
    ),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  ) +
  coord_cartesian(
    clip = "off"
  )
#--------------------------------------------
# RPKM
#--------------------------------------------
RNA_RPKM_plot_data <- bind_rows(
  prepare_RPKM(
    RNA_RPKM,
    sample_cols = c(13, 21, 22),
    cellline = "NCI-H23",
    treatment = "ctrl"
  ),
  prepare_RPKM(
    RNA_RPKM,
    sample_cols = c(23, 14:16),
    cellline = "NCI-H23",
    treatment = "2 Gy"
  ),
  prepare_RPKM(
    RNA_RPKM,
    sample_cols = c(9, 17:19),
    cellline = "HCC44",
    treatment = "ctrl"
  ),
  prepare_RPKM(
    RNA_RPKM,
    sample_cols = c(10:12, 20),
    cellline = "HCC44",
    treatment = "2 Gy"
  )
)

RNA_RPKM_plot_data <- RNA_RPKM_plot_data %>%
  mutate(
    Condition = paste(Cellline, Treatment),
    Condition = factor(
      Condition,
      levels = c(
        "NCI-H23 ctrl",
        "NCI-H23 2 Gy",
        "HCC44 ctrl",
        "HCC44 2 Gy"
      )
    )
  )

RNA_type_counts <- RNA_RPKM %>%
  count(gene_biotype, sort = TRUE)

RNA_top5_RNA_types <- RNA_type_counts %>%
  slice_head(n = 5) %>%
  pull(gene_biotype)

RNA_RPKM_plot_data <- RNA_RPKM_plot_data %>%
  mutate(
    Biotype = ifelse(
      gene_biotype %in% RNA_top5_RNA_types,
      gene_biotype,
      "rest"
    ),
    Biotype = factor(
      Biotype,
      levels = c(RNA_top5_RNA_types, "rest")
    ),
    Condition = factor(
      paste(Cellline, Treatment),
      levels = c(
        "NCI-H23 ctrl",
        "NCI-H23 2 Gy",
        "HCC44 ctrl",
        "HCC44 2 Gy"
      )
    )
  )

names(RPKM_colours) <- c(
  RNA_top5_RNA_types,
  "rest"
)
names(RPKM_colours) <- levels(RNA_RPKM_plot_data$Biotype)

RNA_p_RPKM <- ggplot(
  RNA_RPKM_plot_data,
  aes(
    x = Biotype,
    y = Expression,
    fill = Biotype
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    linewidth = 0.5,
    colour = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 15,
    size = 1.5,
    colour = "black"
  ) +
  facet_wrap(
    ~ Condition,
    nrow = 2
  ) +
  scale_fill_manual(
    values = RPKM_colours,
    name = NULL
  ) +
  scale_y_continuous(
    breaks = c(0, 10, 20),
    minor_breaks = seq(0, 20, by = 5)
  ) +
  coord_cartesian(
    ylim = c(-5, 20)
  ) +
  labs(
    x = NULL,
    y = expression(log[2] * "(RPKM + 1)")
  ) +
  theme_prism_xcell() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c(0.2,0),
    legend.location = "plot",
    legend.text = element_text(size = 8),
    legend.key.width  = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    panel.spacing = grid::unit(5, "mm"),
    plot.margin = margin(
      5, 5, 5, 5
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )
#--------------------------------------------
# TPM
#--------------------------------------------
RNA_TPM_plot_data <- bind_rows(
  prepare_TPM(
    RNA_TPM,
    sample_cols = c(13, 21, 22),
    cellline = "NCI-H23",
    treatment = "ctrl"
  ),
  prepare_TPM(
    RNA_TPM,
    sample_cols = c(23, 14:16),
    cellline = "NCI-H23",
    treatment = "2 Gy"
  ),
  prepare_TPM(
    RNA_TPM,
    sample_cols = c(9, 17:19),
    cellline = "HCC44",
    treatment = "ctrl"
  ),
  prepare_TPM(
    RNA_TPM,
    sample_cols = c(10:12, 20),
    cellline = "HCC44",
    treatment = "2 Gy"
  )
)

RNA_type_counts_TPM <- RNA_TPM %>%
  count(gene_biotype, sort = TRUE)

RNA_top5_RNA_types_TPM <- RNA_type_counts_TPM %>%
  slice_head(n = 5) %>%
  pull(gene_biotype)

RNA_TPM_plot_data <- RNA_TPM_plot_data %>%
  mutate(
    Biotype = ifelse(
      gene_biotype %in% RNA_top5_RNA_types_TPM,
      gene_biotype,
      "rest"
    ),
    Biotype = factor(
      Biotype,
      levels = c(
        RNA_top5_RNA_types_TPM,
        "rest"
      )
    ),
    Condition = factor(
      paste(Cellline, Treatment),
      levels = c(
        "NCI-H23 ctrl",
        "NCI-H23 2 Gy",
        "HCC44 ctrl",
        "HCC44 2 Gy"
      )
    )
  )

names(RPKM_colours) <- c(
  RNA_top5_RNA_types_TPM,
  "rest"
)

RNA_p_TPM <- ggplot(
  RNA_TPM_plot_data,
  aes(
    x = Biotype,
    y = Expression,
    fill = Biotype
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    colour = "black",
    linewidth = 0.5
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 15,
    size = 1.5,
    colour = "black"
  ) +
  facet_wrap(
    ~ Condition,
    nrow = 2
  ) +
  scale_fill_manual(
    values = RPKM_colours,
    name = NULL
  ) +
  scale_y_continuous(
    breaks = c(0, 5, 10, 15, 20),
    minor_breaks = seq(2.5, 17.5, by = 5)
  ) +
  coord_cartesian(
    ylim = c(-1, 20)
  ) +
  labs(
    x = NULL,
    y = expression(log[2] * "(TPM + 1)")
  ) +
  theme_prism_xcell() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c(0.2,0),
    legend.location = "plot",
    legend.text = element_text(size = 8),
    legend.key.width  = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    panel.spacing = grid::unit(5, "mm"),
    plot.margin = margin(
      5, 5, 5, 5
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )
#--------------------------------------------
# Coefficient of Variance
#--------------------------------------------
RNA_cov_RNA <- bind_rows(
  calculate_cov_RNA(
    RNA_RHCC44_24h[, 5:8],
    "HCC44",
    "ctrl",
    "24 h"
  ),
  calculate_cov_RNA(
    RNA_RHCC44_24h[, 9:12],
    "HCC44",
    "2 Gy",
    "24 h"
  ),
  calculate_cov_RNA(
    RNA_RH23_24h[, 5:7],
    "NCI-H23",
    "ctrl",
    "24 h"
  ),
  calculate_cov_RNA(
    RNA_RH23_24h[, 8:11],
    "NCI-H23",
    "2 Gy",
    "24 h"
  )
)

RNA_cov_RNA <- RNA_cov_RNA %>%
  mutate(
    Cellline = factor(
      Cellline,
      levels = c(
        "HCC44",
        "NCI-H23"
      )
    ),
    Treatment = factor(
      Treatment,
      levels = c(
        "ctrl",
        "2 Gy"
      )
    )
  )

RNA_cov_density_RNA <- RNA_cov_RNA %>%
  group_by(
    Cellline,
    Treatment,
    Time
  ) %>%
  group_modify(
    ~ make_density(.x)
  ) %>%
  ungroup() 


RNA_p_cov_RNA <- plot_cov_RNA(
  RNA_cov_density_RNA)
#--------------------------------------------
# saving results
#--------------------------------------------
dir.create(
  file.path(dirs$results, "17_Transcriptome_QualityControl"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "17_Transcriptome_QualityControl"))

#1. saving intensity violine plot
ggsave(
  "17_Transcriptome_measured_intensities.pdf",
  RNA_p_transcriptome_Intensity,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

ggsave(
  "17_Transcriptome_measured_intensities.jpg",
  RNA_p_transcriptome_Intensity,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. Principal component analysis plot
ggsave(
  filename = "17_Transcriptome_PCA_24h.pdf",
  plot = RNA_p_transcriptome_PCA,
  width =21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "17_Transcriptome_PCA_24h.jpg",
  plot = RNA_p_transcriptome_PCA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3. Plotting Venn diagramm of proteome and transcriptome
pdf(
  "17_Transcriptome_VennDiagram.pdf",
  width = 9 / 2.54,
  height = 7 / 2.54,
  useDingbats = FALSE
)

grid.newpage()
grid.draw(RNA_venn_transcriptome)
grid.text(
  "Transcriptome",
  x = 0.5,
  y = 0.75,
  gp = gpar(
    fontsize = 14,
    fontface = "bold",
    fontfamily = "sans"
  )
)
dev.off()

pdf(
  "17_Transcriptome_VennDiagram",
  width = 9 / 2.54,
  height = 7 / 2.54,
  useDingbats = FALSE
)
grid.newpage()
grid.draw(RNA_venn_transcriptome)
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

#4. saving Heatmap
pdf(
  "17_Transcriptome_Heatmap_all_unfiltered.pdf",
  width = 21*1/3,
  height = 3.5
)

grid::grid.newpage()
grid::grid.draw(RNA_p_hm_intensity_all$gtable)

dev.off()

#5. Barplot of identified Transcripts
ggsave(
  "17_Transcriptome_Nr_identified_RNAs.jpeg",
  plot = RNA_p_identifiedRNAs,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "17_Transcriptome_Nr_identified_RNAs.pdf",
  plot = RNA_p_identifiedRNAs,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#6. RPKM violine plot
ggsave(
  filename = "17_Transcriptome_RPKM.pdf",
  plot = RNA_p_RPKM,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "17_Transcriptome_RPKM.jpg",
  plot = RNA_p_RPKM,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

#7. TPM violine plot
ggsave(
  filename = "17_Transcriptome_TPM.pdf",
  plot = RNA_p_identifiedRNAs,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "17_Transcriptome_PM.jpg",
  plot = RNA_p_identifiedRNAs,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

#8. Coefficient of Variation plot
ggsave(
  filename = "17_Transcriptome_COV_all.pdf",
  plot = RNA_p_cov_RNA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "17_Transcriptome_COV_all.jpg",
  plot = RNA_p_cov_RNA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
setwd(dirs$code)
