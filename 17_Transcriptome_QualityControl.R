#-------------------------------------------
# loading or renaming data, preparing colors
#------------------------------------------
load_stat_data_R(
  file.path(dirs$results, "14_Transcriptome_DifferentialExpression"),
  mode = "ready"
)

RH23_24h   <- load_RNA("H23", "24h")
RHCC44_24h <- load_RNA("HCC44", "24h")
mean_pHCC44_24h <- load_prot_mean("HCC44", "24h")
mean_pH23_24h <- load_prot_mean("H23", "24h")

setwd(file.path(dirs$RNA))
RPKM<- data.frame(read.csv("RPKM_Bohnenberger.csv"))
TPM<- data.frame(read.csv("TPM_Bohnenberger.csv"))

RPKM_colours <- c("#872111","#EA4025","#F09837","#5F8F29","#3F8D27","#4294F7")

RNA_colours <- c(
  "HCC44 ctrl"   = "#4294F7",
  "HCC44 2 Gy"    = "#EA4025",
  "NCI-H23 ctrl" = "#5F8F29",
  "NCI-H23 2 Gy"  ="#F09837")

venn_colours <- c(
  "#5F8F29",
  "#F09837",
  "#4294F7",
  "#EA4025"
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
# Transcriptome - Boxplot Intensity
#--------------------------------------------
RNA_box_HCC44.H23_24h <- cbind(
  RH23_24h[,c(1:2,5:11)],
  RHCC44_24h[, 5:12]
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
#--------------------------------------------
# Violin plot of Intensities
#--------------------------------------------
RNA_violin_data <- as.data.frame(
  log10((2^RNA_box_HCC44.H23_24h[, c(9:17, 3:8)]) + 1)
)

colnames(RNA_violin_data) <- c(
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

RNA_violin_long <- RNA_violin_data %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  mutate(
    Sample = factor(
      Sample,
      levels = colnames(RNA_violin_data)
    ),
    Group = case_when(
      grepl("^HCC44_ctrl", Sample) ~ "HCC44 ctrl",
      grepl("^HCC44_2Gy", Sample)  ~ "HCC44 2Gy",
      grepl("^H23_ctrl", Sample)   ~ "NCI-H23 ctrl",
      grepl("^H23_2Gy", Sample)    ~ "NCI-H23 2Gy"
    ),
    Group = factor(
      Group,
      levels = c(
        "HCC44 ctrl",
        "HCC44 2Gy",
        "NCI-H23 ctrl",
        "NCI-H23 2Gy"
      )
    )
  )

p_transcriptome_violin <- ggplot(
  RNA_violin_long,
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
    size = 2.5,
    stroke = 0.5,
    fill = "white",
    colour = "black",
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = RNA_colours,
    breaks = c(
      "HCC44 ctrl",
      "HCC44 2Gy",
      "NCI-H23 ctrl",
      "NCI-H23 2Gy"
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
    y = expression(log[10] * "(Count) + 1")
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
    legend.text = element_text(size = 6.5),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(2.5, "mm"),
    legend.key.spacing.x = grid::unit(1, "mm"),
    legend.spacing.x = grid::unit(0, "mm"),
    legend.box.spacing = grid::unit(2, "mm"),
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
    size = 3
  ) +
  annotate(
    "text",
    x =5,
    y = -0.5,
    label = "HCC44 2Gy",
    angle = 30,
    hjust = 1,
    size = 3
  ) +
  annotate(
    "text",
    x = 9,
    y = -0.5,
    label = "NCI-H23 ctrl",
    angle = 30,
    hjust = 1,
    size = 3
  ) +
  annotate(
    "text",
    x = 12,
    y = -0.5,
    label = "NCI-H23 2Gy",
    angle = 30,
    hjust = 1,
    size =3
  )
  p_transcriptome_violin <- p_transcriptome_violin +
  theme(
    legend.position = "none"
  )
#--------------------------------------------
# Principal Component Analysis
#--------------------------------------------
RNA_base_pca_filt <- data.frame(
log2(RH23_24h[, c(12,14)]),
log2(RHCC44_24h[, c(13,15)]))

RNA_base_pca_filt <- RNA_base_pca_filt %>%
  mutate(
    across(
      everything(),
      ~ replace(.x, !is.finite(.x), NA)
    )
  )

RNA_base_pca_filt <- RNA_base_pca_filt[
  complete.cases(RNA_base_pca_filt),
]

RNA_base_pca_filt <- t(RNA_base_pca_filt)

pc_RNA_base_pca_filt <- prcomp(
  RNA_base_pca_filt,
  center = TRUE,
  scale. = TRUE
)

RNA_pca_variance <- pc_RNA_base_pca_filt$sdev^2 /
  sum(pc_RNA_base_pca_filt$sdev^2) * 100

RNA_PC1_var <- RNA_pca_variance[1]
RNA_PC2_var <- RNA_pca_variance[2]


pca_RNA_data <- as.data.frame(pc_RNA_base_pca_filt$x[, 1:2])

pca_RNA_data$Group <- factor(
  c(
    rep("NCI-H23 ctrl", 1),
    rep("NCI-H23 2 Gy", 1),
    rep("HCC44 ctrl", 1),
    rep("HCC44 2 Gy", 1)
  ),
  levels = c(
    "HCC44 ctrl",
    "HCC44 2 Gy",
    "NCI-H23 ctrl",
    "NCI-H23 2 Gy"
  )
)

p_transcriptome_PCA <- ggplot(
  pca_RNA_data,
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
    expand = expansion(mult = c(0, 0))
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
      b = 25,
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
# Venn Diagram Transcriptome + Proteome
#--------------------------------------------
venn_transcriptome <- venn.diagram(
    x = list(
      "NCI-H23 ctrl" = RH23_24h[complete.cases(RH23_24h[,12]),2],
      "NCI-H23 2 Gy" = RH23_24h[complete.cases(RH23_24h[,14]),2],
      "HCC44 ctrl"   = RHCC44_24h[complete.cases(RHCC44_24h[,13]),2],
      "HCC44 2 Gy"   = RHCC44_24h[complete.cases(RHCC44_24h[,15]),2]
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
    cat.cex = 0.8,
    cat.fontfamily = "sans",
    cat.fontface = "plain",
    cat.col = "black",
    cat.pos = c(347.5, 10, 340, 20),
    cat.dist = c(0.25, 0.25, 0.15, 0.15),
    scaled = FALSE,
    main = "Transcriptome",
    main.cex = 1,
    main.fontfamily = "sans",
    main.fontface = "bold",
    margin = 0.1
  )
  
venn_proteome <- venn.diagram(
    x = list(
      "NCI-H23 ctrl" = mean_pH23_24h[complete.cases(mean_pH23_24h[,5]),1],
      "NCI-H23 2 Gy" = mean_pH23_24h[complete.cases(mean_pH23_24h[,7]),1],
      "HCC44 ctrl"   = mean_pHCC44_24h[complete.cases(mean_pHCC44_24h[,5]),1],
      "HCC44 2 Gy"   = mean_pHCC44_24h[complete.cases(mean_pHCC44_24h[,7]),1]
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
    cat.cex = 0.8,
    cat.fontfamily = "sans",
    cat.fontface = "plain",
    cat.col = "black",
    cat.pos = c(347.5, 10, 340, 20),
    cat.dist = c(0.25, 0.25, 0.15, 0.15),
    scaled = FALSE,
    main = "Proteome",
    main.cex = 1,
    main.fontfamily = "sans",
    main.fontface = "bold",
    margin = 0.1,
  )
#--------------------------------------------
# Heatmap Transcriptome
#--------------------------------------------
intens_RNA_all_24h <- data.frame(
  RH23_24h[, c(1:2,5:11)],
  RHCC44_24h[, 5:12]
)

colnames(intens_RNA_all_24h) <- c(
  "Genes",
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

intens_RNA_all_24h[, 3:17] <- intens_RNA_all_24h[, 3:17] %>%
  mutate(
    across(everything(), as.numeric)
  )

intens_RNA_all_24h[, 3:17] <- intens_RNA_all_24h[, 3:17] %>%
  mutate(
    across(
      everything(),
      ~ replace(.x, !is.finite(.x), NA)
    )
  )

intens_RNA_all_24h <- intens_RNA_all_24h[
  complete.cases(intens_RNA_all_24h[, 3:17]),
]

transcriptome_matrix <- as.matrix(
  intens_RNA_all_24h[, 3:17]
)

rownames(transcriptome_matrix) <- intens_RNA_all_24h$Genes

transcriptome_z <- t(
  scale(
    t(transcriptome_matrix)
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

rownames(annotation_col) <- colnames(transcriptome_matrix)

identical(
  colnames(transcriptome_z),
  rownames(annotation_col)
)

Graph_QC_RNA_2 <- pheatmap(
  transcriptome_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = heatmap_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 9
)
#--------------------------------------------
# Barplots of counts
#--------------------------------------------
n_Rcount <- data.frame(
  H23_24h_c = count(
    RH23_24h[complete.cases(RH23_24h$mean_C), ]
  ),
  H23_24h_2Gy = count(
    RH23_24h[complete.cases(RH23_24h$mean_2Gy), ]
  ),
  HCC44_24h_c = count(
    RHCC44_24h[complete.cases(RHCC44_24h$mean_C), ]
  ),
  HCC44_24h_2Gy = count(
    RHCC44_24h[complete.cases(RHCC44_24h$mean_2Gy), ]
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
    n_Rcount[, c(3, 4, 1, 2)]
  )
)

Transcriptome_identifiedRNAs <- ggplot(
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
    linewidth = 0.75,
    colour = "black"
  ) +
  geom_col(
    width = 0.4,
    colour = "black",
    linewidth = 0.5
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
    vjust = -0.5,
    size = 2.5
  ) +
  annotate(
    "text",
    x = 4.3,
    y = 12251,
    label = "Shared RNAs: 13,751",
    hjust = 0,
    vjust = -0.5,
    size = 2.5
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
      size = 7,
      face = "bold"
    ),
    axis.text.y = element_text(
      size =6,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 6,
      colour = "black",
      angle = 35,
      hjust = 1,
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
RPKM_plot_data <- bind_rows(
  prepare_RPKM(
    RPKM,
    sample_cols = c(13, 21, 22),
    cellline = "NCI-H23",
    treatment = "ctrl"
  ),
  prepare_RPKM(
    RPKM,
    sample_cols = c(23, 14:16),
    cellline = "NCI-H23",
    treatment = "2 Gy"
  ),
  prepare_RPKM(
    RPKM,
    sample_cols = c(9, 17:19),
    cellline = "HCC44",
    treatment = "ctrl"
  ),
  prepare_RPKM(
    RPKM,
    sample_cols = c(10:12, 20),
    cellline = "HCC44",
    treatment = "2 Gy"
  )
)

RPKM_plot_data <- RPKM_plot_data %>%
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

RNA_type_counts <- RPKM %>%
  count(gene_biotype, sort = TRUE)

top5_RNA_types <- RNA_type_counts %>%
  slice_head(n = 5) %>%
  pull(gene_biotype)

RPKM_plot_data <- RPKM_plot_data %>%
  mutate(
    Biotype = ifelse(
      gene_biotype %in% top5_RNA_types,
      gene_biotype,
      "rest"
    ),
    Biotype = factor(
      Biotype,
      levels = c(top5_RNA_types, "rest")
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
  top5_RNA_types,
  "rest"
)
names(RPKM_colours) <- levels(RPKM_plot_data$Biotype)

Transcriptome_RPKM <- ggplot(
  RPKM_plot_data,
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
      size = 10,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7),
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
TPM_plot_data <- bind_rows(
  prepare_TPM(
    TPM,
    sample_cols = c(13, 21, 22),
    cellline = "NCI-H23",
    treatment = "ctrl"
  ),
  prepare_TPM(
    TPM,
    sample_cols = c(23, 14:16),
    cellline = "NCI-H23",
    treatment = "2 Gy"
  ),
  prepare_TPM(
    TPM,
    sample_cols = c(9, 17:19),
    cellline = "HCC44",
    treatment = "ctrl"
  ),
  prepare_TPM(
    TPM,
    sample_cols = c(10:12, 20),
    cellline = "HCC44",
    treatment = "2 Gy"
  )
)

RNA_type_counts_TPM <- TPM %>%
  count(gene_biotype, sort = TRUE)

top5_RNA_types_TPM <- RNA_type_counts_TPM %>%
  slice_head(n = 5) %>%
  pull(gene_biotype)

TPM_plot_data <- TPM_plot_data %>%
  mutate(
    Biotype = ifelse(
      gene_biotype %in% top5_RNA_types_TPM,
      gene_biotype,
      "rest"
    ),
    Biotype = factor(
      Biotype,
      levels = c(
        top5_RNA_types_TPM,
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
  top5_RNA_types_TPM,
  "rest"
)

Transcriptome_TPM <- ggplot(
  TPM_plot_data,
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
      size = 10,
      face = "bold",
      colour = "black"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 7),
    legend.key.width  = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    panel.spacing = grid::unit(5, "mm"),
    
    plot.margin = margin(
      5, 5, 5, 5
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )
#--------------------------------------------
# Coefficient of Variance
#--------------------------------------------
cov_RNA <- bind_rows(
  calculate_cov_RNA(
    RHCC44_24h[, 3:6],
    "HCC44",
    "ctrl",
    "24 h"
  ),
  calculate_cov_RNA(
    RHCC44_24h[, 7:10],
    "HCC44",
    "2 Gy",
    "24 h"
  ),
  calculate_cov_RNA(
    RH23_24h[, 3:5],
    "NCI-H23",
    "ctrl",
    "24 h"
  ),
  calculate_cov_RNA(
    RH23_24h[, 6:9],
    "NCI-H23",
    "2 Gy",
    "24 h"
  )
)

cov_RNA <- cov_RNA %>%
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

cov_density_RNA <- cov_RNA %>%
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
p_cov_RNA <- plot_cov_RNA(
  cov_density_RNA)
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
  "Transcriptome_violineplot.pdf",
  p_transcriptome_violin,
  width = 21*2/3,
  height = 7,
  units = "cm"
)
ggsave(
  "Transcriptome_violineplot.jpg",
  p_transcriptome_violin,
  width = 21*2/3,
  height = 7,
  units = "cm"
)

#2. Principal component analysis plot
ggsave(
  filename = "Transcriptome_qc_PCA.pdf",
  plot = p_transcriptome_PCA,
  width =21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_qc_PCA.jpg",
  plot = p_transcriptome_PCA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#3. Plotting Venn diagramm of proteome and transcriptome
grid.newpage()
grid.draw(venn_transcriptome)

grid.newpage()
grid.draw(venn_proteome)

pushViewport(
  viewport(
    layout.pos.row = 2,
    layout.pos.col = 1
  )
)

grid.draw(venn_transcriptome)
upViewport()
pushViewport(
  viewport(
    layout.pos.row = 2,
    layout.pos.col = 2
  )
)
grid.draw(venn_proteome)
upViewport()

pdf(
  "Venn_Transcriptome.pdf",
  width = 9 / 2.54,
  height = 7 / 2.54,
  useDingbats = FALSE
)

grid.newpage()
grid.draw(venn_transcriptome)

dev.off()

pdf(
  "Venn_Proteome.pdf",
  width = 9 / 2.54,
  height = 7 / 2.54,
  useDingbats = FALSE
)

grid.newpage()
grid.draw(venn_proteome)

dev.off()

#4. saving Heatmap
pdf(
  "Transcriptome_Heatmap.pdf",
  width = 21*1/3,
  height = 3.5
)

grid::grid.newpage()
grid::grid.draw(Graph_QC_RNA_2$gtable)

dev.off()

#5. Barplot of identified Transcripts
ggsave(
  "Transcriptome_identifiedRNAs.jpeg",
  plot = Transcriptome_identifiedRNAs,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "Transcriptome_identifiedRNAs.pdf",
  plot = Transcriptome_identifiedRNAs,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#6. RPKM violine plot
ggsave(
  filename = "Transcriptome_RPKM.pdf",
  plot = Transcriptome_RPKM,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_RPKM.jpg",
  plot = Transcriptome_RPKM,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

#7. TPM violine plot
ggsave(
  filename = "Transcriptome_TPM.pdf",
  plot = Transcriptome_RPKM,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_TPM.jpg",
  plot = Transcriptome_RPKM,
  width = 21*1/2,
  height = 7,
  units = "cm"
)

#8. Coefficient of Variation plot
ggsave(
  filename = "Transcriptome_qc_COV.pdf",
  plot = p_cov_RNA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_qc_COV.jpg",
  plot = p_cov_RNA,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
setwd(dirs$code)
