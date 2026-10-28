#-------------------------------------------
# loading or renaming data
#------------------------------------------
load_stat_data_R(
  file.path(dirs$results, "14_Transcriptome_DifferentialExpression"),
  mode = "ready"
)
#--------------------------------------------
# set path
#--------------------------------------------
dir.create(
  file.path(dirs$results, "18_Transcriptome_Volcano_top10"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "18_Transcriptome_Volcano_top10"))
#--------------------------------------------
# 24h | Volcano HCC44 - Plotting
#--------------------------------------------
p_volcano_RNA_HCC44 <- plot_volcano(
  data = ready_RHCC44_24h,
  title = "Transcriptome, HCC44, 24 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-5, 5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 24h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Volcano_HCC44.pdf",
  plot = p_volcano_RNA_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Volcano_HCC44.jpg",
  plot = p_volcano_RNA_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 HCC44 - filtering
#--------------------------------------------
top10_HCC44_down <- ready_RHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_HCC44_up <- ready_RHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_HCC44 <- bind_rows(
  top10_HCC44_down,
  top10_HCC44_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_HCC44 <- top10_HCC44 %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 24h | Top 10 HCC44 - plotting
#--------------------------------------------
p_RNA_top10_HCC44 <- ggplot(
  top10_HCC44,
  aes(
    x = log2FoldChange,
    y = Gene
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FoldChange,
      y = Gene,
      yend = Gene
    ),
    linewidth = 0.5,
    colour = "grey50"
  ) +
    geom_point(
    aes(fill = neg_log10_p),
    shape = 21,
    size = 2.5,
    stroke = 0.5,
    colour = "black"
  ) +
    geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    colour = "black"
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "#D53E4F",
    name = expression(-log[10] * "(p-value)")
  )+
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Transcriptome differences, HCC44, 24 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
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
  axis.ticks.length = grid::unit(1.5, "mm"),
  axis.text.x = element_text(
      size = 7,
      colour = "black"
    ),
  axis.text.y = element_text(
      size = 6,
      colour = "black"
    ),
  axis.title.x = element_text(
      size = 8,
      face = "bold"
    ),
  plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5
    ),
  legend.position = "right",
  legend.title = element_text(
      size = 6
    ),
  legend.text = element_text(
      size = 5
    ),
  plot.margin = margin(
      10, 10, 5, 10
    )
  )
#--------------------------------------------
# 24h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Top10_HCC44.pdf",
  plot = p_RNA_top10_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Top10_HCC44.jpg",
  plot = p_RNA_top10_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano H23 - Plotting
#--------------------------------------------
p_volcano_RNA_H23 <- plot_volcano(
  data = ready_RH23_24h,
  title = "Transcriptome, NCI-H23, 24 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-5, 5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 24h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Volcano_H23.pdf",
  plot = p_volcano_RNA_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Volcano_H23.jpg",
  plot = p_volcano_RNA_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 H23 - filtering
#--------------------------------------------
top10_H23_down <- ready_RH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_H23_up <- ready_RH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_H23 <- bind_rows(
  top10_H23_down,
  top10_H23_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_H23 <- top10_H23 %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 24h | Top 10 H23 - plotting
#--------------------------------------------
p_RNA_top10_H23 <- ggplot(
  top10_H23,
  aes(
    x = log2FoldChange,
    y = Gene
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FoldChange,
      y = Gene,
      yend = Gene
    ),
    linewidth = 0.5,
    colour = "grey50"
  ) +
  geom_point(
    aes(fill = neg_log10_p),
    shape = 21,
    size = 2.5,
    stroke = 0.5,
    colour = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    colour = "black"
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "#D53E4F",
    name = expression(-log[10] * "(p-value)")
  )+
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Transcriptome differences, NCI-H23, 24 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
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
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 7,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 6,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 8,
      face = "bold"
    ),
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 6
    ),
    legend.text = element_text(
      size = 5
    ),
    plot.margin = margin(
      10, 10, 5, 10
    )
  )
#--------------------------------------------
# 24h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Top10_H23.pdf",
  plot = p_RNA_top10_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Top10_H23.jpg",
  plot = p_RNA_top10_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano Plot Control - Plotting
#--------------------------------------------
p_volcano_RNA_Con <- plot_volcano(
  data = ready_RCon_24h,
  title = "Transcriptome, Control, 24 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 10)
)
#--------------------------------------------
# 24h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Volcano_Con.pdf",
  plot = p_volcano_RNA_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Volcano_Con.jpg",
  plot = p_volcano_RNA_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Control - filtering
#--------------------------------------------
top10_Con_down <- ready_RCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_Con_up <- ready_RCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_Con <- bind_rows(
  top10_Con_down,
  top10_Con_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_Con <- top10_Con %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 24h | Top 10 Control - plotting
#--------------------------------------------
p_RNA_top10_Con <- ggplot(
  top10_Con,
  aes(
    x = log2FoldChange,
    y = Gene
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FoldChange,
      y = Gene,
      yend = Gene
    ),
    linewidth = 0.5,
    colour = "grey50"
  ) +
  geom_point(
    aes(fill = neg_log10_p),
    shape = 21,
    size = 2.5,
    stroke = 0.5,
    colour = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    colour = "black"
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "#D53E4F",
    name = expression(-log[10] * "(p-value)")
  )+
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Transcriptome differences, Control, 24 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
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
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 7,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 6,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 8,
      face = "bold"
    ),
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 6
    ),
    legend.text = element_text(
      size = 5
    ),
    plot.margin = margin(
      10, 10, 5, 10
    )
  )
#--------------------------------------------
# 24h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Top10_Con.pdf",
  plot = p_RNA_top10_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Top10_Con.jpg",
  plot = p_RNA_top10_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Radiation - filtering
#--------------------------------------------
top10_RAD_down <- ready_RRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_RAD_up <- ready_RRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_RAD <- bind_rows(
  top10_RAD_down,
  top10_RAD_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_RAD <- top10_RAD %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 24h | Volcano Plot Radiation - Plotting
#--------------------------------------------
p_volcano_RNA_Rad <- plot_volcano(
  data = ready_RRAD_24h,
  title = "Transcriptome, Radiation, 24 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 10)
)
#--------------------------------------------
# 24h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Volcano_RAD.pdf",
  plot = p_volcano_RNA_Rad,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Volcano_RAD.jpg",
  plot = p_volcano_RNA_Rad,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Radiatiom - filtering
#--------------------------------------------
top10_RAD_down <- ready_RRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_RAD_up <- ready_RRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_RAD <- bind_rows(
  top10_RAD_down,
  top10_RAD_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_RAD <- top10_RAD %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 24h | Top 10 Radiation - plotting
#--------------------------------------------
p_RNA_top10_RAD <- ggplot(
  top10_RAD,
  aes(
    x = log2FoldChange,
    y = Gene
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = log2FoldChange,
      y = Gene,
      yend = Gene
    ),
    linewidth = 0.5,
    colour = "grey50"
  ) +
  geom_point(
    aes(fill = neg_log10_p),
    shape = 21,
    size = 2.5,
    stroke = 0.5,
    colour = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    colour = "black"
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "#D53E4F",
    name = expression(-log[10] * "(p-value)")
  )+
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Transcriptome differences, Radiation, 24 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
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
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 7,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 6,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 8,
      face = "bold"
    ),
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 6
    ),
    legend.text = element_text(
      size = 5
    ),
    plot.margin = margin(
      10, 10, 5, 10
    )
  )
#--------------------------------------------
# 24h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "Transcriptome_Top10_RAD.pdf",
  plot = p_RNA_top10_RAD,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Transcriptome_Top10_RAD.jpg",
  plot = p_RNA_top10_RAD,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#-----------------------------------------------
setwd(dirs$code)
