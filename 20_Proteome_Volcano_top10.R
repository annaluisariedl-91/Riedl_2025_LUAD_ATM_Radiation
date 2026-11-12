#-------------------------------------------
# loading or renaming data
#-------------------------------------------
load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)

filter_proteomes()
#--------------------------------------------
# setting path
#--------------------------------------------
dir.create(
  file.path(dirs$results, "20_Proteome_Volcano_top10"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "20_Proteome_Volcano_top10"))
#--------------------------------------------
# 0h | Volcano HCC44 - Plotting
#--------------------------------------------
Protein_p_volcano_HCC44_0h <- plot_volcano(
  data = ready_pHCC44_0h,
  title = "Proteome, HCC44, 0 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-12.5, 12.5),
  ylim = c(0, 7.5)
  )
#--------------------------------------------
# 0h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_HCC44_0h.pdf",
  plot = Protein_p_volcano_HCC44_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_HCC44_0h.jpg",
  plot = Protein_p_volcano_HCC44_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 0h | Top 10 HCC44 - filtering
#--------------------------------------------
Protein_top10_HCC44_0h_down <- ready_pHCC44_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_HCC44_0h_up <- ready_pHCC44_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_HCC44_0h <- bind_rows(
  Protein_top10_HCC44_0h_down,
  Protein_top10_HCC44_0h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_HCC44_0h <- Protein_top10_HCC44_0h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 0h | Top 10 HCC44 - plotting
#--------------------------------------------
Protein_p_top10_HCC44_0h <- ggplot(
  Protein_top10_HCC44_0h,
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
    linewidth = 0.4,
    colour = "grey50"
  ) +
  geom_point(
    aes(fill = neg_log10_p),
    shape = 21,
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n HCC44, 0 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 8,
      face = "bold"
      ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 0h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_HCC44_0h.pdf",
  plot = Protein_p_top10_HCC44_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_HCC44_0h.jpg",
  plot = Protein_p_top10_HCC44_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 0h | Volcano H23 - Plotting
#--------------------------------------------
p_volcanop_p_H23_0h <- plot_volcano(
  data = ready_pH23_0h,
  title = "Proteome, NCI-H23, 0 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 0h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_H23_0h.pdf",
  plot = p_volcanop_p_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_H23_0h.jpg",
  plot = p_volcanop_p_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 0h | Top 10 H23 - filtering
#--------------------------------------------
Protein_top10_H23_0h_down <- ready_pH23_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_H23_0h_up <- ready_pH23_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_H23_0h <- bind_rows(
  Protein_top10_H23_0h_down,
  Protein_top10_H23_0h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_H23_0h <- Protein_top10_H23_0h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 0h | Top 10 H23 - plotting
#--------------------------------------------
Protein_p_top10_H23_0h <- ggplot(
  Protein_top10_H23_0h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0.2))
  )+
  coord_cartesian(xlim = c(-12.5, 16))+
  labs(
    title = "Proteome differences \n NCI-H23, 0 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  )  +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 0h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_0h.pdf",
  plot = Protein_p_top10_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_0h.jpg",
  plot = Protein_p_top10_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 0h | Volcano Plot Control - Plotting
#--------------------------------------------
Protein_p_volcano_Con_0h <- plot_volcano(
  data = ready_pCon_0h,
  title = "Proteome, Control, 0 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-12.5, 12.5),
  ylim = c(0, 17.5),
  extra_genes = c("CSRP2","SLC2A3","PALLD","MGST1")
)
#--------------------------------------------
# 0h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_Con_0h.pdf",
  plot = Protein_p_volcano_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_Con_0h.jpg",
  plot = Protein_p_volcano_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 0h | Top 10 Control - filtering
#--------------------------------------------
Protein_top10_Con_0h_down <- ready_pCon_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Con_0h_up <- ready_pCon_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Con_0h <- bind_rows(
  Protein_top10_Con_0h_down,
  Protein_top10_Con_0h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Con_0h <- Protein_top10_Con_0h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 0h | Top 10 Control - plotting
#--------------------------------------------
Protein_p_top10_Con_0h <- ggplot(
  Protein_top10_Con_0h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-14, 15))+
  labs(
    title = "Proteome differences \n Control, 0 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 0h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_Con_0h.pdf",
  plot = Protein_p_top10_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_Con_0h.jpg",
  plot = Protein_p_top10_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Volcano HCC44 - Plotting
#--------------------------------------------
Protein_p_volcano_HCC44_2h <- plot_volcano(
  data = ready_pHCC44_2h,
  title = "Proteome, HCC44, 2 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-10, 10),
  ylim = c(0, 5),
  extra_genes = c("CCM2L","L3MBTL3","HDAC6","APT1A2")
)
#--------------------------------------------
# 2h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_HCC44_2h.pdf",
  plot = Protein_p_volcano_HCC44_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_HCC44_2h.jpg",
  plot = Protein_p_volcano_HCC44_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 2h | Top 10 HCC44 - filtering
#--------------------------------------------
Protein_top10_HCC44_2h_down <- ready_pHCC44_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_HCC44_2h_up <- ready_pHCC44_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_HCC44_2h <- bind_rows(
  Protein_top10_HCC44_2h_down,
  Protein_top10_HCC44_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_HCC44_2h <- Protein_top10_HCC44_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 2h | Top 10 HCC44 - plotting
#--------------------------------------------
Protein_p_top10_HCC44_2h <- ggplot(
  Protein_top10_HCC44_2h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-10, 15))+
  labs(
    title = "Proteome differences \n HCC44, 2 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 2h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_HCC44_2h.pdf",
  plot = Protein_p_top10_HCC44_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_HCC44_2h.jpg",
  plot = Protein_p_top10_HCC44_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Volcano H23 - Plotting
#--------------------------------------------
p_volcanop_p_H23_2h <- plot_volcano(
  data = ready_pH23_2h,
  title = "Proteome, NCI-H23, 2 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-12.5, 12.5),
  ylim = c(0, 5),
  extra_genes = c("ZNF682","LRFN3","XK")
  )
#--------------------------------------------
# 2h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_H23_2h.pdf",
  plot = p_volcanop_p_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_H23_2h.jpg",
  plot = p_volcanop_p_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Top 10 H23 - filtering
#--------------------------------------------
Protein_top10_H23_2h_down <- ready_pH23_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_H23_2h_up <- ready_pH23_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_H23_2h <- bind_rows(
  Protein_top10_H23_2h_down,
  Protein_top10_H23_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_H23_2h <- Protein_top10_H23_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 2h | Top 10 H23 - plotting
#--------------------------------------------
Protein_p_top10_H23_2h <- ggplot(
  Protein_top10_H23_2h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n NCI-H23, 2 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 2h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_2h.pdf",
  plot = Protein_p_top10_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_2h.jpg",
  plot = Protein_p_top10_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Volcano Plot Control - Plotting
#--------------------------------------------
Protein_p_volcano_Con_2h <- plot_volcano(
  data = ready_pCon_2h,
  title = "Proteome, Control, 2 h",
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
# 2h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_Con_2h.pdf",
  plot = Protein_p_volcano_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_Con_2h.jpg",
  plot = Protein_p_volcano_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Top 10 Control - filtering
#--------------------------------------------
Protein_top10_Con_2h_down <- ready_pCon_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Con_2h_up <- ready_pCon_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Con_2h <- bind_rows(
  Protein_top10_Con_2h_down,
  Protein_top10_Con_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Con_2h <- Protein_top10_Con_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = as.character(.[[3]])
  ) %>%
  group_by(Gene) %>%
  mutate(
    Gene = if (n() > 1) {
      paste0(Gene, "_", row_number())
    } else {
      Gene
    }
  ) %>%
  ungroup() %>%
  mutate(
    Gene = factor(Gene, levels = Gene)
  )
#--------------------------------------------
# 2h | Top 10 Control - plotting
#--------------------------------------------
Protein_p_top10_Con_2h <- ggplot(
  Protein_top10_Con_2h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Control, 2 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 2h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_Con_2h.pdf",
  plot = Protein_p_top10_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_Con_2h.jpg",
  plot = Protein_p_top10_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Volcano Plot Radiation - Plotting
#--------------------------------------------
Protein_p_volcano_Rad_2h <- plot_volcano(
  data = ready_pRAD_2h,
  title = "Proteome, Radiation, 2 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-12.5, 12.5),
  ylim = c(0, 5),
  extra_genes = c("FAM120C","WNT5A")
)
#--------------------------------------------
# 2h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Proteome_Volcano_RAD_2h.pdf",
  plot = Protein_p_volcano_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_2h.jpg",
  plot = Protein_p_volcano_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Top 10 Radiation - filtering
#--------------------------------------------
Protein_top10_RAD_2h_down <- ready_pRAD_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_RAD_2h_up <- ready_pRAD_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_RAD_2h <- bind_rows(
  Protein_top10_RAD_2h_down,
  Protein_top10_RAD_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_RAD_2h <- Protein_top10_RAD_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 2h | Top 10 Radiation - plotting
#--------------------------------------------
Protein_p_top10_Rad_2h <- ggplot(
  Protein_top10_RAD_2h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Radiation, 2 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 2h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_RAD_2h.pdf",
  plot = Protein_p_top10_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_RAD_2h.jpg",
  plot = Protein_p_top10_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Volcano HCC44 - Plotting
#--------------------------------------------
Protein_p_volcano_HCC44_6h <- plot_volcano(
  data = ready_pHCC44_6h,
  title = "Proteome, HCC44, 6 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-10, 10),
  ylim = c(0, 7.5),
  extra_genes = c("ALD1A1","WDR64","ESYT2","APT1A2")
  )
#--------------------------------------------
# 6h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_HCC44_6h.pdf",
  plot = Protein_p_volcano_HCC44_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_HCC44_6h.jpg",
  plot = Protein_p_volcano_HCC44_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 6h | Top 10 HCC44 - filtering
#--------------------------------------------
Protein_top10_HCC44_6h_down <- ready_pHCC44_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_HCC44_6h_up <- ready_pHCC44_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_HCC44_6h <- bind_rows(
  Protein_top10_HCC44_6h_down,
  Protein_top10_HCC44_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_HCC44_6h <- Protein_top10_HCC44_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 6h | Top 10 HCC44 - plotting
#--------------------------------------------
Protein_p_top10_HCC44_6h <- ggplot(
  Protein_top10_HCC44_6h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n HCC44, 6 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 6h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_HCC44_6h.pdf",
  plot = Protein_p_top10_HCC44_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_HCC44_6h.jpg",
  plot = Protein_p_top10_HCC44_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Volcano H23 - Plotting
#--------------------------------------------
p_volcanop_p_H23_6h <- plot_volcano(
  data = ready_pH23_6h,
  title = "Proteome, NCI-H23, 6 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 5),
  extra_genes = c("TPRA1")
  
)
#--------------------------------------------
# 6h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_H23_6h.pdf",
  plot = p_volcanop_p_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_H23_6h.jpg",
  plot = p_volcanop_p_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 6h | Top 10 H23 - filtering
#--------------------------------------------
Protein_top10_H23_6h_down <- ready_pH23_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_H23_6h_up <- ready_pH23_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_H23_6h <- bind_rows(
  Protein_top10_H23_6h_down,
  Protein_top10_H23_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_H23_6h <- Protein_top10_H23_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 6h | Top 10 H23 - plotting
#--------------------------------------------
Protein_p_top10_H23_6h <- ggplot(
  Protein_top10_H23_6h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n NCI-H23, 6 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 6h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_6h.pdf",
  plot = Protein_p_top10_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_6h.jpg",
  plot = Protein_p_top10_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Volcano Plot Control - Plotting
#--------------------------------------------
Protein_p_volcano_Con_6h <- plot_volcano(
  data = ready_pCon_6h,
  title = "Proteome, Control, 6 h",
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
# 6h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_Con_6h.pdf",
  plot = Protein_p_volcano_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_Con_6h.jpg",
  plot = Protein_p_volcano_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Top 10 Control - filtering
#--------------------------------------------
Protein_top10_Con_6h_down <- ready_pCon_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Con_6h_up <- ready_pCon_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Con_6h <- bind_rows(
  Protein_top10_Con_6h_down,
  Protein_top10_Con_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Con_6h <- Protein_top10_Con_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 6h | Top 10 Control - plotting
#--------------------------------------------
Protein_p_top10_Con_6h <- ggplot(
  Protein_top10_Con_6h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Control, 6 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 6h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_Con_6h.pdf",
  plot = Protein_p_top10_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_Con_6h.jpg",
  plot = Protein_p_top10_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Volcano Plot Radiation - Plotting
#--------------------------------------------
Protein_p_volcano_Rad_6h <- plot_volcano(
  data = ready_pRAD_6h,
  title = "Proteome, Radiation, 6 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-12.5, 12.5),
  ylim = c(0, 10),
  extra_genes = c("CKMT1A","ASPH")
)
#--------------------------------------------
# 6h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Proteome_Volcano_RAD_6h.pdf",
  plot = Protein_p_volcano_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_6h.jpg",
  plot = Protein_p_volcano_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Top 10 Radiation - filtering
#--------------------------------------------
Protein_top10_Rad_6h_down <- ready_pRAD_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Rad_6h_up <- ready_pRAD_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Rad_6h <- bind_rows(
  Protein_top10_Rad_6h_down,
  Protein_top10_Rad_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Rad_6h <- Protein_top10_Rad_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 6h | Top 10 Radiation - plotting
#--------------------------------------------
Protein_p_top10_Rad_6h <- ggplot(
  Protein_top10_Rad_6h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Radiation, 6 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 6h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_RAD_6h.pdf",
  plot = Protein_p_top10_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_RAD_6h.jpg",
  plot = Protein_p_top10_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano HCC44 - Plotting
#--------------------------------------------
Protein_p_volcano_HCC44_24h <- plot_volcano(
  data = ready_pHCC44_24h,
  title = "Proteome, HCC44, 24 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-10, 10),
  ylim = c(0, 7.5),
  extra_genes = c("KDF1","KCNK6")
)
#--------------------------------------------
# 24h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_HCC44_24h.pdf",
  plot = Protein_p_volcano_HCC44_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_HCC44_24h.jpg",
  plot = Protein_p_volcano_HCC44_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 HCC44 - filtering
#--------------------------------------------
Protein_top10_HCC44_24h_down <- ready_pHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_HCC44_24h_up <- ready_pHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_HCC44_24h <- bind_rows(
  Protein_top10_HCC44_24h_down,
  Protein_top10_HCC44_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_HCC44_24h <- Protein_top10_HCC44_24h %>%
  mutate(
    SYMBOL = trimws(as.character(SYMBOL)),
    name_missing = is.na(SYMBOL) |
      SYMBOL %in% c("", "NA"),
    unknown_number = cumsum(name_missing),
    Gene = dplyr::if_else(
      name_missing,
      paste0("Unknown_", unknown_number),
      SYMBOL
    )
  ) |>
  dplyr::select(-name_missing, -unknown_number)

Protein_top10_HCC44_24h <- Protein_top10_HCC44_24h %>%
  arrange(log2FoldChange) %>%
  mutate(
    SYMBOL = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 24h | Top 10 HCC44 - plotting
#--------------------------------------------
Protein_p_top10_HCC44_24h <- ggplot(
  Protein_top10_HCC44_24h,
  aes(
    x = log2FoldChange,
    y = reorder(Gene, log2FoldChange)
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n HCC44, 24 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 24h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_HCC44_24h.pdf",
  plot = Protein_p_top10_HCC44_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_HCC44_24h.jpg",
  plot = Protein_p_top10_HCC44_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano H23 - Plotting
#--------------------------------------------
p_volcanop_p_H23_24h <- plot_volcano(
  data = ready_pH23_24h,
  title = "Proteome, NCI-H23, 24 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-10, 10),
  ylim = c(0,5),
  extra_genes = c("FMN1")
  )
#--------------------------------------------
# 24h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_H23_24h.pdf",
  plot = p_volcanop_p_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_H23_24h.jpg",
  plot = p_volcanop_p_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 H23 - filtering
#--------------------------------------------
Protein_top10_H23_24h_down <- ready_pH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -0.9
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_H23_24h_up <- ready_pH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_H23_24h <- bind_rows(
  Protein_top10_H23_24h_down,
  Protein_top10_H23_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_H23_24h <- Protein_top10_H23_24h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 24h | Top 10 H23 - plotting
#--------------------------------------------
Protein_p_top10_H23_24h <- ggplot(
  Protein_top10_H23_24h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n NCI-H23, 24 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 24h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_24h.pdf",
  plot = Protein_p_top10_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_24h.jpg",
  plot = Protein_p_top10_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano Plot Control - Plotting
#--------------------------------------------
Protein_p_volcano_Con_24h <- plot_volcano_RNA(
  data = ready_pCon_24h,
  title = "Proteome, Control, 24 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 10),
  extra_genes=c("MGST1", "KDF1", "CRIP1" , "MAGEA4")
)
#--------------------------------------------
# 24h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_Con_24h.pdf",
  plot = Protein_p_volcano_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_Con_24h.jpg",
  plot = Protein_p_volcano_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Control - filtering
#--------------------------------------------
Protein_top10_Con_24h_down <- ready_pCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Con_24h_up <- ready_pCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Con_24h <- bind_rows(
  Protein_top10_Con_24h_down,
  Protein_top10_Con_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Con_24h <- Protein_top10_Con_24h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 24h | Top 10 Control - plotting
#--------------------------------------------
Protein_p_top10_Con_24h <- ggplot(
  Protein_top10_Con_24h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Control, 24 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 24h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_Con_24h.pdf",
  plot = Protein_p_top10_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_Con_24h.jpg",
  plot = Protein_p_top10_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano Plot Radiation - Plotting
#--------------------------------------------
Protein_p_volcano_Rad_24h <- plot_volcano(
  data = ready_pRAD_24h,
  title = "Proteome, Radiation, 24 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-12.5, 12.5),
  ylim = c(0, 12.5)
)
#--------------------------------------------
# 24h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Proteome_Volcano_RAD_24h.pdf",
  plot = Protein_p_volcano_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_24h.jpg",
  plot = Protein_p_volcano_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Radiation - filtering
#--------------------------------------------
Protein_top10_Rad_24h_down <- ready_pRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Rad_24h_up <- ready_pRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Rad_24h <- bind_rows(
  Protein_top10_Rad_24h_down,
  Protein_top10_Rad_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Rad_24h <- Protein_top10_Rad_24h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 24h | Top 10 Radiation - plotting
#--------------------------------------------
Protein_p_top10_Rad_24h <- ggplot(
  Protein_top10_Rad_24h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Radiation, 24 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 24h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_RAD_24h.pdf",
  plot = Protein_p_top10_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_RAD_24h.jpg",
  plot = Protein_p_top10_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Volcano HCC44 - Plotting
#--------------------------------------------
Protein_p_volcano_HCC44_48h <- plot_volcano(
  data = ready_pHCC44_48h,
  title = "Proteome, HCC44, 48 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 7.5),
  extra_genes = c("DIAPH2","ZNF580")
)
#--------------------------------------------
# 48h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_HCC44_48h.pdf",
  plot = Protein_p_volcano_HCC44_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_HCC44_48h.jpg",
  plot = Protein_p_volcano_HCC44_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 48h | Top 10 HCC44 - filtering
#--------------------------------------------
Protein_top10_HCC44_48h_down <- ready_pHCC44_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_HCC44_48h_up <- ready_pHCC44_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_HCC44_48h <- bind_rows(
  Protein_top10_HCC44_48h_down,
  Protein_top10_HCC44_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_HCC44_48h <- Protein_top10_HCC44_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 48h | Top 10 HCC44 - plotting
#--------------------------------------------
Protein_p_top10_HCC44_48h <- ggplot(
  Protein_top10_HCC44_48h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n HCC44, 48 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 48h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_HCC44_48h.pdf",
  plot = Protein_p_top10_HCC44_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_HCC44_48h.jpg",
  plot = Protein_p_top10_HCC44_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Volcano H23 - Plotting
#--------------------------------------------
p_volcanop_p_H23_48h <- plot_volcano(
  data = ready_pH23_48h,
  title = "Proteome, NCI-H23, 48 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-10, 10),
  ylim = c(0, 5),
  extra_genes = c("TDRD1","ATP5F1B")
  )
#--------------------------------------------
# 48h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_H23_48h.pdf",
  plot = p_volcanop_p_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_H23_48h.jpg",
  plot = p_volcanop_p_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 48h | Top 10 H23 - filtering
#--------------------------------------------
Protein_top10_H23_48h_down <- ready_pH23_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_H23_48h_up <- ready_pH23_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_H23_48h <- bind_rows(
  Protein_top10_H23_48h_down,
  Protein_top10_H23_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_H23_48h <- Protein_top10_H23_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )

Protein_top10_H23_48h <- Protein_top10_H23_48h %>%
  mutate(
    SYMBOL = trimws(as.character(SYMBOL)),
    name_missing = is.na(SYMBOL) |
      SYMBOL %in% c("", "NA"),
    unknown_number = cumsum(name_missing),
    Gene = dplyr::if_else(
      name_missing,
      paste0("Unknown_", unknown_number),
      SYMBOL
    )
  ) |>
  dplyr::select(-name_missing, -unknown_number)

Protein_top10_H23_48h <- Protein_top10_H23_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    SYMBOL = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 48h | Top 10 H23 - plotting
#--------------------------------------------
Protein_p_top10_H23_48h <- ggplot(
  Protein_top10_H23_48h,
  aes(
    x = log2FoldChange,
    y = reorder(Gene, log2FoldChange)
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n NCI-H23, 48 h",
    x = expression(log[2] * "(2 Gy / Control)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 48h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_48h.pdf",
  plot = Protein_p_top10_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_48h.jpg",
  plot = Protein_p_top10_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Volcano Plot Control - Plotting
#--------------------------------------------
Protein_p_volcano_Con_48h <- plot_volcano(
  data = ready_pCon_48h,
  title = "Proteome, Control, 48 h",
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
# 48h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Volcano_Con_48h.pdf",
  plot = Protein_p_volcano_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Volcano_Con_48h.jpg",
  plot = Protein_p_volcano_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Top 10 Control - filtering
#--------------------------------------------
Protein_top10_Con_48h_down <- ready_pCon_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Con_48h_up <- ready_pCon_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Con_48h <- bind_rows(
  Protein_top10_Con_48h_down,
  Protein_top10_Con_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Con_48h <- Protein_top10_Con_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 48h | Top 10 Control - plotting
#--------------------------------------------
Protein_p_top10_Con_48h <- ggplot(
  Protein_top10_Con_48h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Control, 48 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 48h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_Con_48h.pdf",
  plot = Protein_p_top10_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_Con_48h.jpg",
  plot = Protein_p_top10_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Volcano Plot Radiation - Plotting
#--------------------------------------------
Protein_p_volcano_Rad_48h <- plot_volcano(
  data = ready_pRAD_48h,
  title = "Proteom, Radiation, 48 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 10),
  extra_genes = c("CRIP1","MGST1","OTULINL")
)
#--------------------------------------------
# 48h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Proteome_Volcano_RAD_48h.pdf",
  plot = Protein_p_volcano_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_48h.jpg",
  plot = Protein_p_volcano_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Top 10 Radiation - filtering
#--------------------------------------------
Protein_top10_Rad_48h_down <- ready_pRAD_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

Protein_top10_Rad_48h_up <- ready_pRAD_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

Protein_top10_Rad_48h <- bind_rows(
  Protein_top10_Rad_48h_down,
  Protein_top10_Rad_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

Protein_top10_Rad_48h <- Protein_top10_Rad_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[3]],
      levels = .[[3]]
    )
  )
#--------------------------------------------
# 48h | Top 10 Radiation - plotting
#--------------------------------------------
Protein_p_top10_Rad_48h <- ggplot(
  Protein_top10_Rad_48h,
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
   size =2.5,
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
    breaks = seq(-10, 10, 10),
    minor_breaks = seq(-10, 10,, 5),
    expand = expansion(mult = c(0, 0))
  )+
  coord_cartesian(xlim = c(-12.5, 15))+
  labs(
    title = "Proteome differences \n Radiation, 48 h",
    x = expression(log[2] * "(HCC44 / NCI-H23)"),
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold"
    ),
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = c(0.825,0.25),
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.text = element_text(
      size = 7
    ),
    #legend.justification = c(0.5,0.5),
    plot.margin = margin(
      0, 5, 0, 5
    )
  )+
  guides(
    fill = guide_colourbar(
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(15, "mm")
    )
  )
#--------------------------------------------
# 48h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "20_Proteome_Top10_RAD_48h.pdf",
  plot = Protein_p_top10_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "20_Proteome_Top10_RAD_48h.jpg",
  plot = Protein_p_top10_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#-----------------------------------------------
setwd(dirs$code)
