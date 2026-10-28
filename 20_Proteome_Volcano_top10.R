#-------------------------------------------
# loading or renaming data
#------------------------------------------
load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)
---------------------------------------------
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
p_volcano_p_HCC44_0h <- plot_volcano(
  data = ready_pHCC44_0h,
  title = "Proteome, HCC44, 0 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-10, 10),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 0h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_HCC44_0h.pdf",
  plot = p_volcano_p_HCC44_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_HCC44_0h.jpg",
  plot = p_volcano_p_HCC44_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 0h | Top 10 HCC44 - filtering
#--------------------------------------------
top10_p_HCC44_0h_down <- ready_pHCC44_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_HCC44_0h_up <- ready_pHCC44_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_HCC44_0h <- bind_rows(
  top10_p_HCC44_0h_down,
  top10_p_HCC44_0h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_HCC44_0h <- top10_p_HCC44_0h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 0h | Top 10 HCC44 - plotting
#--------------------------------------------
p_top10_p_HCC44_0h <- ggplot(
  top10_p_HCC44_0h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, HCC44, 0 h",
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
# 0h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_HCC44_0h.pdf",
  plot = p_top10_p_HCC44_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_HCC44_0h.jpg",
  plot = p_top10_p_HCC44_0h,
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
  xlim = c(-7.5, 7.5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 0h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_H23_0h.pdf",
  plot = p_volcanop_p_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_H23_0h.jpg",
  plot = p_volcanop_p_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 0h | Top 10 H23 - filtering
#--------------------------------------------
top10_p_H23_0h_down <- ready_pH23_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_H23_0h_up <- ready_pH23_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_H23_0h <- bind_rows(
  top10_p_H23_0h_down,
  top10_p_H23_0h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_H23_0h <- top10_p_H23_0h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 0h | Top 10 H23 - plotting
#--------------------------------------------
p_top10_p_H23_0h <- ggplot(
  top10_p_H23_0h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, NCI-H23, 0 h",
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
# 0h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_0h.pdf",
  plot = p_top10_p_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_0h.jpg",
  plot = p_top10_p_H23_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 0h | Volcano Plot Control - Plotting
#--------------------------------------------
p_volcano_p_Con_0h <- plot_volcano(
  data = ready_pCon_0h,
  title = "Proteome, Control, 0 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-12.5, 12.5),
  ylim = c(0, 17.5)
)
#--------------------------------------------
# 0h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_Con_0h.pdf",
  plot = p_volcano_p_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_Con_0h.jpg",
  plot = p_volcano_p_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 0h | Top 10 Control - filtering
#--------------------------------------------
top10_p_Con_0h_down <- ready_pCon_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Con_0h_up <- ready_pCon_0h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Con_0h <- bind_rows(
  top10_p_Con_0h_down,
  top10_p_Con_0h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Con_0h <- top10_p_Con_0h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 0h | Top 10 Control - plotting
#--------------------------------------------
p_top10_p_Con_0h <- ggplot(
  top10_p_Con_0h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Control, 0 h",
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
# 0h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_Con_0h.pdf",
  plot = p_top10_p_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_Con_0h.jpg",
  plot = p_top10_p_Con_0h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Volcano HCC44 - Plotting
#--------------------------------------------
p_volcano_p_HCC44_2h <- plot_volcano(
  data = ready_pHCC44_2h,
  title = "Proteome, HCC44, 2 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-7.5, 7.5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 2h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_HCC44_2h.pdf",
  plot = p_volcano_p_HCC44_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_HCC44_2h.jpg",
  plot = p_volcano_p_HCC44_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 2h | Top 10 HCC44 - filtering
#--------------------------------------------
top10_p_HCC44_2h_down <- ready_pHCC44_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_HCC44_2h_up <- ready_pHCC44_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_HCC44_2h <- bind_rows(
  top10_p_HCC44_2h_down,
  top10_p_HCC44_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_HCC44_2h <- top10_p_HCC44_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 2h | Top 10 HCC44 - plotting
#--------------------------------------------
p_top10_p_HCC44_2h <- ggplot(
  top10_p_HCC44_2h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, HCC44, 2 h",
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
# 2h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_HCC44_2h.pdf",
  plot = p_top10_p_HCC44_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_HCC44_2h.jpg",
  plot = p_top10_p_HCC44_2h,
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
  xlim = c(-5, 5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 2h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_H23_2h.pdf",
  plot = p_volcanop_p_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_H23_2h.jpg",
  plot = p_volcanop_p_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Top 10 H23 - filtering
#--------------------------------------------
top10_p_H23_2h_down <- ready_pH23_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_H23_2h_up <- ready_pH23_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_H23_2h <- bind_rows(
  top10_p_H23_2h_down,
  top10_p_H23_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_H23_2h <- top10_p_H23_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 2h | Top 10 H23 - plotting
#--------------------------------------------
p_top10_p_H23_2h <- ggplot(
  top10_p_H23_2h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, NCI-H23, 2 h",
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
# 2h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_2h.pdf",
  plot = p_top10_p_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_2h.jpg",
  plot = p_top10_p_H23_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Volcano Plot Control - Plotting
#--------------------------------------------
p_volcano_p_Con_2h <- plot_volcano(
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
  filename = "Proteom_Volcano_Con_2h.pdf",
  plot = p_volcano_p_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_Con_2h.jpg",
  plot = p_volcano_p_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Top 10 Control - filtering
#--------------------------------------------
top10_p_Con_2h_down <- ready_pCon_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Con_2h_up <- ready_pCon_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Con_2h <- bind_rows(
  top10_p_Con_2h_down,
  top10_p_Con_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Con_2h <- top10_p_Con_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 2h | Top 10 Control - plotting
#--------------------------------------------
p_top10_p_Con_2h <- ggplot(
  top10_p_Con_2h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Control, 2 h",
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
# 2h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_Con_2h.pdf",
  plot = p_top10_p_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_Con_2h.jpg",
  plot = p_top10_p_Con_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Volcano Plot Radiation - Plotting
#--------------------------------------------
p_volcano_p_Rad_2h <- plot_volcano(
  data = ready_pRAD_2h,
  title = "Proteome, Radiation, 2 h",
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
# 2h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Proteome_Volcano_RAD_2h.pdf",
  plot = p_volcano_p_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_2h.jpg",
  plot = p_volcano_p_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 2h | Top 10 Radiation - filtering
#--------------------------------------------
top10_p_RAD_2h_down <- ready_pRAD_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_RAD_2h_up <- ready_pRAD_2h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_RAD_2h <- bind_rows(
  top10_p_RAD_2h_down,
  top10_p_RAD_2h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_RAD_2h <- top10_p_RAD_2h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 2h | Top 10 Radiation - plotting
#--------------------------------------------
p_top10_p_Rad_2h <- ggplot(
  top10_p_RAD_2h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Radiation, 2 h",
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
# 2h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_RAD_2h.pdf",
  plot = p_top10_p_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_RAD_2h.jpg",
  plot = p_top10_p_Rad_2h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Volcano HCC44 - Plotting
#--------------------------------------------
p_volcano_p_HCC44_6h <- plot_volcano(
  data = ready_pHCC44_6h,
  title = "Proteome, HCC44, 6 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-7.5, 7.5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 6h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_HCC44_6h.pdf",
  plot = p_volcano_p_HCC44_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_HCC44_6h.jpg",
  plot = p_volcano_p_HCC44_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 6h | Top 10 HCC44 - filtering
#--------------------------------------------
top10_p_HCC44_6h_down <- ready_pHCC44_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_HCC44_6h_up <- ready_pHCC44_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_HCC44_6h <- bind_rows(
  top10_p_HCC44_6h_down,
  top10_p_HCC44_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_HCC44_6h <- top10_p_HCC44_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 6h | Top 10 HCC44 - plotting
#--------------------------------------------
p_top10_p_HCC44_6h <- ggplot(
  top10_p_HCC44_6h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, HCC44, 6 h",
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
# 6h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_HCC44_6h.pdf",
  plot = p_top10_p_HCC44_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_HCC44_6h.jpg",
  plot = p_top10_p_HCC44_6h,
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
  xlim = c(-5, 5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 6h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_H23_6h.pdf",
  plot = p_volcanop_p_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_H23_6h.jpg",
  plot = p_volcanop_p_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 6h | Top 10 H23 - filtering
#--------------------------------------------
top10_p_H23_6h_down <- ready_pH23_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_H23_6h_up <- ready_pH23_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_H23_6h <- bind_rows(
  top10_p_H23_6h_down,
  top10_p_H23_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_H23_6h <- top10_p_H23_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 6h | Top 10 H23 - plotting
#--------------------------------------------
p_top10_p_H23_6h <- ggplot(
  top10_p_H23_6h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, NCI-H23, 6 h",
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
# 6h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_6h.pdf",
  plot = p_top10_p_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_6h.jpg",
  plot = p_top10_p_H23_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Volcano Plot Control - Plotting
#--------------------------------------------
p_volcano_p_Con_6h <- plot_volcano(
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
  filename = "Proteom_Volcano_Con_6h.pdf",
  plot = p_volcano_p_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_Con_6h.jpg",
  plot = p_volcano_p_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Top 10 Control - filtering
#--------------------------------------------
top10_p_Con_6h_down <- ready_pCon_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Con_6h_up <- ready_pCon_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Con_6h <- bind_rows(
  top10_p_Con_6h_down,
  top10_p_Con_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Con_6h <- top10_p_Con_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 6h | Top 10 Control - plotting
#--------------------------------------------
p_top10_p_Con_6h <- ggplot(
  top10_p_Con_6h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Control, 6 h",
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
# 6h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_Con_6h.pdf",
  plot = p_top10_p_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_Con_6h.jpg",
  plot = p_top10_p_Con_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Volcano Plot Radiation - Plotting
#--------------------------------------------
p_volcano_p_Rad_6h <- plot_volcano(
  data = ready_pRAD_6h,
  title = "Proteome, Radiation, 6 h",
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
# 6h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Proteome_Volcano_RAD_6h.pdf",
  plot = p_volcano_p_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_6h.jpg",
  plot = p_volcano_p_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 6h | Top 10 Radiation - filtering
#--------------------------------------------
top10_p_Rad_6h_down <- ready_pRAD_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Rad_6h_up <- ready_pRAD_6h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Rad_6h <- bind_rows(
  top10_p_Rad_6h_down,
  top10_p_Rad_6h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Rad_6h <- top10_p_Rad_6h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 6h | Top 10 Radiation - plotting
#--------------------------------------------
p_top10_p_Rad_6h <- ggplot(
  top10_p_Rad_6h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Radiation, 6 h",
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
# 6h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_RAD_6h.pdf",
  plot = p_top10_p_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_RAD_6h.jpg",
  plot = p_top10_p_Rad_6h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano HCC44 - Plotting
#--------------------------------------------
p_volcano_p_HCC44_24h <- plot_volcano(
  data = ready_pHCC44_24h,
  title = "Proteome, HCC44, 24 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-10, 10),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 24h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_HCC44_24h.pdf",
  plot = p_volcano_p_HCC44_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_HCC44_24h.jpg",
  plot = p_volcano_p_HCC44_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 HCC44 - filtering
#--------------------------------------------
top10_p_HCC44_24h_down <- ready_pHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_HCC44_24h_up <- ready_pHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_HCC44_24h <- bind_rows(
  top10_p_HCC44_24h_down,
  top10_p_HCC44_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_HCC44_24h <- top10_p_HCC44_24h %>%
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

top10_p_HCC44_24h <- top10_p_HCC44_24h %>%
  arrange(log2FoldChange) %>%
  mutate(
    SYMBOL = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 24h | Top 10 HCC44 - plotting
#--------------------------------------------
p_top10_p_HCC44_24h <- ggplot(
  top10_p_HCC44_24h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, HCC44, 24 h",
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
  filename = "Proteom_Top10_HCC44_24h.pdf",
  plot = p_top10_p_HCC44_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_HCC44_24h.jpg",
  plot = p_top10_p_HCC44_24h,
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
  xlim = c(-5, 5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 24h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_H23_24h.pdf",
  plot = p_volcanop_p_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_H23_24h.jpg",
  plot = p_volcanop_p_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 H23 - filtering
#--------------------------------------------
top10_p_H23_24h_down <- ready_pH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_H23_24h_up <- ready_pH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_H23_24h <- bind_rows(
  top10_p_H23_24h_down,
  top10_p_H23_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_H23_24h <- top10_p_H23_24h %>%
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
p_top10_p_H23_24h <- ggplot(
  top10_p_H23_24h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, NCI-H23, 24 h",
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
  filename = "Protein_Top10_H23_24h.pdf",
  plot = p_top10_p_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_24h.jpg",
  plot = p_top10_p_H23_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano Plot Control - Plotting
#--------------------------------------------
p_volcano_p_Con_24h <- plot_volcano(
  data = ready_pCon_24h,
  title = "Proteome, Control, 24 h",
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
  filename = "Proteom_Volcano_Con_24h.pdf",
  plot = p_volcano_p_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_Con_24h.jpg",
  plot = p_volcano_p_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Control - filtering
#--------------------------------------------
top10_p_Con_24h_down <- ready_pCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Con_24h_up <- ready_pCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Con_24h <- bind_rows(
  top10_p_Con_24h_down,
  top10_p_Con_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Con_24h <- top10_p_Con_24h %>%
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
p_top10_p_Con_24h <- ggplot(
  top10_p_Con_24h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Control, 24 h",
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
  filename = "Proteom_Top10_Con_24h.pdf",
  plot = p_top10_p_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_Con_24h.jpg",
  plot = p_top10_p_Con_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano Plot Radiation - Plotting
#--------------------------------------------
p_volcano_p_Rad_24h <- plot_volcano(
  data = ready_pRAD_24h,
  title = "Proteome, Radiation, 24 h",
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
  filename = "Proteome_Volcano_RAD_24h.pdf",
  plot = p_volcano_p_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_24h.jpg",
  plot = p_volcano_p_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Radiation - filtering
#--------------------------------------------
top10_p_Rad_24h_down <- ready_pRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Rad_24h_up <- ready_pRAD_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Rad_24h <- bind_rows(
  top10_p_Rad_24h_down,
  top10_p_Rad_24h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Rad_24h <- top10_p_Rad_24h %>%
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
p_top10_p_Rad_24h <- ggplot(
  top10_p_Rad_24h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Radiation, 24 h",
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
  filename = "Proteom_Top10_RAD_24h.pdf",
  plot = p_top10_p_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_RAD_24h.jpg",
  plot = p_top10_p_Rad_24h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Volcano HCC44 - Plotting
#--------------------------------------------
p_volcano_p_HCC44_48h <- plot_volcano(
  data = ready_pHCC44_48h,
  title = "Proteome, HCC44, 48 h",
  x_title = expression(
    log[2] * "(2 Gy / Control)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-7.5, 7.5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 48h | Volcano HCC44 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_HCC44_48h.pdf",
  plot = p_volcano_p_HCC44_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_HCC44_48h.jpg",
  plot = p_volcano_p_HCC44_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 48h | Top 10 HCC44 - filtering
#--------------------------------------------
top10_p_HCC44_48h_down <- ready_pHCC44_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_HCC44_48h_up <- ready_pHCC44_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_HCC44_48h <- bind_rows(
  top10_p_HCC44_48h_down,
  top10_p_HCC44_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_HCC44_48h <- top10_p_HCC44_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 48h | Top 10 HCC44 - plotting
#--------------------------------------------
p_top10_p_HCC44_48h <- ggplot(
  top10_p_HCC44_48h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, HCC44, 48 h",
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
# 48h | Top 10 HCC44 - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_HCC44_48h.pdf",
  plot = p_top10_p_HCC44_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_HCC44_48h.jpg",
  plot = p_top10_p_HCC44_48h,
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
  xlim = c(-5, 5),
  ylim = c(0, 7.5)
)
#--------------------------------------------
# 48h | Volcano H23 - Saving
#--------------------------------------------
ggsave(
  filename = "Proteom_Volcano_H23_48h.pdf",
  plot = p_volcanop_p_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_H23_48h.jpg",
  plot = p_volcanop_p_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 48h | Top 10 H23 - filtering
#--------------------------------------------
top10_p_H23_48h_down <- ready_pH23_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_H23_48h_up <- ready_pH23_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_H23_48h <- bind_rows(
  top10_p_H23_48h_down,
  top10_p_H23_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_H23_48h <- top10_p_H23_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )

top10_p_H23_48h <- top10_p_H23_48h %>%
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

top10_p_H23_48h <- top10_p_H23_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    SYMBOL = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 48h | Top 10 H23 - plotting
#--------------------------------------------
p_top10_p_H23_48h <- ggplot(
  top10_p_H23_48h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, NCI-H23, 48 h",
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
# 48h | Top 10 H23 - saving results
#--------------------------------------------
ggsave(
  filename = "Protein_Top10_H23_48h.pdf",
  plot = p_top10_p_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Protein_Top10_H23_48h.jpg",
  plot = p_top10_p_H23_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Volcano Plot Control - Plotting
#--------------------------------------------
p_volcano_p_Con_48h <- plot_volcano(
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
  filename = "Proteom_Volcano_Con_48h.pdf",
  plot = p_volcano_p_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Volcano_Con_48h.jpg",
  plot = p_volcano_p_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Top 10 Control - filtering
#--------------------------------------------
top10_p_Con_48h_down <- ready_pCon_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Con_48h_up <- ready_pCon_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Con_48h <- bind_rows(
  top10_p_Con_48h_down,
  top10_p_Con_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Con_48h <- top10_p_Con_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 48h | Top 10 Control - plotting
#--------------------------------------------
p_top10_p_Con_48h <- ggplot(
  top10_p_Con_48h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Control, 48 h",
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
# 48h | Top 10 Control - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_Con_48h.pdf",
  plot = p_top10_p_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_Con_48h.jpg",
  plot = p_top10_p_Con_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Volcano Plot Radiation - Plotting
#--------------------------------------------
p_volcano_p_Rad_48h <- plot_volcano(
  data = ready_pRAD_48h,
  title = "Proteom, Radiation, 48 h",
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
# 48h | Volcano Plot Radiation - Saving
#--------------------------------------------
ggsave(
  filename = "Proteome_Volcano_RAD_48h.pdf",
  plot = p_volcano_p_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteome_Volcano_RAD_48h.jpg",
  plot = p_volcano_p_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 48h | Top 10 Radiation - filtering
#--------------------------------------------
top10_p_Rad_48h_down <- ready_pRAD_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

top10_p_Rad_48h_up <- ready_pRAD_48h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

top10_p_Rad_48h <- bind_rows(
  top10_p_Rad_48h_down,
  top10_p_Rad_48h_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

top10_p_Rad_48h <- top10_p_Rad_48h %>%
  arrange(log2FoldChange) %>%
  mutate(
    Gene = factor(
      .[[1]],
      levels = .[[1]]
    )
  )
#--------------------------------------------
# 48h | Top 10 Radiation - plotting
#--------------------------------------------
p_top10_p_Rad_48h <- ggplot(
  top10_p_Rad_48h,
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
    breaks = scales::breaks_pretty(n = 5),
    minor_breaks = scales::minor_breaks_n(1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Proteome differences, Radiation, 48 h",
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
# 48h | Top 10 Radiation - saving results
#--------------------------------------------
ggsave(
  filename = "Proteom_Top10_RAD_48h.pdf",
  plot = p_top10_p_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "Proteom_Top10_RAD_48h.jpg",
  plot = p_top10_p_Rad_48h,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#-----------------------------------------------
setwd(dirs$code)
