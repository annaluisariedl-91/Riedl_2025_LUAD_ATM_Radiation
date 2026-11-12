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
RNA_p_volcano_HCC44 <- plot_volcano_RNA(
  data = ready_RHCC44_24h,
  title = "Transcriptome \nHCC44, 24 h",
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
  filename = "18_Transcriptome_Volcano_HCC44.pdf",
  plot = RNA_p_volcano_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Volcano_HCC44.jpg",
  plot = RNA_p_volcano_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 HCC44 - filtering
#--------------------------------------------
RNA_top10_HCC44_down <- ready_RHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

RNA_top10_HCC44_up <- ready_RHCC44_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

RNA_top10_HCC44 <- bind_rows(
  RNA_top10_HCC44_down,
  RNA_top10_HCC44_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

RNA_top10_HCC44 <- RNA_top10_HCC44 %>%
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
RNA_p_top10_HCC44 <- ggplot(
  RNA_top10_HCC44,
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
  filename = "18_Transcriptome_Top10_HCC44.pdf",
  plot = RNA_p_top10_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Top10_HCC44.jpg",
  plot = RNA_p_top10_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano H23 - Plotting
#--------------------------------------------
RNA_p_volcano_H23 <- plot_volcano_RNA(
  data = ready_RH23_24h,
  title = "Transcriptome \nNCI-H23, 24 h",
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
  filename = "18_Transcriptome_Volcano_H23.pdf",
  plot = RNA_p_volcano_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Volcano_H23.jpg",
  plot = RNA_p_volcano_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

#--------------------------------------------
# 24h | Top 10 H23 - filtering
#--------------------------------------------
RNA_top10_H23_down <- ready_RH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

RNA_top10_H23_up <- ready_RH23_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

RNA_top10_H23 <- bind_rows(
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

RNA_top10_H23 <- RNA_top10_H23 %>%
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
RNA_p_top10_H23 <- ggplot(
  RNA_top10_H23,
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
  filename = "18_Transcriptome_Top10_H23.pdf",
  plot = RNA_p_top10_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Top10_H23.jpg",
  plot = RNA_p_top10_H23,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Volcano Plot Control - Plotting
#--------------------------------------------
RNA_p_volcano_Con <- plot_volcano_RNA(
  data = ready_RCon_24h,
  title = "Transcriptome \nControl, 24 h",
  x_title = expression(
    log[2] * "(HCC44 / NCI-H23)"
  ),
  FC_cutoff = 1,
  p_cutoff = 0.05,
  n_labels = 2,
  xlim = c(-15, 15),
  ylim = c(0, 10),
  extra_genes=c( "SPOCK2","CCND2","TNS4","KRT81"))
#--------------------------------------------
# 24h | Volcano Plot Control - Saving
#--------------------------------------------
ggsave(
  filename = "18_Transcriptome_Volcano_Con.pdf",
  plot = RNA_p_volcano_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Volcano_Con.jpg",
  plot = RNA_p_volcano_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Control - filtering
#--------------------------------------------
RNA_top10_Con_down <- ready_RCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

RNA_top10_Con_up <- ready_RCon_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

RNA_top10_Con <- bind_rows(
  RNA_top10_Con_down,
  RNA_top10_Con_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

RNA_top10_Con <- RNA_top10_Con %>%
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
RNA_p_top10_Con <- ggplot(
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
    title = "Transcriptome differences \nControl, 24 h",
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
  filename = "18_Transcriptome_Top10_Con.pdf",
  plot = RNA_p_top10_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Top10_Con.jpg",
  plot = RNA_p_top10_Con,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Radiation - filtering
#--------------------------------------------
RNA_top10_Rad_down <- ready_RRad_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

RNA_top10_Rad_up <- ready_RRad_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

RNA_top10_Rad <- bind_rows(
  RNA_top10_Rad_down,
  RNA_top10_Rad_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

RNA_top10_Rad <- RNA_top10_Rad %>%
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
RNA_p_volcano_Rad <- plot_volcano_RNA(
  data = ready_RRad_24h,
  title = "Transcriptome \nRadiation, 24 h",
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
  filename = "18_Transcriptome_Volcano_Rad.pdf",
  plot = RNA_p_volcano_Rad,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Volcano_Rad.jpg",
  plot = RNA_p_volcano_Rad,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#--------------------------------------------
# 24h | Top 10 Radiatiom - filtering
#--------------------------------------------
RNA_top10_Rad_down <- ready_RRad_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange < -1
  ) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

RNA_top10_Rad_up <- ready_RRad_24h %>%
  filter(
    pval < 0.05,
    log2FoldChange > 1
  ) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10)

RNA_top10_Rad <- bind_rows(
  RNA_top10_Rad_down,
  RNA_top10_Rad_up
) %>%
  mutate(
    neg_log10_p = -log10(pval),
    Regulation = if_else(
      log2FoldChange > 0,
      "Upregulated",
      "Downregulated"
    )
  )

RNA_top10_Rad <- RNA_top10_Rad %>%
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
RNA_p_top10_Rad <- ggplot(
  top10_Rad,
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
  filename = "18_Transcriptome_Top10_Rad.pdf",
  plot = RNA_p_top10_Rad,
  width = 21*1/3,
  height = 7,
  units = "cm"
)

ggsave(
  filename = "18_Transcriptome_Top10_Rad.jpg",
  plot = RNA_p_top10_Rad,
  width = 21*1/3,
  height = 7,
  units = "cm"
)
#-----------------------------------------------
setwd(dirs$code)
