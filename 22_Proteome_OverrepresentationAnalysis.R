#--------------------------------------------
# loading spline regression results
#--------------------------------------------
setwd(file.path(dirs$results, "21_Proteome_SplineRegression"))

HCC44_clusters <- read.csv(
  "HCC44_Splineregression_Clusters.csv")

HCC44_clusters_info <- read.csv(
    "HCC44_Splineregression_significantProteins.csv"
  ) %>%
    dplyr::select(1, 9:12)

HCC44_clusters <- HCC44_clusters_info %>%
    left_join(
      HCC44_clusters,
      by = "ProteinKey"
    )
HCC44_clusters<-HCC44_clusters[complete.cases(HCC44_clusters$ProteinKey),]
HCC44_cluster3 <- HCC44_clusters %>%
    filter(HCC44_clusters$Cluster == 3)  

H23_clusters <- read.csv(
  "NCI-H23_Splineregression_Clusters.csv")

H23_clusters_info <- read.csv(
  "NCI-H23_Splineregression_significantProteins.csv"
) %>%
  dplyr::select(1, 9:12)

H23_clusters <- H23_clusters_info %>%
  left_join(
    H23_clusters,
    by = "ProteinKey"
  )
H23_cluster2 <- H23_clusters %>%
  filter(H23_clusters$Cluster == 2)  
#--------------------------------------------
# setting path
#--------------------------------------------
setwd(file.path(dirs$results, "15_Proteome_DifferentialExpression"))

combinations <- expand.grid(
  replicate = paste0("Rep_", 1:4),
  treatment = c("c", "2Gy"),
  time = c("0h", "2h", "6h", "24h", "48h"),
  cellline = c("H23", "HCC44")
)

column_names <- with(
  combinations,
  paste(cellline, time, treatment, replicate, sep = "_")
)

combinations_mean <- expand.grid(
  treatment = c("c", "2Gy"),
  time = c("0h", "2h", "6h", "24h", "48h"),
  cellline = c("H23", "HCC44")
)

column_mean_names <- with(
  combinations_mean,
  paste("Mean",cellline, time, treatment, sep = "_")
)
#--------------------------------------------
# Binding data
#--------------------------------------------
HCC44_timeData <- make_timeData(
    ready_pHCC44_0h,
    ready_pHCC44_2h,
    ready_pHCC44_6h,
    ready_pHCC44_24h,
    ready_pHCC44_48h
  )

HCC44_timeData <- HCC44_timeData[
  complete.cases(HCC44_timeData[, 3:4]),
]

HCC44_timeData_all <- make_timeData_all("HCC44")
colnames(HCC44_timeData_all)<-c("ENTREZID","SYMBOL","ProteinKey","ENSEMBL",column_names[41:80])

HCC44_Cluster_data <- HCC44_timeData_all %>%
  inner_join(HCC44_clusters, by =c("ProteinKey"))

res_sig_HCC44_time <- results_HCC44_time %>%
  filter(adj.P.Val < 0.05)

H23_timeData <- make_timeData(
  ready_pH23_0h,
  ready_pH23_2h,
  ready_pH23_6h,
  ready_pH23_24h,
  ready_pH23_48h
)

H23_timeData_all <- make_timeData_all("H23")

res_sig_H23_time <- results_H23_time %>%
  filter(adj.P.Val < 0.05)
colnames(H23_timeData_all)<-c("ENTREZID","SYMBOL","ProteinKey","ENSEMBL",column_names[1:40])

H23_Cluster <- H23_timeData_all %>%
  inner_join(H23_clusters, by =c("ProteinKey"))

Data_all <- merge(
  H23_timeData_all,
  HCC44_timeData_all[, c(3, 5:44)],
  by = "ProteinKey"
)

setwd(file.path(dirs$results, "15_Proteome_DifferentialExpression"))
Prot_timeData_mean <- cbind(
    make_meanData("H23", first = TRUE),
    make_meanData("HCC44")
  )

colnames(Prot_timeData_mean)<-c("ProteinKey","ENSEMBL","ENTREZID","SYMBOL",column_mean_names)
#--------------------------------------------
# HCC44 | ORA
#--------------------------------------------
GO_HCC44_results <- make_GO_ORA(
  HCC44_clusters  %>% dplyr::select(ProteinKey, SYMBOL, ENSEMBL),
  HCC44_timeData %>% dplyr::select(ProteinKey, ENTREZID),
  top_n = 10
)

res_pHCC44_Cluster <- GO_HCC44_results$enrichment
go_pHCC44_Cluster_plot <- GO_HCC44_results$plot_data
#--------------------------------------------
# HCC44 | GO enrichment plot
#--------------------------------------------
GO_ORA_HCC44 <- ggplot(
  go_pHCC44_Cluster_plot,
  aes(
    x = gene_ratio_numeric,
    y = GO_position
  )
) +
  geom_point(
    aes(
      size = Count,
      fill = p.adjust_log
    ),
    shape = 21,
    colour = "black",
    stroke = 0.5,
    alpha = 0.85
  ) +
  scale_fill_distiller(
    palette = "Greens",
    breaks = c(3, 5, 7, 9),
    limits = c(2, 10),
    direction = 1,
    name = expression(-log[10](p))
  ) +
  scale_size(
    range = c(3, 8),
    name = "# Proteins"
  ) +
  scale_x_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = go_pHCC44_Cluster_plot$GO_position,
    labels = go_pHCC44_Cluster_plot$Description,
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title = "HCC44",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    plot.title = element_text(
      size = 9,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(
      size = 8,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 9,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 7,
      colour = "black",
      hjust = 1
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    legend.text = element_text(
      size = 6
    ),
    legend.key = element_blank(),
    plot.margin = margin(
      t =10,
      r = 5,
      b = 5,
      l = 10
    )
  ) +
  guides(
    colour = guide_colourbar(
      order = 1,
      barwidth  = grid::unit(3, "mm"),
      barheight = grid::unit(20, "mm")
    ),
    size = guide_legend(
      order = 2,
      title = "# Proteins"
    )
  ) 

#--------------------------------------------
# NCI-H23 | ORA
#--------------------------------------------
GO_H23_results <- make_GO_ORA(
  H23_clusters  %>%  dplyr::select(ProteinKey, SYMBOL, ENSEMBL),
  H23_timeData %>% dplyr::select(ProteinKey, ENTREZID),
  top_n = 8
)

res_pH23_cluster <- GO_H23_results$enrichment
go_pH23_cluster_plot <- GO_H23_results$plot_data
#--------------------------------------------
# NCI-H23 | GO enrichment plot
#--------------------------------------------
GO_ORA_H23 <- ggplot(
  go_pH23_cluster_plot,
  aes(
    x = gene_ratio_numeric,
    y = GO_position
  )
) +
  geom_point(
    aes(
      size = Count,
      fill = p.adjust_log
    ),
    shape = 21,
    colour = "black",
    stroke = 0.5,
    alpha = 0.85
  ) +
  scale_fill_distiller(
    palette = "Greens",
    breaks = c(1, 2, 4, 6),
    limits = c(1, 7),
    direction = 1,
    name = expression(-log[10](p[adj]))
  ) +
  scale_size(
    range = c(3, 8),
    name = "# Proteins"
  ) +
  scale_x_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = go_pH23_cluster_plot$GO_position,
    labels = go_pH23_cluster_plot$Description,
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title = "NCI-H23",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    # Titel
    plot.title = element_text(
      size = 9,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(
      size = 8,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 9,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 7,
      colour = "black",
      hjust = 1
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    legend.text = element_text(
      size = 6
    ),
    legend.key = element_blank(),
    plot.margin = margin(
      t =10,
      r = 5,
      b = 5,
      l = 10
    )
  ) +
  guides(
    colour = guide_colourbar(
      order = 1,
      barwidth  = grid::unit(3, "mm"),
      barheight = grid::unit(20, "mm")
    ),
    size = guide_legend(
      order = 2,
      title = "# Proteins"
    )
  ) 
#--------------------------------------------
# Analysis of ATM content time course
#--------------------------------------------
ATM_data <- make_protein_timecourse(
  Prot_timeData_mean,
  gene = "ATM"
)
#--------------------------------------------
# ATM content time course Plot
#--------------------------------------------
Cluster3_ATM <- ggplot(
  ATM_data,
  aes(
    x = Time,
    y = Intensity,
    colour = Cellline,
    group = Cellline
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey60"
  ) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    aes(fill = Cellline),
    shape = 21,
    size = 2.5,
    stroke = 0.5,
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(
      "HCC44 2 Gy" =  "#EA4025",
      "NCI-H23 2 Gy" = "#5F8F29"
    )
  ) +
  scale_fill_manual(
    values = c(
      "HCC44 2 Gy" =  "#EA4025",
      "NCI-H23 2 Gy" ="#5F8F29"
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 2, 6, 24, 48),
    limits = c(0, 48)
  ) +
  scale_y_continuous(
    breaks = seq(0.95, 1.05, by = 0.05),
    minor_breaks = seq(0.975, 1.025, by = 0.025),
    limits = c(0.95, 1.05)
  ) +
  labs(
    title = "ATM",
    x = "Time [h]",
    y = expression(log[2]~"(Intensity), normalized")
  ) +
  theme_prism_xcell() +
  theme(
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "bottom"
  ) +
  guides(
    colour = "none",
    fill = guide_legend(title = NULL)
  )

#--------------------------------------------
# Preparing Analysis of Cell cycle proteins
#--------------------------------------------
time_points <- c(0,2,6,24,48)

Cluster3_cell_cycle_data <- bind_rows(
  make_cluster_data(
    Prot_timeData_mean,
    cols = c(16, 18, 20, 22, 24),
    cellline = "HCC44 2 Gy"
  ),
  make_cluster_data(
    Prot_timeData_mean,
    cols = c(6, 8, 10, 12, 14),
    cellline = "NCI-H23 2 Gy"
  )
) %>%
  mutate(
    Cellline = factor(
      Cellline,
      levels = c("HCC44 2 Gy", "NCI-H23 2 Gy")
    )
  )

# excluding ATM and extracting cell cylce proteins
Cluster3_cell_cycle_data <- Cluster3_cell_cycle_data %>%
  dplyr::filter(
    !is.na(Protein),
    Protein != "",
    Protein != "ATM",
    Protein != "ANLN",
    Protein %in% HCC44_cluster3[[5]]
  )
#--------------------------------------------
# Protein Content of Cluster 3 (Cell cylce)
#--------------------------------------------
Cluster3_other <- ggplot(
  Cluster3_cell_cycle_data,
  aes(
    x = Time,
    y = Intensity,
    colour = Cellline,
    group = Cellline
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dotted",
    linewidth = 0.6,
    colour = "grey40"
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    aes(fill = Cellline),
    shape = 21,
    size = 2,
    stroke = 0.4,
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(
      "HCC44 2 Gy" = "#EA4025",
      "NCI-H23 2 Gy" = "#5F8F29"
    )
  ) +
  scale_fill_manual(
    values = c(
      "HCC44 2 Gy" = "#EA4025",
      "NCI-H23 2 Gy" = "#5F8F29"
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 2, 6, 24, 48),
    limits = c(0, 48)
  ) +
  facet_wrap(
    ~ Protein,
    ncol = 4
  ) +
  labs(
    x = "Time [h]",
    y = expression(log[2]~"(Intensity), normalized")
  ) +
  theme_prism_xcell() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 9,
      face = "bold",
      colour = "black"
    ),
    legend.position = "bottom"
  ) +
  guides(
    colour = "none",
    fill = guide_legend(
      title = NULL,
      nrow = 1
    )
  )    
#--------------------------------------------
# saving results
#--------------------------------------------
dir.create(
  file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis"))

#1 Saving ORA HCC44 plots and top10 paths
ggsave(
  "GO_ORA_HCC44.jpeg",
  plot = GO_ORA_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "GO_ORA_HCC44.pdf",
  plot = GO_ORA_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)  

write.csv(
  go_pHCC44_Cluster_plot,
  file = file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis", "go_pHCC44_Cluster_top10.csv"),
  row.names = FALSE
)

#2 Saving ORA NCI-H23 plots and top10 paths
ggsave(
  "GO_ORA_H23.jpeg",
  plot = GO_ORA_H23,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "GO_ORA_H23.pdf",
  plot = GO_ORA_H23,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)  

write.csv(
  go_pHCC44_Cluster_plot,
  file = file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis", "go_pH23_Cluster_top10.csv"),
  row.names = FALSE
)

#4 Saving plot of ATM content depending on time
ggsave(
  "Cluster3_ATM.jpeg",
  Cluster3_ATM,
  width = 21*1/2,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "Cluster3_ATM.pdf",
  Cluster3_ATM,
  width = 21*1/2,
  height = 7,
  units = "cm",
  dpi = 300
)    

write.csv(
  ATM_intensity_proteome,
  file = file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis", "ATM_intensity_proteome.csv"),
  row.names = FALSE
)

#5 Saving plot of Cluster 3 protein content depending on time
ggsave(
  "Cluster3_other.jpeg",
  Cluster3_other,
  width = 21*3/3,
  height = 14,
  units = "cm",
  dpi = 300
)

ggsave(
  "Cluster3_other.pdf",
  Cluster3_other,
  width = 21*3/3,
  height = 14,
  units = "cm",
  dpi = 300
)   

write.csv(
  ATM_intensity_proteome,
  file = file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis", "ATM_intensity_proteome.csv"),
  row.names = FALSE
)
#--------------------------------------------
setwd(dirs$code)
