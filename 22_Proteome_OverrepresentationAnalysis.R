#--------------------------------------------
# loading spline regression results and proteome data
#--------------------------------------------
setwd(file.path(dirs$results, "21_Proteome_SplineRegression"))

HCC44_clusters <- read.csv(
  "HCC44_Splineregression_Clusters_label.csv")%>%
  mutate(
    ProteinKey=as.character(ProteinKey)
  )

HCC44_cluster3 <- HCC44_clusters %>%
    filter(HCC44_clusters$Cluster == 5)  

H23_clusters <- read.csv(
  "NCI_H23_Splineregression_Clusters_label.csv")%>%
  mutate(
    ProteinKey=as.character(ProteinKey)
  )

H23_cluster2 <- H23_clusters %>%
  filter(H23_clusters$Cluster == 2) 

load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)
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
#--------------------------------------------
# Binding data
#--------------------------------------------
ORA_HCC44_timeData_all <- make_timeData_all("HCC44")
colnames(ORA_HCC44_timeData_all)<-c("ENTREZID","SYMBOL","ProteinKey","ENSEMBL",column_names[41:80])
ORA_HCC44_timeData_all <- ORA_HCC44_timeData_all %>%
  mutate(
    ProteinKey=as.character(ProteinKey)
  )

ORA_H23_timeData_all <- make_timeData_all("H23")
colnames(ORA_H23_timeData_all)<-c("ENTREZID","SYMBOL","ProteinKey","ENSEMBL",column_names[41:80])
ORA_H23_timeData_all <- ORA_H23_timeData_all %>%
  mutate(
    ProteinKey=as.character(ProteinKey)
  )

setwd(file.path(dirs$results, "15_Proteome_DifferentialExpression"))
#--------------------------------------------
# HCC44 | ORA , names, ENSEMBL
#--------------------------------------------
ORA_GO_HCC44_results <- make_GO_ORA(
  HCC44_clusters  %>% dplyr::select(ProteinKey),
  ORA_HCC44_timeData_all %>% dplyr::select(ProteinKey, SYMBOL),
  top_n = 10
)

ORA_GO_res_pHCC44_Cluster <- ORA_GO_HCC44_results$enrichment
ORA_GO_pHCC44_Cluster_plot <- ORA_GO_HCC44_results$plot_data
ORA_GO_pHCC44_Cluster_plot <- ORA_GO_pHCC44_Cluster_plot %>%
  dplyr::arrange(dplyr::desc(gene_ratio_numeric)) %>%
  dplyr::mutate(
    GO_position = dplyr::row_number()
  )
#--------------------------------------------
# HCC44 | GO enrichment plot
#--------------------------------------------
ORA_GO_p_HCC44 <- ggplot(
  ORA_GO_pHCC44_Cluster_plot,
  aes(
    x = rev(gene_ratio_numeric),
    y = (GO_position)
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
    breaks = c(2,4,6,8,10),
    limits = c(2, 10),
    direction = 1,
    name =expression(-log[10](p[adj]))
  ) +
  scale_size(
    range = c(3, 7),
    name = "# Proteins",
    guide = guide_legend(
      position = "inside"
    )
  )  +
  scale_x_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    expand = expansion(mult = c(0, 0.02))
  ) +  
  coord_cartesian(
    xlim = c(0, 0.25),
    clip = "off"
  )+
  scale_y_continuous(
    breaks = ORA_GO_pHCC44_Cluster_plot$GO_position,
    labels = rev(ORA_GO_pHCC44_Cluster_plot$Description),
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
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black",
      hjust = 1
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    legend.text = element_text(
      size = 8
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
    fill = guide_colourbar(
      order = 2,
      barwidth = grid::unit(3, "mm"),
      barheight = grid::unit(12.5, "mm")
    ),
    size = guide_legend(
      order = 1,
      title = "# Proteins",
      override.aes = list(
        fill = "white",
        colour = "black",
        stroke = 0.3,
        alpha = 1
      ),
      keyheight = grid::unit(2, "mm"),
      keywidth = grid::unit(3, "mm")
    )
  )+
  theme(
    legend.position.inside = c(1.1, 0.8)
  )
#--------------------------------------------
# NCI-H23 | ORA
#--------------------------------------------
ORA_GO_H23_results <- make_GO_ORA(
  H23_clusters  %>%  dplyr::select(ProteinKey),
  ORA_H23_timeData_all %>% dplyr::select( SYMBOL,ProteinKey),
  top_n = 8
)

ORA_GO_res_pH23_Cluster <- ORA_GO_H23_results$enrichment
ORA_GO_pH23_Cluster_plot <- ORA_GO_H23_results$plot_data
ORA_GO_pH23_Cluster_plot <- ORA_GO_pH23_Cluster_plot %>%
  dplyr::arrange(dplyr::desc(gene_ratio_numeric)) %>%
  dplyr::mutate(
    GO_position = dplyr::row_number()
  )
#--------------------------------------------
# NCI-H23 | GO enrichment plot
#--------------------------------------------
ORA_GO_p_H23 <- ggplot(
  ORA_GO_pH23_Cluster_plot,
  aes(
    x = (gene_ratio_numeric),
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
    breaks = c(5,10,15),
    limits = c(1, 15),
    direction = 1,
    name = expression(-log[10](p[adj]))
  ) +
  scale_size(
    range = c(3, 7),
    name = "# Proteins",
    guide = guide_legend(
      position = "inside"
    )
  )  +
  scale_x_continuous(
    limits = c(0, 0.25),
    breaks = seq(0, 0.25, by = 0.05),
    labels = seq(0, 0.25, by = 0.05),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(
    xlim = c(0, 0.25)
  )+
  scale_y_continuous(
    breaks = (ORA_GO_pH23_Cluster_plot$GO_position),
    labels = rev(ORA_GO_pH23_Cluster_plot$Description),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title = "NCI-H23",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_prism_xcell() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black",
      hjust = 1
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    legend.text = element_text(
      size = 8
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
    fill = guide_colourbar(
      order = 2,
      barwidth = grid::unit(3, "mm"),
      barheight = grid::unit(20, "mm")
    ),
    size = guide_legend(
      order = 1,
      title = "# Proteins",
      override.aes = list(
        fill = "white",
        colour = "black",
        stroke = 0.3,
        alpha = 1
      ),
      keyheight = grid::unit(2, "mm"),
      keywidth = grid::unit(3, "mm")
    )
  )+
  theme(
    legend.position.inside = c(1.1, 0.8)
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
  "22_Proteome_GO_ORA_HCC44.jpeg",
  plot = ORA_GO_p_HCC44,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "22_Proteome_GO_ORA_HCC44.pdf",
  plot = ORA_GO_p_HCC44,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)  

write.csv(
  ORA_GO_pHCC44_Cluster_plot,
  file = file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis", "go_pHCC44_Cluster_Top10.csv"),
  row.names = FALSE
)

#2 Saving ORA NCI-H23 plots and top10 paths
ggsave(
  "22_Proteome_GO_ORA_H23.jpeg",
  plot = ORA_GO_p_H23,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "22_Proteome_GO_ORA_H23.pdf",
  plot = ORA_GO_p_H23,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)  

write.csv(
  ORA_GO_pH23_cluster_plot,
  file = file.path(dirs$results, "22_Proteome_OverrepresentationAnalysis", "go_pH23_Cluster_top10.csv"),
  row.names = FALSE
)
#--------------------------------------------
setwd(dirs$code)
