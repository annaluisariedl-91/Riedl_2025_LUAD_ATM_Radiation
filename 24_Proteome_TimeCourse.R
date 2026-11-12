#--------------------------------------------
# Loading data
#--------------------------------------------
setwd(file.path(dirs$results, "23_protein_ProteinProteinInteraction"))
HCC44_cluster3_network <- read.csv(
  "23_PIP_HCC44_Cluster3_STRING_network_statistics.csv")

setwd(file.path(dirs$results, "21_Proteome_SplineRegression"))
HCC44_labels_row <- read.csv(
  "HCC44_Splineregression_Clusters_label.csv")

#load_proteome_data(file.path(dirs$results, "15_Proteome_DifferentialExpression"))
Prot_timeData_mean <- cbind(
  make_meanData("H23", first = TRUE,file.path(dirs$results, "15_Proteome_DifferentialExpression")),
  make_meanData("HCC44",file.path(dirs$results, "15_Proteome_DifferentialExpression"))
)

times<-c(0,2,6,24,48)
Cellline<-c("NCI-H23","HCC44")
treatments<-c("c","2Gy")
sample_names <- unlist(lapply(Cellline, function(cellline) {
  unlist(lapply(times, function(time) {
    unlist(lapply(treatments, function(treatment) {
      paste0(cellline, "_", time, "_", treatment)
    }))
  }))
}))

TimeCourse_timeData_mean<-Prot_timeData_mean[,c(1:14,19:28)]
colnames(TimeCourse_timeData_mean)<-c("ProteinKey","ENSEMBL","ENTREZID","SYMBOL",sample_names)
#--------------------------------------------
# Analysis of ATM content time course
#--------------------------------------------
TimeCourse_ATM_data <- TimeCourse_timeData_mean %>%
  dplyr::filter(SYMBOL == "ATM")

TimeCourse_ATM_data_long <- TimeCourse_ATM_data %>%
  tidyr::pivot_longer(
    cols = -c(ProteinKey, ENSEMBL, ENTREZID, SYMBOL),
    names_to = c("Cellline", "Time", "Treatment"),
    names_sep = "_",
    values_to = "Intensity"
  ) %>%
  dplyr::mutate(
    Time = as.numeric(gsub("h", "", Time)),
    Treatment= ifelse(Treatment=="2Gy", "2 Gy",Treatment)
  )

TimeCourse_ATM_data_long <- TimeCourse_ATM_data_long %>%
  group_by(Cellline, Treatment) %>%
  mutate(
    rel_Intensity = Intensity / Intensity[Time == 0],
    group = paste(Cellline, Treatment)
  ) %>%
  ungroup()

TimeCourse_ATM_data_long<-TimeCourse_ATM_data_long %>%
  filter(
    Treatment == "2 Gy"
  )
#--------------------------------------------
# ATM content time course Plot
#--------------------------------------------
TimeCourse_p_ATM <- ggplot(
  TimeCourse_ATM_data_long,
  aes(
    x = Time,
    y = rel_Intensity,
    colour = group,
    group = group
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
    aes(fill = group),
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
    title = "Abundance - ATM",
    x = "Time [h]",
    y = expression(log[2]~"(Intensity), normalized")
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
    legend.position = "bottom",
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
    colour = "none",
    fill = guide_legend(title = NULL), 
    
  )
#--------------------------------------------
# Preparing Analysis of Cell cycle proteins
#--------------------------------------------
TimeCourse_cell_cycle_data <- data.frame(
  make_cluster_data(
    TimeCourse_timeData_mean
)) %>%
  mutate(
    Cellline = factor(
      Cellline,
      levels = c("HCC44", "NCI-H23")
    ),
    ProteinKey = as.character(ProteinKey)
  )

# excluding ATM and extracting cell cylce proteins
TimeCourse_cell_cycle_data <- TimeCourse_cell_cycle_data %>%
  dplyr::inner_join(
    HCC44_cluster3_network %>%
      mutate(SYMBOL = as.character(Gene))%>%
      dplyr::select(SYMBOL),
    by ="SYMBOL"
  ) %>%
  filter(
    Treatment=="2 Gy"
  )

TimeCourse_cell_cycle_data <- TimeCourse_cell_cycle_data %>%
  group_by(Cellline, Treatment, ProteinKey) %>%
  mutate(
    rel_Intensity =Intensity/Intensity[Time == 0],
    group = paste(Cellline, Treatment)
  ) %>%
  ungroup()
#--------------------------------------------
# Protein Content of Cluster 3 (Cell cylce)
#--------------------------------------------
TimeCourse_p_cell_cycle <- ggplot(
  TimeCourse_cell_cycle_data,
  aes(
    x = Time,
    y = rel_Intensity,
    colour = group,
    group = group
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
    ~ SYMBOL,
    ncol = 4
  ) +
  labs(
    x = "Time [h]",
    y = expression(log[2]~"(Intensity), normalized")
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
    legend.position =c(0.9,0.15),
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    legend.key.width = grid::unit(3, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    legend.spacing.x = grid::unit(1, "mm"),
    legend.text = element_text(
      size = 10
    ),
    legend.key = element_blank(),
    plot.margin = margin(
      t =10,
      r = 5,
      b = 5,
      l = 10
    )
  ) 
#--------------------------------------------
# Saving results
#--------------------------------------------
dir.create(
  file.path(dirs$results, "24_Proteome_TimeCourse"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "24_Proteome_TimeCourse"))

#1 Saving plot of ATM content depending on time
ggsave(
  "24_Proteome_TimeCourse_ATM.jpeg",
  TimeCourse_p_ATM,
  width = 21*1/2,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  "24_Proteome_TimeCourse_ATM.pdf",
  TimeCourse_p_ATM,
  width = 21*1/2,
  height = 7,
  units = "cm",
  dpi = 300
)    

#2 Saving plot of Cluster 3 protein content depending on time
ggsave(
  "24_Proteome_TimeCourse_Cluster3.jpeg",
  TimeCourse_p_cell_cycle,
  width = 21*3/3,
  height = 14,
  units = "cm",
  dpi = 300
)

ggsave(
  "24_Proteome_TimeCourse_Cluster3.pdf",
  TimeCourse_p_cell_cycle,
  width = 21*3/3,
  height = 14,
  units = "cm",
  dpi = 300
)   
#-------------------------------------------
setwd(dirs$code)
