#--------------------------------------------
# loading spline regression results
#--------------------------------------------
setwd(file.path(dirs$results, "21_Proteome_SplineRegression"))

times <- c("0h", "2h", "6h", "24h", "48h")

load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)

string_db <- STRINGdb$new(
  version = "12.5",
  species = 9606,
  score_threshold = 700,
  input_directory = ""
)

HCC44_clusters <- read.csv(
  "HCC44_Splineregression_Clusters_label.csv")
HCC44_spline_res<- read.csv(
  "HCC44_Splineregression_SignifcantProteins.csv"
) 

HCC44_timeData <- make_timeData(
  ready_pHCC44_0h,
  ready_pHCC44_2h,
  ready_pHCC44_6h,
  ready_pHCC44_24h,
  ready_pHCC44_48h
)

colnames(HCC44_timeData) <- c(
  "ProteinKey", "SYMBOL", "ENTREZID", "ENSEMBL",
  as.vector(rbind(
    paste0(times, "_log2FoldChange"),
    paste0(times, "_pval")
  ))
)

PIP_HCC44_fc_data <- HCC44_timeData %>%
  tidyr::pivot_longer(
    cols = matches("_(log2FoldChange|pval)$"),
    names_to = c("Time", ".value"),
    names_pattern = "^(.*)_(log2FoldChange|pval)$"
  ) %>%
  mutate(
    Time = factor(Time, levels = times)
  )
#--------------------------------------------
# HCC44 | Extracting Cluster 1 
#--------------------------------------------
PIP_HCC44_cluster1_ppi_input <- HCC44_clusters %>%
  filter(Cluster == 1) %>%
  rename(
    Gene = SYMBOL 
  )

PIP_HCC44_cluster1_ppi_input <- PIP_HCC44_cluster1_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

PIP_HCC44_cluster1_ppi_mapped <- string_db$map(
  PIP_HCC44_cluster1_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(PIP_HCC44_cluster1_ppi_input),
  "\nMapped proteins:",
  nrow(PIP_HCC44_cluster1_ppi_mapped),
  "\n"
)

PIP_HCC44_cluster1_unmapped <- PIP_HCC44_cluster1_ppi_input %>%
  filter(
    !Gene %in% PIP_HCC44_cluster1_ppi_mapped$Gene
  )

PIP_HCC44_cluster1_ppi_interactions <- string_db$get_interactions(
  PIP_HCC44_cluster1_ppi_mapped$STRING_id
)

PIP_HCC44_cluster1_id_to_gene <- PIP_HCC44_cluster1_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

PIP_HCC44_cluster1_ppi_edges <- PIP_HCC44_cluster1_ppi_interactions %>%
  left_join(
    PIP_HCC44_cluster1_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    PIP_HCC44_cluster1_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

PIP_HCC44_cluster1_ppi_edges_clean <- PIP_HCC44_cluster1_ppi_edges %>%
  filter(
    !is.na(from_gene),
    !is.na(to_gene),
    from_gene != to_gene
  ) %>%
  distinct(
    from_gene,
    to_gene,
    .keep_all = TRUE
  )

PIP_HCC44_cluster1_ppi_graph <- graph_from_data_frame(
  PIP_HCC44_cluster1_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

PIP_HCC44_cluster1_network_stats <- data.frame(
  Gene = V(PIP_HCC44_cluster1_ppi_graph)$name,
  Degree = degree(
    PIP_HCC44_cluster1_ppi_graph
  ),
  Betweenness = betweenness(
    PIP_HCC44_cluster1_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    PIP_HCC44_cluster1_ppi_graph,
    normalized = TRUE
  )
)

PIP_HCC44_cluster1_network_stats <- PIP_HCC44_cluster1_network_stats %>%
  arrange(
    desc(Degree)
  )

PIP_HCC44_cluster1_degree_cutoff <- quantile(
  PIP_HCC44_cluster1_network_stats$Degree,
  0.75
)

PIP_HCC44_cluster1_network_stats <- PIP_HCC44_cluster1_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= PIP_HCC44_cluster1_degree_cutoff,
      "Hub",
      "Other"
    )
  )

PIP_HCC44_cluster1_network_stats %>%
  filter(Hub == "Hub")

PIP_HCC44_cluster1_ppi_graph_plot <- as_tbl_graph(
  PIP_HCC44_cluster1_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    PIP_HCC44_cluster1_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# HCC44 | Cluster 1 - plotting protein-protein interaction
#--------------------------------------------
PIP_p_HCC44_cluster1_PPI <- ggraph(
  PIP_HCC44_cluster1_ppi_graph_plot,
  layout = "fr"
) +
  geom_edge_link(
    colour = "grey60",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_node_point(
    aes(
      size = Degree,
      fill = Hub
    ),
    shape = 21,
    colour = "black",
    stroke = 0.6
  ) +
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 5
  ) +
  scale_size_continuous(
    range = c(3, 8)
  ) +
  scale_fill_manual(
    values = c(
      "Hub" = "#EA4025",
      "Other" = "#D9D9D9"
    )
  ) +
  labs(
    size = "Degree",
    fill = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(
      size = 10
    )
  )
#--------------------------------------------
# HCC44 | Extracting Cluster 2 
#--------------------------------------------
PIP_HCC44_cluster2_ppi_input <- HCC44_clusters %>%
  filter(Cluster == 2) %>%
  rename(
    Gene = SYMBOL 
  )

PIP_HCC44_cluster2_ppi_input <- PIP_HCC44_cluster2_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

PIP_HCC44_cluster2_ppi_mapped <- string_db$map(
  PIP_HCC44_cluster2_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(PIP_HCC44_cluster2_ppi_input),
  "\nMapped proteins:",
  nrow(PIP_HCC44_cluster2_ppi_mapped),
  "\n"
)

PIP_HCC44_cluster2_unmapped <- PIP_HCC44_cluster2_ppi_input %>%
  filter(
    !Gene %in% PIP_HCC44_cluster2_ppi_mapped$Gene
  )

PIP_HCC44_cluster2_ppi_interactions <- string_db$get_interactions(
  PIP_HCC44_cluster2_ppi_mapped$STRING_id
)

PIP_HCC44_cluster2_id_to_gene <- PIP_HCC44_cluster2_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

PIP_HCC44_cluster2_ppi_edges <- PIP_HCC44_cluster2_ppi_interactions %>%
  left_join(
    PIP_HCC44_cluster2_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    PIP_HCC44_cluster2_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

PIP_HCC44_cluster2_ppi_edges_clean <- PIP_HCC44_cluster2_ppi_edges %>%
  filter(
    !is.na(from_gene),
    !is.na(to_gene),
    from_gene != to_gene
  ) %>%
  distinct(
    from_gene,
    to_gene,
    .keep_all = TRUE
  )

PIP_HCC44_cluster2_ppi_graph <- graph_from_data_frame(
  PIP_HCC44_cluster2_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

PIP_HCC44_cluster2_network_stats <- data.frame(
  Gene = V(PIP_HCC44_cluster2_ppi_graph)$name,
  Degree = degree(
    PIP_HCC44_cluster2_ppi_graph
  ),
  Betweenness = betweenness(
    PIP_HCC44_cluster2_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    PIP_HCC44_cluster2_ppi_graph,
    normalized = TRUE
  )
)

PIP_HCC44_cluster2_network_stats <- PIP_HCC44_cluster2_network_stats %>%
  arrange(
    desc(Degree)
  )

PIP_HCC44_cluster2_degree_cutoff <- quantile(
  PIP_HCC44_cluster2_network_stats$Degree,
  0.75
)

PIP_HCC44_cluster2_network_stats <- PIP_HCC44_cluster2_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= PIP_HCC44_cluster2_degree_cutoff,
      "Hub",
      "Other"
    )
  )

PIP_HCC44_cluster2_network_stats %>%
  filter(Hub == "Hub")

PIP_HCC44_cluster2_ppi_graph_plot <- as_tbl_graph(
  PIP_HCC44_cluster2_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    PIP_HCC44_cluster2_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# HCC44 | Cluster 2 - plotting protein-protein interaction
#--------------------------------------------
PIP_p_HCC44_cluster2_PPI <- ggraph(
  PIP_HCC44_cluster2_ppi_graph_plot,
  layout = "fr"
) +
  geom_edge_link(
    colour = "grey60",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_node_point(
    aes(
      size = Degree,
      fill = Hub
    ),
    shape = 21,
    colour = "black",
    stroke = 0.6
  ) +
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 5
  ) +
  scale_size_continuous(
    range = c(3, 8)
  ) +
  scale_fill_manual(
    values = c(
      "Hub" = "#EA4025",
      "Other" = "#D9D9D9"
    )
  ) +
  labs(
    size = "Degree",
    fill = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(
      size = 10
    )
  )
#--------------------------------------------
# HCC44 | Extracting Cluster 3 
#--------------------------------------------
PIP_HCC44_cluster3_ppi_input <- HCC44_clusters %>%
  filter(Cluster == 4) %>%
  rename(
    Gene =SYMBOL 
  )


PIP_HCC44_cluster3_ppi_input <- PIP_HCC44_cluster3_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

PIP_HCC44_cluster3_ppi_mapped <- string_db$map(
  PIP_HCC44_cluster3_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(PIP_HCC44_cluster3_ppi_input),
  "\nMapped proteins:",
  nrow(PIP_HCC44_cluster3_ppi_mapped),
  "\n"
)

PIP_HCC44_cluster3_unmapped <- PIP_HCC44_cluster3_ppi_input %>%
  filter(
    !Gene %in% PIP_HCC44_cluster3_ppi_mapped$Gene
  )

PIP_HCC44_cluster3_ppi_interactions <- string_db$get_interactions(
  PIP_HCC44_cluster3_ppi_mapped$STRING_id
)

PIP_HCC44_cluster3_id_to_gene <- PIP_HCC44_cluster3_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

PIP_HCC44_cluster3_ppi_edges <- PIP_HCC44_cluster3_ppi_interactions %>%
  left_join(
    PIP_HCC44_cluster3_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    PIP_HCC44_cluster3_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

PIP_HCC44_cluster3_ppi_edges_clean <- PIP_HCC44_cluster3_ppi_edges %>%
  filter(
    !is.na(from_gene),
    !is.na(to_gene),
    from_gene != to_gene
  ) %>%
  distinct(
    from_gene,
    to_gene,
    .keep_all = TRUE
  )

PIP_HCC44_cluster3_ppi_graph <- graph_from_data_frame(
  PIP_HCC44_cluster3_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

PIP_HCC44_cluster3_network_stats <- data.frame(
  Gene = V(PIP_HCC44_cluster3_ppi_graph)$name,
  Degree = degree(
    PIP_HCC44_cluster3_ppi_graph
  ),
  Betweenness = betweenness(
    PIP_HCC44_cluster3_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    PIP_HCC44_cluster3_ppi_graph,
    normalized = TRUE
  )
)

PIP_HCC44_cluster3_network_stats <- PIP_HCC44_cluster3_network_stats %>%
  arrange(
    desc(Degree)
  )

PIP_HCC44_cluster3_degree_cutoff <- quantile(
  PIP_HCC44_cluster3_network_stats$Degree,
  0.75
)

PIP_HCC44_cluster3_network_stats <- PIP_HCC44_cluster3_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= PIP_HCC44_cluster3_degree_cutoff,
      "Hub",
      "Other"
    )
  )

PIP_HCC44_cluster3_network_stats %>%
  filter(Hub == "Hub")

PIP_HCC44_cluster3_ppi_graph_plot <- as_tbl_graph(
  PIP_HCC44_cluster3_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    PIP_HCC44_cluster3_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# HCC44 | Cluster 3 - plotting protein-protein interaction
#--------------------------------------------
set.seed(123)

PIP_p_layout_HCC44_cluster3 <- igraph::layout_with_fr(
  PIP_HCC44_cluster3_ppi_graph_plot,
  niter = 2000,
  weights = NA,
  grid = "nogrid"
)

rownames(PIP_p_layout_HCC44_cluster3) <-
  igraph::V(PIP_HCC44_cluster3_ppi_graph_plot)$name

PIP_p_layout_HCC44_cluster3 <- PIP_p_layout_HCC44_cluster3 * 2.5

PIP_p_HCC44_cluster3_PPI <- ggraph(
  PIP_HCC44_cluster3_ppi_graph_plot,
  layout = "manual",
  x = PIP_p_layout_HCC44_cluster3[, 1],
  y = PIP_p_layout_HCC44_cluster3[, 2]
) +
  geom_edge_link(
    colour = "grey40",
    linewidth = 0.8,
    alpha = 0.7
  ) +
  geom_node_point(
    aes(
      size = Degree,
      fill = Hub
    ),
    shape = 21,
    colour = "black",
    stroke = 0.6
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 5,
    max.overlaps = Inf,
    fontface="bold"
  ) +
  scale_size_continuous(
    range = c(1, 10)
  ) +
  scale_fill_manual(
    values = c(
      "Hub" = "#EA4025",
      "Other" = "#D9D9D9"
    )
  ) +
  labs(
    size = "Degree",
    fill = NULL
  ) +
  coord_equal() +
  theme_void() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.text = element_text(
      size = 12
    )
  )

#--------------------------------------------
# NCI-23 | preparing data
#--------------------------------------------
H23_clusters <- read.csv(
  "NCI_H23_Splineregression_Clusters_label.csv")
H23_spline_res <- read.csv(
  "NCI-H23_Splineregression_SignifcantProteins.csv"
) 

load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)

H23_timeData <- make_timeData(
  ready_pH23_0h,
  ready_pH23_2h,
  ready_pH23_6h,
  ready_pH23_24h,
  ready_pH23_48h
)

colnames(H23_timeData) <- c(
  "ProteinKey", "SYMBOL", "ENTREZID", "ENSEMBL",
  as.vector(rbind(
    paste0(times, "_log2FoldChange"),
    paste0(times, "_pval")
  ))
)

H23_fc_data <- H23_timeData %>%
  tidyr::pivot_longer(
    cols = matches("_(log2FoldChange|pval)$"),
    names_to = c("Time", ".value"),
    names_pattern = "^(.*)_(log2FoldChange|pval)$"
  ) %>%
  mutate(
    Time = factor(Time, levels = times)
  )
#--------------------------------------------
# NCI-H23 | Extracting Cluster 1 
#--------------------------------------------
PIP_H23_cluster1_ppi_input <- H23_clusters %>%
  filter(Cluster == 1) %>%
  rename(
    Gene = SYMBOL
  )

PIP_H23_cluster1_ppi_input <- PIP_H23_cluster1_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

PIP_H23_cluster1_ppi_mapped <- string_db$map(
  PIP_H23_cluster1_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(PIP_H23_cluster1_ppi_input),
  "\nMapped proteins:",
  nrow(PIP_H23_cluster1_ppi_mapped),
  "\n"
)

PIP_H23_cluster1_unmapped <- PIP_H23_cluster1_ppi_input %>%
  filter(
    !Gene %in% PIP_H23_cluster1_ppi_mapped$Gene
  )

PIP_H23_cluster1_ppi_interactions <- string_db$get_interactions(
  PIP_H23_cluster1_ppi_mapped$STRING_id
)

PIP_H23_cluster1_id_to_gene <- PIP_H23_cluster1_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

PIP_H23_cluster1_ppi_edges <- PIP_H23_cluster1_ppi_interactions %>%
  left_join(
    PIP_H23_cluster1_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    PIP_H23_cluster1_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

PIP_H23_cluster1_ppi_edges_clean <- PIP_H23_cluster1_ppi_edges %>%
  filter(
    !is.na(from_gene),
    !is.na(to_gene),
    from_gene != to_gene
  ) %>%
  distinct(
    from_gene,
    to_gene,
    .keep_all = TRUE
  )

PIP_H23_cluster1_ppi_graph <- graph_from_data_frame(
  PIP_H23_cluster1_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

PIP_H23_cluster1_network_stats <- data.frame(
  Gene = V(PIP_H23_cluster1_ppi_graph)$name,
  Degree = degree(
    PIP_H23_cluster1_ppi_graph
  ),
  Betweenness = betweenness(
    PIP_H23_cluster1_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    PIP_H23_cluster1_ppi_graph,
    normalized = TRUE
  )
)

PIP_H23_cluster1_network_stats <- PIP_H23_cluster1_network_stats %>%
  arrange(
    desc(Degree)
  )

PIP_H23_cluster1_degree_cutoff <- quantile(
  PIP_H23_cluster1_network_stats$Degree,
  0.75
)

PIP_H23_cluster1_network_stats <- PIP_H23_cluster1_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= PIP_H23_cluster1_degree_cutoff,
      "Hub",
      "Other"
    )
  )

PIP_H23_cluster1_network_stats %>%
  filter(Hub == "Hub")

PIP_H23_cluster1_ppi_graph_plot <- as_tbl_graph(
  PIP_H23_cluster1_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    PIP_H23_cluster1_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# NCI-H23 | Cluster 1 - plotting protein-protein interaction
#--------------------------------------------
PIP_p_H23_cluster1_PPI <- ggraph(
  PIP_H23_cluster1_ppi_graph_plot,
  layout = "fr"
) +
  geom_edge_link(
    colour = "grey60",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_node_point(
    aes(
      size = Degree,
      fill = Hub
    ),
    shape = 21,
    colour = "black",
    stroke = 0.6
  ) +
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 5
  ) +
  scale_size_continuous(
    range = c(3, 8)
  ) +
  scale_fill_manual(
    values = c(
      "Hub" = "#EA4025",
      "Other" = "#D9D9D9"
    )
  ) +
  labs(
    size = "Degree",
    fill = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(
      size = 10
    )
  )
#--------------------------------------------
# NCI-H23 | Extracting Cluster 2 
#--------------------------------------------
PIP_H23_cluster2_ppi_input <- H23_clusters %>%
  filter(Cluster == 2) %>%
  rename(
    Gene = SYMBOL
  )

PIP_H23_cluster2_ppi_input <- PIP_H23_cluster2_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

PIP_H23_cluster2_ppi_mapped <- string_db$map(
  PIP_H23_cluster2_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(PIP_H23_cluster2_ppi_input),
  "\nMapped proteins:",
  nrow(PIP_H23_cluster2_ppi_mapped),
  "\n"
)

PIP_H23_cluster2_unmapped <- PIP_H23_cluster2_ppi_input %>%
  filter(
    !Gene %in% PIP_H23_cluster2_ppi_mapped$Gene
  )

PIP_H23_cluster2_ppi_interactions <- string_db$get_interactions(
  PIP_H23_cluster2_ppi_mapped$STRING_id
)

PIP_H23_cluster2_id_to_gene <- PIP_H23_cluster2_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

PIP_H23_cluster2_ppi_edges <- PIP_H23_cluster2_ppi_interactions %>%
  left_join(
    PIP_H23_cluster2_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    PIP_H23_cluster2_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

PIP_H23_cluster2_ppi_edges_clean <- PIP_H23_cluster2_ppi_edges %>%
  filter(
    !is.na(from_gene),
    !is.na(to_gene),
    from_gene != to_gene
  ) %>%
  distinct(
    from_gene,
    to_gene,
    .keep_all = TRUE
  )

PIP_H23_cluster2_ppi_graph <- graph_from_data_frame(
  PIP_H23_cluster2_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

PIP_H23_cluster2_network_stats <- data.frame(
  Gene = V(PIP_H23_cluster2_ppi_graph)$name,
  Degree = degree(
    PIP_H23_cluster2_ppi_graph
  ),
  Betweenness = betweenness(
    PIP_H23_cluster2_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    PIP_H23_cluster2_ppi_graph,
    normalized = TRUE
  )
)

PIP_H23_cluster2_network_stats <- PIP_H23_cluster2_network_stats %>%
  arrange(
    desc(Degree)
  )

PIP_H23_cluster2_degree_cutoff <- quantile(
  PIP_H23_cluster2_network_stats$Degree,
  0.75
)

PIP_H23_cluster2_network_stats <- PIP_H23_cluster2_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= PIP_H23_cluster2_degree_cutoff,
      "Hub",
      "Other"
    )
  )

PIP_H23_cluster2_network_stats %>%
  filter(Hub == "Hub")

PIP_H23_cluster2_ppi_graph_plot <- as_tbl_graph(
  PIP_H23_cluster2_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    PIP_H23_cluster2_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# NCI-H23 | Cluster 2 - plotting protein-protein interaction
#--------------------------------------------
PIP_p_H23_cluster2_PPI <- ggraph(
  PIP_H23_cluster2_ppi_graph_plot,
  layout = "fr"
) +
  geom_edge_link(
    colour = "grey60",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_node_point(
    aes(
      size = Degree,
      fill = Hub
    ),
    shape = 21,
    colour = "black",
    stroke = 0.6
  ) +
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 5
  ) +
  scale_size_continuous(
    range = c(3, 8)
  ) +
  scale_fill_manual(
    values = c(
      "Hub" = "#EA4025",
      "Other" = "#D9D9D9"
    )
  ) +
  labs(
    size = "Degree",
    fill = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(
      size = 10
    )
  )

#--------------------------------------------
# Saving results
#--------------------------------------------
dir.create(
  file.path(dirs$results, "23_protein_ProteinProteinInteraction"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "23_protein_ProteinProteinInteraction"))

#1 saving PIP HCC44 Cluster 1
ggsave(
  "23_PIP_HCC44_Cluster1_STRING_PPI.pdf",
  plot = PIP_p_HCC44_cluster1_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "23_PIP_HCC44_Cluster1_STRING_PPI.jpg",
  plot = PIP_p_HCC44_cluster1_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  PIP_HCC44_cluster1_network_stats,
  "23_PIP_HCC44_Cluster1_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  PIP_HCC44_cluster1_ppi_edges_clean,
  "23_PIP_HCC44_Cluster1_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  PIP_HCC44_cluster1_unmapped,
  "23_PIP_HCC44_Cluster1_STRING_unmapped.csv",
  row.names = FALSE
)

#2 saving PIP HCC44 Cluster 2
ggsave(
  "23_PIP_HCC44_Cluster2_STRING_PPI.pdf",
  plot = PIP_p_HCC44_cluster2_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "23_PIP_HCC44_Cluster2_STRING_PPI.jpg",
  plot = PIP_p_HCC44_cluster2_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  PIP_HCC44_cluster2_network_stats,
  "23_PIP_HCC44_Cluster2_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  PIP_HCC44_cluster2_ppi_edges_clean,
  "23_PIP_HCC44_Cluster2_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  PIP_HCC44_cluster2_unmapped,
  "23_PIP_HCC44_Cluster2_STRING_unmapped.csv",
  row.names = FALSE
)

#3 saving PIP HCC44 Cluster 3
ggsave(
  "23_PIP_HCC44_Cluster3_STRING_PPI.pdf",
  plot = PIP_p_HCC44_cluster3_PPI,
  width = 21*2/3,
  height = 14,
  units = "cm"
)

ggsave(
  "23_PIP_HCC44_Cluster3_STRING_PPI.jpg",
  plot = PIP_p_HCC44_cluster3_PPI,
  width = 21*2/3,
  height = 14,
  units = "cm",
  dpi = 300
)

write.csv(
  PIP_HCC44_cluster3_network_stats,
  "23_PIP_HCC44_Cluster3_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  PIP_HCC44_cluster3_ppi_edges_clean,
  "23_PIP_HCC44_Cluster3_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  PIP_HCC44_cluster3_unmapped,
  "23_PIP_HCC44_Cluster3_STRING_unmapped.csv",
  row.names = FALSE
)

#4 saving PIP NCI-H23 Cluster 1
ggsave(
  "23_PIP_H23_Cluster1_STRING_PPI.pdf",
  plot = PIP_p_H23_cluster1_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "23_PIP_H23_Cluster1_STRING_PPI.jpg",
  plot = PIP_p_H23_cluster1_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  PIP_H23_cluster1_network_stats,
  "23_PIP_H23_Cluster1_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  PIP_H23_cluster1_ppi_edges_clean,
  "23_PIP_H23_Cluster1_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  PIP_H23_cluster1_unmapped,
  "23_PIP_H23_Cluster1_STRING_unmapped.csv",
  row.names = FALSE
)

#5 saving PIP NCI-H23 Cluster 2
ggsave(
  "23_PIP_H23_Cluster2_STRING_PPI.pdf",
  plot = PIP_p_H23_cluster2_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "23_PIP_H23_Cluster2_STRING_PPI.jpg",
  plot = PIP_p_H23_cluster2_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  PIP_H23_cluster2_network_stats,
  "23_PIP_H23_Cluster2_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  PIP_H23_cluster2_ppi_edges_clean,
  "23_PIP_H23_Cluster2_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  PIP_H23_cluster2_unmapped,
  "23_PIP_H23_Cluster2_STRING_unmapped.csv",
  row.names = FALSE
)
#-------------------------------------------
setwd(dirs$code)
