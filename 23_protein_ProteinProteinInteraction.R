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
#--------------------------------------------
# HCC44 - preparing data
#--------------------------------------------
HCC44_clusters <- read.csv(
  "HCC44_Splineregression_Clusters.csv")
HCC44_spline_res<- read.csv(
  "HCC44_Splineregression_significantProteins.csv"
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

HCC44_fc_data <- HCC44_timeData %>%
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
HCC44_cluster1_ppi_input <- HCC44_clusters %>%
  filter(Cluster == 1) %>%
  rename(
    Gene = names
  )

HCC44_cluster1_ppi_input <- HCC44_cluster1_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

HCC44_cluster1_ppi_mapped <- string_db$map(
  HCC44_cluster1_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(HCC44_cluster1_ppi_input),
  "\nMapped proteins:",
  nrow(HCC44_cluster1_ppi_mapped),
  "\n"
)

HCC44_cluster1_unmapped <- HCC44_cluster1_ppi_input %>%
  filter(
    !Gene %in% HCC44_cluster1_ppi_mapped$Gene
  )

HCC44_cluster1_ppi_interactions <- string_db$get_interactions(
  HCC44_cluster1_ppi_mapped$STRING_id
)

HCC44_cluster1_id_to_gene <- HCC44_cluster1_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

HCC44_cluster1_ppi_edges <- HCC44_cluster1_ppi_interactions %>%
  left_join(
    HCC44_cluster1_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    HCC44_cluster1_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

HCC44_cluster1_ppi_edges_clean <- HCC44_cluster1_ppi_edges %>%
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

HCC44_cluster1_ppi_graph <- graph_from_data_frame(
  HCC44_cluster1_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

HCC44_cluster1_network_stats <- data.frame(
  Gene = V(HCC44_cluster1_ppi_graph)$name,
  Degree = degree(
    HCC44_cluster1_ppi_graph
  ),
  Betweenness = betweenness(
    HCC44_cluster1_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    HCC44_cluster1_ppi_graph,
    normalized = TRUE
  )
)

HCC44_cluster1_network_stats <- HCC44_cluster1_network_stats %>%
  arrange(
    desc(Degree)
  )

HCC44_cluster1_degree_cutoff <- quantile(
  HCC44_cluster1_network_stats$Degree,
  0.75
)

HCC44_cluster1_network_stats <- HCC44_cluster1_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= HCC44_cluster1_degree_cutoff,
      "Hub",
      "Other"
    )
  )

HCC44_cluster1_network_stats %>%
  filter(Hub == "Hub")

HCC44_cluster1_ppi_graph_plot <- as_tbl_graph(
  HCC44_cluster1_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    HCC44_cluster1_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# HCC44 | Cluster 1 - plotting protein-protein interaction
#--------------------------------------------
HCC44_cluster1_p_PPI <- ggraph(
  HCC44_cluster1_ppi_graph_plot,
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
HCC44_cluster2_ppi_input <- HCC44_clusters %>%
  filter(Cluster == 2) %>%
  rename(
    Gene = names
  )

HCC44_cluster2_ppi_input <- HCC44_cluster2_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

HCC44_cluster2_ppi_mapped <- string_db$map(
  HCC44_cluster2_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(HCC44_cluster2_ppi_input),
  "\nMapped proteins:",
  nrow(HCC44_cluster2_ppi_mapped),
  "\n"
)

HCC44_cluster2_unmapped <- HCC44_cluster2_ppi_input %>%
  filter(
    !Gene %in% HCC44_cluster2_ppi_mapped$Gene
  )

HCC44_cluster2_ppi_interactions <- string_db$get_interactions(
  HCC44_cluster2_ppi_mapped$STRING_id
)

HCC44_cluster2_id_to_gene <- HCC44_cluster2_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

HCC44_cluster2_ppi_edges <- HCC44_cluster2_ppi_interactions %>%
  left_join(
    HCC44_cluster2_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    HCC44_cluster2_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

HCC44_cluster2_ppi_edges_clean <- HCC44_cluster2_ppi_edges %>%
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

HCC44_cluster2_ppi_graph <- graph_from_data_frame(
  HCC44_cluster2_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

HCC44_cluster2_network_stats <- data.frame(
  Gene = V(HCC44_cluster2_ppi_graph)$name,
  Degree = degree(
    HCC44_cluster2_ppi_graph
  ),
  Betweenness = betweenness(
    HCC44_cluster2_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    HCC44_cluster2_ppi_graph,
    normalized = TRUE
  )
)

HCC44_cluster2_network_stats <- HCC44_cluster2_network_stats %>%
  arrange(
    desc(Degree)
  )

HCC44_cluster2_degree_cutoff <- quantile(
  HCC44_cluster2_network_stats$Degree,
  0.75
)

HCC44_cluster2_network_stats <- HCC44_cluster2_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= HCC44_cluster2_degree_cutoff,
      "Hub",
      "Other"
    )
  )

HCC44_cluster2_network_stats %>%
  filter(Hub == "Hub")

HCC44_cluster2_ppi_graph_plot <- as_tbl_graph(
  HCC44_cluster2_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    HCC44_cluster2_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# HCC44 | Cluster 2 - plotting protein-protein interaction
#--------------------------------------------
HCC44_cluster2_p_PPI <- ggraph(
  HCC44_cluster2_ppi_graph_plot,
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
HCC44_cluster3_ppi_input <- HCC44_clusters %>%
  filter(Cluster == 3) %>%
  rename(
    Gene = names
  )

HCC44_cluster3_ppi_input <- HCC44_cluster3_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

HCC44_cluster3_ppi_mapped <- string_db$map(
  HCC44_cluster3_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(HCC44_cluster3_ppi_input),
  "\nMapped proteins:",
  nrow(HCC44_cluster3_ppi_mapped),
  "\n"
)

HCC44_cluster3_unmapped <- HCC44_cluster3_ppi_input %>%
  filter(
    !Gene %in% HCC44_cluster3_ppi_mapped$Gene
  )

HCC44_cluster3_ppi_interactions <- string_db$get_interactions(
  HCC44_cluster3_ppi_mapped$STRING_id
)

HCC44_cluster3_id_to_gene <- HCC44_cluster3_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

HCC44_cluster3_ppi_edges <- HCC44_cluster3_ppi_interactions %>%
  left_join(
    HCC44_cluster3_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    HCC44_cluster3_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

HCC44_cluster3_ppi_edges_clean <- HCC44_cluster3_ppi_edges %>%
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

HCC44_cluster3_ppi_graph <- graph_from_data_frame(
  HCC44_cluster3_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

HCC44_cluster3_network_stats <- data.frame(
  Gene = V(HCC44_cluster3_ppi_graph)$name,
  Degree = degree(
    HCC44_cluster3_ppi_graph
  ),
  Betweenness = betweenness(
    HCC44_cluster3_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    HCC44_cluster3_ppi_graph,
    normalized = TRUE
  )
)

HCC44_cluster3_network_stats <- HCC44_cluster3_network_stats %>%
  arrange(
    desc(Degree)
  )

HCC44_cluster3_degree_cutoff <- quantile(
  HCC44_cluster3_network_stats$Degree,
  0.75
)

HCC44_cluster3_network_stats <- HCC44_cluster3_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= HCC44_cluster3_degree_cutoff,
      "Hub",
      "Other"
    )
  )

HCC44_cluster3_network_stats %>%
  filter(Hub == "Hub")

HCC44_cluster3_ppi_graph_plot <- as_tbl_graph(
  HCC44_cluster3_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    HCC44_cluster3_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# HCC44 | Cluster 3 - plotting protein-protein interaction
#--------------------------------------------
HCC44_cluster3_p_PPI <- ggraph(
  HCC44_cluster3_ppi_graph_plot,
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
# HCC44 | Extracting all spline significant proteins 
#--------------------------------------------
HCC44_clusterAll_ppi_input <- data.frame(
  Gene =HCC44_spline_res[,12]
)

HCC44_clusterAll_ppi_mapped <- string_db$map(
  HCC44_clusterAll_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(HCC44_clusterAll_ppi_input),
  "\nMapped proteins:",
  nrow(HCC44_clusterAll_ppi_mapped),
  "\n"
)

HCC44_clusterAll_unmapped <- HCC44_clusterAll_ppi_input %>%
  filter(
    !Gene %in% HCC44_clusterAll_ppi_mapped$Gene
  )

HCC44_clusterAll_ppi_interactions <- string_db$get_interactions(
  HCC44_clusterAll_ppi_mapped$STRING_id
)

colnames(HCC44_clusterAll_ppi_interactions)
HCC44_clusterAll_id_to_gene <- HCC44_clusterAll_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

HCC44_clusterAll_ppi_edges <- HCC44_clusterAll_ppi_interactions %>%
  left_join(
    HCC44_clusterAll_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    HCC44_clusterAll_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

HCC44_clusterAll_ppi_edges_clean <- HCC44_clusterAll_ppi_edges %>%
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

HCC44_clusterAll_ppi_graph <- graph_from_data_frame(
  HCC44_clusterAll_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

HCC44_clusterAll_network_stats <- data.frame(
  Gene = V(HCC44_clusterAll_ppi_graph)$name,
  Degree = degree(
    HCC44_clusterAll_ppi_graph
  ),
  Betweenness = betweenness(
    HCC44_clusterAll_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    HCC44_clusterAll_ppi_graph,
    normalized = TRUE
  )
)

HCC44_clusterAll_network_stats <- HCC44_clusterAll_network_stats %>%
  arrange(
    desc(Degree)
  )

HCC44_clusterAll_network_stats <- HCC44_clusterAll_network_stats %>%
  dplyr::mutate(
    Hub = ifelse(
      Degree >= quantile(Degree, 0.75, na.rm = TRUE),
      "Hub",
      "Other"
    )
  )

HCC44_clusterAll_network_stats %>%
  filter(Hub == "Hub")

HCC44_clusterAll_ppi_graph_plot <- as_tbl_graph(
  HCC44_clusterAll_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    HCC44_clusterAll_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# HCC44 | Plotting weighted PIP network
#--------------------------------------------
HCC44_clusterAll_ppi_input <- HCC44_clusterAll_ppi_input %>%
  filter(
    !is.na(Gene ),
    Gene  != ""
  )%>%
  distinct(Gene ) %>%
  rename(
    SYMBOL = Gene 
  )

HCC44_spline_PPI <- make_temporal_PPI(
  protein_data = HCC44_clusterAll_ppi_input,
  fc_data = HCC44_fc_data,
  string_db = string_db
)

p_HCC44_spline_PPI <- plot_temporal_PPI(
  HCC44_spline_PPI
)
#--------------------------------------------
# HCC44 | all Cluster - plotting protein-protein interaction
#--------------------------------------------
HCC44_clusterAll_p_PPI <- ggraph(
  HCC44_clusterAll_ppi_graph_plot,
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
# HCC44 | select significant genes 
#--------------------------------------------
HCC44_significant_proteins <- HCC44_fc_data %>%
  filter(
    !is.na(pval),
    pval < 0.05,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > 1
  ) %>%
  distinct(
    ProteinKey, SYMBOL, ENTREZID, ENSEMBL,
    Time, log2FoldChange, pval
  )

HCC44_significant_PPI_input <- HCC44_significant_proteins %>%
  filter(
    !is.na(SYMBOL),
    trimws(SYMBOL) != ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

HCC44_significant_PPI_excluded <- HCC44_significant_proteins %>%
  filter(
    is.na(SYMBOL) | trimws(SYMBOL) == ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

HCC44_PPI_input <- HCC44_significant_proteins %>%
  dplyr::filter(
    !is.na(SYMBOL),
    trimws(SYMBOL) != ""
  ) %>%
  dplyr::transmute(Gene = SYMBOL) %>%
  dplyr::distinct()

gene_groups <- split(
  HCC44_PPI_input$Gene,
  ceiling(seq_along(HCC44_PPI_input$Gene) / 50)
)

HCC44_significant_PPI_mapped <- lapply(names(gene_groups), function(group) {
  string_db$map(
    data.frame(Gene = gene_groups[[group]]),
    "Gene",
    removeUnmappedRows = TRUE
  ) %>%
    dplyr::mutate(Group = group)
  
}) %>%
  dplyr::bind_rows() %>%
  dplyr::distinct()
dplyr::n_distinct(HCC44_significant_PPI_mapped$Gene)
HCC44_significant_PPI_unmapped <- HCC44_PPI_input %>%
  dplyr::anti_join(
    HCC44_significant_PPI_mapped,
    by = "Gene"
  )

HCC44_PPI_significant <- make_temporal_PPI(
  protein_data = HCC44_significant_PPI_input,
  fc_data = HCC44_significant_proteins,
  string_db = string_db
)
  
p_HCC44_PPI_significant <- plot_temporal_PPI(
  HCC44_PPI_significant
)
#--------------------------------------------
# NCI-23 | preparing data
#--------------------------------------------
H23_clusters <- read.csv(
  "NCI-H23_Splineregression_Clusters.csv")
H23_spline_res <- read.csv(
  "NCI-H23_Splineregression_significantProteins.csv"
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
H23_cluster1_ppi_input <- H23_clusters %>%
  filter(Cluster == 1) %>%
  rename(
    Gene = names
  )

H23_cluster1_ppi_input <- H23_cluster1_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

H23_cluster1_ppi_mapped <- string_db$map(
  H23_cluster1_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(H23_cluster1_ppi_input),
  "\nMapped proteins:",
  nrow(H23_cluster1_ppi_mapped),
  "\n"
)

H23_cluster1_unmapped <- H23_cluster1_ppi_input %>%
  filter(
    !Gene %in% H23_cluster1_ppi_mapped$Gene
  )

H23_cluster1_ppi_interactions <- string_db$get_interactions(
  H23_cluster1_ppi_mapped$STRING_id
)

H23_cluster1_id_to_gene <- H23_cluster1_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

H23_cluster1_ppi_edges <- H23_cluster1_ppi_interactions %>%
  left_join(
    H23_cluster1_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    H23_cluster1_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

H23_cluster1_ppi_edges_clean <- H23_cluster1_ppi_edges %>%
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

H23_cluster1_ppi_graph <- graph_from_data_frame(
  H23_cluster1_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

H23_cluster1_network_stats <- data.frame(
  Gene = V(H23_cluster1_ppi_graph)$name,
  Degree = degree(
    H23_cluster1_ppi_graph
  ),
  Betweenness = betweenness(
    H23_cluster1_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    H23_cluster1_ppi_graph,
    normalized = TRUE
  )
)

H23_cluster1_network_stats <- H23_cluster1_network_stats %>%
  arrange(
    desc(Degree)
  )

H23_cluster1_degree_cutoff <- quantile(
  H23_cluster1_network_stats$Degree,
  0.75
)

H23_cluster1_network_stats <- H23_cluster1_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= H23_cluster1_degree_cutoff,
      "Hub",
      "Other"
    )
  )

H23_cluster1_network_stats %>%
  filter(Hub == "Hub")

H23_cluster1_ppi_graph_plot <- as_tbl_graph(
  H23_cluster1_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    H23_cluster1_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# NCI-H23 | Cluster 1 - plotting protein-protein interaction
#--------------------------------------------
H23_cluster1_p_PPI <- ggraph(
  H23_cluster1_ppi_graph_plot,
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
H23_cluster2_ppi_input <- H23_clusters %>%
  filter(Cluster == 2) %>%
  rename(
    Gene = names
  )

H23_cluster2_ppi_input <- H23_cluster2_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  distinct(Gene)

H23_cluster2_ppi_mapped <- string_db$map(
  H23_cluster2_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(H23_cluster2_ppi_input),
  "\nMapped proteins:",
  nrow(H23_cluster2_ppi_mapped),
  "\n"
)

H23_cluster2_unmapped <- H23_cluster2_ppi_input %>%
  filter(
    !Gene %in% H23_cluster2_ppi_mapped$Gene
  )

H23_cluster2_ppi_interactions <- string_db$get_interactions(
  H23_cluster2_ppi_mapped$STRING_id
)

H23_cluster2_id_to_gene <- H23_cluster2_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

H23_cluster2_ppi_edges <- H23_cluster2_ppi_interactions %>%
  left_join(
    H23_cluster2_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    H23_cluster2_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

H23_cluster2_ppi_edges_clean <- H23_cluster2_ppi_edges %>%
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

H23_cluster2_ppi_graph <- graph_from_data_frame(
  H23_cluster2_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

H23_cluster2_network_stats <- data.frame(
  Gene = V(H23_cluster2_ppi_graph)$name,
  Degree = degree(
    H23_cluster2_ppi_graph
  ),
  Betweenness = betweenness(
    H23_cluster2_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    H23_cluster2_ppi_graph,
    normalized = TRUE
  )
)

H23_cluster2_network_stats <- H23_cluster2_network_stats %>%
  arrange(
    desc(Degree)
  )

H23_cluster2_degree_cutoff <- quantile(
  H23_cluster2_network_stats$Degree,
  0.75
)

H23_cluster2_network_stats <- H23_cluster2_network_stats %>%
  mutate(
    Hub = ifelse(
      Degree >= H23_cluster2_degree_cutoff,
      "Hub",
      "Other"
    )
  )

H23_cluster2_network_stats %>%
  filter(Hub == "Hub")

H23_cluster2_ppi_graph_plot <- as_tbl_graph(
  H23_cluster2_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    H23_cluster2_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# NCI-H23 | Cluster 2 - plotting protein-protein interaction
#--------------------------------------------
H23_cluster2_p_PPI <- ggraph(
  H23_cluster2_ppi_graph_plot,
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
# NCI-H23 | Extracting all spline significant proteins 
#--------------------------------------------
H23_clusterAll_ppi_input <- data.frame(
  Gene =H23_spline_res[,12]
)

H23_clusterAll_ppi_mapped <- string_db$map(
  H23_clusterAll_ppi_input,
  "Gene",
  removeUnmappedRows = TRUE
)

cat(
  "Input proteins:",
  nrow(H23_clusterAll_ppi_input),
  "\nMapped proteins:",
  nrow(H23_clusterAll_ppi_mapped),
  "\n"
)

H23_clusterAll_unmapped <- H23_clusterAll_ppi_input %>%
  filter(
    !Gene %in% H23_clusterAll_ppi_mapped$Gene
  )

H23_clusterAll_ppi_interactions <- string_db$get_interactions(
  H23_clusterAll_ppi_mapped$STRING_id
)

colnames(H23_clusterAll_ppi_interactions)
H23_clusterAll_id_to_gene <- H23_clusterAll_ppi_mapped %>%
  select(
    STRING_id,
    Gene
  ) %>%
  distinct()

H23_clusterAll_ppi_edges <- H23_clusterAll_ppi_interactions %>%
  left_join(
    H23_clusterAll_id_to_gene,
    by = c("from" = "STRING_id")
  ) %>%
  rename(
    from_gene = Gene
  ) %>%
  left_join(
    H23_clusterAll_id_to_gene,
    by = c("to" = "STRING_id")
  ) %>%
  rename(
    to_gene = Gene
  )

H23_clusterAll_ppi_edges_clean <- H23_clusterAll_ppi_edges %>%
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

H23_clusterAll_ppi_graph <- graph_from_data_frame(
  H23_clusterAll_ppi_edges_clean %>%
    select(
      from_gene,
      to_gene
    ),
  directed = FALSE
)

H23_clusterAll_network_stats <- data.frame(
  Gene = V(H23_clusterAll_ppi_graph)$name,
  Degree = degree(
    H23_clusterAll_ppi_graph
  ),
  Betweenness = betweenness(
    H23_clusterAll_ppi_graph,
    normalized = TRUE
  ),
  Closeness = closeness(
    H23_clusterAll_ppi_graph,
    normalized = TRUE
  )
)

H23_clusterAll_network_stats <- H23_clusterAll_network_stats %>%
  arrange(
    desc(Degree)
  )

H23_clusterAll_network_stats <- H23_clusterAll_network_stats %>%
  dplyr::mutate(
    Hub = ifelse(
      Degree >= quantile(Degree, 0.75, na.rm = TRUE),
      "Hub",
      "Other"
    )
  )

H23_clusterAll_network_stats %>%
  filter(Hub == "Hub")

H23_clusterAll_ppi_graph_plot <- as_tbl_graph(
  H23_clusterAll_ppi_graph
) %>%
  activate(nodes) %>%
  left_join(
    H23_clusterAll_network_stats,
    by = c("name" = "Gene")
  )
#--------------------------------------------
# NCI-H23 | Plotting weighted PIP network
#--------------------------------------------
H23_clusterAll_ppi_input <- H23_clusterAll_ppi_input %>%
  filter(
    !is.na(Gene),
    Gene != ""
  )%>%
  distinct(Gene) %>%
  rename(
    SYMBOL = Gene
  )

H23_fc_data <-H23_fc_data%>%
  rename(
#    Gene = SYMBOL
  )

H23_spline_PPI <- make_temporal_PPI(
  protein_data = H23_clusterAll_ppi_input,
  fc_data = H23_fc_data,
  string_db = string_db
)

p_H23_spline_PPI <- plot_temporal_PPI(
  H23_spline_PPI
)
#--------------------------------------------
# NCI-H23 | all Cluster - plotting protein-protein interaction
#--------------------------------------------
H23_clusterAll_p_PPI <- ggraph(
  H23_clusterAll_ppi_graph_plot,
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
# NCI-H23 | select significant genes 
#--------------------------------------------
H23_significant_proteins <- H23_fc_data %>%
  filter(
    !is.na(pval),
    pval < 0.05,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > 1
  ) %>%
  distinct(
    ProteinKey, SYMBOL, ENTREZID, ENSEMBL,
    Time, log2FoldChange, pval
  )

H23_significant_PPI_input <- H23_significant_proteins %>%
  filter(
    !is.na(SYMBOL),
    trimws(SYMBOL) != ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

H23_significant_PPI_excluded <- H23_significant_proteins %>%
  filter(
    is.na(SYMBOL) | trimws(SYMBOL) == ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

H23_significant_PPI_mapped <- lapply(gene_groups, function(genes) {
  string_db$map(
    data.frame(Gene = genes),
    "Gene",
    removeUnmappedRows = TRUE
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::distinct()

dplyr::n_distinct(H23_significant_PPI_mapped$Gene)
H23_significant_PPI_unmapped <- H23_significant_PPI_input %>%
  dplyr::anti_join(
    H23_significant_PPI_mapped,
    by = c("SYMBOL"="Gene")
  )

H23_PPI_significant <- make_temporal_PPI(
  protein_data = H23_significant_PPI_input,
  fc_data = H23_significant_proteins,
  string_db = string_db
)

p_H23_PPI_significant <- plot_temporal_PPI(
  H23_PPI_significant
)
#--------------------------------------------
# Control | preparing data
#--------------------------------------------
Con_timeData <- make_timeData(
  ready_pCon_0h,
  ready_pCon_2h,
  ready_pCon_6h,
  ready_pCon_24h,
  ready_pCon_48h
)

colnames(Con_timeData) <- c(
  "ProteinKey", "SYMBOL", "ENTREZID", "ENSEMBL",
  as.vector(rbind(
    paste0(times, "_log2FoldChange"),
    paste0(times, "_pval")
  ))
)

Con_fc_data <- Con_timeData %>%
  tidyr::pivot_longer(
    cols = matches("_(log2FoldChange|pval)$"),
    names_to = c("Time", ".value"),
    names_pattern = "^(.*)_(log2FoldChange|pval)$"
  ) %>%
  mutate(
    Time = factor(Time, levels = times)
  )
#--------------------------------------------
# Control | select significant genes 
#--------------------------------------------
Con_significant_proteins <- Con_fc_data %>%
  filter(
    !is.na(pval),
    pval < 0.05,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > 5
  ) %>%
  distinct(
    ProteinKey, SYMBOL, ENTREZID, ENSEMBL,
    Time, log2FoldChange, pval
  )

Con_significant_PPI_input <- Con_significant_proteins %>%
  filter(
    !is.na(SYMBOL),
    trimws(SYMBOL) != ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

Con_significant_PPI_excluded <- Con_significant_proteins %>%
  filter(
    is.na(SYMBOL) | trimws(SYMBOL) == ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

Con_significant_PPI_mapped <- lapply(gene_groups, function(genes) {
  string_db$map(
    data.frame(Gene = genes),
    "Gene",
    removeUnmappedRows = TRUE
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::distinct()

dplyr::n_distinct(Con_significant_PPI_mapped$Gene)
Con_significant_PPI_unmapped <- Con_significant_PPI_input %>%
  dplyr::anti_join(
    Con_significant_PPI_mapped,
    by = c("SYMBOL"="Gene")
  )

Con_PPI_significant <- make_temporal_PPI(
  protein_data = Con_significant_PPI_input,
  fc_data = Con_significant_proteins,
  string_db = string_db
)

p_Con_PPI_significant <- plot_temporal_PPI(
  Con_PPI_significant
)
#--------------------------------------------
# Radiation | preparing Data
#--------------------------------------------
Rad_timeData <- make_timeData_Rad(
  ready_pRAD_2h,
  ready_pRAD_6h,
  ready_pRAD_24h,
  ready_pRAD_48h
)

colnames(Rad_timeData) <- c(
  "ProteinKey", "SYMBOL", "ENTREZID", "ENSEMBL",
  as.vector(rbind(
    paste0(times[2:5], "_log2FoldChange"),
    paste0(times[2:5], "_pval")
  ))
)

Rad_fc_data <- Rad_timeData %>%
  tidyr::pivot_longer(
    cols = matches("_(log2FoldChange|pval)$"),
    names_to = c("Time", ".value"),
    names_pattern = "^(.*)_(log2FoldChange|pval)$"
  ) %>%
  mutate(
    Time = factor(Time, levels = times)
  )
#--------------------------------------------
# Radiation | select significant genes 
#--------------------------------------------
Rad_significant_proteins <- Rad_fc_data %>%
  filter(
    !is.na(pval),
    pval < 0.05,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > 5
  ) %>%
  distinct(
    ProteinKey, SYMBOL, ENTREZID, ENSEMBL,
    Time, log2FoldChange, pval
  )

Rad_significant_PPI_input <- Rad_significant_proteins %>%
  filter(
    !is.na(SYMBOL),
    trimws(SYMBOL) != ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

Rad_significant_PPI_excluded <- Rad_significant_proteins %>%
  filter(
    is.na(SYMBOL) | trimws(SYMBOL) == ""
  ) %>%
  distinct(ProteinKey, SYMBOL)

Rad_significant_PPI_mapped <- lapply(gene_groups, function(genes) {
  string_db$map(
    data.frame(Gene = genes),
    "Gene",
    removeUnmappedRows = TRUE
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::distinct()

dplyr::n_distinct(Rad_significant_PPI_mapped$Gene)
Rad_significant_PPI_unmapped <- Rad_significant_PPI_input %>%
  dplyr::anti_join(
    Rad_significant_PPI_mapped,
    by = c("SYMBOL"="Gene")
  )

Rad_PPI_significant <- make_temporal_PPI_Rad(
  protein_data = Rad_significant_PPI_input,
  fc_data = Rad_significant_proteins,
  string_db = string_db
)

p_Rad_PPI_significant <- plot_temporal_PPI(
  Rad_PPI_significant
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
  "HCC44_Cluster1_STRING_PPI.pdf",
  plot = HCC44_cluster1_p_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "HCC44_Cluster1_STRING_PPI.jpg",
  plot = HCC44_cluster1_p_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  HCC44_cluster1_network_stats,
  "HCC44_Cluster1_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  HCC44_cluster1_ppi_edges_clean,
  "HCC44_Cluster1_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  HCC44_cluster1_unmapped,
  "HCC44_Cluster1_STRING_unmapped.csv",
  row.names = FALSE
)

#2 saving PIP HCC44 Cluster 2
ggsave(
  "HCC44_Cluster2_STRING_PPI.pdf",
  plot = HCC44_cluster2_p_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "HCC44_Cluster2_STRING_PPI.jpg",
  plot = HCC44_cluster2_p_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  HCC44_cluster2_network_stats,
  "HCC44_Cluster2_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  HCC44_cluster2_ppi_edges_clean,
  "HCC44_Cluster2_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  HCC44_cluster2_unmapped,
  "HCC44_Cluster2_STRING_unmapped.csv",
  row.names = FALSE
)

#3 saving PIP HCC44 Cluster 3
ggsave(
  "HCC44_Cluster3_STRING_PPI.pdf",
  plot = HCC44_cluster3_p_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "HCC44_Cluster3_STRING_PPI.jpg",
  plot = HCC44_cluster3_p_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  HCC44_cluster3_network_stats,
  "HCC44_Cluster3_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  HCC44_cluster3_ppi_edges_clean,
  "HCC44_Cluster3_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  HCC44_cluster3_unmapped,
  "HCC44_Cluster3_STRING_unmapped.csv",
  row.names = FALSE
)

#4 saving PIP all spline regression significant proteins in HCC44
ggsave(
  "HCC44_ClusterAll_STRING_PPI.pdf",
  plot = HCC44_clusterAll_p_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "HCC44_ClusterAll_STRING_PPI.jpg",
  plot = HCC44_clusterAll_p_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  HCC44_clusterAll_network_stats,
  "HCC44_clusterAll_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  HCC44_clusterAll_ppi_edges_clean,
  "HCC44_clusterAll_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  HCC44_clusterAll_unmapped,
  "HCC44_clusterAll_STRING_unmapped.csv",
  row.names = FALSE
)

#5 saving PIP all spline regression significant proteins in HCC44, weighted
ggsave(
  "HCC44_Cluster1_STRING_PPI.pdf",
  plot = p_HCC44_spline_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "HCC44_splineWeight_STRING_PPI.jpg",
  plot = p_HCC44_spline_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#6 saving PIP with all significant proteins in HCC44
ggsave(
  "HCC44_significant_STRING_PPI.pdf",
  plot = p_HCC44_PPI_significant,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "HCC44_significant_STRING_PPI.jpg",
  plot = p_HCC44_PPI_significant,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#7 saving PIP NCI-H23 Cluster 1
ggsave(
  "H23_Cluster1_STRING_PPI.pdf",
  plot = H23_cluster1_p_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "H23_Cluster1_STRING_PPI.jpg",
  plot = H23_cluster1_p_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  H23_cluster1_network_stats,
  "H23_Cluster1_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  H23_cluster1_ppi_edges_clean,
  "H23_Cluster1_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  H23_cluster1_unmapped,
  "H23_Cluster1_STRING_unmapped.csv",
  row.names = FALSE
)

#8 saving PIP NCI-H23 Cluster 2
ggsave(
  "H23_Cluster2_STRING_PPI.pdf",
  plot = H23_cluster2_p_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "H23_Cluster2_STRING_PPI.jpg",
  plot = H23_cluster2_p_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  H23_cluster2_network_stats,
  "H23_Cluster2_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  H23_cluster2_ppi_edges_clean,
  "H23_Cluster2_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  H23_cluster2_unmapped,
  "H23_Cluster2_STRING_unmapped.csv",
  row.names = FALSE
)

#9 saving PIP all spline regression significant proteins in NCI-H23
ggsave(
  "H23_ClusterAll_STRING_PPI.pdf",
  plot = H23_clusterAll_p_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "H23_ClusterAll_STRING_PPI.jpg",
  plot = H23_clusterAll_p_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

write.csv(
  H23_clusterAll_network_stats,
  "H23_clusterAll_STRING_network_statistics.csv",
  row.names = FALSE
)

write.csv(
  H23_clusterAll_ppi_edges_clean,
  "H23_clusterAll_STRING_interactions.csv",
  row.names = FALSE
)

write.csv(
  H23_clusterAll_unmapped,
  "H23_clusterAll_STRING_unmapped.csv",
  row.names = FALSE
)

#10 saving PIP all spline regression significant proteins in NCI-H23, weighted
ggsave(
  "H23_Cluster1_STRING_PPI.pdf",
  plot = p_H23_spline_PPI,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "H23_splineWeight_STRING_PPI.jpg",
  plot = p_H23_spline_PPI,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#11 saving PIP with all significant proteins in NCI-H23
ggsave(
  "H23_significant_STRING_PPI.pdf",
  plot = p_H23_PPI_significant,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "H23_significant_STRING_PPI.jpg",
  plot = p_H23_PPI_significant,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#12 saving PIP with all significant proteins in control
ggsave(
  "Control_significant_STRING_PPI.pdf",
  plot = p_Con_PPI_significant,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "Control_significant_STRING_PPI.jpg",
  plot = p_Con_PPI_significant,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#12 saving PIP with all significant proteins in control
ggsave(
  "Radiation_significant_STRING_PPI.pdf",
  plot = p_Rad_PPI_significant,
  width = 42*2/3,
  height = 28,
  units = "cm"
)

ggsave(
  "Radiation_significant_STRING_PPI.jpg",
  plot = p_Rad_PPI_significant,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)
#-------------------------------------------
setwd(dirs$code)
