#--------------------------------------------
# load and prepare data
#--------------------------------------------
load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)

names <- list()
for (cellline in c("H23", "HCC44")) {
  name <- paste0("spline_", cellline)
  names[[cellline]] <- name
  assign(
    name,
    dplyr::bind_rows(
      prepare_spline_data(cellline, "c"),
      prepare_spline_data(cellline, "2Gy")
    )
  )
}


annotation_col <- data.frame(
  Time = factor(
    rep(c("0h", "2h", "6h", "24h", "48h"), each = 8),
    levels = c("0h", "2h", "6h", "24h", "48h")
  ),
  Treatment = factor(
    rep(c(rep("0 Gy", 4), rep("2 Gy", 4)), 5),
    levels = c("0 Gy", "2 Gy")
  )
)

annotation_rad <- data.frame(
  Cellline = factor(
    rep("HCC44", 20)
  ),
  Treatment = factor(
    rep("2Gy", 20)
  ),
  Time = factor(
    rep(
      c("0h", "2h", "6h", "24h", "48h"),
      each = 4
    ),
    levels = c("0h", "2h", "6h", "24h", "48h")
  )
)

annotation_all <- data.frame(
  Cellline = factor(
    rep("HCC44", 20)
  ),
  Treatment = factor(
    rep(c("c","2Gy"), each=4)
  ),
  Time = factor(
    rep(
      c("0h", "2h", "6h", "24h", "48h"),
      each = 8
    ),
    levels = c("0h", "2h", "6h", "24h", "48h")
  )
)

heatmap_colors <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))

rad_cols <- unlist(lapply(
  c("0h", "2h", "6h", "24h", "48h"),
  function(t) paste0(t, "_2Gy_Rep_", 1:4)
))

genes_to_label <- c(
  "ANLN", "AURKA", "BUB1", "CCNB1", "CCNB2", "CDC20",
  "CENPE", "CKAP2", "GTSE1", "HMMR", "PLK1", "TROAP"
)

annotation_colors <- list(
  Time = c(
    "0h" = "#4294F7",
    "2h" = "#5F8F29",
    "6h" = "#F09837",
    "24h" = "#EA4025",
    "48h" = "#872111"
  )
)

times <- c("0h", "2h", "6h", "24h", "48h")
treatments <- c("c", "2Gy")

names <- unlist(lapply(times, function(time) {
  paste0("raw_HCC44_", time, "_", treatments)
}))

timepoints <- c(0, 2, 6, 24, 48)
#--------------------------------------------
# HCC44 | preparing HCC44 data for spline regression
#--------------------------------------------
spline_HCC44_2Gy <- cbind(raw_HCC44_0h_2Gy[,c(2:3,5:8)],
                        raw_HCC44_2h_2Gy[,c(5:8)],
                        raw_HCC44_6h_2Gy[,c(5:8)],
                        raw_HCC44_24h_2Gy[,c(5:8)],
                        raw_HCC44_48h_2Gy[,c(5:8)])

spline_HCC44_2Gy[, 1:2] <- lapply(
  spline_HCC44_2Gy[, 1:2],
  as.character
)

colnames(spline_HCC44_2Gy) <- c("ID","SYMBOL",rad_cols)

names_HCC44 <- unlist(lapply(times, function(time) {
  paste0("raw_HCC44_", time, "_", treatments)
}))

spline_HCC44_all <- purrr::reduce(
  lapply(names_HCC44, function(name) {
    data <- get(name, envir = .GlobalEnv)
    data %>%
      dplyr::select(
        ProteinKey,
        dplyr::all_of(paste0("Rep_", 1:4))
      ) %>%
      dplyr::rename_with(
        ~ paste0(sub("^raw_HCC44_", "", name), "_", .x),
        dplyr::starts_with("Rep_")
      )
  }),
  dplyr::full_join,
  by = "ProteinKey"
)%>%
  mutate(
    ProteinKey=as.character(ProteinKey)    
  )

sample_names <- unlist(lapply(times, function(time) {
  unlist(lapply(treatments, function(treatment) {
    paste0(time, "_", treatment, "_Rep_", 1:4)
  }))
}))

expr_matrix_HCC44 <- as.matrix(spline_HCC44_2Gy[, -(1:2)])
rownames(expr_matrix_HCC44) <- spline_HCC44_2Gy$ProteinKey
#--------------------------------------------
# HCC44 | Regression model
#--------------------------------------------
targets_HCC44_time <- data.frame(
  Sample = colnames(expr_matrix_HCC44),
  Time = rep(c(0, 2, 6, 24, 48), each = 4)
)

design_HCC44 <- model.matrix(
  ~ ns(Time, df = 3),
  data = targets_HCC44_time
)

colnames(design_HCC44) <- c(
  "Intercept",
  "Spline1",
  "Spline2",
  "Spline3"
)

fit_HCC44_time <- lmFit(
  expr_matrix_HCC44,
  design_HCC44
)

fit2_HCC44_time <- eBayes(fit_HCC44_time)

results_HCC44_time <- topTable(
  fit2_HCC44_time,
  coef = 2:4,
  adjust.method = "BH",
  number = Inf
) %>%
  rownames_to_column("ProteinKey") %>%
  left_join(annotation, by = "ProteinKey")

#--------------------------------------------
# HCC44 | selecting significant proteins
#--------------------------------------------
res_sig_HCC44_time <- results_HCC44_time %>%
  filter(adj.P.Val < 0.05)%>%
  mutate(
    ProteinKey=as.character(ProteinKey),
    SYMBOL = as.character(PG.Genes)
  )

spline_sig_HCC44 <- spline_HCC44_all %>%
  inner_join(res_sig_HCC44_time, by ="ProteinKey")

mat_Spline_HCC44 <- spline_sig_HCC44[complete.cases(spline_sig_HCC44[,2:42]),]
mat_Spline_HCC44_z <- as.matrix(mat_Spline_HCC44[,2:41])
mat_Spline_HCC44_z <- t(scale(t(mat_Spline_HCC44_z)))

rownames(annotation_col) <- colnames(mat_Spline_HCC44_z)
#--------------------------------------------
# HCC44 | generating names
#-------------------------------------------
heatmap_matrix_HCC44_rad <- as.matrix(
  mat_Spline_HCC44_z[, rad_cols]
)

rownames(heatmap_matrix_HCC44_rad)<-mat_Spline_HCC44$ProteinKey

heatmap_matrix_HCC44_2 <- mat_Spline_HCC44_z[
  rownames(mat_Spline_HCC44_z),
  ,
  drop = FALSE
]

hm_cluster_HCC44<-pheatmap(
  heatmap_matrix_HCC44_rad,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = c(4, 8, 12, 16),
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  width = 21/3,
  height = 14
)

HCC44_clusters <- data.frame(
  ProteinKey = as.character(names(cutree(hm_cluster_HCC44$tree_row, k = 3))),
  Cluster = as.integer(cutree(hm_cluster_HCC44$tree_row, k = 3))
) %>%
  dplyr::inner_join(
    res_sig_HCC44_time[, c(1,9:12)],
    by = "ProteinKey"
  )

HCC44_protein_order <- rownames(heatmap_matrix_HCC44_rad)[
  hm_cluster_HCC44$tree_row$order
]

HCC44_cluster_order <- unique(
  HCC44_clusters$Cluster[
    match(HCC44_protein_order, HCC44_clusters$ProteinKey)
  ]
)
HCC44_missing_names <- is.na(HCC44_protein_order) | HCC44_protein_order == ""
HCC44_protein_names <- mat_Spline_HCC44$SYMBOL
HCC44_protein_names[HCC44_missing_names] <- paste0(
  "NA_",
  1:1
)

HCC44_protein_order <- make.unique(HCC44_protein_names)
HCC44_target_cluster <- HCC44_cluster_order[3]
HCC44_labels_row <- data.frame(
  Label = rep("", length(HCC44_protein_order)),
  names = as.character(HCC44_protein_order),
  ProteinKey=as.character(HCC44_clusters$ProteinKey),
  Cluster=as.character(HCC44_clusters$Cluster)
)

HCC44_genes_in_cluster <- HCC44_clusters %>%
  dplyr::filter(Cluster == HCC44_target_cluster) %>%
  dplyr::pull(SYMBOL)

HCC44_labels_row$Label <- ifelse(
  HCC44_labels_row$names %in% HCC44_genes_in_cluster,
  HCC44_labels_row$names,
  ""
)
#-------------------------------------------
# HCC44 2Gy | hierarchical clustering and mapping
#-------------------------------------------
rownames(annotation_rad) <- colnames(heatmap_matrix_HCC44_rad)

pheatmap(
  heatmap_matrix_HCC44_rad,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = hm_cluster_HCC44$tree_row,
  cluster_cols = FALSE,
  gaps_col = c(4, 8, 12, 16),
  annotation_col = annotation_rad,
  annotation_colors = annotation_colors,
  labels_row = HCC44_labels_row[,1],
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  filename =   file.path(dirs$results, "21_Proteome_SplineRegression",
  "HCC44_Spline_2Gy_heatmap.pdf"),
  width = 21/3,
  height = 14
)
#-------------------------------------------
# HCC44 all | hierarchical clustering and mapping
#-------------------------------------------
mat_Spline_HCC44_2 <- as.matrix(
  mat_Spline_HCC44_z)

rownames(annotation_all) <- colnames(mat_Spline_HCC44_2)

pheatmap(
  mat_Spline_HCC44_2,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = c( 8, 16, 24, 32),
  annotation_col = annotation_all,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  filename =   file.path(dirs$results, "21_Proteome_SplineRegression",
    "HCC44_Spline_all_heatmap.pdf"),
  width = 21/3,
  height = 14
)
#--------------------------------------------
# NCI-H23 | preparing HCC44 data for spline regression
#--------------------------------------------
spline_H23_2Gy <- cbind(raw_H23_0h_2Gy[,c(2:3,5:8)],
                          raw_H23_2h_2Gy[,c(5:8)],
                          raw_H23_6h_2Gy[,c(5:8)],
                          raw_H23_24h_2Gy[,c(5:8)],
                          raw_H23_48h_2Gy[,c(5:8)])

colnames(spline_H23_2Gy) <- c("ID","SYMBOL",rad_cols)

names_H23 <- unlist(lapply(times, function(time) {
  paste0("raw_H23_", time, "_", treatments)
}))

spline_H23_all <- purrr::reduce(
  lapply(names_H23, function(name) {
    data <- get(name, envir = .GlobalEnv)
    data %>%
      dplyr::select(
        ProteinKey,
        dplyr::all_of(paste0("Rep_", 1:4))
      ) %>%
      dplyr::rename_with(
        ~ paste0(sub("^raw_H23_", "", name), "_", .x),
        dplyr::starts_with("Rep_")
      )
  }),
  dplyr::full_join,
  by = "ProteinKey"
)%>%
  mutate(
    ProteinKey=as.character(ProteinKey)    
  )

sample_names <- unlist(lapply(times, function(time) {
  unlist(lapply(treatments, function(treatment) {
    paste0(time, "_", treatment, "_Rep_", 1:4)
  }))
}))

expr_matrix_H23 <- as.matrix(spline_H23_2Gy[, -(1:2)])

rownames(expr_matrix_H23) <- spline_H23_2Gy$ProteinKey
#--------------------------------------------
# NCI-H23 | Regression model
#--------------------------------------------
targets_H23_time <- data.frame(
  Sample = colnames(expr_matrix_H23),
  Time = rep(c(0, 2, 6, 24, 48), each = 4)
)

design_H23 <- model.matrix(
  ~ ns(Time, df = 3),
  data = targets_H23_time
)

colnames(design_H23) <- c(
  "Intercept",
  "Spline1",
  "Spline2",
  "Spline3"
)

fit_H23_time <- lmFit(
  expr_matrix_H23,
  design_H23
)

fit2_H23_time <- eBayes(fit_H23_time)

results_H23_time <- topTable(
  fit2_H23_time,
  coef = 2:4,
  adjust.method = "BH",
  number = Inf
) %>%
  rownames_to_column("ProteinKey") %>%
  left_join(annotation, by = "ProteinKey")
#--------------------------------------------
# NCI-H23 | selecting significant proteins
#--------------------------------------------
res_sig_H23_time <- results_H23_time %>%
  filter(adj.P.Val < 0.05)%>%
  mutate(
    ProteinKey=as.character(ProteinKey),
    SYMBOL = as.character(PG.Genes)
  )

spline_sig_H23 <- spline_H23_all %>%
  inner_join(res_sig_H23_time, by =c("ProteinKey"))

mat_Spline_H23 <- spline_sig_H23[complete.cases(spline_sig_H23[,2:42]),]
mat_Spline_H23_z <- as.matrix(mat_Spline_H23[,2:41])
mat_Spline_H23_z <- t(scale(t(mat_Spline_H23_z)))

rownames(annotation_col) <- colnames(mat_Spline_H23_z)
#--------------------------------------------
# NCI-H23 | generating names
#-------------------------------------------
heatmap_matrix_H23_rad <- as.matrix(
  mat_Spline_H23_z[, rad_cols]
)

rownames(heatmap_matrix_H23_rad)<-mat_Spline_H23$ProteinKey

heatmap_matrix_H23_rad_2 <- mat_Spline_H23_z[
  rownames(mat_Spline_H23_z),
  ,
  drop = FALSE
]

hm_cluster_H23<-pheatmap(
  heatmap_matrix_H23_rad,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = c(4, 8, 12, 16),
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  width = 21/3,
  height = 14
)

H23_clusters <- data.frame(
  ProteinKey = names(cutree(hm_cluster_H23$tree_row, k = 2)),
  Cluster = as.integer(cutree(hm_cluster_H23$tree_row, k = 2))
) %>%
  dplyr::inner_join(
    res_sig_H23_time[, c(1,9:12)],
    by = "ProteinKey"
  )

H23_protein_order <- rownames(heatmap_matrix_H23_rad)[
  hm_cluster_HCC44$tree_row$order
]

H23_cluster_order <- unique(
  H23_clusters$Cluster[
    match(H23_protein_order, H23_clusters$ProteinKey)
  ]
)
H23_missing_names <- is.na(H23_protein_order) | H23_protein_order == ""
H23_protein_names <- mat_Spline_H23$SYMBOL

H23_protein_names[H23_missing_names] <- paste0(
  "NA_",
  1:1
)

H23_protein_order <- make.unique(H23_protein_names)
H23_target_cluster <- H23_cluster_order[1]
H23_labels_row <- data.frame(
  Label = rep("", length(H23_protein_order)),
  names = as.character(H23_protein_order),
  ProteinKey=as.character(H23_clusters$ProteinKey),
  Cluster=as.character(H23_clusters$Cluster)
)

H23_genes_in_cluster <- H23_clusters %>%
  dplyr::filter(
    Cluster == H23_target_cluster,
    grepl("^MRPL", SYMBOL)
  ) %>%
  dplyr::pull(SYMBOL)

H23_labels_row$Label <- ifelse(
  H23_labels_row$names %in% H23_genes_in_cluster,
  H23_labels_row$names,
  ""
)
#-------------------------------------------
# NCI-H23 2Gy | hierarchical clustering and mapping
#-------------------------------------------
rownames(annotation_rad) <- colnames(heatmap_matrix_H23_rad)

pheatmap(
  heatmap_matrix_H23_rad,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = hm_cluster_H23$tree_row,
  cluster_cols = FALSE,
  gaps_col = c(4, 8, 12, 16),
  annotation_col = annotation_rad,
  annotation_colors = annotation_colors,
  labels_row = H23_labels_row[,1],
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  filename =   file.path(dirs$results, "21_Proteome_SplineRegression",
                         "H23_Spline_2Gy_heatmap.pdf"),
  width = 21/3,
  height = 14
)
#-------------------------------------------
# NCI-H23 all | hierarchical clustering and mapping
#-------------------------------------------
mat_Spline_H23_2 <- as.matrix(
  mat_Spline_H23_z)

rownames(annotation_all) <- colnames(mat_Spline_H23_2)

pheatmap(
  mat_Spline_H23_2,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = c( 8, 16, 24, 32),
  annotation_col = annotation_all,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 8,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  filename =   file.path(dirs$results, "21_Proteome_SplineRegression",
                         "H23_Spline_all_heatmap.pdf"),
  width = 21/3,
  height = 14
)
#-------------------------------------------
# saving results
#-------------------------------------------
dir.create(
  file.path(dirs$results, "21_Proteome_SplineRegression"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "21_Proteome_SplineRegression"))

write.csv(
  H23_labels_row,
  file = file.path(dirs$results, "21_Proteome_SplineRegression", "NCI-H23_Splineregression_Clusters.csv"),
  row.names = FALSE
)

write.csv(
  res_sig_H23_time,
  file = file.path(dirs$results, "21_Proteome_SplineRegression", "NCI-H23_Splineregression_significantProteins.csv"),
  row.names = FALSE
)

write.csv(
  res_sig_HCC44_time,
  file = file.path(dirs$results, "21_Proteome_SplineRegression", "HCC44_Splineregression_significantProteins.csv"),
  row.names = FALSE
)

write.csv(
  HCC44_labels_row,
  file = file.path(dirs$results, "21_Proteome_SplineRegression", "HCC44_Splineregression_Clusters.csv"),
  row.names = FALSE
)
#-------------------------------------------
setwd(dirs$code)
