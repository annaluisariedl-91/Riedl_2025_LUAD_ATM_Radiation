spline_targets_HCC44_time <- data.frame(
  Sample = colnames(spline_expr_matrix_HCC44),
  Time = rep(c(0, 2, 6, 24, 48), each = 4)
)

spline_design_HCC44 <- model.matrix(
  ~ ns(Time, df = 3),
  data = spline_targets_HCC44_time
)

colnames(spline_design_HCC44) <- c(
  "Intercept",
  "spline1",
  "spline2",
  "spline3"# ,
  #"spline4"
  
)

spline_fit_HCC44_time <- lmFit(
  spline_expr_matrix_HCC44,
  spline_design_HCC44
)

spline_fit2_HCC44_time <- eBayes(spline_fit_HCC44_time)

spline_results_HCC44_time <- topTable(
  spline_fit2_HCC44_time,
  coef = 2:4,
  adjust.method = "BH",
  number = Inf
) %>%
  rownames_to_column("SYMBOL")
#--------------------------------------------
# HCC44 | selecting significant proteins
#--------------------------------------------
spline_sig_HCC44 <- spline_results_HCC44_time %>%
  dplyr::mutate(ProteinKey = as.character(.[[1]])) %>%
  dplyr::inner_join(spline_HCC44_2Gy[-1], by = "ProteinKey") %>%
  dplyr::filter(adj.P.Val < 0.05)

spline_mat_HCC44 <- spline_sig_HCC44 %>%
  dplyr::select(dplyr::all_of(rad_cols)) %>%
  as.matrix()

rownames(spline_mat_HCC44) <- make.unique(
  as.character(spline_sig_HCC44[[2]])
)

spline_mat_HCC44 <- spline_mat_HCC44[
  complete.cases(spline_mat_HCC44),
  ,
  drop = FALSE
]

spline_mat_HCC444_z <- t(scale(t(spline_mat_HCC44)))
#--------------------------------------------
# HCC44 | generating names
#-------------------------------------------
spline_heatmap_matrix_HCC44_rad <- as.matrix(
  spline_mat_HCC444_z[, rad_cols]
)

spline_hm_cluster_HCC44<-pheatmap(
  spline_heatmap_matrix_HCC44_rad[,rad_cols],
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
  filename =   file.path(dirs$results, "21_Proteome_splineRegression",
                         "21_Proteome_HCC44_spline_2Gy_names.pdf"),
  annotation_legend = TRUE,
  width = 21/3,
  height = 14
)

spline_HCC44_cluster_assignment <- cutree(
  spline_hm_cluster_HCC44$tree_row,
  k = 9
)

spline_HCC44_clusters <- data.frame(
  SYMBOL = names(spline_HCC44_cluster_assignment),
  Cluster = as.integer(spline_HCC44_cluster_assignment)
) %>%
  dplyr::left_join(
    spline_sig_HCC44 %>%
      dplyr::select(ProteinKey,ID) %>%
      dplyr::distinct(ID, .keep_all = TRUE),
    by = c("SYMBOL"="ID")
  )

spline_HCC44_protein_order <- rownames(spline_heatmap_matrix_HCC44_rad)[
  spline_hm_cluster_HCC44$tree_row$order
]

spline_HCC44_clusters_ordered <- spline_HCC44_clusters[
  match(
    spline_HCC44_protein_order,
    spline_HCC44_clusters$SYMBOL
  ),
]

spline_cluster_sizes_HCC44 <- rle(spline_HCC44_clusters_ordered$Cluster)$lengths

spline_HCC44_cluster_order <- rle(
  spline_HCC44_clusters_ordered$Cluster
)$values

spline_HCC44_gaps_row <- head(cumsum(spline_cluster_sizes_HCC44), -1)

spline_HCC44_target_cluster <- spline_HCC44_cluster_order[4]

spline_HCC44_genes_in_cluster <- spline_HCC44_clusters %>%
  dplyr::filter(Cluster == spline_HCC44_target_cluster) %>%
  dplyr::pull(ProteinKey ,SYMBOL)

spline_HCC44_labels <- spline_HCC44_clusters %>%
   dplyr::mutate(
    Label = ifelse(
      ProteinKey  %in% spline_HCC44_genes_in_cluster,
      SYMBOL ,
      ""
    ),
    Label = dplyr::coalesce(Label, "")
  )

rownames(annotation_rad) <- colnames(spline_heatmap_matrix_HCC44_rad)
#-------------------------------------------
# HCC44 2Gy | hierarchical clustering and mapping
#-------------------------------------------
spline_HCC44_hm_Cluster_mark<- pheatmap(
  spline_heatmap_matrix_HCC44_rad,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = spline_hm_cluster_HCC44$tree_row,
  cluster_cols = FALSE,
  gaps_col = c(4, 8, 12, 16),
  annotation_col = annotation_rad,
  annotation_colors = annotation_colors,
  labels_row = spline_HCC44_labels$Label,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 16,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  filename =   file.path(dirs$results, "21_Proteome_splineRegression",
  "21_Proteome_splineRegression_splineregression_Heatmap_HCC44.pdf"),
  width = 21/3,
  height = 14
)
#-------------------------------------------
# HCC44 2Gy |marking relevant clusters
#-------------------------------------------
spline_HCC44_tree_pos <- which(spline_HCC44_hm_Cluster_mark$gtable$layout$name == "row_tree")
spline_HCC44_hm_Cluster_mark$gtable <- gtable_add_cols(
  spline_HCC44_hm_Cluster_mark$gtable,
  unit(12, "mm"),
  pos = max(spline_HCC44_hm_Cluster_mark$gtable$layout$l[spline_HCC44_tree_pos])-1
)

spline_HCC44_tree_pos <- which(spline_HCC44_hm_Cluster_mark$gtable$layout$name == "row_tree")
spline_HCC44_matrix_pos <- which(spline_HCC44_hm_Cluster_mark$gtable$layout$name == "matrix")
spline_HCC44_n <- nrow(spline_HCC44_hm_Cluster_mark)

spline_cluster_labels <- c("A", "B", "C")

spline_HCC44_cluster_ranges <- data.frame(
  Cluster = c("A", "B", "C"),
  Start = c(1,(sum(spline_cluster_sizes_HCC44[c(1:3)])+1),(sum(spline_cluster_sizes_HCC44[c(1:3)])+1)),
  End = c(sum(spline_cluster_sizes_HCC44[c(1:3)]), sum(spline_cluster_sizes_HCC44),  sum(spline_cluster_sizes_HCC44[c(1:4)])))

spline_HCC44_brackets <- grobTree(
  children = do.call(
    gList,
    lapply(seq_len(nrow(spline_HCC44_cluster_ranges)), function(i) {
      n<- max(spline_HCC44_cluster_ranges$End)
      cluster <- spline_HCC44_cluster_ranges$Cluster[i]
      y_top <- 1 - (spline_HCC44_cluster_ranges$Start[i] - 1) / n
      y_bottom <- 1 - spline_HCC44_cluster_ranges$End[i] / n
      y_mid <- (y_top + y_bottom) / 2
      x <- if (cluster == "C") 0.35 else 0.75
      grobTree(
        segmentsGrob(
          x0 = x, x1 = x,
          y0 = y_bottom, y1 = y_top,
          gp = gpar(lwd = 2.5)
        ),
        segmentsGrob(
          x0 = c(x, x),
          x1 = c(x + 0.20, x + 0.20),
          y0 = c(y_top, y_bottom),
          y1 = c(y_top, y_bottom),
          gp = gpar(lwd = 2.5)
        ),
        textGrob(
          cluster,
          x = x - 0.08,
          y = y_mid,
          just = "right",
          gp = gpar(
            fontsize = 14,
            fontface = "bold"
          )
        )
      )
    })
  )
)

spline_HCC44_hm_Cluster_mark$gtable <- gtable_add_grob(
  spline_HCC44_hm_Cluster_mark$gtable,
  grobs = spline_HCC44_brackets,
  t = spline_HCC44_hm_Cluster_mark$gtable$layout$t[spline_HCC44_matrix_pos],
  b = spline_HCC44_hm_Cluster_mark$gtable$layout$b[spline_HCC44_matrix_pos],
  l = spline_HCC44_hm_Cluster_mark$gtable$layout$l[spline_HCC44_tree_pos]-1,
  r = spline_HCC44_hm_Cluster_mark$gtable$layout$r[spline_HCC44_tree_pos]-1,
  name = "cluster_brackets"
)
#-------------------------------------------
# HCC44 2Gy | saving results
#-------------------------------------------
pdf(
  file.path(
    dirs$results,
    "21_Proteome_splineRegression",
    "21_Proteome_splineRegression_splineregression_Heatmap_HCC44.pdf"
  ),
  width = 10.5 / 2.54,
  height = 14 / 2.54
)

grid.newpage()
grid.draw(spline_HCC44_hm_Cluster_mark$gtable)

dev.off()
#--------------------------------------------
# NCI-H23 | preparing HCC44 data for spline regression
#--------------------------------------------
spline_H23_2Gy <- cbind(raw_pH23_0h_2Gy[,c(2:3,5:8)],
                          raw_pH23_2h_2Gy[,c(5:8)],
                          raw_pH23_6h_2Gy[,c(5:8)],
                          raw_pH23_24h_2Gy[,c(5:8)],
                          raw_pH23_48h_2Gy[,c(5:8)])

spline_H23_2Gy[, 1:2] <- lapply(
  spline_H23_2Gy[, 1:2],
  as.character
)

colnames(spline_H23_2Gy) <- c("SYMBOL","ProteinKey",rad_cols)

spline_H23_all <- spline_H23_2Gy %>%
  tidyr::pivot_longer(
    cols = -c(SYMBOL, ProteinKey),
    names_to = c("Time", "Treatment", "Replicate"),
    names_pattern = "(\\d+)h_(\\d+Gy)_Rep_(\\d+)",
    values_to = "Intensity"
  ) %>%
  dplyr::mutate(
    Time = as.numeric(Time),
    Replicate = as.integer(Replicate),
    Sample = paste0(Time, "h_", Treatment, "_", Replicate)
  ) %>%
  dplyr::select(
    ProteinKey,
    Sample,
    Intensity
  ) %>%
  tidyr::pivot_wider(
    names_from = Sample,
    values_from = Intensity
  ) %>%
  dplyr::mutate(
    ProteinKey = as.character(ProteinKey)
  )

spline_expr_matrix_H23 <- as.matrix(spline_H23_2Gy[,3:22])
rownames(spline_expr_matrix_H23) <- makeUnique(spline_H23_2Gy$SYMBOL)
#--------------------------------------------
# NCI-H23 | Regression model
#--------------------------------------------
spline_targets_H23_time <- data.frame(
  Sample = colnames(spline_expr_matrix_H23),
  Time = rep(c(0, 2, 6, 24, 48), each = 4)
)

spline_design_H23 <- model.matrix(
  ~ splines::ns(Time, df = 3),
  data = spline_targets_H23_time
)

colnames(spline_design_H23) <- c(
  "Intercept",
  "spline1",
  "spline2",
  "spline3"
)

spline_fit_H23_time <- lmFit(
  spline_expr_matrix_H23,
  spline_design_H23
)

spline_fit2_H23_time <- eBayes(spline_fit_H23_time)

spline_results_H23_time <- topTable(
  spline_fit2_H23_time,
  coef = 2:4,
  adjust.method = "BH",
  number = Inf
) %>%
  rename(
    SYMBOL=ID) 
#--------------------------------------------
# NCI-H23 | selecting significant proteins
#--------------------------------------------
spline_sig_H23 <- spline_results_H23_time %>%
  dplyr::mutate(ProteinKey = as.character(rownames(spline_results_H23_time)))  %>%
  dplyr::inner_join(spline_H23_2Gy[-1], by = "ProteinKey") %>%
  dplyr::filter(adj.P.Val < 0.05)

spline_mat_H23 <- spline_sig_H23 %>%
  dplyr::select(dplyr::all_of(rad_cols)) %>%
  as.matrix()

rownames(spline_mat_H23) <- make.unique(
  as.character(spline_sig_H23[[1]])
)

spline_mat_H23 <- spline_mat_H23[
  complete.cases(spline_mat_H23),
  ,
  drop = FALSE
]

spline_mat_H23_z <- t(scale(t(spline_mat_H23)))
#--------------------------------------------
# NCI-H23 | generating names
#-------------------------------------------
spline_heatmap_matrix_H23_rad <- as.matrix(
  spline_mat_H23_z[, rad_cols]
)

spline_hm_cluster_H23<-pheatmap(
  spline_heatmap_matrix_H23_rad[,rad_cols],
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
  filename =   file.path(dirs$results, "21_Proteome_splineRegression",
                         "21_Proteome_H23_spline_2Gy_names.pdf"),
  annotation_legend = TRUE,
  width = 21/3,
  height = 14
)

spline_H23_cluster_assignment <- cutree(
  spline_hm_cluster_H23$tree_row,
  k = 2
)

spline_H23_clusters <- data.frame(
  SYMBOL = names(spline_H23_cluster_assignment),
  Cluster = as.integer(spline_H23_cluster_assignment)
) %>%
  dplyr::left_join(
    spline_sig_H23 %>%
      dplyr::select(SYMBOL,ProteinKey) %>%
      dplyr::distinct(SYMBOL,ProteinKey, .keep_all = TRUE),
    by = c("SYMBOL")
  )

spline_H23_protein_order <- rownames(spline_heatmap_matrix_H23_rad)[
  spline_hm_cluster_H23$tree_row$order
]

spline_H23_clusters_ordered <- spline_H23_clusters[
  match(
    spline_H23_protein_order,
    spline_H23_clusters$SYMBOL
  ),
]

spline_cluster_sizes_H23 <- rle(spline_H23_clusters_ordered$Cluster)$lengths

spline_H23_cluster_order <- rle(
  spline_H23_clusters_ordered$Cluster
)$values

spline_H23_genes_in_cluster <- spline_H23_clusters %>%
  dplyr::filter(
    grepl("MRPL", SYMBOL, ignore.case = TRUE)
  ) %>%
  dplyr::select(SYMBOL)

spline_H23_labels <-  spline_H23_clusters %>%
  dplyr::mutate(
    Label = ifelse(
      SYMBOL %in% spline_H23_genes_in_cluster$SYMBOL,
      SYMBOL,
      ""
    ),
    Label = dplyr::coalesce(Label, "")
  )

rownames(annotation_rad) <- colnames(spline_heatmap_matrix_H23_rad)

spline_H23_labels_row = spline_H23_labels$Label[
  match(
    rownames(spline_heatmap_matrix_H23_rad),
    spline_H23_labels$SYMBOL
  )
]
#-------------------------------------------
# NCI-H23 2Gy | hierarchical clustering and mapping
#-------------------------------------------
spline_H23_hm_Cluster_mark<- pheatmap(
  spline_heatmap_matrix_H23_rad,
  color = heatmap_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = spline_hm_cluster_H23$tree_row,
  cluster_cols = FALSE,
  gaps_col = c(4, 8, 12, 16),
  annotation_col = annotation_rad,
  annotation_colors = annotation_colors,
  labels_row = spline_H23_labels_row,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 16,
  fontsize_row = 5,
  fontsize_col = 8,
  border_color = NA,
  cellwidth = 12,
  cellheight = 4,
  legend = TRUE,
  annotation_legend = TRUE,
  filename =   file.path(dirs$results, "21_Proteome_splineRegression",
                         "21_Proteome_splineRegression_splineregression_Heatmap_H23.pdf"),
  width = 21/3,
  height = 14
)
#-------------------------------------------
# NCI-H23 2Gy |marking relevant clusters
#-------------------------------------------
spline_H23_clusters <- cutree(
  spline_H23_hm_Cluster_mark$tree_row,
  k = 2
)

spline_H23_row_order <- spline_H23_hm_Cluster_mark$tree_row$order

spline_H23_clusters_ordered <- spline_H23_clusters[spline_H23_row_order]

spline_H23_gaps_row <- which(
  diff(spline_H23_clusters_ordered) != 0
)

spline_H23_cluster_sizes<-table(as.numeric(spline_H23_clusters_ordered))

spline_H23_cluster_ids <- rle(spline_H23_clusters_ordered)$values

spline_H23_cluster_labels <- c("A", "B")
spline_H23_cluster_ranges <- data.frame(
  Cluster = c("A", "B"),
  Start = c(1,(sum(spline_cluster_sizes_H23[c(1)])+1)),
  End = c(sum(spline_H23_cluster_sizes[c(1)]),sum(spline_H23_cluster_sizes[c(1:2)])))

spline_H23_tree_pos <- which(spline_H23_hm_Cluster_mark$gtable$layout$name == "row_tree")

spline_H23_hm_Cluster_mark$gtable <- gtable_add_cols(
  spline_H23_hm_Cluster_mark$gtable,
  unit(12, "mm"),
  pos = max(spline_H23_hm_Cluster_mark$gtable$layout$l[spline_H23_tree_pos])-1
)

spline_H23_tree_pos <- which(spline_H23_hm_Cluster_mark$gtable$layout$name == "row_tree")
spline_H23_matrix_pos <- which(spline_H23_hm_Cluster_mark$gtable$layout$name == "matrix")
spline_H23_brackets <- grobTree(
  children = do.call(
    gList,
    lapply(seq_along(spline_H23_cluster_sizes), function(i) {
      n <- sum(spline_H23_cluster_sizes)
      y_top <- 1 - sum(spline_H23_cluster_sizes[seq_len(i - 1)]) / n
      y_bottom <- 1 - sum(spline_H23_cluster_sizes[seq_len(i)]) / n
      y_mid <- (y_top + y_bottom) / 2
      grobTree(
       segmentsGrob(
          x0 = unit(0.6, "npc"),
          x1 = unit(0.6, "npc"),
          y0 = unit(y_bottom, "npc"),
          y1 = unit(y_top, "npc"),
          gp = gpar(lwd = 2.5)
        ),
        segmentsGrob(
          x0 = unit(0.6, "npc"),
          x1 = unit(0.9, "npc"),
          y0 = unit(y_top, "npc"),
          y1 = unit(y_top, "npc"),
          gp = gpar(lwd = 2.5)
        ),
        segmentsGrob(
          x0 = unit(0.6, "npc"),
          x1 = unit(0.9, "npc"),
          y0 = unit(y_bottom, "npc"),
          y1 = unit(y_bottom, "npc"),
          gp = gpar(lwd = 2.5)
        ),
        textGrob(
          spline_H23_cluster_labels[i],
          x = unit(0.3, "npc"),
          y = unit(y_mid, "npc"),
          just = "right",
          gp = gpar(
            fontsize = 14,
            fontface = "bold"
          )
        )
      )
    })
  ))
  
spline_H23_hm_Cluster_mark$gtable <- gtable_add_grob(
  spline_H23_hm_Cluster_mark$gtable,
    grobs = spline_H23_brackets,
    t = spline_H23_hm_Cluster_mark$gtable$layout$t[spline_H23_matrix_pos],
    b = spline_H23_hm_Cluster_mark$gtable$layout$b[spline_H23_matrix_pos],
    l = spline_H23_hm_Cluster_mark$gtable$layout$l[spline_H23_tree_pos]-1,
    r = spline_H23_hm_Cluster_mark$gtable$layout$r[spline_H23_tree_pos]-1,
    name = "cluster_brackets"
  )
#-------------------------------------------
# NCI-H23 2Gy | saving results
#-------------------------------------------
pdf(
  file.path(
    dirs$results,
    "21_Proteome_splineRegression",
    "21_Proteome_H23_spline_2Gy_heatmap_2.pdf"
  ),
  width = 21*3/6,
  height = 14
)

grid.newpage()
grid.draw(spline_H23_hm_Cluster_mark$gtable)

dev.off()
#-------------------------------------------
# saving results
#-------------------------------------------
setwd(file.path(dirs$results, "21_Proteome_splineRegression"))

write.csv(
  spline_H23_labels,
  file = file.path(dirs$results, "21_Proteome_splineRegression", "NCI_H23_splineregression_Clusters_label.csv"),
  row.names = FALSE
)

write.csv(
  spline_results_H23_time,
  file = file.path(dirs$results, "21_Proteome_splineRegression", "NCI-H23_splineregression_SignifcantProteins.csv"),
  row.names = FALSE
)

write.csv(
  spline_HCC44_labels,
  file = file.path(dirs$results, "21_Proteome_splineRegression", "HCC44_splineregression_Clusters_label.csv"),
  row.names = FALSE
)

write.csv(
  spline_results_HCC44_time ,
  file = file.path(dirs$results, "21_Proteome_splineRegression", "HCC44_splineregression_SignifcantProteins.csv"),
  row.names = FALSE
) 

spline_HCC44_cluster3_assignment <- cutree(
  spline_hm_cluster_HCC44$tree_row,
  k = 3
)

spline_HCC44_clusters3 <- data.frame(
  SYMBOL = names(spline_HCC44_cluster3_assignment),
  Cluster = as.integer(spline_HCC44_cluster3_assignment)
) %>%
  dplyr::left_join(
    spline_sig_HCC44 %>%
      dplyr::select(ProteinKey,ID) %>%
      dplyr::distinct(ID, .keep_all = TRUE),
    by = c("SYMBOL"="ID")
  )

write.csv(
  spline_HCC44_clusters3 ,
  file = file.path(dirs$results, "21_Proteome_splineRegression", "HCC44_splineregression_SignifcantProteins_3Clusters.csv"),
  row.names = FALSE)
  
#-------------------------------------------
setwd(dirs$code)
