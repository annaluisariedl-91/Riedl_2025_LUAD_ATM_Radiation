#-------------------------------------------
# loading or renaming data
#------------------------------------------
load_stat_data_R(
  file.path(dirs$results, "14_Transcriptome_DifferentialExpression"),
  mode = "ready"
)

load_stat_data_p(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  mode = "ready"
)
#-------------------------------------------
# NCI-H23 - binding proteome and transcriptome data
#-------------------------------------------
corr_H23 <- merge(
  x = ready_pH23_24h[, c(1,2, 3, 6,9)] %>%#[, c(2, 6, 9)] %>%
    dplyr::rename(
      Proteome = log2FoldChange,#
      p_proteome = pval
    ),
  y = ready_RH23_24h[, c(1, 3, 6)] %>%#[, c(1, 3, 6)] %>%
    dplyr::rename(
      Transcriptome = log2FoldChange,
      p_transcriptome = pval
    ),
  by.x = "SYMBOL",
  by.y = "GeneSymbol"
) %>%
  # filtering for genes with symbols and values
  filter(
    complete.cases(.),
    PG.Genes != "",
    !duplicated(ProteinKey)
  ) %>%
# sort for variance
  mutate(
    significance = case_when(
      p_proteome < 0.05 &
        p_transcriptome < 0.05 ~
        "significant in both",
      p_proteome < 0.05 ~
        "significant in proteome",
      p_transcriptome < 0.05 ~
        "significant in transcriptome",
      TRUE ~
        "not significant"
    ),
    significance = factor(
      significance,
      levels = c(
        "not significant",
        "significant in proteome",
        "significant in transcriptome",
        "significant in both"
      )
    )
  )   
#-------------------------------------------
# HCC44 - binding proteome and transcriptome data
#-------------------------------------------
corr_HCC44 <- merge(
  x = ready_pHCC44_24h[,c(2, 6, 9)] %>%
    dplyr::rename(
      Proteome = log2FoldChange,
      p_proteome = pval
    ),
  y = ready_RHCC44_24h[, c(1, 3, 6)] %>%
    dplyr::rename(
      Transcriptome = log2FoldChange,
      p_transcriptome = pval
    ),
  by.x = "SYMBOL",
  by.y = "GeneSymbol"
) %>%
  # filtering for genes with symbols and values
  filter(
    complete.cases(.),
    SYMBOL != "",
    !duplicated(Proteome)
  ) %>%
# sort for variance
  mutate(
    significance = case_when(
      p_proteome < 0.05 &
        p_transcriptome < 0.05 ~
        "significant in both",
      p_proteome < 0.05 ~
        "significant in proteome",
      p_transcriptome < 0.05 ~
        "significant in transcriptome",
      TRUE ~
        "not significant"
    ),
    significance = factor(
      significance,
      levels = c(
        "not significant",
        "significant in proteome",
        "significant in transcriptome",
        "significant in both"
      )
    )
  )
#-------------------------------------------
# Correlation
#-------------------------------------------
corr_data_H23 <- corr_H23 |>
  dplyr::filter(
    is.finite(Transcriptome),
    is.finite(Proteome)
  )

corr_n_pairs_H23 <- nrow(corr_data_H23)

corr_H23_r <- if (
  corr_n_pairs_H23 >= 3 &&
  sd(corr_data_H23$Transcriptome, na.rm = TRUE) > 0 &&
  sd(corr_data_H23$Proteome, na.rm = TRUE) > 0
) {
  cor(
    corr_data_H23$Transcriptome,
    corr_data_H23$Proteome,
    method = "pearson",
    use = "complete.obs"
  )
} else {
  NA_real_
}

corr_label_H23 <- paste0(
  "Pearson r = ", if (is.na(corr_H23_r)) "NA" else sprintf("%.2f", corr_H23_r),
  "\nn = ", corr_n_pairs_H23
)

corr_data_HCC44 <- corr_HCC44 |>
  dplyr::filter(
    is.finite(Transcriptome),
    is.finite(Proteome)
  )

corr_n_pairs_HCC44 <- nrow(corr_data_HCC44)

corr_HCC44_r <- if (
  corr_n_pairs_HCC44 >= 3 &&
  sd(corr_data_HCC44$Transcriptome) > 0 &&
  sd(corr_data_HCC44$Proteome) > 0
) {
  cor(
    corr_data_HCC44$Transcriptome,
    corr_data_HCC44$Proteome,
    method = "pearson"
  )
} else {
  NA_real_
}

corr_label_HCC44 <- paste0(
  "Pearson r = ", if (is.na(corr_HCC44_r)) "NA" else sprintf("%.2f", corr_HCC44_r),
  "\nn = ", corr_n_pairs_HCC44
)
#-------------------------------------------
# preparing for plotting
#-------------------------------------------
corr_colours <- c(
  "not significant" = "black",
  "significant in proteome" = "#5F8F29",
  "significant in transcriptome" = "#4294F7",
  "significant in both" ="#E84125"
)
#-------------------------------------------
# Plotting correlation
#-------------------------------------------
corr_p_proteome_transcriptome_HCC44 <- plot_correlation(
  data = corr_HCC44,
  title = "Comparison Proteome vs Transcriptome\nHCC44"
)

corr_p_proteome_transcriptome_HCC44 <- corr_p_proteome_transcriptome_HCC44 +
  annotate(
    "text",
    x = -7.5,
    y = 0.5,
    label = corr_label_HCC44,
    hjust = 0,
    vjust = 0,
    size = 3
  )

corr_p_proteome_transcriptome_H23 <- plot_correlation(
  data = corr_data_H23,
  title = "Comparison Proteome vs Transcriptome\nNCI-H23"
)

corr_p_proteome_transcriptome_H23 <- corr_p_proteome_transcriptome_H23 +
  annotate(
    "text",
    x = -7.5,
    y = 0.5,
    label = corr_label_H23,
    hjust = 0,
    vjust = 0,
    size = 3
  )  

#-------------------------------------------
# Combining plots
#-------------------------------------------
corr_p_proteome_transcriptome_final <-
  corr_p_proteome_transcriptome_HCC44 /
  corr_p_proteome_transcriptome_H23 +
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )
#-------------------------------------------
# Saving results
#-------------------------------------------
dir.create(
  file.path(dirs$results, "16_Proteom_Transcriptome_Correlation"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "16_Proteom_Transcriptome_Correlation"))

#1. saving combined results
ggsave(
  filename = "16_Proteom_Transcriptome_Correlation_plot.jpeg",
  plot = corr_p_proteome_transcriptome_final,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = "16_Proteom_Transcriptome_Correlation_plot.pdf",
  plot = corr_p_proteome_transcriptome_final,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#2. saving HCC44 results
ggsave(
  filename = "16_Proteom_Transcriptome_Correlation_plot_HCC44.jpeg",
  plot = corr_p_proteome_transcriptome_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = "16_Proteom_Transcriptome_Correlation_plot_HCC44.pdf",
  plot = corr_p_proteome_transcriptome_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#3. saving NCI-H23 results
ggsave(
  filename = "16_Proteom_Transcriptome_Correlation_plot_NCI-H23.jpeg",
  plot = corr_p_proteome_transcriptome_H23,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = "16_Proteom_Transcriptome_Correlation_plot_NCI-H23.pdf",
  plot = corr_p_proteome_transcriptome_H23,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#4. saving NCI-H23 results table
write.csv(
  corr_H23,
  file = file.path(dirs$results, "16_Proteom_Transcriptome_Correlation", "16_Proteom_Transcriptome_Correlation_plot_NCI_H23_results.csv"),
  row.names = FALSE
)


#5. saving HCC44 results table
write.csv(
  corr_HCC44,
  file = file.path(dirs$results, "16_Proteom_Transcriptome_Correlation", "16_Proteom_Transcriptome_Correlation_plot_HCC44_results.csv"),
  row.names = FALSE
)
#-----------------------------------------------
setwd(dirs$code)
