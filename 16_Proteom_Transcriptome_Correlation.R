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
compare_H23 <- merge(
  x = ready_pH23_24h[, c(2, 6, 9)] %>%
    dplyr::rename(
      Proteome = log2FoldChange,
      p_proteome = pval
    ),
  y = ready_RH23_24h[, c(1, 3, 6)] %>%
    dplyr::rename(
      Transcriptome = log2FoldChange,
      p_transcriptome = pval
    ),
  by.x = "PG.Genes",
  by.y = "GeneSymbol"
) %>%
  # filtering for genes with symbols and values
  filter(
    complete.cases(.),
    PG.Genes != "",
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
# HCC44 - binding proteome and transcriptome data
#-------------------------------------------
compare_HCC44 <- merge(
  x = ready_pHCC44_24h[, c(2, 6, 9)] %>%
    dplyr::rename(
      Proteome = log2FoldChange,
      p_proteome = pval
    ),
  y = ready_RHCC44_24h[, c(1, 3, 6)] %>%
    dplyr::rename(
      Transcriptome = log2FoldChange,
      p_transcriptome = pval
    ),
  by.x = "PG.Genes",
  by.y = "GeneSymbol"
) %>%
  # filtering for genes with symbols and values
  filter(
    complete.cases(.),
    PG.Genes != "",
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
corr_data_H23 <- compare_H23 |>
  dplyr::filter(
    is.finite(Transcriptome),
    is.finite(Proteome)
  )

n_pairs <- nrow(corr_data_H23)

r <- if (
  n_pairs >= 3 &&
  sd(corr_data_H23$Transcriptome) > 0 &&
  sd(corr_data_H23$Proteome) > 0
) {
  cor(
    corr_data_H23$Transcriptome,
    corr_data_H23$Proteome,
    method = "pearson"
  )
} else {
  NA_real_
}

cor_label_H23 <- paste0(
  "Pearson r = ", if (is.na(r)) "NA" else sprintf("%.2f", r),
  "\nn = ", n_pairs
)

corr_data_HCC44 <- compare_HCC44 |>
  dplyr::filter(
    is.finite(Transcriptome),
    is.finite(Proteome)
  )

n_pairs <- nrow(corr_data_HCC44)

r <- if (
  n_pairs >= 3 &&
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

cor_label_HCC44 <- paste0(
  "Pearson r = ", if (is.na(r)) "NA" else sprintf("%.2f", r),
  "\nn = ", n_pairs
)
#-------------------------------------------
# preparing for plotting
#-------------------------------------------
significance_colours <- c(
  "not significant" = "black",
  "significant in proteome" = "#5F8F29",
  "significant in transcriptome" = "#4294F7",
  "significant in both" ="#E84125"
)

#-------------------------------------------
# Plotting correlation
#-------------------------------------------
p_cor_proteome_transcriptome_HCC44 <- plot_correlation(
  data = compare_HCC44,
  title = "Comparison Proteome vs Transcriptome\nHCC44"
)

p_cor_proteome_transcriptome_HCC44 <- p_cor_proteome_transcriptome_HCC44 +
  annotate(
    "text",
    x = -5,
    y = Inf,
    label = cor_label_HCC44,
    hjust = 0.5,
    vjust = 1.2,
    size = 3
  )

p_cor_proteome_transcriptome_H23 <- plot_correlation(
  data = corr_data_H23,
  title = "Comparison Proteome vs Transcriptome\nNCI-H23"
)

p_cor_proteome_transcriptome_H23 <- p_cor_proteome_transcriptome_H23 +
  annotate(
    "text",
    x = -5,
    y = Inf,
    label = cor_label_H23,
    hjust = 0.5,
    vjust = 1.2,
    size = 3
  )  

#-------------------------------------------
# Combining plots
#-------------------------------------------
Correlation_Proteom_Transcriptome <-
  p_cor_proteome_transcriptome_HCC44 /
  p_cor_proteome_transcriptome_H23 +
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
  filename = "Correlation_Proteom_Transcriptome.jpeg",
  plot = Correlation_Proteom_Transcriptome,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = "Correlation_Proteom_Transcriptome.pdf",
  plot = Correlation_Proteom_Transcriptome,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#1. saving HCC44 results
ggsave(
  filename = "Correlation_Proteom_Transcriptome_HCC44.jpeg",
  plot = p_cor_proteome_transcriptome_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = "Correlation_Proteom_Transcriptome_HCC44.pdf",
  plot = p_cor_proteome_transcriptome_HCC44,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

#1. saving NCI-H23 results
ggsave(
  filename = "Correlation_Proteom_Transcriptome_NCI-H23.jpeg",
  plot = p_cor_proteome_transcriptome_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = "Correlation_Proteom_Transcriptome_NCI-H23.pdf",
  plot = p_cor_proteome_transcriptome_HCC44,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)
#-----------------------------------------------
setwd(dirs$code)
