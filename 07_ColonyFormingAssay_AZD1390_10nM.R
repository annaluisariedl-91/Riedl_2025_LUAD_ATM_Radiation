#-------------------------------------------
# set directory
#-------------------------------------------
setwd(file.path(base_dir,"functional_assays","CFA","Exp4_ATMinhibition","raw")) #replace with path on your own PC
#dose <-c(0,2,4,8)
CFA_AZD1930_10nM_A549<- read.csv("CFA_AZD1390_A549.csv",sep=",")
CFA_AZD1930_10nM_HCC827<- read.csv("CFA_AZD1390_HCC827.csv",sep=",")
CFA_AZD1930_10nM_HCC78<- read.csv("CFA_AZD1390_HCC78.csv",sep=",")
#-------------------------------------------
# importing A549
#-------------------------------------------
rep_cols <- list(
  2:4,
  7:9,
  12:14,
  17:19
)

CFA_AZD1930_10nM_A549_DMSO_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      CFA_AZD1930_10nM_A549,
      rows = 19:30,
      cols = c(cols)
    )
  }
)

CFA_AZD1930_10nM_A549_10nM_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      CFA_AZD1930_10nM_A549,
      rows = 3:14,
      cols = cols
    )
  }
)
#-------------------------------------------
# importing HCC827
#-------------------------------------------
rep_cols_HCC827 <- list(
  2:4,
  7:9,
  12:14
)

CFA_AZD1930_10nM_HCC827_DMSO_rep <- lapply(
  rep_cols_HCC827,
  function(cols) {
    analyse_CFA_rep(
      CFA_AZD1930_10nM_HCC827,
      rows = 19:30,
      cols = cols
    )
  }
)

CFA_AZD1930_10nM_HCC827_10nM_rep <- lapply(
  rep_cols_HCC827,
  function(cols) {
    analyse_CFA_rep(
      CFA_AZD1930_10nM_HCC827,
      rows = 3:14,
      cols = cols
    )
  }
)
#-------------------------------------------
# importing HCC78
#-------------------------------------------
CFA_AZD1930_10nM_HCC78_DMSO_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      CFA_AZD1930_10nM_HCC78,
      rows = 19:30,
      cols = cols
    )
  }
)

CFA_AZD1930_10nM_HCC78_10nM_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      CFA_AZD1930_10nM_HCC78,
      rows = 3:14,
      cols = cols
    )
  }
)
#-------------------------------------------
# Binding data
#-------------------------------------------
CFA_AZD1930_10nM_data <- bind_rows(
  # A549
  make_CFA_data(
    CFA_AZD1930_10nM_A549_DMSO_rep,
    "A549",
    "DMSO"
  ),
  make_CFA_data(
    CFA_AZD1930_10nM_A549_10nM_rep,
    "A549",
    "10 nM AZD1390"
  ),
  # HCC827
  make_CFA_data(
    CFA_AZD1930_10nM_HCC827_DMSO_rep,
    "HCC827",
    "DMSO"
  ),
  make_CFA_data(
    CFA_AZD1930_10nM_HCC827_10nM_rep,
    "HCC827",
    "10 nM AZD1390"
  ),
  # HCC78
  make_CFA_data(
    CFA_AZD1930_10nM_HCC78_DMSO_rep,
    "HCC78",
    "DMSO"
  ),
  make_CFA_data(
    CFA_AZD1930_10nM_HCC78_10nM_rep,
    "HCC78",
    "10 nM AZD1390"
  )
)

CFA_AZD1930_10nM_data <- CFA_AZD1930_10nM_data %>%
  mutate(
    Cellline = factor(Cellline),
      Treatment = factor(
      Treatment,
      levels = c(
        "DMSO",
        "10 nM AZD1390"
      )
    ),
    Experiment = factor(experiment),
    adj_sf = pmax(SF, 1e-6),
    adj_log_sf = log(adj_sf)
  )

CFA_AZD1930_10nM_mean_stat <- CFA_AZD1930_10nM_data %>%
  summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    n       = sum(!is.na(SF)),
    adj_log_sf = mean(adj_sf, na.rm = TRUE),
    .by = c(Cellline, Treatment, Dose)
  )
#-------------------------------------------
# global models
#-------------------------------------------
CFA_AZD1930_10nM_model_global <- lm(
  adj_log_sf ~ Experiment  +
    factor(Dose) * Treatment * Cellline,
  data = CFA_AZD1930_10nM_data
)

CFA_AZD1930_10nM_res_aov_global<-anova(CFA_AZD1930_10nM_model_global)

CFA_AZD1930_10nM_emm_global<- emmeans(
  CFA_AZD1930_10nM_model_global,
  ~ Treatment | Dose + Cellline
)

CFA_AZD1930_10nM_posthoc_global  <- pairs(
  CFA_AZD1930_10nM_emm_global,
  adjust = "tukey"
)
#-------------------------------------------
# linear-quadratic model
#-------------------------------------------
CFA_AZD1930_10nM_stat_lq <- CFA_AZD1930_10nM_mean_stat %>%
  mutate(
    adj_sf = pmax(SF_mean, 1e-6),
    adj_log_sf = log(adj_sf)
  )

CFA_AZD1930_10nM_LQ_A549<- fit_LQ_treatments(
  CFA_AZD1930_10nM_stat_lq , "A549"
)

CFA_AZD1930_10nM_LQ_HCC78 <- fit_LQ_treatments(
  CFA_AZD1930_10nM_stat_lq , "HCC78"
)

CFA_AZD1930_10nM_LQ_HCC827  <- fit_LQ_treatments(
  CFA_AZD1930_10nM_stat_lq , "HCC827"
)

CFA_AZD1930_10nM_coefficients_all  <- bind_rows(A549 = as.data.frame(CFA_AZD1930_10nM_LQ_A549 $coefficients) %>%
                         tibble::rownames_to_column("coefficient"),
                       HCC78  = as.data.frame(CFA_AZD1930_10nM_LQ_HCC78 $coefficients) %>%
                         tibble::rownames_to_column("coefficient"),
                       HCC827 = as.data.frame(CFA_AZD1930_10nM_LQ_HCC827 $coefficients) %>%
                         tibble::rownames_to_column("coefficient"),
                       .id = "Cellline")
#-------------------------------------------
# Dose enhancement factor (DEF)
#-------------------------------------------
CFA_AZD1930_10nM_coefficients_A549   <- CFA_AZD1930_10nM_LQ_A549$coefficients
CFA_AZD1930_10nM_coefficients_HCC78  <- CFA_AZD1930_10nM_LQ_HCC78$coefficients
CFA_AZD1930_10nM_coefficients_HCC827 <- CFA_AZD1930_10nM_LQ_HCC827$coefficients

CFA_AZD1930_10nM_DEF_results_A549 <- calc_DEF(
  CFA_AZD1930_10nM_coefficients_A549
)

CFA_AZD1930_10nM_DEF_label_A549 <- paste0(
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    CFA_AZD1930_10nM_DEF_results_A549$DEF_10nM[1],
    1
  )
)

CFA_AZD1930_10nM_DEF_results_HCC78 <- calc_DEF(
  CFA_AZD1930_10nM_coefficients_HCC78
)

CFA_AZD1930_10nM_DEF_label_HCC78<- paste0(
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    CFA_AZD1930_10nM_DEF_results_HCC78$DEF_10nM[1],
    1
  )
)

CFA_AZD1930_10nM_DEF_results_HCC827 <- calc_DEF(
  CFA_AZD1930_10nM_coefficients_HCC827
)

CFA_AZD1930_10nM_DEF_label_HCC827 <- paste0(
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    CFA_AZD1930_10nM_DEF_results_HCC827$DEF_10nM[1],
    1
  )
)

CFA_AZD1930_10nM_DEF_all <- bind_rows(  A549   = as.data.frame(CFA_AZD1930_10nM_DEF_results_A549) %>%
                         tibble::rownames_to_column("Effect"),
                       HCC78  = as.data.frame(CFA_AZD1930_10nM_DEF_results_HCC78) %>%
                         tibble::rownames_to_column("Effect"),
                       HCC827 = as.data.frame(CFA_AZD1930_10nM_DEF_results_HCC827) %>%
                         tibble::rownames_to_column("Effect"),
                       .id = "Cellline")
#-------------------------------------------
# Anova
#-------------------------------------------
CFA_AZD1930_10nM_A549 <- dplyr::bind_rows(
  lapply(seq_along(CFA_AZD1930_10nM_A549_DMSO_rep), function(i) {
    CFA_AZD1930_10nM_A549_DMSO_rep[[i]]$summary %>%
      dplyr::mutate(
        Experiment = i,
        Treatment = "DMSO",
        SF = N_mean / (N0 * CFA_AZD1930_10nM_A549_DMSO_rep[[i]]$PE)
      )
  }),
  lapply(seq_along(CFA_AZD1930_10nM_A549_10nM_rep), function(i) {
    CFA_AZD1930_10nM_A549_10nM_rep[[i]]$summary %>%
      dplyr::mutate(
        Experiment = i,
        Treatment = "10 nM AZD1390",
        SF = N_mean / (N0 * CFA_AZD1930_10nM_A549_10nM_rep[[i]]$PE)
      )
  })
)

CFA_AZD1930_10nM_HCC78 <- dplyr::bind_rows(
  lapply(seq_along(CFA_AZD1930_10nM_HCC78_DMSO_rep), function(i) {
    CFA_AZD1930_10nM_HCC78_DMSO_rep[[i]]$summary %>%
      dplyr::mutate(
        Experiment = i,
        Treatment = "DMSO",
        SF = N_mean / (N0 * CFA_AZD1930_10nM_HCC78_DMSO_rep[[i]]$PE)
      )
  }),
  lapply(seq_along(CFA_AZD1930_10nM_HCC78_10nM_rep), function(i) {
    CFA_AZD1930_10nM_HCC78_10nM_rep[[i]]$summary %>%
      dplyr::mutate(
        Experiment = i,
        Treatment = "10 nM AZD1390",
        SF = N_mean / (N0 * CFA_AZD1930_10nM_HCC78_10nM_rep[[i]]$PE)
      )
  })
)

CFA_AZD1930_10nM_HCC827 <- dplyr::bind_rows(
  lapply(seq_along(CFA_AZD1930_10nM_HCC827_DMSO_rep), function(i) {
    CFA_AZD1930_10nM_HCC827_DMSO_rep[[i]]$summary %>%
      dplyr::mutate(
        Experiment = i,
        Treatment = "DMSO",
        SF = N_mean / (N0 * CFA_AZD1930_10nM_HCC827_DMSO_rep[[i]]$PE)
      )
  }),
  lapply(seq_along(CFA_AZD1930_10nM_HCC827_10nM_rep), function(i) {
    CFA_AZD1930_10nM_HCC827_10nM_rep[[i]]$summary %>%
      dplyr::mutate(
        Experiment = i,
        Treatment = "10 nM AZD1390",
        SF = N_mean / (N0 * CFA_AZD1930_10nM_HCC827_10nM_rep[[i]]$PE)
      )
  })
)

CFA_AZD1930_10nM_ANOVA_summary <- dplyr::bind_rows(
  lapply(
    list(
      ANOVA_A549 = aov(SF ~ factor(Dose) * Treatment, data=CFA_AZD1930_10nM_A549),
      ANOVA_HCC827 = aov(SF ~ factor(Dose) * Treatment, data=CFA_AZD1930_10nM_HCC827),
      ANOVA_HCC78 = aov(SF ~ factor(Dose) * Treatment, data=CFA_AZD1930_10nM_HCC78)
    ),
    function(model) {
      result <- as.data.frame(anova(model))
      result$Effect <- rownames(result)
      rownames(result) <- NULL
      result
    }
  ),
  .id = "Model"
) %>%
  dplyr::mutate(
    Cellline = gsub("^ANOVA_", "", Model)
  ) %>%
  dplyr::relocate(Cellline, Model, Effect)
#--------------------------------------------
# prepare plotting
#--------------------------------------------
CFA_AZD_colours <- c(
  "DMSO"          = "#000000",
  "10 nM AZD1390" = "#F09938"
)

CFA_AZD_shapes <- c(
  "DMSO"          = 21,
  "10 nM AZD1390" = 23
)

CFA_AZD1930_10nM_ANOVA_significance <- CFA_AZD1930_10nM_ANOVA_summary %>%
  dplyr::filter(Effect == "factor(Dose):Treatment") %>%
  dplyr::mutate(
    significance = dplyr::case_when(
      `Pr(>F)` < 0.0001 ~ "****",
      `Pr(>F)` < 0.001  ~ "***",
      `Pr(>F)` < 0.01   ~ "**",
      `Pr(>F)` < 0.05   ~ "*",
      TRUE              ~ "ns"
    )
  ) %>%
  dplyr::select(Model, Cellline, `Pr(>F)`, significance)%>%
  dplyr::left_join(
    data.frame(
      Cellline = c("A549", "HCC78","HCC827"),
      label = "2Gy:10nM AZD1390"
    ),
    by = c("Cellline")
  ) %>%
  dplyr::select(Cellline, Model, `Pr(>F)`,significance, label)
#--------------------------------------------
# Plotting figures
#--------------------------------------------
CFA_AZD1930_10nM_p_A549_10nM <- plot_CFA_AZD(
  plot_data = CFA_AZD1930_10nM_mean_stat,
  coefficients = CFA_AZD1930_10nM_coefficients_A549,
  DEF_label = data.frame(CFA_AZD1930_10nM_DEF_label_A549),
  selected_cellline = "A549",
  xmax = 8
)+
  ggplot2::geom_text(
    data = CFA_AZD1930_10nM_ANOVA_significance %>%
      dplyr::filter(Cellline == "A549", label != "ns"),
    ggplot2::aes(x = 2.5, y = 9, label = paste(label,significance)),
    inherit.aes = FALSE,
    size = 2.5
  )

CFA_AZD1930_10nM_p_HCC78_10nM <- plot_CFA_AZD(
  plot_data = CFA_AZD1930_10nM_mean_stat,
  coefficients = CFA_AZD1930_10nM_coefficients_HCC78,
  DEF_label =  data.frame(CFA_AZD1930_10nM_DEF_label_HCC78),
  selected_cellline = "HCC78",
  xmax = 8
)+
  ggplot2::geom_text(
    data = CFA_AZD1930_10nM_ANOVA_significance %>%
      dplyr::filter(Cellline == "HCC78", label != "ns"),
    ggplot2::aes(x = 2.5, y = 9, label = paste(label,significance)),
    inherit.aes = FALSE,
    size = 2.5
  )

CFA_AZD1930_10nM_p_HCC827_10nM <- plot_CFA_AZD(
  plot_data = CFA_AZD1930_10nM_mean_stat,
  coefficients = CFA_AZD1930_10nM_coefficients_HCC827,
  DEF_label =  data.frame(CFA_AZD1930_10nM_DEF_label_HCC827),
  selected_cellline = "HCC827",
  xmax = 8
)+
  ggplot2::geom_text(
    data = CFA_AZD1930_10nM_ANOVA_significance %>%
      dplyr::filter(Cellline == "HCC827", label != "ns"),
    ggplot2::aes(x = 2.5, y = 9, label = paste(label,significance)),
    inherit.aes = FALSE,
    size = 2.5
  )

CFA_AZD1930_10nM_p_CFA_final <- (
  CFA_AZD1930_10nM_p_A549_10nM |
    CFA_AZD1930_10nM_p_HCC78_10nM |
    CFA_AZD1930_10nM_p_HCC827_10nM
) +
  patchwork::plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )
#--------------------------------------------
# saving results
#--------------------------------------------
dir.create(
  file.path(dirs$results, "07_ColonyFormingAssay_AZD1390_10nM_AZD1390"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "07_ColonyFormingAssay_AZD1390_10nM_AZD1390"))
         
# 1. Mean, SD and SEM of CFA
write.csv(
  CFA_AZD1930_10nM_mean_stat,
  file = file.path(dirs$results, "07_ColonyFormingAssay_AZD1390_10nM_AZD1390", "CFA_AZD1390_summary.csv"),
  row.names = FALSE
)

# 2. LQ coefficient 
write.csv(
  CFA_AZD1930_10nM_coefficients_all,
  file = file.path(dirs$results, "07_ColonyFormingAssay_AZD1390_10nM", "CFA_AZD1390_LQ_all.csv"),
  row.names = FALSE
)

# 3. DEF - Dose enhancement factor
write.csv(
  CFA_AZD1930_10nM_DEF_all,
  file = file.path(dirs$results, "07_ColonyFormingAssay_AZD1390_10nM", "CFA_AZD1390_DEF_all.csv"),
  row.names = FALSE
)

# 4. ANOVA results
write.csv(
  CFA_AZD1930_10nM_ANOVA_summary,
  file = file.path(dirs$results, "07_ColonyFormingAssay_AZD1390_10nM", "CFA_AZD1390_ANOVA.csv"),
  row.names = FALSE
)

# 5. global ANOVA results
write.csv(
  CFA_AZD1930_10nM_posthoc_global,
  file = file.path(dirs$results, "07_ColonyFormingAssay_AZD1390_10nM", "CFA_AZD1390_ANOVA_global.csv"),
  row.names = FALSE
)

# 6. saving combined Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "CFA_AZD1390_10nM_all.jpeg"
  ),
  plot = CFA_AZD1930_10nM_p_CFA_final,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "CFA_AZD1390_10nM_all.pdf"
  ),
  plot = CFA_AZD1930_10nM_p_CFA_final,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 7. saving A549 Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "07_CFA_AZD1390_10nM_A549.jpeg"
  ),
  plot = CFA_AZD1930_10nM_p_A549_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "07_CFA_AZD1390_10nM_A549.pdf"
  ),
  plot = CFA_AZD1930_10nM_p_A549_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)
 
# 9. saving HCC78 Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "CFA_AZD1390_10nM_HCC78.jpeg"
  ),
  plot = CFA_AZD1930_10nM_p_HCC78_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "CFA_AZD1390_10nM_HCC78.pdf"
  ),
  plot = CFA_AZD1930_10nM_p_HCC78_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 10. saving HCC827 Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "CFA_AZD1390_10nM_HCC827.jpeg"
  ),
  plot = CFA_AZD1930_10nM_p_HCC827_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_ColonyFormingAssay_AZD1390_10nM",
    "CFA_AZD1390_10nM_HCC827.pdf"
  ),
  plot = CFA_AZD1930_10nM_p_HCC827_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)
#------------------------------------------
setwd(dirs$code)
