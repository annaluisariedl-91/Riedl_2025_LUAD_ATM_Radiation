#------------------------------------------
# #set directory
#------------------------------------------
setwd(file.path(dirs$CFA,"Exp3_ATMinhibition","raw","HCC44")) #replace with path on your own PC

CFA_AZD1390_HCC44_DMSO_rep <- lapply(
  paste0("CFA_HCC44_DMSO_Rep", 1:2, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

CFA_AZD1390_HCC44_1nM_rep <- lapply(
  paste0("CFA_HCC44_1nM_AZD1390_Rep", 1:2, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

CFA_AZD1390_HCC44_10nM_rep <- lapply(
  paste0("CFA_HCC44_10nM_AZD1390_Rep", 1:2, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

setwd(file.path(base_dir,"functional_assays","CFA","Exp3_ATMinhibition","raw","NCI-H23")) #replace with path on your own PC

CFA_AZD1390_H23_DMSO_rep <- lapply(
  paste0("CFA_H23_DMSO_Rep", 1:3, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

CFA_AZD1390_H23_1nM_rep <- lapply(
  paste0("CFA_H23_1nM_AZD1390_Rep", 1:3, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

CFA_AZD1390_H23_10nM_rep <- lapply(
  paste0("CFA_H23_10nM_AZD1390_Rep", 1:3, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)
#------------------------------------------
# binding data and calculation
#------------------------------------------
CFA_AZD1390_data <- bind_rows(
  
  make_CFA_data(CFA_AZD1390_HCC44_DMSO_rep, "HCC44", "DMSO"),
  make_CFA_data(CFA_AZD1390_HCC44_1nM_rep,  "HCC44", "1 nM AZD1390"),
  make_CFA_data(CFA_AZD1390_HCC44_10nM_rep, "HCC44", "10 nM AZD1390"),
  
  make_CFA_data(CFA_AZD1390_H23_DMSO_rep, "NCI-H23", "DMSO"),
  make_CFA_data(CFA_AZD1390_H23_1nM_rep,  "NCI-H23", "1 nM AZD1390"),
  make_CFA_data(CFA_AZD1390_H23_10nM_rep, "NCI-H23", "10 nM AZD1390")
)

CFA_AZD1390_summary <- CFA_AZD1390_data %>%
  mutate(
    Dose = as.numeric(as.character(Dose))
  ) %>%
  group_by(Cellline, Treatment, Dose) %>%
  summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    adj_sf = pmax(SF_mean, 1e-6),
    adj_log_sf = log(adj_sf),
    n       = sum(!is.na(SF)),
    .groups = "drop"
  )
#------------------------------------------
# binding for statistics
#------------------------------------------
CFA_AZD1390_data_stat_all <- CFA_AZD1390_data %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("DMSO", "1 nM AZD1390", "10 nM AZD1390")
    ),
    Cellline = factor(Cellline),
    Experiment = factor(experiment),
    adj_sf = pmax(SF, 1e-6),
    adj_log_sf = log(adj_sf),
    Dose = as.numeric(Dose)
  )
#------------------------------------------
# linear-quadratic modell
#------------------------------------------
CFA_AZD1390_LQ_H23 <- fit_LQ_treatments(
  CFA_AZD1390_data_stat_all,
  "NCI-H23"
)

CFA_AZD1390_coefficients_H23 <- CFA_AZD1390_LQ_H23$coefficients

CFA_AZD1390_LQ_HCC44 <- fit_LQ_treatments(
  CFA_AZD1390_data_stat_all,
  "HCC44"
)

CFA_AZD1390_coefficients_HCC44 <- CFA_AZD1390_LQ_HCC44$coefficients
#-------------------------------------------
# ANOVA
#-------------------------------------------
# HCC44
# DMSO as Reference
CFA_AZD1390_data_stat_all$Treatment <- relevel(
  factor(CFA_AZD1390_data_stat_all$Treatment),
  ref = "DMSO",
  factor(CFA_AZD1390_data_stat_all$Dose),
)

# Model
CFA_AZD1390_model_HCC44_lm <- lm(
  adj_log_sf ~  Dose * Treatment,
  data = CFA_AZD1390_data_stat_all[CFA_AZD1390_data_stat_all$Cellline=="HCC44",]
)

CFA_AZD1390_emm_HCC44 <- emmeans(
  CFA_AZD1390_model_HCC44_lm,
  ~ Treatment | Dose 
)

# 1 nM and 10 nM vs DMSO
CFA_AZD1390_posthoc_HCC44 <- contrast(
  CFA_AZD1390_emm_HCC44,
  method = "trt.vs.ctrl",
  ref = "DMSO",
  adjust = "dunnett"
)

# Result HCC44
CFA_AZD1390_posthoc_HCC44_df <- as.data.frame(CFA_AZD1390_posthoc_HCC44)

# NCI-H23
# DMSO as Reference
CFA_AZD1390_model_H23_lm <- lm(
  adj_log_sf ~ Experiment + Dose * Treatment,
  data = CFA_AZD1390_data_stat_all[CFA_AZD1390_data_stat_all$Cellline=="NCI-H23",]
)

# Model
CFA_AZD1390_emm_H23 <- emmeans(
  CFA_AZD1390_model_H23_lm,
  ~ Treatment | Dose 
)

# 1 nM and 10 nM vs DMSO
CFA_AZD1390_posthoc_H23 <- contrast(
  CFA_AZD1390_emm_H23,
  method = "trt.vs.ctrl",
  ref = "DMSO",
  adjust = "dunnett"
)

CFA_AZD1390_posthoc_H23_df <- as.data.frame(CFA_AZD1390_posthoc_H23)

# Results
CFA_AZD1390_posthoc<-rbind(CFA_AZD1390_posthoc_HCC44_df,CFA_AZD1390_posthoc_H23_df)
CFA_AZD1390_posthoc$Cellline <- c(rep("HCC44",2),rep("NCI-H23",2))
#------------------------------------------
# Dose Enhancement Factor DEF
#------------------------------------------
CFA_AZD1390_DEF_results_H23 <- calculate_DEF(
  CFA_AZD1390_coefficients_H23
)

CFA_AZD1390_DEF_label_H23<- paste0(
    "1 nM AZD1390: DEF[0.1] = ",
    round(
      CFA_AZD1390_DEF_results_H23$DEF_1nM[1],
      1
    ),
    "\n",
    "10 nM AZD1390: DEF[0.1] = ",
    round(
      CFA_AZD1390_DEF_results_H23$DEF_10nM[1],
      1
  )
)

CFA_AZD1390_DEF_results_HCC44 <- calculate_DEF(
  CFA_AZD1390_coefficients_HCC44
)

CFA_AZD1390_DEF_label_HCC44 <- paste0(
  "1 nM AZD1390: DEF[0.1] = ",
  round(
    CFA_AZD1390_DEF_results_HCC44$DEF_1nM[1],
    1
  ),
  "\n",
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    CFA_AZD1390_DEF_results_HCC44$DEF_10nM[1],
    1
  )
)
#------------------------------------------
# Preparing graphs
#------------------------------------------
CFA_AZD_colours <- c(
  "DMSO"          = "#000000",
  "1 nM AZD1390"  = "#4294F7",
  "10 nM AZD1390" = "#F09938"
)

CFA_AZD_shapes <- c(
  "DMSO"          = 21,
  "1 nM AZD1390"  = 23,
  "10 nM AZD1390" = 22
)

CFA_AZD1390_significance <- as.data.frame(CFA_AZD1390_posthoc) %>%
  dplyr::mutate(
    label = dplyr::case_when(
      p.value < 0.0001 ~ "****",
      p.value < 0.001  ~ "***",
      p.value < 0.01   ~ "**",
      p.value < 0.05   ~ "*",
      TRUE             ~ "ns"
    ),
    Treatment = trimws(sub("- .*$", "", contrast))
  ) %>%
  dplyr::select(
    Cellline,
    Treatment,
    contrast,
    p.value,
    label
  ) %>%
  dplyr::left_join(
    CFA_AZD1390_summary %>%
      dplyr::filter(Dose == 2 & Treatment!="DMSO") %>%
      dplyr::select(
        Cellline,
        Treatment,
        y.pos = SF_mean
      ),
    by = c("Cellline", "Treatment")
  )
#------------------------------------------
# Plotting figure
#------------------------------------------
CFA_AZD1390_p_HCC44 <- plot_CFA_AZD(
  plot_data = CFA_AZD1390_summary,
  coefficients = CFA_AZD1390_coefficients_HCC44,
  DEF_label = data.frame(CFA_AZD1390_DEF_label_HCC44),
  selected_cellline = "HCC44",
  xmax = 8
)+
  ggplot2::geom_text(
    data = CFA_AZD1390_significance %>%
      dplyr::filter(Cellline == "HCC44"),
    ggplot2::aes(x =0.2, y =c(20,5), label =paste(contrast, label),hjust = 0, vjust = 0.5,),
    inherit.aes = FALSE,
    size = 2.25
  )

CFA_AZD1390_p_NCI_H23 <- plot_CFA_AZD(
  plot_data = CFA_AZD1390_summary,
  coefficients = CFA_AZD1390_coefficients_H23,
  DEF_label = data.frame(CFA_AZD1390_DEF_label_H23),
  selected_cellline = "NCI-H23",
  xmax = 4
)+
  ggplot2::geom_text(
    data = CFA_AZD1390_significance %>%
      dplyr::filter(Cellline == "NCI-H23"),
    ggplot2::aes(x =0.2, y =c(20,5), label =paste(contrast, label)),
    inherit.aes = FALSE,
    size = 2.25
  )

CFA_AZD1390_p_HCC44_H23_final <- (
  CFA_AZD1390_p_HCC44 |
    CFA_AZD1390_p_NCI_H23
) +
  patchwork::plot_layout(
    guides = "collect"
  )&
  theme(
    legend.position = "bottom"
  )
#------------------------------------------
# saving results
#------------------------------------------
dir.create(
  file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23"))

# 1. Mean, SD and SEM of CFA
write.csv(
  CFA_AZD1390_summary,
  file = file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23", "06_CFA_AZD1390_summary.csv"),
  row.names = FALSE
)

# 2. LQ coefficient HCC44
CFA_AZD1390_LQ_results_HCC44 <- do.call(
  rbind,
  lapply(names(CFA_AZD1390_coefficients_HCC44), function(treatment) {
    data.frame(
      Treatment = treatment,
      alpha = unname(CFA_AZD1390_coefficients_HCC44[[treatment]]["alpha"]),
      beta  = unname(CFA_AZD1390_coefficients_HCC44[[treatment]]["beta"])
    )
  })
)

write.csv(
  CFA_AZD1390_LQ_results_HCC44,
  file = file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23", "06_CFA_AZD1390_LQ_HCC44.csv"),
  row.names = FALSE
)

# 3. LQ coefficient NCI-H23
CFA_AZD1390_LQ_results_H23<- do.call(
  rbind,
  lapply(names(CFA_AZD1390_coefficients_H23), function(treatment) {
    data.frame(
      Treatment = treatment,
      alpha = unname(CFA_AZD1390_coefficients_H23[[treatment]]["alpha"]),
      beta  = unname(CFA_AZD1390_coefficients_H23[[treatment]]["beta"])
    )
  })
)

write.csv(
  CFA_AZD1390_LQ_results_H23,
  file = file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23", "06_CFA_AZD1390_LQ_H23.csv"),
  row.names = FALSE
)

# 4. DEF - Dose enhancement factor
write.csv(
  CFA_AZD1390_DEF_results_HCC44,
  file = file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23", "06_CFA_AZD1390_DEF_HCC44.csv"),
  row.names = FALSE
)

write.csv(
  CFA_AZD1390_DEF_results_H23,
  file = file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23", "06_CFA_AZD1390_DEF_H23.csv"),
  row.names = FALSE
)

# 5. ANOVA results
write.csv(
  CFA_AZD1390_posthoc,
  file = file.path(dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23", "06_CFA_AZD1390_ANOVA.csv"),
  row.names = FALSE
)

# 6. saving combined Plot 
ggsave(
  filename = file.path(
    dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23",
    "06_CFA_AZD1390_HCC44_H23_1_10nM.jpeg"
  ),
  plot = CFA_AZD1390_p_HCC44_H23_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23",
    "06_CFA_AZD1390_HCC44_H23_1_10nM.pdf"
  ),
  plot = CFA_AZD1390_p_HCC44_H23_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 7. saving HCC44 Plot 
ggsave(
  filename = file.path(
    dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23",
    "06_CFA_AZD1390_HCC44_1_10nM.jpeg"
  ),
  plot = CFA_AZD1390_p_HCC44,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23",
    "06_CFA_AZD1390_HCC44_1_10nM.pdf"
  ),
  plot = CFA_AZD1390_p_HCC44,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300
)

# 8. saving NCI-H23 Plot 
ggsave(
  filename = file.path(
    dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23",
    "06_CFA_AZD1390_NCI-H23_1_10nM.jpeg"
  ),
  plot = CFA_AZD1390_p_NCI_H23,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "06_ColonyFormingAssay_AZD1390_HCC44_H23",
    "06_CFA_AZD1390_NCI-H23_1_10nM.pdf"
  ),
  plot =CFA_AZD1390_p_NCI_H23,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300)

#------------------------------------------
setwd(dirs$code)
