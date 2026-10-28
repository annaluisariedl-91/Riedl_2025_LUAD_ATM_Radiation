#------------------------------------------
# #set directory
#------------------------------------------
dose_HCC44 <-c(0,1,2,4,8)
dose_H23<-c(0,0.5,1,2,4)
#------------------------------------------
# loading raw data
#------------------------------------------
setwd(file.path(base_dir,"functional_assays","CFA","Exp3_ATMinhibition","raw","HCC44")) #replace with path on your own PC

HCC44_DMSO_rep <- lapply(
  paste0("CFA_HCC44_DMSO_Rep", 1:2, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

HCC44_1nM_rep <- lapply(
  paste0("CFA_HCC44_1nM_AZD1390_Rep", 1:2, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

HCC44_10nM_rep <- lapply(
  paste0("CFA_HCC44_10nM_AZD1390_Rep", 1:2, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

setwd(file.path(base_dir,"functional_assays","CFA","Exp3_ATMinhibition","raw","NCI-H23")) #replace with path on your own PC

H23_DMSO_rep <- lapply(
  paste0("CFA_H23_DMSO_Rep", 1:3, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

H23_1nM_rep <- lapply(
  paste0("CFA_H23_1nM_AZD1390_Rep", 1:3, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)

H23_10nM_rep <- lapply(
  paste0("CFA_H23_10nM_AZD1390_Rep", 1:3, ".csv"),
  \(file) analyse_CFA_rep(
    read.csv(file),
    remove_first_row = TRUE
  )
)
#------------------------------------------
# binding data and calculation
#------------------------------------------
CFA_ATMi_data <- bind_rows(
  
  make_CFA_data(HCC44_DMSO_rep, "HCC44", "DMSO"),
  make_CFA_data(HCC44_1nM_rep,  "HCC44", "1 nM AZD1390"),
  make_CFA_data(HCC44_10nM_rep, "HCC44", "10 nM AZD1390"),
  
  make_CFA_data(H23_DMSO_rep, "NCI-H23", "DMSO"),
  make_CFA_data(H23_1nM_rep,  "NCI-H23", "1 nM AZD1390"),
  make_CFA_data(H23_10nM_rep, "NCI-H23", "10 nM AZD1390")
)

CFA_ATMi_summary <- CFA_ATMi_data %>%
  group_by(cellline, Treatment, Dose) %>%
  summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    n       = sum(!is.na(SF)),
    .groups = "drop"
  )
#------------------------------------------
# binding for statistics
#------------------------------------------
data_stat_all_ATMi <- CFA_ATMi_data %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c("DMSO", "1 nM AZD1390", "10 nM AZD1390")
    ),
    Cellline = factor(cellline),
    Experiment = factor(experiment),
    adj_sf = pmax(SF, 1e-6),
    adj_log_sf = log(adj_sf)
  )
#-------------------------------------------
# ANOVA
#------------------------------------------
model_all_ATMi_lm <- lm(
  adj_log_sf ~ Experiment + Dose * Treatment * Cellline,
  data = data_stat_all_ATMi
)

summary(model_all_ATMi_lm)
aov_all_ATMi_lm<-anova(model_all_ATMi_lm)

emm_ATMi <- emmeans(
    model_all_ATMi_lm,
    ~ Treatment | Dose + Cellline
  )

posthoc_ATMi <- pairs(
  emm_ATMi,
  adjust = "tukey"
)
#------------------------------------------
# linear-quadrativ modell
#------------------------------------------
LQ_H23_ATMi <- fit_LQ_treatments(
  data_stat_all_ATMi,
  "NCI-H23"
)

coefficients_H23_ATMi <- LQ_H23_ATMi$coefficients

LQ_HCC44_ATMi <- fit_LQ_treatments(
  data_stat_all_ATMi,
  "HCC44"
)

coefficients_HCC44_ATMi <- LQ_HCC44_ATMi$coefficients
#------------------------------------------
# Dose Enhancement Factor DEF
#------------------------------------------
DEF_results_H23_ATMi <- calculate_DEF(
  coefficients_H23_ATMi
)

DEF_label_H23_ATMi <- paste0(
    "1 nM AZD1390: DEF[0.1] = ",
    round(
      DEF_results_H23_ATMi$DEF_1nM[1],
      1
    ),
    "\n",
    "10 nM AZD1390: DEF[0.1] = ",
    round(
      DEF_results_H23_ATMi$DEF_10nM[1],
      1
  )
)

DEF_results_HCC44_ATMi <- calculate_DEF(
  coefficients_HCC44_ATMi
)

DEF_label_HCC44_ATMi <- paste0(
  "1 nM AZD1390: DEF[0.1] = ",
  round(
    DEF_results_HCC44_ATMi$DEF_1nM[1],
    1
  ),
  "\n",
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    DEF_results_HCC44_ATMi$DEF_10nM[1],
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
#------------------------------------------
# Plotting figure
#------------------------------------------
p_CFA_HCC44_AZD <- plot_CFA_AZD(
  plot_data = CFA_ATMi_summary,
  coefficients = coefficients_HCC44_ATMi,
  DEF_label = data.frame(DEF_label_HCC44_ATMi),
  selected_cellline = "HCC44",
  xmax = 8
)

p_CFA_NCI_H23_AZD <- plot_CFA_AZD(
  plot_data = CFA_ATMi_summary,
  coefficients = coefficients_H23_ATMi,
  DEF_label = data.frame(DEF_label_H23_ATMi),
  selected_cellline = "NCI-H23",
  xmax = 4
)

CFA_AZD_final <- (
  p_CFA_HCC44_AZD |
    p_CFA_NCI_H23_AZD
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
  file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23"))

# 1. Mean, SD and SEM of CFA
write.csv(
  CFA_ATMi_summary,
  file = file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23", "CFA_AZD1390_summary.csv"),
  row.names = FALSE
)

# 2. LQ coefficient HCC44
LQ_results_ATMi_HCC44 <- do.call(
  rbind,
  lapply(names(coefficients_HCC44_ATMi), function(treatment) {
    data.frame(
      Treatment = treatment,
      alpha = unname(coefficients_HCC44_ATMi[[treatment]]["alpha"]),
      beta  = unname(coefficients_HCC44_ATMi[[treatment]]["beta"])
    )
  })
)

write.csv(
  LQ_results_ATMi_HCC44,
  file = file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23", "CFA_AZD1390_LQ_HCC44.csv"),
  row.names = FALSE
)

# 3. LQ coefficient NCI-H23
LQ_results_H23_ATMi <- do.call(
  rbind,
  lapply(names(coefficients_H23_ATMi), function(treatment) {
    data.frame(
      Treatment = treatment,
      alpha = unname(coefficients_H23_ATMi[[treatment]]["alpha"]),
      beta  = unname(coefficients_H23_ATMi[[treatment]]["beta"])
    )
  })
)

write.csv(
  LQ_results_H23_ATMi,
  file = file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23", "CFA_AZD1390_LQ_H23.csv"),
  row.names = FALSE
)

# 4. DEF - Dose enhancement factor
write.csv(
  DEF_results_HCC44_ATMi,
  file = file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23", "CFA_AZD1390_DEF_HCC44.csv"),
  row.names = FALSE
)

write.csv(
  DEF_results_H23_ATMi,
  file = file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23", "CFA_AZD1390_DEF_H23.csv"),
  row.names = FALSE
)

# 5. ANOVA results
ANOVA_results_CFA_AZD <-  as.data.frame(
  anova(model_all_ATMi_lm)
)

ANOVA_results_CFA_AZD$Effect <- rownames(
  ANOVA_results_CFA_AZD
)

rownames(ANOVA_results_CFA_AZD) <- NULL

ANOVA_results_CFA_AZD <- ANOVA_results_CFA_AZD %>%
  dplyr::select(
    Effect,
    dplyr::everything()
  )

write.csv(
  ANOVA_results_CFA_AZD,
  file = file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23", "CFA_AZD1390_ANOVA.csv"),
  row.names = FALSE
)

# 6. Post-hoc results
posthoc_results_CFA_AZD <- as.data.frame(
  summary(posthoc_ATMi)
)

write.csv(
  posthoc_results_CFA_AZD,
  file = file.path(dirs$results, "06_CFA_ATMi_HCC44_NCI_H23", "CFA_AZD1390_posthoc.csv"),
  row.names = FALSE
)

# 7. saving combined Plot 
ggsave(
  filename = file.path(
    dirs$results, "06_CFA_ATMi_HCC44_NCI_H23",
    "CFA_AZD1390_1_10nM.jpeg"
  ),
  plot = CFA_AZD_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "06_CFA_ATMi_HCC44_NCI_H23",
    "CFA_AZD1390_1_10nM.pdf"
  ),
  plot = CFA_AZD_final,
  width = 21*2/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 8. saving HCC44 Plot 
ggsave(
  filename = file.path(
    dirs$results, "06_CFA_ATMi_HCC44_NCI_H23",
    "CFA_AZD1390_HCC44_1_10nM.jpeg"
  ),
  plot = p_CFA_HCC44_AZD,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "06_CFA_ATMi_HCC44_NCI_H23",
    "CFA_AZD1390_HCC44_1_10nM.pdf"
  ),
  plot = p_CFA_HCC44_AZD,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300
)

# 8. saving NCI-H23 Plot 
ggsave(
  filename = file.path(
    dirs$results, "06_CFA_ATMi_HCC44_NCI_H23",
    "CFA_AZD1390_NCI-H23_1_10nM.jpeg"
  ),
  plot = p_CFA_NCI_H23_AZD,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "06_CFA_ATMi_HCC44_NCI_H23",
    "CFA_AZD1390_NCI-H23_1_10nM.pdf"
  ),
  plot =p_CFA_NCI_H23_AZD,
  width = 21*4/9,
  height = 7,
  units = "cm",
  dpi = 300
)
#------------------------------------------
setwd(dirs$code)
