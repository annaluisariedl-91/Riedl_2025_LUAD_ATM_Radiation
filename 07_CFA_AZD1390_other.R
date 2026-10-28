#-------------------------------------------
# set directory
#-------------------------------------------
setwd(file.path(base_dir,"functional_assays","CFA","Exp4_ATMinhibition","raw")) #replace with path on your own PC
dose <-c(0,2,4,8)
A549<- read.csv("CFA_AZD1390_A549.csv",sep=",")
HCC827<- read.csv("CFA_AZD1390_HCC827.csv",sep=",")
HCC78<- read.csv("CFA_AZD1390_HCC78.csv",sep=",")
#-------------------------------------------
# importing A549
#-------------------------------------------
rep_cols <- list(
  2:4,
  7:9,
  12:14,
  17:19
)

A549_DMSO_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      A549,
      rows = 19:30,
      cols = c(cols)
    )
  }
)

A549_10nM_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      A549,
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

HCC827_DMSO_rep <- lapply(
  rep_cols_HCC827,
  function(cols) {
    analyse_CFA_rep(
      HCC827,
      rows = 19:30,
      cols = cols
    )
  }
)

HCC827_10nM_rep <- lapply(
  rep_cols_HCC827,
  function(cols) {
    analyse_CFA_rep(
      HCC827,
      rows = 3:14,
      cols = cols
    )
  }
)
#-------------------------------------------
# importing HCC78
#-------------------------------------------
HCC78_DMSO_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      HCC78,
      rows = 19:30,
      cols = cols
    )
  }
)

HCC78_10nM_rep <- lapply(
  rep_cols,
  function(cols) {
    analyse_CFA_rep(
      HCC78,
      rows = 3:14,
      cols = cols
    )
  }
)
#-------------------------------------------
# Binding data
#-------------------------------------------
CFA_ATMi_3cell_10nM_data <- bind_rows(
  # A549
  make_CFA_data(
    A549_DMSO_rep,
    "A549",
    "DMSO"
  ),
  make_CFA_data(
    A549_10nM_rep,
    "A549",
    "10 nM AZD1390"
  ),
  # HCC827
  make_CFA_data(
    HCC827_DMSO_rep,
    "HCC827",
    "DMSO"
  ),
  make_CFA_data(
    HCC827_10nM_rep,
    "HCC827",
    "10 nM AZD1390"
  ),
  # HCC78
  make_CFA_data(
    HCC78_DMSO_rep,
    "HCC78",
    "DMSO"
  ),
  make_CFA_data(
    HCC78_10nM_rep,
    "HCC78",
    "10 nM AZD1390"
  )
)

stat_3cell_ATMi_10nM <- CFA_ATMi_3cell_10nM_data %>%
  mutate(
    Cellline = factor(cellline),
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

mean_3cell_ATMi_10nM <- stat_3cell_ATMi_10nM %>%
  group_by(Cellline, Treatment, Dose) %>%
  summarise(
    SF_mean = mean(SF, na.rm = TRUE),
    SF_SD   = sd(SF, na.rm = TRUE),
    SF_SEM  = SF_SD / sqrt(sum(!is.na(SF))),
    n       = sum(!is.na(SF)),
    .groups = "drop"
  )
#-------------------------------------------
# Statistics
#-------------------------------------------
results_A549_ATMi <- analyse_CFA_model(
  stat_3cell_ATMi_10nM,
  "A549"
)

results_HCC78_ATMi <- analyse_CFA_model(
  stat_3cell_ATMi_10nM,
  "HCC78"
)

results_HCC827_ATMi <- analyse_CFA_model(
  stat_3cell_ATMi_10nM,
  "HCC827"
)

results_CFA_3cell_ATMi_postHoc <- bind_rows(
  A549   = as.data.frame(results_A549_ATMi$posthoc),
  HCC78  = as.data.frame(results_HCC78_ATMi$posthoc),
  HCC827 = as.data.frame(results_HCC827_ATMi$posthoc),
  .id = "Cellline"
)

results_CFA_3cell_ATM_AOV <- bind_rows(
  A549   = as.data.frame(results_A549_ATMi$anova) %>%
    tibble::rownames_to_column("Effect"),
  HCC78  = as.data.frame(results_HCC78_ATMi$anova) %>%
    tibble::rownames_to_column("Effect"),
  HCC827 = as.data.frame(results_HCC827_ATMi$anova) %>%
    tibble::rownames_to_column("Effect"),
  .id = "Cellline"
)
#-------------------------------------------
# global models
#-------------------------------------------
model_global_3cell_ATM <- lm(
  adj_log_sf ~ Experiment  +
    factor(Dose) * Treatment * Cellline,
  data = stat_3cell_ATMi_10nM
)

res_aov_global_3cell_ATM<-anova(model_global_3cell_ATM)

emm_global_3cell_ATM <- emmeans(
  model_global_3cell_ATM,
  ~ Treatment | Dose + Cellline
)

posthoc_global_3cell_ATM <- pairs(
  emm_global_3cell_ATM,
  adjust = "tukey"
)
#-------------------------------------------
# linear-quadratic model
#-------------------------------------------
stat_3cell_ATMi_10nM <- mean_3cell_ATMi_10nM %>%
  mutate(
    adj_sf = pmax(SF_mean, 1e-6),
    adj_log_sf = log(adj_sf)
  )

LQ_A549_ATMi_10nM <- fit_LQ_treatments(
  stat_3cell_ATMi_10nM, "A549"
)

LQ_HCC78_ATMi_10nM <- fit_LQ_treatments(
  stat_3cell_ATMi_10nM, "HCC78"
)

LQ_HCC827_ATMi_10nM <- fit_LQ_treatments(
  stat_3cell_ATMi_10nM, "HCC827"
)

coefficients_A549_ATMi_10nM   <- LQ_A549_ATMi_10nM$coefficients
coefficients_HCC78_ATMi_10nM  <- LQ_HCC78_ATMi_10nM$coefficients
coefficients_HCC827_ATMi_10nM <- LQ_HCC827_ATMi_10nM$coefficients

coefficients_all_ATMi_10nM <- bind_rows(A549 = as.data.frame(LQ_A549_ATMi_10nM$coefficients) %>%
                         tibble::rownames_to_column("coefficient"),
                       HCC78  = as.data.frame(LQ_HCC78_ATMi_10nM$coefficients) %>%
                         tibble::rownames_to_column("coefficient"),
                       HCC827 = as.data.frame(LQ_HCC827_ATMi_10nM$coefficients) %>%
                         tibble::rownames_to_column("coefficient"),
                       .id = "Cellline")
#-------------------------------------------
# Dose enhancement factor (DEF)
#-------------------------------------------
DEF_results_A549_ATMi_10nM <- calc_DEF(
  coefficients_A549_ATMi_10nM
)

DEF_label_A549_ATMi_10nM <- paste0(
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    DEF_results_A549_ATMi_10nM$DEF_10nM[1],
    1
  )
)

DEF_results_HCC78_ATMi_10nM <- calc_DEF(
  coefficients_HCC78_ATMi_10nM
)

DEF_label_HCC78_ATMi_10nM <- paste0(
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    DEF_results_HCC78_ATMi_10nM$DEF_10nM[1],
    1
  )
)

DEF_results_HCC827_ATMi_10nM <- calc_DEF(
  coefficients_HCC827_ATMi_10nM
)

DEF_label_HCC827_ATMi_10nM <- paste0(
  "10 nM AZD1390: DEF[0.1] = ",
  round(
    DEF_results_HCC827_ATMi_10nM$DEF_10nM[1],
    1
  )
)

DEF_all_ATMi_10nM <- bind_rows(  A549   = as.data.frame(DEF_results_A549_ATMi_10nM) %>%
                         tibble::rownames_to_column("Effect"),
                       HCC78  = as.data.frame(DEF_results_HCC78_ATMi_10nM) %>%
                         tibble::rownames_to_column("Effect"),
                       HCC827 = as.data.frame(DEF_results_HCC827_ATMi_10nM) %>%
                         tibble::rownames_to_column("Effect"),
                       .id = "Cellline")
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

CFA_significance <- data.frame(
  Cellline = c(
    "A549", "A549", "A549",
    "HCC78", "HCC78",
    "HCC827", "HCC827"
  ),
  Dose = c(
    2, 4, 8,
    2, 4,
    2, 4
  ),
  label = c(
    "*", "***", "*",
    "*", "***",
    "***", "***"
  ),
  y = c(
    2, 0.15, 0.002,
    1, 0.02,
    1, 0.02
  )
)
#--------------------------------------------
# Plotting figures
#--------------------------------------------
mean_3cell_ATMi_10nM<- mean_3cell_ATMi_10nM %>%
  mutate(
    SF = SF_mean,
    cellline=Cellline 
  )

p_CFA_A549_ATMi_10nM <- plot_CFA_AZD(
  plot_data = mean_3cell_ATMi_10nM,
  coefficients = coefficients_A549_ATMi_10nM,
  DEF_label = data.frame(DEF_label_A549_ATMi_10nM),
  selected_cellline = "A549",
  xmax = 8
)

p_CFA_HCC78_ATMi_10nM <- plot_CFA_AZD(
  plot_data = mean_3cell_ATMi_10nM,
  coefficients = coefficients_HCC78_ATMi_10nM,
  DEF_label =  data.frame(DEF_label_HCC78_ATMi_10nM),
  selected_cellline = "HCC78",
  xmax = 8
)

p_CFA_HCC827_ATMi_10nM <- plot_CFA_AZD(
  plot_data = mean_3cell_ATMi_10nM,
  coefficients = coefficients_HCC827_ATMi_10nM,
  DEF_label =  data.frame(DEF_label_HCC827_ATMi_10nM),
  selected_cellline = "HCC827",
  xmax = 8
)

CFA_AZD_3Cell_ATMi_10nM_final <- (
  p_CFA_A549_ATMi_10nM |
    p_CFA_HCC78_ATMi_10nM |
    p_CFA_HCC827_ATMi_10nM
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
  file.path(dirs$results, "07_CFA_AZD1390_other"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "07_CFA_AZD1390_other"))
         
# 1. Mean, SD and SEM of CFA
write.csv(
  mean_3cell_ATMi_10nM,
  file = file.path(dirs$results, "07_CFA_AZD1390_other", "CFA_AZD1390_summary.csv"),
  row.names = FALSE
)

# 2. LQ coefficient 
write.csv(
  coefficients_all_ATMi_10nM,
  file = file.path(dirs$results, "07_CFA_AZD1390_other", "CFA_AZD1390_LQ_all.csv"),
  row.names = FALSE
)

# 3. DEF - Dose enhancement factor
write.csv(
  DEF_all_ATMi_10nM,
  file = file.path(dirs$results, "07_CFA_AZD1390_other", "CFA_AZD1390_DEF_all.csv"),
  row.names = FALSE
)

# 4. ANOVA results
write.csv(
  results_CFA_3cell_ATM_AOV,
  file = file.path(dirs$results, "07_CFA_AZD1390_other", "CFA_AZD1390_ANOVA.csv"),
  row.names = FALSE
)

# 5. global ANOVA results
write.csv(
  res_aov_global_3cell_ATM,
  file = file.path(dirs$results, "07_CFA_AZD1390_other", "CFA_AZD1390_ANOVA_global.csv"),
  row.names = FALSE
)

# 6. Post-hoc results
write.csv(
  results_CFA_3cell_ATMi_postHoc,
  file = file.path(dirs$results, "07_CFA_AZD1390_other", "CFA_AZD1390_posthoc.csv"),
  row.names = FALSE
)

# 7. Post-hoc results
write.csv(
  posthoc_global_3cell_ATM,
  file = file.path(dirs$results, "07_CFA_AZD1390_other", "CFA_AZD1390_posthoc_global.csv"),
  row.names = FALSE
)
 
# 8. saving combined Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_all.jpeg"
  ),
  plot = CFA_AZD_3Cell_ATMi_10nM_final,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_all.pdf"
  ),
  plot = CFA_AZD_3Cell_ATMi_10nM_final,
  width = 21*3/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 9. saving A549 Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_A549.jpeg"
  ),
  plot = p_CFA_A549_ATMi_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_A549.pdf"
  ),
  plot = p_CFA_A549_ATMi_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)
 
# 10. saving HCC78 Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_HCC78.jpeg"
  ),
  plot = p_CFA_HCC78_ATMi_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_HCC78.pdf"
  ),
  plot = p_CFA_HCC78_ATMi_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

# 11. saving HCC827 Plot 
ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_HCC827.jpeg"
  ),
  plot = p_CFA_HCC827_ATMi_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)

ggsave(
  filename = file.path(
    dirs$results, "07_CFA_AZD1390_other",
    "CFA_AZD1390_10nM_HCC827.pdf"
  ),
  plot = p_CFA_HCC827_ATMi_10nM,
  width = 21*1/3,
  height = 7,
  units = "cm",
  dpi = 300
)
#------------------------------------------
setwd(dirs$code)
