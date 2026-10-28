#------------------------------------------
# generating directories # chose path in your own directory
#------------------------------------------
base_dir <- "C:/Users/annar/Paper_Thesis/Codes+Data"

dirs <- list(
  code       = file.path(base_dir, "r codes"),
  RNA        = file.path(base_dir, "GEO"),
  proteome   = file.path(base_dir, "Pride", "results"),
  data       = file.path(base_dir, "data"),
  results    = file.path(base_dir, "results"),
  plots      = file.path(base_dir, "results", "plots")
)

dir.create(dirs$results, recursive = TRUE, showWarnings = FALSE)
dir.create(dirs$plots, recursive = TRUE, showWarnings = FALSE)

setwd(dirs$code)
#------------------------------------------
# loading parts of analysis
#------------------------------------------
#1. Setup
source("01_importing_packages.R")
source("02_defining functions.R")
source("03_importing_data.R")

#2. Analysis of functional assays
source("04_CFA_8Gy.R")
source("05_CFA_2Gy.R")
source("06_CFA_AZD1390_HCC44_H23.R")
source("07_CFA_AZD1390_other.R")

source("08_xcelligence_growth.R")
source("09_xcelligence_ATMinhibition.R")

source("10_MTS_DoseResponse.R")
source("11_MTS_ATMinhibition.R")
source("12_MTS_IC50_AZD1390.R")

source("13_WesternBlot_ATM.R")

#3. Calculation of differetial expression
source("14_Transcriptome_DifferentialExpression.R")
source("15_Proteome_DifferentialExpression.R")

#4. Transcriptome vs Proteome
source("16_Proteom_Transcriptome_Correlation.R")

#5. Quality analysis of Transcriptome data
source("17_Transcriptome_QualityControl.R")
source("18_Transcriptome_Volcano_top10.R")

#6. Analysis of Proteome data
source("19_Proteome_QualityControl.R")
source("20_Proteome_Volcano_top10.R")

#7. Temporal Expression Pattern
source("21_Proteome_SplineRegression.R")
source("22_Proteome_OverrepresentationAnalysis.R")
source("23_protein_ProteinProteinInteraction.R")
