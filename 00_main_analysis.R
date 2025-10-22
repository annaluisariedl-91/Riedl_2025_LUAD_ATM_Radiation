#saving directory
save_wd <- paste("C:/Users/annar/Desktop/Codes+Data/r codes/test")
#loading codes
source_wd <- paste("C:/Users/annar/Desktop/Codes+Data/r codes")
#loading functional assays
functionalAssays_wd <- paste("C:/Users/annar/Desktop/Codes+Data/r codes/test")
#loading proteome data
proteome_wd <- paste("C:/Users/annar/Desktop/Codes+Data/r codes")
## main_analysis.R
setwd("C:/Users/annar/Desktop/Codes+Data/r codes")
# 1. Setup
source("01_importing_packages.R")
source("02_importing_data.R")

# 2. Analysis of functional assays
source("03_CFA_8Gy.R")
source("04_CFA_2Gy.R")
source("05_CFA_AZD1390_HCC44_H23.R")
source("06_CFA_AZD1390_other.R")
source("07_xcelligence_growth.R")
source("08_xcelligence_ATMinhibition.R")
source("09_MTS_DoseResponse.R")
source("10_MTS_ATMinhibition.R")
source("11_MTS_IC50_AZD1390.R")
source("12_wb_ATM.R")

# 3. Analysis of Transcriptome data
source("13_RNA_calc_DE.R")
source("14_RNA_qc.R")
source("15_RNA_Vulcano_Top10.R")

# 4. Analysis of Proteome data
source("16_protein_calc_DE.R")
source("17_protein_qc.R")
source("18_protein_control_Vulcano_Top10.R")
source("19_protein_Radiation_Vulcano_Top10.R")
source("20_protein_HCC44_Vulcano_Top10.R")
source("21_protein_H23_Vulcano_Top10.R")

# 5. Transcriptome vs Proteome
source("22_correlation_plot.R")

# 6. Temporal Expression Pattern
source("23_protein_splineRegression.R")
setwd(source_wd)
source("24_protein_ORA.R")
