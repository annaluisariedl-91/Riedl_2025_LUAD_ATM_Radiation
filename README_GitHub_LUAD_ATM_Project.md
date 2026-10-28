# LUAD ATM Project – Functional, Transcriptomic and Proteomic Analysis

This repository contains R scripts and data processing workflows used to analyze 
functional assays, transcriptome, and proteome data from lung adenocarcinoma (LUAD) 
cell lines with different ATM status. The aim is to explore temporal radiation responses,
identifing determinatns and role of different ATM status for radiosensitivity across multiple experimental layers.

---

## 📁 Directory Structure

```
Codes+Data/
├── r codes/
│   ├── 00_main_analysis.R
│   ├── 01_importing_packages.R
│   ├── 02_defining functions.R
│   ├── 03_importing_data.R
│   ├── 04_CFA_8Gy.R
│   ├── 05_CFA_2Gy.R
│   ├── 06_CFA_AZD1390_HCC44_H23.R
│   ├── 07_CFA_AZD1390_other.R
│   ├── 08_xcelligence_growth.R
│   ├── 09_xcelligence_ATMinhibition.R
│   ├── 10_MTS_DoseResponse.R
│   ├── 11_MTS_ATMinhibition.R
│   ├── 12_MTS_IC50_AZD1390.R
│   ├── 13_WesternBlot_ATM.R
│   ├── 14_Transcriptome_DifferentialExpression.R
│   ├── 15_Proteome_DifferentialExpression.R
│   ├── 16_Proteom_Transcriptome_Correlation.R
│   ├── 17_Transcriptome_QualityControl.R
│   ├── 18_Transcriptome_Volcano_top10.R
│   ├── 19_Proteome_QualityControl.R
│   ├── 20_Proteome_Volcano_top10.R
│   ├── 21_Proteome_SplineRegression.R
│   ├── 22_Proteome_OverrepresentationAnalysis.R
│   ├── 23_protein_ProteinProteinInteraction.R
│   └── 24_clinical data_cohort_analysis
|       ├── requirements
|       ├── README
|       ├── msk_met_atm_analysis
|       ├── data
|       |   ├── ATM_MSK-MET_lung_extract.csv
|       |   └── ATM_MSK-MET_variants.csv
|       └── results  
└── test/                # example outputs or intermediate results
```

---

## ⚙️ Setup

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

## 🧠 Analysis Workflow

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

## 📦 Required Packages
cran_packages <- c(
  "broom",
  "corrplot",
  "data.table",
  "drc",
  "dplyr",
  "emmeans",
  "factoextra",
  "forcats",
  "fpc",
  "ggplot2",
  "ggprism",
  "ggraph",
  "gprofiler2",
  "gridExtra",
  "HGNChelper",
  "Hmisc",
  "httr",
  "igraph",
  "jsonlite",
  "lme4",
  "lmerTest",
  "Matrix",
  "minpack.lm",
  "nlme",
  "patchwork",
  "pheatmap",
  "pracma",
  "purrr",
  "RColorBrewer",
  "readr",
  "rrvgo",
  "stringr",
  "survminer",
  "tibble",
  "tidygraph",
  "tidyr",
  "VennDiagram",
  "viridis"
)

# Bioconductor-Packages
bioc_packages <- c(
  "AnnotationDbi",
  "clusterProfiler",
  "DOSE",
  "EnhancedVolcano",
  "genefilter",
  "Homo.sapiens",
  "limma",
  "msigdbr",
  "protti",
  "STRINGdb"
)

missing_cran <- setdiff(
  cran_packages,
  rownames(installed.packages())
)
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
missing_bioc <- setdiff(
  bioc_packages,
  rownames(installed.packages())
)
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE)
}
#------------------------------------------
# Loading packages
#------------------------------------------
library(AnnotationDbi)
library(data.table)
library(dplyr)
library(forcats)
library(HGNChelper)
library(Homo.sapiens)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# Statistical analysis
library(broom)
library(drc)
library(emmeans)
library(limma)
library(lme4)
library(lmerTest)
library(minpack.lm)
library(nlme)
library(splines)

# Proteomics and enrichment
library(clusterProfiler)
library(DOSE)
library(EnhancedVolcano)
library(gprofiler2)
library(msigdbr)
library(protti)
library(rrvgo)

# Clustering and networks
library(factoextra)
library(fpc)
library(ggraph)
library(igraph)
library(pheatmap)
library(STRINGdb)
library(tidygraph)

# Visualization
library(corrplot)
library(ggplot2)
library(ggprism)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(VennDiagram)
library(viridis)

# Survival analysis
library(survival)
library(survminer)
## 🧾 Output
All generated plots, processed data, and results are saved under:
```
save_wd <- "C:/Users/annar/Desktop/Codes+Data/r codes/test"
```

Output types:
- `.png`/`.tiff` plots (CFA, MTS, xCELLigence, Volcano, Correlation)
- `.tsv`/`.txt` processed result tables (DE, normalized data, metrics)
- `.RData` intermediate files (optional)

---

## 🧪 Data References
Raw and processed datasets (RNA, Proteome, Functional Assays) correspond to 
datasets described in GEO, PRIDE, and institutional repositories.  
Design matrices and README files for each experiment are stored in `metadata/`.

---

## 🧑‍💻 Author
**Hanibal Bohnenberger**  
Institute of Pathology  
University Medical Center Göttingen  
📧 Hanibal.bohnenberger@med.uni-goettingen.de

---

## 📜 License
This repository is released under the MIT License.  
You are free to use, modify, and distribute with attribution.

---

## 🧭 Citation
If you use this code or derived data in publications, please cite:

> Riedl et al., *Toward biomarker-guided radiotherapy: Multi-omics analysis links ATM status to radioresistance in lung adenocarcinoma*, 2025.
