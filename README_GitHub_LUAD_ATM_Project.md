# LUAD ATM Project – Functional, Transcriptomic and Proteomic Analysis

This repository contains the R and Python scripts used for the statistical analysis and visualization of experimental and clinical data associated with the manuscript "Toward biomarker-guided radiotherapy: Multi-omics analysis links ATM status to radioresistance in lung adenocarcinoma".

The analyses include functional assays, transcriptomic and proteomic profiling, and clinical cohort analyses to investigate the association between ATM status and radiosensitivity in lung adenocarcinoma (LUAD).

The repository provides the analytical workflows, supporting datasets, and results to facilitate transparency and reproducibility.
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
│       ├── requirements.txt
│       ├── README.md
│       ├── msk_met_atm_analysis.py
|       ├── data
|       |   ├── ATM_MSK-MET_lung_extract.csv
|       |   └── ATM_MSK-MET_variants.csv
|       └── results  
└── test/                # example outputs or intermediate results
```

---

## ⚙️ Setup

base_dir <- getwd()

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

The R scripts can be executed either as a complete analysis pipeline or individually to reproduce specific results.

Option 1: Complete Analysis
---
To reproduce all analyses, run the main script:

source("00_main_analysis.R")

The main script initializes the analysis environment, loads the required packages, functions, and datasets, and executes the individual analysis scripts in the predefined order.

Option 2: Individual Analyses
---
Individual scripts can also be executed independently to reproduce specific analyses without running the entire pipeline.

For example, to reproduce the clonogenic survival analysis:

source("04_CFA_8Gy.R")

Each analysis script initializes the required environment if it has not already been loaded.

Analysis Structure

The workflow is organized into the following components:

Scripts	Description
---
00–03	Main workflow, package loading, function definitions, and data import
04–07	Clonogenic survival assays and ATM inhibitor experiments
08–09	Impedance-based cell growth and radiation response
10–12	MTS assays and AZD1390 dose-response analysis
13	Western blot quantification
14–15	Differential expression analysis
16	Transcriptome–proteome correlation
17–18	Transcriptomic quality control and visualization
19–20	Proteomic quality control and visualization
21–23	Temporal expression analysis, pathway enrichment, and protein–protein interaction networks

The clinical MSK-MET cohort analysis is implemented separately in Python. Instructions for reproducing this analysis are provided in its corresponding README.

## 📦 Required Packages
The R analyses require several packages from CRAN and Bioconductor.

All required packages are listed in 01_importing_packages.R, which automatically checks for missing dependencies and installs them if necessary.

To install and load the required packages, run:

source("01_importing_packages.R")

The analyses were performed using R [VERSION]. Additional information about the computational environment and package versions can be obtained using sessionInfo().
---

## 🧪 Data References
Raw and processed datasets (RNA, Proteome, Functional Assays) correspond to 
datasets described in GEO, PRIDE, and institutional repositories.  
Design matrices and README files for each experiment are stored in `metadata/`.

## Output
---
The analysis scripts generate statistical results, processed datasets, and publication-quality figures.

Generated files are organized in the results/ directory, with figures stored in results/plots/.

Output formats include:

.png, .pdf, and .tiff files for figures and visualizations.
.csv, .tsv, and .txt files for statistical results and processed datasets.
.RData files for intermediate analysis objects, where applicable.

The Python-based MSK-MET analysis generates its results in the corresponding results/ subdirectory.

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
