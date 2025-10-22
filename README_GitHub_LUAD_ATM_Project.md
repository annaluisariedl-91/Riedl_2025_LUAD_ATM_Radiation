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
│   ├── 02_importing_data.R
│   ├── 03_CFA_8Gy.R
│   ├── 04_CFA_2Gy.R
│   ├── 05_CFA_AZD1390_HCC44_H23.R
│   ├── 06_CFA_AZD1390_other.R
│   ├── 07_xcelligence_growth.R
│   ├── 08_xcelligence_ATMinhibition.R
│   ├── 09_MTS_DoseResponse.R
│   ├── 10_MTS_ATMinhibition.R
│   ├── 11_MTS_IC50_AZD1390.R
│   ├── 12_wb_ATM.R
│   ├── 13_RNA_calc_DE.R
│   ├── 14_RNA_qc.R
│   ├── 15_RNA_Vulcano_Top10.R
│   ├── 16_protein_calc_DE.R
│   ├── 17_protein_qc.R
│   ├── 18_protein_control_Vulcano_Top10.R
│   ├── 19_protein_Radiation_Vulcano_Top10.R
│   ├── 20_protein_HCC44_Vulcano_Top10.R
│   ├── 21_protein_H23_Vulcano_Top10.R
│   ├── 22_correlation_plot.R
│   ├── 23_protein_splineRegression.R
│   └── 24_protein_ORA.R
└── test/                # example outputs or intermediate results
```

---

## ⚙️ Setup

**Working directories:**
```r
# set save / source paths
save_wd      <- "C:/Users/annar/Desktop/Codes+Data/r codes/test"
source_wd    <- "C:/Users/annar/Desktop/Codes+Data/r codes"
```

All file paths are currently local (`C:/Users/...`).  
For public repositories, please use relative paths or environment variables 
(e.g. `here::here()` or `Sys.getenv("PROJECT_DIR")`).

---

## 🧠 Analysis Workflow

### **1. Setup and Data Import**
```r
source("01_importing_packages.R")  # installs/loads all required R packages
source("02_importing_data.R")      # imports metadata, design matrices and datasets
```

### **2. Functional Assays**
```r
source("03_CFA_8Gy.R")             # clonogenic assay – 8 Gy dose response
source("04_CFA_2Gy.R")             # clonogenic assay – 2 Gy response
source("05_CFA_AZD1390_HCC44_H23.R") 
source("06_CFA_AZD1390_other.R")   # AZD1390 effects across other LUAD lines
source("07_xcelligence_growth.R")  # impedance-based growth curves
source("08_xcelligence_ATMinhibition.R")  # ATM inhibitor ± irradiation
source("09_MTS_DoseResponse.R")    # metabolic dose response (MTS)
source("10_MTS_ATMinhibition.R")   # MTS ± AZD1390
source("11_MTS_IC50_AZD1390.R")    # IC50 determination for AZD1390
source("12_wb_ATM.R")              # Western blot ATM quantification
```

### **3. Transcriptome (RNA-seq)**
```r
source("13_RNA_calc_DE.R")         # differential expression analysis
source("14_RNA_qc.R")              # quality control
source("15_RNA_Vulcano_Top10.R")   # volcano plots and top genes
```

### **4. Proteome**
```r
source("16_protein_calc_DE.R")     # differential protein expression (Protti)
source("17_protein_qc.R")          # QC of proteomics data
source("18_protein_control_Vulcano_Top10.R")
source("19_protein_Radiation_Vulcano_Top10.R")
source("20_protein_HCC44_Vulcano_Top10.R")
source("21_protein_H23_Vulcano_Top10.R")
```

### **5. Correlation and Integration**
```r
source("22_correlation_plot.R")    # transcriptome vs proteome correlation
```

### **6. Temporal Expression and ORA**
```r
source("23_protein_splineRegression.R")  # spline modeling of protein time courses
source("24_protein_ORA.R")               # over-representation analysis (pathways)
```

---

## 📦 Required Packages
Typical dependencies (specified in `01_importing_packages.R`):
```
tidyverse, ggplot2, DESeq2, clusterProfiler, limma, edgeR, cowplot,
pheatmap, ggpubr, readxl, openxlsx, RColorBrewer, here, data.table,
msigdbr, splines, nlme, ggpmisc
```

Install them via:
```r
source("01_importing_packages.R")
```

---

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
