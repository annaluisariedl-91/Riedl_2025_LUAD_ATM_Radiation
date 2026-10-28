# ATM alterations in the MSK-MET lung cancer cohort

This directory contains the Python-based clinical cohort analysis accompanying the manuscript "Toward biomarker-guided radiotherapy: Multi-omics analysis links ATM status to radioresistance in lung adenocarcinoma". The analysis is part of the LUAD ATM Project repository
and used for figure 7.

Code and data for Figure 7 and the Results section *"ATM alterations in lung adenocarcinoma: frequency, genomic context and outcome in the MSK-MET cohort"* of

> Riedl AL, Reiß JJ, Niemeier T, et al. Toward biomarker-guided radiotherapy: Multi-omics analysis links ATM status to radioresistance in lung adenocarcinoma. *Cancer Cell International* (revised manuscript, 2026).

## Contents

```
cohort_analysis/
├── msk_met_atm_analysis.py        analysis script (statistics + Figure 7)
├── requirements.txt               Python dependencies
├── data/
│   ├── ATM_MSK-MET_lung_extract.csv   one row per LUAD/LUSC sample (n = 4,601)
│   └── ATM_MSK-MET_variants.csv       all ATM variants of the lung samples (n = 380)
└── results/                       output of the script as used in the manuscript
    ├── cohort_statistics.txt      all numbers reported in the text
    ├── coalterations.csv          odds ratios, 95 % CI and Fisher p-values (Figure 7C)
    ├── Figure_7_170mm.pdf
    └── Figure_7_170mm.png
```

## Data provenance

* Source: MSK-MET (Nguyen B, et al. Genomic characterization of metastatic patterns from prospective clinical sequencing of 25,000 patients. *Cell* 2022;185:563-575.e11), cBioPortal study `msk_met_2021`, <https://www.cbioportal.org/study/summary?id=msk_met_2021>, accessed 26 October 2026. Tumour-normal matched MSK-IMPACT sequencing (341/410/468-gene panels); all variants are somatic. Public data without access restriction.
* Extract: all samples with OncoTree code `LUAD` (n = 4,064) or `LUSC` (n = 537), one sample per patient.
* `ATM_any` = 1 if the sample carries any non-synonymous ATM variant; `ATM_truncating` = 1 if at least one nonsense, frameshift or splice-site variant. Missense variants are predominantly of unknown significance; analyses of "ATM loss" use `ATM_truncating`.
* Co-alteration flags (`KRAS_mut`, `EGFR_mut`, `TP53_mut`, `STK11_mut`, `KEAP1_mut`, `BRAF_mut`, `MET_mut` = MET exon 14 skipping, `PIK3CA_mut`): 1 for oncogenic or likely oncogenic variants by hotspot and truncation rules (activating hotspots for KRAS, EGFR, BRAF, PIK3CA; truncating variants and known pathogenic missense variants for TP53, STK11 and KEAP1). `MET_amp` and `CDKN2A_del` from the copy-number calls of cBioPortal.
* Other columns: `SAMPLE_TYPE`, `SITE_GROUP`, `METASTATIC_SITE`, `MET_SITE_COUNT`, `MET_COUNT`, `TMB_NONSYNONYMOUS` (non-synonymous mutations per Mb), `FGA`, `MSI_SCORE`, `MSI_TYPE`, `GENE_PANEL`, `AGE_AT_SEQUENCING`, `SEX`, `RACE`, `OS_MONTHS` and `OS_STATUS` (overall survival from the time of sequencing). Smoking status and radiotherapy are not recorded in MSK-MET.

## How to run
Run the following commands from the 24_clinical_data_cohort_analysis directory.

```bash
pip install -r requirements.txt
python msk_met_atm_analysis.py            # reads data/ATM_MSK-MET_lung_extract.csv, writes results/
```

Tested with Python 3.12.3, pandas 2.3.3, SciPy 1.17.1, lifelines 0.30.3, matplotlib 3.10.8.

## Statistics

* Frequencies compared by Fisher's exact test (LUAD vs. LUSC; primary vs. metastasis).
* Co-alterations: ATM-truncated (n = 120) and ATM-altered (n = 311) versus ATM wild-type LUAD (n = 3,753); odds ratios with Woolf 95 % confidence intervals (Haldane correction for zero cells) and Fisher's exact test.
* Tumour mutational burden: Mann-Whitney U test.
* Overall survival: Kaplan-Meier estimates and log-rank test versus ATM wild-type (4,041 LUAD with follow-up).

Key results (reproduced by the script): ATM-altered 7.7 % and ATM-truncated 3.0 % of LUAD versus 4.1 % and 0.6 % of LUSC; KRAS 57.9 % vs. 30.7 % (OR 3.10), STK11 27.0 % vs. 14.8 % (OR 2.13), EGFR 8.4 % vs. 27.8 % (OR 0.24), TP53 34.4 % vs. 48.4 % (OR 0.56) in ATM-altered vs. wild-type LUAD; median TMB 9.5 vs. 4.9 mutations/Mb; median overall survival 40.8 (wild-type), 38.1 (missense/in-frame only) and 29.5 months (truncated), log-rank p = 0.10 (any alteration) and 0.13 (truncated).
