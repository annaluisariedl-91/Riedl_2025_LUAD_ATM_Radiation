#------------------------------------------
# installing packages, if needed
#------------------------------------------
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

#-----------------------------------------------
setwd(dirs$code)
