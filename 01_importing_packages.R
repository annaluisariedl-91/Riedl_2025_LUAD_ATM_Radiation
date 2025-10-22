####loading packages####
library(splines)
library(limma)
library(Matrix)
library(factoextra)
library(dplyr)
library(HGNChelper)
library(protti)
library(corrplot)
library(RColorBrewer)
library(pheatmap)
library(vioplot)
library(VennDiagram)
library(pathfindR)
library(Hmisc)
library(tidyverse)
library(ggplot2)
library(viridis)
library(DescTools)
library(stringr)
library(drc)
library(minpack.lm)
library(hgnc)
library(org.Dr.eg.db)
library(Homo.sapiens)
#library(rlang)
library(dtw)
library(reshape2)
library(proxy)
library(clusterSim)
library(tidyr)
library(cluster)     
library(fpc)          
library(tibble)
library(genefilter)
library(forcats)
library(msigdbr)
library(purrr)
library(patchwork)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(mixOmics)
library(ReactomeGSA)
library(DOSE)
library(rrvgo)
library(STRINGdb)
library(igraph)
library(gprofiler2)
library(data.table)
#library(lme4)
library(readr)
library(AnnotationDbi) 
#####
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

save_pheatmap_png <- function(x, filename, width=2400, height=2000, res = 250) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
####
cols_Con=c( "#E5F5E0" , "#A1D99B", "#74C476" , "#238B45", "#00441B")
pal_Con = colorRampPalette(cols_Con)
cols_HCC44=c( "#DEEBF7" , "#9ECAE1", "#6BAED6", "#2171B5", "#08306B")
pal_HCC44 = colorRampPalette(cols_HCC44)
cols_RAD=c( "#FEE0D2" , "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D")
pal_RAD = colorRampPalette(cols_RAD)

color<-data.frame(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"))
colnames(color)<- "color"
ann_colors <- list(
  Radiation = c("control" = "grey", "2 Gy" = "#EF3B2C"),
  Cell_line = c(HCC44 =  "#6BAED6", "NCI-H23" = "#41AB5D"),
  Time = c("0h" ="#FB8072" ,"2h"= "#FDB462" ,"6h"= "#8DD3C7" , "24h"="#80B1D3", "48h"="#BEBADA")
) 

my_colors <-rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))

names_prot <-c("HCC44 c 0h.1","HCC44 c 0h.2","HCC44 c 0h.3","HCC44 c 0h.4","HCC44 2Gy 0h.1","HCC44 2Gy 0h.2","HCC44 2Gy 0h.3","HCC44 2Gy 0h.4",
                    "HCC44 c 2h.1","HCC44 c 2h.2","HCC44 c 2h.3","HCC44 c 2h.4","HCC44 2Gy 2h.1","HCC44 2Gy 2h.2","HCC44 2Gy 2h.3","HCC44 2Gy 2h.4",
                    "HCC44 c 6h.1","HCC44 c 6h.2","HCC44 c 6h.3","HCC44 c 6h.4","HCC44 2Gy 6h.1","HCC44 2Gy 6h.2","HCC44 2Gy 6h.3","HCC44 2Gy 6h.4",
                    "HCC44 c 24h.1","HCC44 c 24h.2","HCC44 c 24h.3","HCC44 c 24h.4","HCC44 2Gy 24h.1","HCC44 2Gy 24h.2","HCC44 2Gy 24h.3","HCC44 2Gy 24h.4",
                    "HCC44 c 48h.1","HCC44 c 48h.2","HCC44 c 48h.3","HCC44 c 48h.4","HCC44 2Gy 48h.1","HCC44 2Gy 48h.2","HCC44 2Gy 48h.3","HCC44 2Gy 48h.4",
                    "H23 c 0h.1","H23 c 0h.2","H23 c 0h.3","H23 c 0h.4","H23 2Gy 0h.1","H23 2Gy 0h.2","H23 2Gy 0h.3","H23 2Gy 0h.4",
                    "H23 c 2h.1","H23 c 2h.2","H23 c 2h.3","H23 c 2h.4","H23 2Gy 2h.1","H23 2Gy 2h.2","H23 2Gy 2h.3","H23 2Gy 2h.4",
                    "H23 c 6h.1","H23 c 6h.2","H23 c 6h.3","H23 c 6h.4","H23 2Gy 6h.1","H23 2Gy 6h.2","H23 2Gy 6h.3","H23 2Gy 6h.4",
                    "H23 c 24h.1","H23 c 24h.2","H23 c 24h.3","H23 c 24h.4","H23 2Gy 24h.1","H23 2Gy 24h.2","H23 2Gy 24h.3","H23 2Gy 24h.4",
                    "H23 c 48h.1","H23 c 48h.2","H23 c 48h.3","H23 c 48h.4","H23 2Gy 48h.1","H23 2Gy 48h.2","H23 2Gy 48h.3","H23 2Gy 48h.4")

annotation_prot <- data.frame(
  Radiation = factor(rep(c((rep("control",4)),rep("2 Gy",4)),10)),
  Radiation = factor(rep(c((rep("0h",8)),rep("2h",8),(rep("6h",8)),rep("24h",8),(rep("48h",8))),2)),
  Cell_line = factor(c(rep("HCC44",40),rep("NCI-H23",40)))
)
rownames(annotation_prot)<- names_prot
#####
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)

cls <- columns(Homo.sapiens)
cls <- cls[c(10,25,26,40)]
kts <- keytypes(Homo.sapiens)
kt <- kts[31]
ks <- keys(Homo.sapiens, keytype=kt)
res <-AnnotationDbi::select(Homo.sapiens, keys = ks, columns = cls, keytype = kt)
