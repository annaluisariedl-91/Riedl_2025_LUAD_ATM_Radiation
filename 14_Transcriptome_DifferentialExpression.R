#------------------------------------------
# preparing HCC44
#------------------------------------------
RNA_HCC44 <- Rdata %>%
  transmute(
    GeneSymbol = .[[2]],
    GeneID = .[[24]],
    ENSEMBL = .[[1]],
    gene_biotype = .[[8]],
    
    C1 = .[[9]],
    C2 = .[[10]],
    C3 = .[[11]],
    C4 = .[[12]],
    
    R1 = .[[13]],
    R2 = .[[14]],
    R3 = .[[15]],
    R4 = .[[16]]
  ) %>%
  filter(gene_biotype == "protein_coding")
#------------------------------------------
# preparing NCI-H23
#------------------------------------------
RNA_H23 <- Rdata %>%
  transmute(
    GeneSymbol = .[[2]],
    GeneID = .[[24]],
    ENSEMBL = .[[1]],
    gene_biotype = .[[8]],
    
    C1 = .[[17]],
    C2 = .[[18]],
    C3 = .[[19]],
    
    R1 = .[[20]],
    R2 = .[[21]],
    R3 = .[[22]],
    R4 = .[[23]]
  ) %>%
  filter(gene_biotype == "protein_coding")
#------------------------------------------
# binding NCI-H23 and HCC44
#-------------------------------------------
RNA_combined <- RNA_HCC44 %>%
  rename_with(
    ~ paste0("HCC44_", .x),
    matches("^[CR][1-4]$")
  ) %>%
  inner_join(
    RNA_H23 %>%
      dplyr::select(GeneID, matches("^[CR][1-4]$")) %>%
      rename_with(
        ~ paste0("H23_", .x),
        matches("^[CR][1-4]$")
      ),
    by = "GeneID"
  )
#-------------------------------------------
# HCC44 2Gy vs. control
#-------------------------------------------
ready_RHCC44_24h <- bind_cols(
  RNA_HCC44 %>%
    dplyr::select(GeneSymbol, GeneID),
  make_DE_protti(
    RNA_HCC44,
    group1 = paste0("R", 1:4),
    group2 = paste0("C", 1:4)
  ),
  RNA_HCC44 %>%
    dplyr::select(ENSEMBL, gene_biotype)
)
#-------------------------------------------
# NCI-H23 2 Gy vs. control
#-------------------------------------------
ready_RH23_24h <- bind_cols(
  RNA_H23 %>%
    dplyr::select(GeneSymbol, GeneID),
  make_DE_protti(
    RNA_H23,
    group1 = paste0("R", 1:4),
    group2 = paste0("C", 1:3)
  ),
  RNA_H23 %>%
    dplyr::select(ENSEMBL, gene_biotype)
)
#-------------------------------------------
# Radiation HCC44 vs. NCI-H23
#-------------------------------------------
ready_RRAD_24h <- bind_cols(
  RNA_combined %>%
    dplyr::select(GeneSymbol, GeneID),
  make_DE_protti(
    RNA_combined,
    group1 = paste0("HCC44_R", 1:4),
    group2 = paste0("H23_R", 1:4)
  ),
  RNA_combined %>%
    dplyr::select(ENSEMBL, gene_biotype)
)
#-------------------------------------------
# Control HCC44 vs. NCI-H23
#-------------------------------------------
ready_RCon_24h <- bind_cols(
  RNA_combined %>%
    dplyr::select(GeneSymbol, GeneID),
  make_DE_protti(
    RNA_combined,
    group1 = paste0("HCC44_C", 1:4),
    group2 = paste0("H23_C", 1:3)
  ),
  RNA_combined %>%
    dplyr::select(ENSEMBL, gene_biotype)
)
#-------------------------------------------
# saving results
#-------------------------------------------
dir.create(
  file.path(dirs$results, "14_Transcriptome_DifferentialExpression"),
  showWarnings = FALSE,
  recursive = TRUE
)
setwd(file.path(dirs$results, "14_Transcriptome_DifferentialExpression"))

DE_results <- list(
  "DE_HCC44_2Gy_vs_0Gy.csv" = ready_RHCC44_24h,
  "DE_H23_2Gy_vs_0Gy.csv" = ready_RH23_24h,
  "DE_HCC44_vs_H23_2Gy.csv" = ready_RRAD_24h,
  "DE_HCC44_vs_H23_0Gy.csv" = ready_RCon_24h
)

purrr::iwalk(
  DE_results,
  ~ write.csv(
    .x,
    file = file.path(dirs$results, "14_Transcriptome_DifferentialExpression", .y),
    row.names = FALSE,
    quote = FALSE
  )
)

RNA_HCC44_summary <- make_RNA_summary(RNA_HCC44)
write.csv(
  RNA_HCC44_summary,
  file.path(dirs$results, "14_Transcriptome_DifferentialExpression", "RNA_HCC44_24h_intensities.csv"),
  row.names = FALSE
)

RNA_H23_summary <- make_RNA_summary(RNA_H23)
write.csv(
  RNA_H23_summary,
  file.path(dirs$results, "14_Transcriptome_DifferentialExpression", "RNA_H23_24h_intensities.csv"),
  row.names = FALSE
)
#------------------------------------------
setwd(dirs$code)
