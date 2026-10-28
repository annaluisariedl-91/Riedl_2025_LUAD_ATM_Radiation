#--------------------------------------------
# defining columns
#--------------------------------------------
index <- list(
  HCC44 = list(
    `0h`  = list(c = c(18,41,68,57), rad = c(17,40,67,56)),
    `2h`  = list(c = c(22,32,70,59), rad = c(21,31,69,58)),
    `6h`  = list(c = c(35,42,79,61), rad = c(34,33,71,60)),
    `24h` = list(c = c(37,44,81,63), rad = c(36,43,80,62)),
    `48h` = list(c = c(39,46,83,65), rad = c(38,45,82,64))
  ),
  H23 = list(
    `0h`  = list(c = c(8,20,48,73), rad = c(7,19,47,72)),
    `2h`  = list(c = c(10,24,50,75), rad = c(9,23,49,74)),
    `6h`  = list(c = c(12,26,90,85), rad = c(11,25,84,76)),
    `24h` = list(c = c(14,28,66,87), rad = c(13,27,55,86)),
    `48h` = list(c = c(16,30,78,89), rad = c(15,29,77,88))
  )
)

annotation <- tibble::tibble(
  ProteinKey = as.character(seq_len(nrow(proteome_data))),
  ID = as.character(proteome_data[[91]]),
  ENSEMBL = as.character(proteome_data[[92]]),
  PG.Genes = as.character(proteome_data[[2]]),
)

#--------------------------------------------
# saving results
#--------------------------------------------
dir.create(
  file.path(dirs$results, "15_Proteome_DifferentialExpression"),
  showWarnings = FALSE,
  recursive = TRUE
)
#--------------------------------------------
# binding data
#--------------------------------------------
sample_map <- tibble::tibble(
  source_column = names(proteome_data)
) %>%
  dplyr::filter(
    stringr::str_detect(
      source_column,
      "H23_|HCC44_"
    )
  ) %>%
  dplyr::mutate(
    sample = stringr::str_remove(
      source_column,
      "^X\\.\\d+\\.+"
    )
  ) %>%
  tidyr::extract(
    sample,
    into = c("Cellline", "Treatment", "Time", "Replicate"),
    regex = "^(H23|HCC44)_(0Gy|2Gy)_(0h|2h|6h|24h|48h)_rep([1-4])$",
    remove = FALSE
  ) %>%
  dplyr::mutate(
    Replicate = as.integer(Replicate),
    Treatment = factor(Treatment, levels = c("0Gy", "2Gy")),
    Cellline = factor(Cellline, levels = c("H23", "HCC44")),
    Time = factor(Time, levels = c("0h", "2h", "6h", "24h", "48h"))
  ) %>%
  dplyr::arrange(Cellline, Time, Treatment, Replicate)

intensity_mat <- vapply(
  sample_map$source_column,
  function(col) {
    x <- proteome_data[[col]]
    if (is.factor(x)) {
      x <- as.character(x)
    }
    intensity <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & is.na(intensity))) {
      stop(
        "non-numerical intensity in column: ",
        col
      )
    }
    intensity
  },
  numeric(nrow(proteome_data))
)

proteome_long <- proteome_data %>%
  dplyr::select(
    ProteinKey,
    ENSEMBL,
    SYMBOL,
    ENTREZID,
    dplyr::all_of(sample_map$source_column)
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(ProteinKey, ENSEMBL, SYMBOL, ENTREZID),
      as.character
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_map$source_column),
    names_to = "source_column",
    values_to = "intensity",
    values_transform = list(intensity = as.numeric)
  ) %>%
  dplyr::left_join(
    sample_map,
    by = "source_column"
  ) %>%
  dplyr::mutate(
    log2_intensity = log2(
      dplyr::if_else(
        is.finite(intensity) & intensity > 0,
        intensity,
        NA_real_
      )
    )
  )

#test completeness
stopifnot(
  nrow(sample_map) == 80,
  !anyNA(sample_map$Cellline),
  !anyNA(sample_map$Treatment),
  !anyNA(sample_map$Time),
  !anyNA(sample_map$Replicate),
  !anyDuplicated(sample_map$sample)
)
#--------------------------------------------
# calculation
#--------------------------------------------
# comparison between radiadted and controll cells
for (cellline in c("H23", "HCC44")) {
  for (time in c("0h", "2h", "6h", "24h", "48h")) {
    wide <- get_wide(cellline, time)
    for (treatment in c("0Gy", "2Gy")) {
      cols <- paste0(treatment, "_", 1:4)
      raw <- wide %>%
        dplyr::select(
          ENTREZID,
          SYMBOL,
          ProteinKey,
          ENSEMBL,
          dplyr::all_of(cols)
        ) %>%
        dplyr::mutate(
          mean = safe_row_mean(
            as.matrix(dplyr::pick(dplyr::all_of(cols)))
          ),
          sd = safe_row_sd(
            as.matrix(dplyr::pick(dplyr::all_of(cols)))
          )
        ) %>%
        dplyr::rename_with(
          ~ paste0("Rep_", 1:4),
          dplyr::all_of(cols)
        )
      treatment_file <- ifelse(treatment == "0Gy", "c", "2Gy")
      write_output(
        raw,
        paste0(
          "raw_p", cellline, "_", time, "_",
          treatment_file, ".csv"
        )
      )
    }
    cvals <- as.matrix(wide[, paste0("0Gy_", 1:4)])
    rvals <- as.matrix(wide[, paste0("2Gy_", 1:4)])
    means <- wide %>%
      dplyr::select(ProteinKey, ENSEMBL,ENTREZID,SYMBOL)
    means$mean_c <- safe_row_mean(cvals)
    means$SD_c <- safe_row_sd(cvals)
    means$mean_2Gy <- safe_row_mean(rvals)
    means$SD_2Gy <- safe_row_sd(rvals)
    means$ENSEMBL <- as.character(means$ENSEMBL)
    means$SYMBOL <- as.character(means$SYMBOL)
    write_output(
      means,
      paste0("mean_p", cellline, "_", time, ".csv")
    )
    comparison_dat_intraCell <- proteome_long %>%
      dplyr::filter(
        .data$Cellline == .env$cellline,
        .data$Time == .env$time,
        Treatment %in% c("0Gy", "2Gy")
      ) %>%
      dplyr::mutate(
        condition = if_else(
          Treatment == "2Gy",
          "Rad",
          "Con"
        )
      )
    stat_results_intraCell <- run_comparison(
      comparison_dat_intraCell,
      "Rad",
      "Con"
    )
    write_output(
      stat_results_intraCell,
      paste0("stat_p", cellline, "_", time, ".csv")
    )
  }
}

# inter cell line comparison: 2/6/24/48h je 2Gy und 0Gy.
for (time in c("2h", "6h", "24h", "48h")) {
  for (treatment in c("2Gy", "0Gy")) {
    comparison_dat_inter <- proteome_long %>%
      dplyr::filter(
        .data$Time == .env$time,
        .data$Treatment == .env$treatment
      ) %>%
      dplyr::mutate(
        condition = as.character(Cellline)
      )
    prefix <- ifelse(
      treatment == "2Gy",
      "stat_pRAD_",
      "stat_pCon_"
    )
    write_output(
      run_comparison(
        comparison_dat_inter,
        "HCC44",
        "H23"
      ),
      paste0(prefix, time, ".csv")
    )
  }
}

# calculation of base value with all 0h values
comparison_0h <- proteome_long |>
  dplyr::filter(Time=='0h') |>
  dplyr::mutate(condition=Cellline)
write_output(run_comparison(comparison_0h,'HCC44','H23'),'stat_pCon_0h.csv')
#-----------------------------------------------
setwd(dirs$code)
