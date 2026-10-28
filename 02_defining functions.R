#------------------------------------------
# General functions
#------------------------------------------
#defining scaling function
range01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] - rng[1] == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

# defining base design without package
theme_prism_xcell <- function() {
  theme_classic(base_size = 9) +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      axis.ticks.length = unit(2, "mm"),
      axis.text = element_text(
        colour = "black",
        size = 9
      ),
      axis.title = element_text(
        colour = "black",
        size = 10,
        face = "bold"
      ),
      axis.title.x = element_text(
        margin = margin(t = 6)
      ),
      axis.title.y = element_text(
        margin = margin(r = 6)
      ),
      plot.title = element_text(
        colour = "black",
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(
        colour = "black",
        size = 8
      ),
      legend.key = element_blank(),
      legend.key.height = unit(4, "mm"),
      legend.spacing.y = unit(1, "mm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# additional design
theme_prism_compact <- function(
    base_size = 9,
    legend_position = "bottom"
) {
  theme_prism_custom(
    base_size = base_size,
    legend_position = legend_position
  ) +
    theme(
      legend.text = element_text(
        size = 6.5
      ),
      legend.key.width = grid::unit(
        1,
        "mm"
      ),
      legend.key.height = grid::unit(
        2.5,
        "mm"
      ),
      legend.spacing.x = grid::unit(
        1,
        "mm"
      ),
      legend.box.spacing = grid::unit(
        2,
        "mm"
      ),
      legend.direction = "horizontal",
      plot.margin = margin(
        5, 5, 5, 5
      )
    )
}
#------------------------------------------
# Functional assays | Colony forming assays (CFA)
#------------------------------------------
# binding counting data and calculation of platting efficiency
analyse_CFA_rep <- function(data,
                            rows = NULL,
                            cols = NULL,
                            remove_first_row = FALSE) {
  if (remove_first_row) {
    data <- data[-1, , drop = FALSE]
  }
  if (!is.null(rows)) {
    data <- data[rows, , drop = FALSE]
  }
  if (!is.null(cols)) {
    data <- data[, cols, drop = FALSE]
  }
  dat <- data[, 1:3, drop = FALSE]
  colnames(dat) <- c("Dose", "N0", "N")
  dat <- dat %>%
    dplyr::mutate(
      Dose = ifelse(Dose == "Control", "0", Dose),
      Dose = as.numeric(gsub(" Gy", "", Dose)),
      N0   = as.numeric(N0),
      N    = as.numeric(N)
    )
  PE <- dat %>%
    dplyr::filter(Dose == 0) %>%
    dplyr::summarise(
      PE = mean(N / N0, na.rm = TRUE)
    ) %>%
    pull(PE)
  dat <- dat %>%
    dplyr::mutate(
      SF = N / (PE * N0)
    )
  summary <- dat %>%
    group_by(Dose) %>%
    dplyr::summarise(
      N0     = first(N0),
      N_mean = mean(N, na.rm = TRUE),
      N_SD   = sd(N, na.rm = TRUE),
      N_SEM  = N_SD / sqrt(sum(!is.na(N))),
      SF     = mean(SF, na.rm = TRUE),
      .groups = "drop"
    )
  list(
    raw = dat,
    summary = summary,
    PE = PE
  )
}

# statistics calculation for counting data
make_stat_CFA <- function(data, cellline) {
  purrr::map_dfr(
    data,
    ~ .x$summary[, c("Dose", "SF")],
    .id = "experiment"
  ) %>%
    mutate(Cellline = cellline)
}

# linear quadratic model function
fit_lq_model <- function(df) {
  fit <- lm(
    mean_log_sf ~ 0 + Dose + I(Dose^2),
    data = df
  )
  c(
    alpha = -unname(coef(fit)["Dose"]),
    beta  = -unname(coef(fit)["I(Dose^2)"])
  )
}

# binding for 2 Gy data of HCC44 and NCI-H23
make_CFA_data <- function(data,
                          cellline,
                          treatment = NULL) {
  result <- purrr::map_dfr(
    data,
    ~ .x$summary[, c("Dose", "SF")],
    .id = "experiment"
  ) %>%
    mutate(
      cellline = cellline,
      experiment = factor(experiment)
    )
  if (!is.null(treatment)) {
    result <- result %>%
      mutate(
        Treatment = treatment
      )
  }
  result
}
#------------------------------------------
# function assays | CFA with AZD1390 (HCC44,NCI-H23)
#------------------------------------------
# linear quadratic model for AZD1390 and DMSO treated samples
fit_LQ_treatments <- function(data, selected_cellline) {
  lq_data <- data %>%
    filter(Cellline == selected_cellline) %>%
    group_by(Treatment, Dose) %>%
    summarise(
      mean_log_sf = mean(adj_log_sf, na.rm = TRUE),
      .groups = "drop"
    )
  coefficients <- lq_data %>%
    split(.$Treatment) %>%
    lapply(fit_lq_model)
  list(
    data = lq_data,
    coefficients = coefficients
  )
}

# calculating dose enhancement factor (DEF) of 1nM and 10nM AZD1390 for HCC44 and NCI-H23
calculate_DEF <- function(coefficients, treat_10nm, treat_1nm,sf = 0.1) {
  control <- coefficients[["DMSO"]]
  treat_10nM   <- coefficients[["10 nM AZD1390"]]
  treat_1nM   <- coefficients[["1 nM AZD1390"]]
  d_control <- dose_at_sf(
    control["alpha"],
    control["beta"],
    sf
  )
  d_treat_10nM  <- dose_at_sf(
    treat_10nM["alpha"],
    treat_10nM["beta"],
    sf
  )
  d_treat_1nM <- dose_at_sf(
    treat_1nM["alpha"],
    treat_1nM["beta"],
    sf
  )
  data.frame(
    Treatment_10nM = treat_10nM,
    Treatment_1nM = treat_1nM,
    SF = sf,
    Dose_control = d_control,
    Dose_treatment_10nM = d_treat_10nM,
    DEF_10nM = d_control / d_treat_10nM,
    Dose_treatment_1nM = d_treat_1nM,
    DEF_1nM = d_control / d_treat_1nM
  )
}

# avoiding 0 during plotting
dose_at_sf <- function(alpha, beta, sf) {
  disc <- alpha^2 - 4 * beta * log(sf)
  if (disc < 0) return(NA)
  if (sf <= 0 || sf >= 1) return(NA)
  if (abs(beta) < 1e-8) {
    if (alpha <= 0) return(NA)
    return(-log(sf) / alpha)
  }
  if (abs(alpha) < 1e-8) {
    if (beta <= 0) return(NA)
    return(sqrt(-log(sf) / beta))
  }
  dose <- (-alpha + sqrt(disc)) / (2 * beta)
  return(dose)
}

# plotting cell line with DMSO an AZD1390
plot_CFA_AZD <- function(plot_data,
                         coefficients,
                         DEF_label,
                         selected_cellline,
                         xmax,
                         significance = NULL) {
  dat <- plot_data %>%
    filter(cellline == selected_cellline)
  dseq <- seq(0, xmax, length.out = 200)
  lq_dat <- purrr::imap_dfr(
    coefficients,
    ~ data.frame(
      Treatment = .y,
      Dose = dseq,
      SF = exp(
        -(.x["alpha"] * dseq +
            .x["beta"] * dseq^2)
      )
    )
  )
    p <- ggplot(
    dat,
    aes(
      x = Dose,
      y = SF_mean,
      fill = Treatment,
      shape = Treatment
    )
  ) +
    geom_line(
      data = lq_dat,
      aes(
        x = Dose,
        y = SF,
        group = Treatment
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.6,
      linetype = "dashed"
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(SF_mean - SF_SEM, 1e-4),
        ymax = SF_mean + SF_SEM
      ),
      width = 0.08,
      linewidth = 0.5,
      colour = "black"
    ) +
    geom_point(
      size = 2.5,
      stroke = 0.6,
      colour = "black"
    ) +
    annotate(
      "text",
      x = xmax * 0.25,
      y = 100,
      label = DEF_label,
      hjust = 0,
      vjust = 1,
      size = 2.5
    ) +
    scale_fill_manual(
      values = CFA_AZD_colours
    ) +
    scale_shape_manual(
      values = CFA_AZD_shapes
    ) +
    scale_y_log10(
      limits = c(1e-4, 1000),
      breaks = c(
        1e-3,
        1
      ),
      labels = c(
        expression(10^-3),
        "1"
      )
    ) +
     scale_x_continuous(
      limits = c(0, xmax),
      breaks = c(0, 4, 8)[c(0, 4, 8) <= xmax],
      minor_breaks = seq(0, xmax, 2),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = selected_cellline,
      x = "Dose [Gy]",
      y = "Surviving Fraction",
      fill = NULL,
      shape = NULL
    ) +
    ggprism::theme_prism(
      base_size = 8
    ) +
    theme(
      axis.title = element_text(
        size = 8,
        face = "bold"
      ),
      axis.text = element_text(
        size = 7,
        colour = "black"
      ),
      plot.title = element_text(
        size = 9,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(
        size = 6
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          size = 2,
          stroke = 0.5
        )
      ),
      shape = guide_legend(
        override.aes = list(
          size = 2,
          stroke = 0.5
        )
      )
    )
    if (!is.null(significance)) {
      sig_dat <- significance %>%
        filter(Cellline == selected_cellline)
      p <- p +
        geom_text(
          data = sig_dat,
          aes(
            x = Dose,
            y = y,
            label = label
          ),
          inherit.aes = FALSE,
          size = 3
        )
    }
    return(p)
                         }

# inserting label with DEF in plot
make_DEF_label <- function(def_data, sf = 0.10) {
  dat <- def_data %>%
    filter(SF == sf)
  paste0(
    dat$Treatment,
    ": DEF",
    format(sf, trim = TRUE),
    " = ",
    round(dat$DEF, 2),
    collapse = "\n"
  )
}
#------------------------------------------
# function assays | CFA with AZD1390 (A529, HCC78, HCC827)
#------------------------------------------
# statistics for AZD1390 treatment of A529, HCC78 and HCC827 
analyse_CFA_model <- function(data, cellline) {
  dat <- data %>%
    filter(Cellline == cellline) %>%
    droplevels()
  model <- lm(
    adj_log_sf ~ Experiment + factor(Dose) * Treatment,
    data = dat
  )
  emm <- emmeans(
    model,
    ~ Treatment | Dose
  )
  posthoc <- pairs(
    emm,
    adjust = "tukey"
  )
  list(
    data = dat,
    model = model,
    anova = anova(model),
    emm = emm,
    posthoc = posthoc
  )
}

# calculating dose enhancement factor (DEF) of 10nM AZD1390 for HCC44 and NCI-H23
calc_DEF <- function(coefficients, sf = 0.1) {
  control <- coefficients[["DMSO"]]
  treat_10nM   <- coefficients[["10 nM AZD1390"]]
  d_control <- dose_at_sf(
    control["alpha"],
    control["beta"],
    sf
  )
  d_treat_10nM  <- dose_at_sf(
    treat_10nM["alpha"],
    treat_10nM["beta"],
    sf
  )
  data.frame(
    Treatment_10nM = treat_10nM,
    SF = sf,
    Dose_control = d_control,
    Dose_treatment_10nM = d_treat_10nM,
    DEF_10nM = d_control / d_treat_10nM
  )
}
#------------------------------------------
# function assays | xCelligence impedance measurment
#-----------------------------------------
# binding measurments for analysis
prepare_xCell <- function(file) {
  dat <- read.csv(file)
  dat <- dat[-(1:37), 2:8]
  names(dat) <- c(
    "time",
    "2 Gy",
    "1.5 Gy",
    "1 Gy",
    "0.5 Gy",
    "0 Gy",
    "Medium"
  )
  zeit <- stringr::str_pad(
    as.character(dat$time),
    9,
    side = "left",
    pad = "0"
  )
  dat$Time <- (
    as.numeric(substr(zeit, 1, 3)) * 3600 +
      as.numeric(substr(zeit, 5, 6)) * 60 +
      as.numeric(substr(zeit, 8, 9))
  ) / 3600
  return(dat)
}

# plotting growth impedance data
plot_xcell <- function(data,
                       title,
                       colours,
                       breaks = NULL,
                       xmax,
                       ymax,
                       ymin,
                       x_breaks = NULL,
                       x_minor = NULL,
                       y_breaks = NULL,
                       y_minor = NULL) {
  ggplot(
    data,
    aes(
      x = Time,
      y = CellIndex,
      colour = Group
    )
  ) +
    geom_point(
      size = 0.35,
      alpha = 0.9
    ) +
    scale_colour_manual(
      values = colours,
      breaks = breaks
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      minor_breaks = x_minor,
      limits = c(0, xmax),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      minor_breaks = y_minor,
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      ylim = c(ymin, ymax)
    ) +
    labs(
      title = title,
      x = "Time [h]",
      y = "Cell Index [a.u.]"
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 2.5,
          alpha = 1
        )
      )
    ) +
    theme_prism_xcell()
}
#-----------------------------------------
# function assays | xCelligence with AZD1390 and DMSO
#-----------------------------------------
# prepare impedance data of treatment with AZD1390 and DMSO combined with radiation
prepare_xcell_AZD <- function(file,
                              cellline,
                              replicate,
                              cols,
                              doses) {
  dat <- read.csv(file, sep = "\t")
  dat <- dat[-(1:15), cols]
  names(dat) <- c(
    "time",
    paste0(doses[1], "_AZD1390"),
    paste0(doses[1], "_DMSO"),
    paste0(doses[2], "_AZD1390"),
    paste0(doses[2], "_DMSO"),
    paste0(doses[3], "_AZD1390"),
    paste0(doses[3], "_DMSO"),
    "Medium"
  )
  dat$Timepoint <- seq_len(nrow(dat))
  zeit <- stringr::str_pad(
    as.character(dat$time),
    9,
    side = "left",
    pad = "0"
  )
  dat$Time <- (
    as.numeric(substr(zeit, 1, 3)) * 3600 +
      as.numeric(substr(zeit, 5, 6)) * 60 +
      as.numeric(substr(zeit, 8, 9))
  ) / 3600
  dat <- dat[, !names(dat) %in% c("time", "Medium")]
  dat <- tidyr::pivot_longer(
    dat,
    cols = -c(Time, Timepoint),
    names_to = c("Dose", "Treatment"),
    names_pattern = "(.+)_(AZD1390|DMSO)",
    values_to = "CellIndex",
    values_transform = list(
      CellIndex = as.numeric
    )
  )
  dat %>%
    dplyr::mutate(
      Cellline = cellline,
      Replicate = replicate,
      Dose = gsub("Gy", " Gy", Dose),
      Dose = gsub("^05", "0.5", Dose)
    )
}
#------------------------------
# function assays | MTS Viability Assay
#-----------------------------
# loading MTS raw data
prepare_growth_data <- function(data, cellline, dose, cols) {
  times <- c(0, 24, 48, 72, 96, 120)
  rows  <- c(1:10, 13:22, 25:34)
  purrr::map_dfr(
    seq_along(data),
    function(rep) {
      dat <- data[[rep]]
      purrr::map2_dfr(
        cols,
        times,
        function(col, time) {
          tibble(
            Intensity = as.numeric(dat[rows, col]),
            Time = time,
            Plate = rep(1:3, each = 10),
            BioReplicate = rep,
            Cellline = cellline,
            Dose = dose
          )
        }
      )
    }
  )
}

# design for plotting
theme_MTS_growth <- function() {
  theme_classic(base_size = 7) +
    theme(
      axis.line = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      axis.ticks.length = unit(2, "mm"),
      axis.text = element_text(
        size = 7,
        colour = "black"
      ),
      axis.title = element_text(
        size = 8,
        face = "bold",
        margin = margin(8, 4, 8, 4)
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      plot.margin = margin(10, 10, 2, 10)
    )
}
#------------------------------
# function assays | ATM inhibition
#-----------------------------
# reading, sorting and binding of raw data
read_replicates_MTS_AZD1390 <- function(cellline) {
  lapply(
    paste0("ATMinhibition_", cellline, "_rep", 1:4, ".csv"),
    function(file) {
      dat <- read.csv(file)
      data.frame(lapply(dat, as.numeric))
    }
  )
}

# calculation
prepare_MTS_ATMi <- function(data, cellline, treatment, dose) {
  times <- c(0, 24, 48, 72, 96, 120)
  rows <- if (dose == "2 Gy") {
      3:7
    } else if (dose == "0 Gy") {
      11:15
    } else {
      stop("Dose muss '0 Gy' oder '2 Gy' sein.")
    }
  cols <- if (treatment == "DMSO") {
      3:8
    } else if (treatment == "AZD1390") {
      11:16
    } else {
      stop("Treatment muss 'DMSO' oder 'AZD1390' sein.")
    }
  purrr::map_dfr(
    seq_along(data),
    function(rep) {
      dat <- data[[rep]]
      purrr::map2_dfr(
        cols,
        times,
        function(col, time) {
          tibble(
            Intensity = as.numeric(dat[[col]][rows]),
            Time = time,
            Treatment = treatment,
            Dose = dose,
            BioReplicate = rep,
            Cellline = cellline
        )
      }
    )
  }
) %>%
  filter(!is.na(Intensity))
}
#------------------------------
# function assays | MTS IC50 Assay
#-----------------------------------------
# binding data of inhibition data, unradiated 
read_IC50 <- function(cellline) {
  dplyr::bind_rows(
    lapply(1:3, function(rep) {
      file <- paste0(
        "IC50_AZD1390_", cellline, "_rep", rep, ".csv"
      )
      dat <- read.delim(
        file,
        header = FALSE,
        dec = ",",
        fileEncoding = "Windows-1252",
        skip = 0,
        blank.lines.skip = FALSE
      )
      block_rows <- list(
        3:8,
        12:17,
        21:26
      )
      dplyr::bind_rows(
        lapply(seq_along(block_rows), function(b) {
          rows <- block_rows[[b]]
          block_data <- dat[rows, 1:11]
          names(block_data) <- c(
            "Concentration",
            paste0("Measurement_", 1:10)
          )
          block_data %>%
            tidyr::pivot_longer(
              cols = starts_with("Measurement_"),
              names_to = "TechnicalReplicate",
              values_to = "Intensity"
            ) %>%
            dplyr::mutate(
              Cellline = cellline,
              BioReplicate = rep,
              Block = b,
              Concentration = as.numeric(Concentration),
              Intensity = as.numeric(Intensity),
              TechnicalReplicate = as.integer(
                sub("Measurement_", "", TechnicalReplicate)
              )
            )
        })
      )
    })
  )
}

# plotting of IC50 calculations
MTS_IC50_plotting <- function(
    global_IC_plotting,
    IC50_predictions,
    ic50_HCC44 = NULL,
    ic50_H23 = NULL,
    colours = c(
      "HCC44" = "#EA4025",
      "NCI-H23" = "#5F8F29"
    ),
    ymax = 120
) {
  p <- ggplot(
    global_IC_plotting,
    aes(
      x = Dose_log,
      y = mean_rel_Intensity_m,
      fill = Cellline
    )
  ) +
    # reference line at 50 %
    geom_hline(
      yintercept = 50,
      linetype = "dashed",
      colour = "grey50",
      linewidth = 0.5
    ) +
    # modellig of dosis-reaction curves
    geom_line(
      data = IC50_predictions,
      aes(
        x = Dose_log,
        y = Prediction,
        colour = Cellline,
        group = Cellline
      ),
      inherit.aes = FALSE,
      linewidth = 0.8
    ) +
    # error bars
    geom_errorbar(
      aes(
        ymin = mean_rel_Intensity_m - sem_rel,
        ymax = mean_rel_Intensity_m + sem_rel
      ),
      width = 0.1,
      linewidth = 0.45,
      colour = "black"
    ) +
    geom_point(
      shape = 21,
      size = 2,
      stroke = 0.5,
      colour = "black"
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_colour_manual(
      values = colours
    ) +
    scale_x_continuous(
      name = expression(log[10] * "(AZD1390 concentration [µM])")
    ) +
    scale_y_continuous(
      name = "Relative cell viability (%)",
      limits = c(0, ymax),
      breaks = seq(0, ymax, by = 20)
    ) +
    ggprism::theme_prism(
      base_size = 8
    ) +
    theme(
      axis.title = element_text(
        size = 8,
        face = "bold"
      ),
      axis.text = element_text(
        size = 7,
        colour = "black"
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 7)
    ) +
    guides(
      colour = "none",
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          size = 2,
          colour = "black"
        )
      )
    )
  # IC50 HCC44
  if (!is.null(ic50_HCC44) && is.finite(ic50_HCC44)) {
    p <- p +
      geom_vline(
        xintercept = ic50_HCC44,
        linetype = "dotted",
        colour = colours["HCC44"],
        linewidth = 0.6
      )
  }
  # IC50 NCI-H23
  if (!is.null(ic50_H23) && is.finite(ic50_H23)) {
    p <- p +
      geom_vline(
        xintercept = ic50_H23,
        linetype = "dotted",
        colour = colours["NCI-H23"],
        linewidth = 0.6
      )
  }
  return(p)
}
#-----------------------------------------
# Transcriptome | Differential Expression, Transcriptome
#-----------------------------------------
# calculation of differential expression and statistics is based on protti package
make_DE_protti <- function(data, group1, group2) {
  x1 <- log2(as.matrix(data[, group1]))
  x2 <- log2(as.matrix(data[, group2]))
  x1[!is.finite(x1)] <- NA_real_
  x2[!is.finite(x2)] <- NA_real_
  mean1 <- rowMeans(x1, na.rm = TRUE)
  mean2 <- rowMeans(x2, na.rm = TRUE)
  sd1 <- apply(x1, 1, sd, na.rm = TRUE)
  sd2 <- apply(x2, 1, sd, na.rm = TRUE)
  n1 <- rowSums(is.finite(x1))
  n2 <- rowSums(is.finite(x2))
  res <- ttest_protti(
    mean1 = mean1,
    mean2 = mean2,
    sd1 = sd1,
    sd2 = sd2,
    n1 = n1,
    n2 = n2,
    log_values = TRUE
  )
  res$adj_pVal <- p.adjust(
    res$pval,
    method = "BH"
  )
  res$mean_group1 <- mean1
  res$mean_group2 <- mean2
  colnames(res)[colnames(res) == "diff"] <- "log2FoldChange"
  return(res)
}

# generating summary of transcriptome raw data
make_RNA_summary <- function(data) {
  control_cols <- grep("^C[1-4]$", names(data), value = TRUE)
  rad_cols <- grep("^R[1-4]$", names(data), value = TRUE)
  data %>%
    mutate(
      mean_C = rowMeans(pick(all_of(control_cols)), na.rm = TRUE),
      SD_C = apply(pick(all_of(control_cols)), 1, sd, na.rm = TRUE),
      mean_2Gy = rowMeans(pick(all_of(rad_cols)), na.rm = TRUE),
      SD_2Gy = apply(pick(all_of(rad_cols)), 1, sd, na.rm = TRUE)
  )
}
#-----------------------------------------
# Proteome | Differential Expression, Proteome
#-----------------------------------------
# binding data for calculation
get_wide <- function(cellline, time) {
  proteome_long %>%
    dplyr::filter(
      Cellline == cellline,
      Time == time
    ) %>%
    dplyr::mutate(
      Condition = paste0(Treatment, "_", Replicate)
    ) %>%
    dplyr::select(
      ProteinKey,
      SYMBOL,
      ENTREZID,
      ENSEMBL,
      Condition,
      log2_intensity
    ) %>%
    tidyr::pivot_wider(
      names_from = Condition,
      values_from = log2_intensity
    ) %>%
    dplyr::left_join(
      annotation %>%
        dplyr::select(
          ProteinKey,
          ID
        ),
      by = "ProteinKey"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(ENSEMBL, SYMBOL, ENTREZID),
        as.character
      )
    )
}

# calculation of mean and standard deviation for raw summary data
safe_row_mean <- function(x) { 
  z <- rowMeans(x,na.rm=TRUE); z[is.nan(z)] <- NA_real_; z }
safe_row_sd <- function(x) apply(x,1,function(v) if(sum(is.finite(v))>=2) stats::sd(v[is.finite(v)]) else NA_real_)

# saving sorted raw and mean data
write_output <- function(x, filename) {
  utils::write.csv(x, file=file.path(dirs$results, "15_Proteome_DifferentialExpression",filename), row.names=FALSE, na='')
}

# calculation of differential expression and statistics is based on protti package
run_comparison <- function(dat, first, second) {
  dat <- dat |>
    dplyr::filter(condition %in% c(first, second)) |>
    dplyr::mutate(
      ProteinKey = as.character(ProteinKey),
      condition = factor(condition, levels = c(second, first))
    )
  miss <- protti::assign_missingness(
    data = dat,
    sample = sample,
    condition = condition,
    grouping = ProteinKey,
    intensity = log2_intensity,
    ref_condition = second,
    retain_columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
  )
  ans <- protti::calculate_diff_abundance(
    data = miss,
    sample = sample,
    condition = condition,
    grouping = ProteinKey,
    intensity_log2 = log2_intensity,
    missingness = missingness,
    comparison = comparison,
    method = "moderated_t-test",
    p_adj_method = "BH",
    filter_NA_missingness = FALSE,
    retain_columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
  )
  if (!all(c('ProteinKey','diff','pval','adj_pval') %in% names(ans))) {
    stop('unexpected Protti result: ',paste(names(ans),collapse=', '))
  }
  ans <- ans |> dplyr::filter(as.character(comparison)==paste0(first,'_vs_',second))
  if(anyDuplicated(ans$ProteinKey)) stop('multiple Protti results per ProteinKey')
  annotation |>
    dplyr::select(
      ProteinKey,
      PG.Genes
    ) |>
    dplyr::left_join(
      ans,
      by = "ProteinKey"
    ) |>
    dplyr::mutate(
      dplyr::across(
        c(SYMBOL, ENTREZID, ENSEMBL),
        as.character
      )
    ) |>
    dplyr::rename(
      log2FoldChange = diff,
      adj_pVal = adj_pval
    )  |>
  dplyr::select(
    -comparison,
    -missingness,
    -CI_2.5,
    -CI_97.5,
    -B,
    -n_obs,
  )
}
#------------------------------------------
# Omics data | correlation between Proteome and Transcriptome data
#------------------------------------------
# loading differential protein expression data and renaming data frame
load_stat_data_p <- function(path, mode = "both") {
  mode <- match.arg(mode, c("stat", "ready", "both"))
  files <- list.files(
    path = path,
    pattern = "^stat_p.*\\.csv$",
    full.names = TRUE
  )
  for (file in files) {
    name <- tools::file_path_sans_ext(basename(file))
    data <- read.csv(
      file,
      check.names = FALSE
    )
    if (mode %in% c("stat", "both")) {
      assign(name, data, envir = .GlobalEnv)
    }
    if (mode %in% c("ready", "both")) {
      ready_name <- sub("^stat_", "ready_", name)
      assign(ready_name, data, envir = .GlobalEnv)
    }
  }
  invisible(NULL)
}

# loading transcriptome differential expression data and renaming data frame
load_stat_data_R <- function(path, mode = "both") {
  mode <- match.arg(mode, c("stat", "ready", "both"))
  files <- list.files(
    path = path,
    pattern = "^stat_.*\\.csv$",
    full.names = TRUE
  )
  for (file in files) {
    name <- tools::file_path_sans_ext(basename(file))
    data <- read.csv(
      file,
      check.names = FALSE
    )
    if (mode %in% c("stat", "both")) {
      assign(name, data, envir = .GlobalEnv)
    }
    if (mode %in% c("ready", "both")) {
      ready_name <- sub("^stat_", "ready_", name)
      assign(ready_name, data, envir = .GlobalEnv)
    }
    message("Geladen: ", name)
  }
  invisible(NULL)
}

# plotting correlation between Proteome and Transcriptome with Pearson correlation model
plot_correlation <- function(data, title) {
  # linear model
  model <- lm(
    Proteome ~ Transcriptome,
    data = data
  )
  # R² of Pearson correlation
  r2 <- summary(model)$r.squared
  ggplot(
    data,
    aes(
      x = Proteome,
      y = Transcriptome,
      colour = significance
    )
  ) +
    #base line at 0
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.75,
    colour = "grey40"
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      linewidth = 0.75,
      colour = "grey40"
    ) +
  # mean intensities Proteome vs Transcriptome
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  # linear regression
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = "#D53E4F",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 2,
    y = 4,
    label = paste0(
      "R² = ",
      round(r2, 5)
    ),
    hjust = 0,
    size = 3
  ) +
  # color coding
  scale_colour_manual(
    values = significance_colours,
    breaks = c(
      "significant in proteome",
      "significant in transcriptome",
      "significant in both"
    ),
    drop = FALSE,
    name = NULL
  ) +
  # axis
  scale_x_continuous(
    breaks = c(
      -7.5,
      0,
      7.5
    ),
    minor_breaks = seq(
      -5,
      5,
      by = 2.5
    )
  ) +
    scale_y_continuous(
      breaks = c(
        -2.5,
        0,
        2.5
      ),
      minor_breaks = seq(
        -2.5,
        2.5,
        by =1
      )
    ) +
    coord_cartesian(
      xlim = c(-7.5, 7.5),
      ylim = c(-2.5, 2.5)
    ) +
  # annotations
  labs(
    title = title,
    x = expression(
      Transcriptome~log[2]~"(2 Gy / Control)"
    ),
    y = expression(
      Proteome~log[2]~"(2 Gy / Control)"
    )
  ) +
  theme_prism_xcell() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(
        size = 8
      ),
      legend.key.width = grid::unit(
        4,
        "mm"
      ),
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(
        size = 10,
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        size = 9,
        colour = "black"
      )
    )
}
#------------------------------------------
# Omics data | Volcano plots
#------------------------------------------
# preparing data for plotting, volcano plot based on EnhancedVolcano
plot_volcano <- function(
    data,
    title,
    x_title,
    FC_cutoff = 1,
    p_cutoff = 0.05,
    n_labels = 2,
    xlim = c(-5, 5),
    ylim = c(0, 7.5)
) {
  # prepare data
  data <- data %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      pval = as.numeric(pval)
    ) %>%
    filter(
      is.finite(log2FoldChange),
      is.finite(pval),
      pval > 0
    )
  # identifying significant changes 
  sig_data <- data %>%
    filter(
      pval < p_cutoff,
      abs(log2FoldChange) >= FC_cutoff
    )
  # identifying 10 most differential regulated hits
    top_up <- sig_data %>%
    filter(log2FoldChange > 0) %>%
    slice_max(
      order_by = log2FoldChange,
      n = n_labels,
      with_ties = FALSE
    )
  top_down <- sig_data %>%
    filter(log2FoldChange < 0) %>%
    slice_min(
      order_by = log2FoldChange,
      n = n_labels,
      with_ties = FALSE
    )
  # Automatic annotation, column 1 contains names
  selected_genes <- c(
    top_up[[1]],
    top_down[[1]]
  )
  # plotting
  p <- EnhancedVolcano(
    data,
    lab = data[[1]],
    x = "log2FoldChange",
    y = "pval",
    FCcutoff = FC_cutoff,
    pCutoff = p_cutoff,
    selectLab = selected_genes,
    title = title,
    subtitle = NULL,
    caption = NULL,
    xlab = x_title,
    ylab = expression(-log[10] * "(p-value)"),
    pointSize = 1.5,
    labSize = 3,
    # NS / just FC / just p / FC + p
    col = c(
      "grey40",
      "grey40",
      "grey40",
      "#D53E4F"
    ),
    colAlpha = 0.8,
    # Cutoff lines
    cutoffLineType = "dotted",
    cutoffLineWidth = 0.75,
    cutoffLineCol = "black",
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    colConnectors = "grey40",
    border = "full",
    borderWidth = 0.5,
    borderColour = "black",
    legendPosition = "none"
  )
  p <- p +
    coord_cartesian(
      xlim = xlim,
      ylim = ylim,
      clip = "off"
    ) +
    theme_prism_xcell() +
    theme(
      axis.line = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),
      axis.ticks.length =
        grid::unit(1.5, "mm"),
      axis.text = element_text(
        size = 9,
        colour = "black"
      ),
      axis.title = element_text(
        size = 10,
        face = "bold",
        colour = "black"
      ),
      plot.title = element_text(
        size = 11,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "none",
      plot.margin = margin(
        10, 10, 5, 10
      )
    )
  return(p)
}
#------------------------------------------
# Transcriptome | Quality Control
#------------------------------------------
# loading and binding RNA data
load_RNA <- function(cellline, time) {
  file <- file.path(
    dirs$results,
    "14_Transcriptome_DifferentialExpression",
    paste0("RNA_", cellline, "_", time, "_intensities.csv")
  )
  read.csv(file) %>%
    dplyr::filter(gene_biotype == "protein_coding")
}

# prepare RPKM data for plotting
prepare_RPKM <- function(data, sample_cols, cellline, treatment) {
  data.frame(
    gene_name = data$gene_name,
    gene_biotype = data$gene_biotype,
    Expression = log2(
      rowMeans(
        data[, sample_cols],
        na.rm = TRUE
      ) + 1
    ),
    Cellline = cellline,
    Treatment = treatment
  ) %>%
    filter(
      is.finite(Expression),
      !is.na(gene_biotype)
    )
}

# prepare TPM data for plotting
prepare_TPM <- function(data, sample_cols, cellline, treatment) {
  data.frame(
    gene_name = data$gene_name,
    gene_biotype = data$gene_biotype,
    Expression = log2(
      rowMeans(
        data[, sample_cols],
        na.rm = TRUE
      ) + 1
    ),
    Cellline = cellline,
    Treatment = treatment
  ) %>%
    filter(
      is.finite(Expression),
      !is.na(gene_biotype)
    )
}

# calculate coefficient of variation for transcriptome data
calculate_cov_RNA <- function(data, cellline, treatment, time) {
  data <- as.data.frame(data)
  data <- data %>%
    mutate(
      across(
        everything(),
        ~ as.numeric(as.character(.x))
      )
    )
  # log2 values back to linear scale
  data_linear <- 2^data
  mean_value <- rowMeans(
    data_linear,
    na.rm = TRUE
  )
  sd_value <- apply(
    data_linear,
    1,
    sd,
    na.rm = TRUE
  )
  cv_value <- (sd_value / mean_value) * 100
  data.frame(
    mean = mean_value,
    SD = sd_value,
    CV = cv_value,
    Cellline = cellline,
    Treatment = treatment,
    Time = time
  ) %>%
    filter(is.finite(CV))
}

# calculate density for COV data
make_density <- function(data) {
  dens <- density(
    data$CV,
    na.rm = TRUE
  )
  data.frame(
    CV = dens$x,
    Density = dens$y * nrow(data)
  )
}

# plotting COV data 
plot_cov_RNA <- function(data) {
  data <- data %>%
    mutate(
      Group = paste(Cellline, Treatment),
      Group = factor(
        Group,
        levels = c(
          "HCC44 ctrl",
          "HCC44 2 Gy",
          "NCI-H23 ctrl",
          "NCI-H23 2 Gy"
        )
      )
    )
  RNA_colours <- c(
    "HCC44 ctrl"   = "#4294F7",
    "HCC44 2 Gy"   = "#EA4025",
    "NCI-H23 ctrl" = "#5F8F29",
    "NCI-H23 2 Gy" = "#F09837"
  )
  ggplot(
    data,
    aes(
      x = CV,
      colour = Group,
      fill = Group
    )
  ) +
    geom_density(
      alpha = 0.15,
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = RNA_colours
    ) +
    scale_fill_manual(
      values = RNA_colours
    ) +
    scale_x_continuous(
      breaks = seq(0, 500, by = 100),
      minor_breaks = seq(50, 450, by = 100),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      xlim = c(-20, 200),
      clip = "off"
    ) +
    labs(
      x = "Coefficient of variation [%]",
      y = "Density",
      colour = NULL,
      fill = NULL
    ) +
    theme_prism_xcell() +
    theme(
      axis.title = element_text(
        size = 7,
        face = "bold"
      ),
      axis.text = element_text(
        size = 7,
        colour = "black"
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 4),
      legend.key.width = grid::unit(3, "mm"),
      legend.key.height = grid::unit(2, "mm"),
      plot.margin = margin(
        5, 10, 5, 5
      )
    ) +
    guides(
      fill = "none",
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          linewidth = 0.4,
          alpha = 1
        )
      )
    )}

#------------------------------------------
# Proteome | Quality Control
#------------------------------------------
# loading proteome mean data
load_prot_mean <- function(cellline, time) {
  file <- file.path(
    dirs$results,
    "15_Proteome_DifferentialExpression",
    paste0("mean_p", cellline, "_", time, ".csv")
  )
  read.csv(file, check.names = FALSE)
}

# loading proteome raw data
load_prot_raw <- function(cellline, time, Treatment) {
  file <- file.path(
    dirs$results,
    "15_Proteome_DifferentialExpression",
    paste0("raw_p", cellline, "_", time, "_", Treatment, ".csv")
  )
  read.csv(file, check.names = FALSE)
}

# binding data for box plot
make_box_all <- function(cellline) {
  times <- c("0h", "24h")
  treatments <- c("c", "2Gy")
  names <- as.vector(outer(
    times, treatments,
    function(time, treatment)
      paste0("raw_", cellline, "_", time, "_", treatment)
  ))
  box <- do.call(
    cbind,
    lapply(names, function(x) get(x)[, c(2:3,5:8)])
  )
  labels <- sub("^raw_", "", names)
  labels <- gsub("_", " ", labels)
  list(box = box, labels = labels)
}

# plotting of coefficient of calculation for proteome data
plot_cov <- function(
    data,
    cellline,
    treatment,
    title
) {
  plot_data <- data %>%
    filter(
      Cellline == cellline,
      Treatment == treatment
    )
  ggplot(
    plot_data,
    aes(
      x = log_CV,
      y = Density,
      colour = Time,
      group = Time
    )
  ) +
    geom_line(
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = log10(20 + 1),
      linetype = "dashed",
      linewidth = 0.4,
      colour = "grey40"
    )+
    scale_colour_manual(
      values = time_colours,
      breaks = c(
        "0 h",
        "2 h",
        "6 h",
        "24 h",
        "48 h"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = c(-2.5, 0, 2.5),
      minor_breaks = c(-1.25, 1.25),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = c(
        0,
        150,
        300,
        450
      ),
      minor_breaks = c(
        75,
        225,
        375
      ),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      xlim = c(-2.5, 2.5),
      ylim = c(0, 450),
      clip = "off"
    ) +
    labs(
      title = title,
      x = expression(log[10] * "(CV + 1)"),
      y = "Density"
    ) +
    theme_prism_xcell() +
    theme(
      plot.title = element_text(
        size = 9,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(
        size = 8,
        face = "bold"
      ),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.4
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.25
      ),
      axis.text = element_text(
        size = 7,
        colour = "black"
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box.just = "left",
      legend.text = element_text(
        size =8
      ),
      legend.key.width =
        grid::unit(4, "mm"),
      legend.key.height =
        grid::unit(2.5, "mm"),
      legend.key.spacing.x =
        grid::unit(1, "mm"),
      legend.spacing.x =
        grid::unit(0, "mm"),
      legend.box.spacing =
        grid::unit(2, "mm"),
      plot.margin = margin(
        5, 10, 5, 10
      )
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        override.aes = list(
          linewidth = 1
        )
      )
    )
}

# calculation Coefficient of Variation
calculate_cov <- function(data, cellline, treatment, time) {
  values <- data[, 4:7]
  values <- values %>%
    mutate(
      across(
        everything(),
        ~ as.numeric(as.character(.x))
      )
    )
  mean_value <- rowMeans(
    values,
    na.rm = TRUE
  )
  sd_value <- apply(
    values,
    1,
    sd,
    na.rm = TRUE
  )
  cv_value <- (
    sd_value / mean_value
  ) * 100
  cov_data <- data.frame(
    mean = mean_value,
    SD = sd_value,
    CV = cv_value,
    Cellline = cellline,
    Treatment = treatment,
    Time = time
  )
  cov_data <- cov_data %>%
    filter(
      !is.na(CV),
      is.finite(CV)
    )
  return(cov_data)
}

# calculating density of CoV
make_density <- function(data) {
  dens <- density(
    data$CV,
    na.rm = TRUE
  )
  data.frame(
    CV = dens$x,
    Density = dens$y * nrow(data)
  )
}
#------------------------------------------
# Proteome | Spline Regression
#------------------------------------------
# binding data for spline regression
prepare_spline_data <- function(cellline, treatment) {
  times <- c("0h", "2h", "6h", "24h", "48h")
  spline_data <- lapply(times, function(time) {
    # load dataframe
    name <- paste0(
      "raw_", cellline, "_", time, "_", treatment
    )
    data <- get(name, envir = .GlobalEnv)
    data %>%
      dplyr::select(
        ProteinKey,
        SYMBOL,
        ENTREZID,
        ENSEMBL,
        dplyr::all_of(paste0("Rep_", 1:4))
      ) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("Rep_"),
        names_to = "Replicate",
        values_to = "Intensity"
      ) %>%
      dplyr::mutate(
        Cellline = cellline,
        Treatment = treatment,
        Time = as.numeric(sub("h", "", time))
      )
  }) %>%
    dplyr::bind_rows()
  return(spline_data)
}
#------------------------------------------
# Proteome | Over-Representation Analysis (ORA)
#------------------------------------------
# binding data with annotation, log2FoldChange and pValue
make_timeData <- function(data_0h, data_2h, data_6h, data_24h, data_48h) {
  cbind(
    data_0h[, c(1, 3, 4, 5, 6,9)],
    data_2h[, c(6,9)],
    data_6h[, c(6,9)],
    data_24h[, c(6,9)],
    data_48h[, c(6,9)]
  )
}

# binding data with annotation, log2FoldChange and pValue minus 0h time point
make_timeData_Rad <- function(data_2h, data_6h, data_24h, data_48h) {
  cbind(
    data_2h[, c(1, 3, 4, 5, 6,9)],
    data_6h[, c(6,9)],
    data_24h[, c(6,9)],
    data_48h[, c(6,9)]
  )
}

# loading and bind raw data
make_timeData_all <- function(cellline) {
  times <- c("0h", "2h", "6h", "24h", "48h")
  treatments <- c("c", "2Gy")
  datasets <- lapply(times, function(time) {
    lapply(treatments, function(treatment) {
      get(paste0(
        "raw_", cellline, "_", time, "_", treatment
      ))
    })
  })
  datasets <- unlist(datasets, recursive = FALSE)
  annotation <- datasets[[1]][, 1:4]
  intensities <- lapply(datasets, function(x) {
    log2(x[, 5:8])
  })
  cbind(annotation, do.call(cbind, intensities))
}

# loading and binding mean data
cellline <- c("HCC44", "NCI-H23")
make_meanData <- function(cellline, first = FALSE) {
  times <- c("0h", "2h", "6h", "24h", "48h")
  # import CSV Data 
  datasets <- lapply(times, function(time) {
    name <- paste0("mean_p", cellline, "_", time)
    data <- read.csv(paste0(name, ".csv"))
    colnames(data)[5] <- paste0(name, "_c")
    colnames(data)[7] <- paste0(name, "_2Gy")
    return(data)
  })
  # Annotation
  if (first) {
    first_data <- datasets[[1]][, c(1:4, 5, 7)]
  } else {
    first_data <- datasets[[1]][, c(5, 7)]
  }
  remaining <- lapply(datasets[-1], function(x) {
    x[, c(5, 7)]
  })
  meanData <- cbind(first_data, do.call(cbind, remaining))
  return(meanData)
}

# binding data with ENTREZID and overrepressentation analysis
make_GO_ORA <- function(cluster_data, time_data, top_n = 10) {
  # connect Cluster with annotation 
  GO_data <- merge(
    cluster_data,
    time_data[, 1:2],
    by = "ProteinKey"
  )
  # transform SYMBOL in ENTREZIDs
  entrez_ids <- bitr(
    GO_data[, 2],
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  # GO-Enrichment
  res_GO <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  # cleansing results
  go_data <- as.data.frame(res_GO) %>%
    distinct(geneID, .keep_all = TRUE) %>%
    arrange(desc(Count)) %>%
    slice_head(n = top_n) %>%
    mutate(
      p.adjust_log = -log10(p.adjust),
      gene_ratio_numeric = as.numeric(
        sub("/.*", "", GeneRatio)
      ) / as.numeric(
        sub(".*/", "", GeneRatio)
      ),
      GO_position = rev(seq_len(n()))
    )
  return(list(
    GO_data = GO_data,
    entrez_ids = entrez_ids,
    enrichment = res_GO,
    plot_data = go_data
  ))
}
#------------------------------------------
# Proteome | Analysis of Cell cycle proteins
#------------------------------------------
# binding proteins for time course analysis
make_protein_timecourse <- function(data, gene) {
  # select proteins
  Symbol <- data[data[[4]] == gene, ]
  if (nrow(Symbol) != 1) {
    stop("Gene was not identified as unique")
  }
  H23   <- as.numeric(unlist(Symbol[1, seq(6, 14, 2)]))
  HCC44 <- as.numeric(unlist(Symbol[1, seq(16, 24, 2)]))
  data.frame(
    Time = rep(c(0, 2, 6, 24, 48), 2),
    Intensity = c(
      HCC44 / HCC44[1],
      H23 / H23[1]
    ),
    Cellline = factor(
      rep(c("HCC44 2 Gy", "NCI-H23 2 Gy"), each = 5),
      levels = c("HCC44 2 Gy", "NCI-H23 2 Gy")
    )
  )
}

# preparing data
make_cluster_data <- function(data, cols, cellline) {
  data %>%
    transmute(
      Protein = .[[4]],
      across(all_of(cols))
    ) %>%
    pivot_longer(
      cols = -Protein,
      names_to = "Sample",
      values_to = "Intensity"
    ) %>%
    group_by(Protein) %>%
    mutate(
      Time = rep(time_points, length.out = n()),
      Intensity = as.numeric(Intensity),
      Intensity = Intensity / Intensity[Time == 0],
      Cellline = cellline
    ) %>%
    ungroup()
}
#------------------------------------------
# Proteome | Protein-Protein Interaction networks
#------------------------------------------
# bind data for PIP analysis and PIP analysis
make_temporal_PPI <- function(
    protein_data,
    fc_data,
    string_db,
    times = c("0h", "2h", "6h", "24h", "48h")
) {
  # select proteins
  ppi_input <- protein_data %>%
    filter(
      !is.na(SYMBOL),
      trimws(SYMBOL) != ""
    ) %>%
    transmute(Gene = SYMBOL) %>%
    distinct()
  # STRING-Mapping
  gene_groups <- split(
    ppi_input$Gene,
    ceiling(seq_along(ppi_input$Gene) / 50)
  )
  ppi_mapped <- lapply(gene_groups, function(genes) {
    string_db$map(
      data.frame(Gene = genes),
      "Gene",
      removeUnmappedRows = TRUE
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct()
  # call Interactions 
  ppi_interactions <- string_db$get_interactions(
    ppi_mapped$STRING_id
  ) %>%
    filter(
      from %in% ppi_mapped$STRING_id,
      to %in% ppi_mapped$STRING_id,
      from != to
    )
  # convert STRING-IDs to gene SYMBOL
  id_to_gene <- ppi_mapped %>%
    select(STRING_id, Gene) %>%
    distinct()
  ppi_edges <- ppi_interactions %>%
    left_join(
      id_to_gene,
      by = c("from" = "STRING_id")
    ) %>%
    rename(from_gene = Gene) %>%
    left_join(
      id_to_gene,
      by = c("to" = "STRING_id")
    ) %>%
    rename(to_gene = Gene) %>%
    filter(
      !is.na(from_gene),
      !is.na(to_gene),
      from_gene != to_gene
    ) %>%
    distinct(from_gene, to_gene, .keep_all = TRUE)
  # making network
  ppi_graph <- graph_from_data_frame(
    ppi_edges %>%
      select(from_gene, to_gene),
    directed = FALSE,
    vertices = data.frame(name = ppi_input$Gene)
  )
  # statistics
  network_stats <- data.frame(
    Gene = V(ppi_graph)$name,
    Degree = degree(ppi_graph),
    Betweenness = betweenness(
      ppi_graph,
      normalized = TRUE
    )
  )
  # add own experimental values
  plot_data <- tidyr::expand_grid(
    Gene = V(ppi_graph)$name,
    Time = times
  ) %>%
    left_join(
      fc_data %>%
        transmute(
          Gene = SYMBOL,
          Time = as.character(Time),
          p.adjust = pval,
          log2FoldChange = as.numeric(log2FoldChange)
        ),
      by = c("Gene", "Time")
    ) %>%
    left_join(
      network_stats,
      by = "Gene"
    ) %>%
    mutate(
      Time = factor(Time, levels = times)
    )
  return(list(
    graph = ppi_graph,
    stats = network_stats,
    plot_data = plot_data,
    mapped = ppi_mapped
  ))
}

# bind data for PIP analysis and PIP analysis for radiated data
make_temporal_PPI_Rad <- function(
    protein_data,
    fc_data,
    string_db,
    times = c("2h", "6h", "24h", "48h")
) {
  # select proteins
  ppi_input <- protein_data %>%
    filter(
      !is.na(SYMBOL),
      trimws(SYMBOL) != ""
    ) %>%
    transmute(Gene = SYMBOL) %>%
    distinct()
  gene_groups <- split(
    ppi_input$Gene,
    ceiling(seq_along(ppi_input$Gene) / 50)
  )
  ppi_mapped <- lapply(gene_groups, function(genes) {
    string_db$map(
      data.frame(Gene = genes),
      "Gene",
      removeUnmappedRows = TRUE
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct()
  # convert STRING-IDs to gene SYMBOL
  ppi_interactions <- string_db$get_interactions(
    ppi_mapped$STRING_id
  ) %>%
    filter(
      from %in% ppi_mapped$STRING_id,
      to %in% ppi_mapped$STRING_id,
      from != to
    )
  # convert STRING-IDs to gene SYMBOL
  id_to_gene <- ppi_mapped %>%
    select(STRING_id, Gene) %>%
    distinct()
  ppi_edges <- ppi_interactions %>%
    left_join(
      id_to_gene,
      by = c("from" = "STRING_id")
    ) %>%
    rename(from_gene = Gene) %>%
    left_join(
      id_to_gene,
      by = c("to" = "STRING_id")
    ) %>%
    rename(to_gene = Gene) %>%
    filter(
      !is.na(from_gene),
      !is.na(to_gene),
      from_gene != to_gene
    ) %>%
    distinct(from_gene, to_gene, .keep_all = TRUE)
  # making network
  ppi_graph <- graph_from_data_frame(
    ppi_edges %>%
      select(from_gene, to_gene),
    directed = FALSE,
    vertices = data.frame(name = ppi_input$Gene)
  )
  # statistics
  network_stats <- data.frame(
    Gene = V(ppi_graph)$name,
    Degree = degree(ppi_graph),
    Betweenness = betweenness(
      ppi_graph,
      normalized = TRUE
    )
  )
  # add own experimental values
  plot_data <- tidyr::expand_grid(
    Gene = V(ppi_graph)$name,
    Time = times
  ) %>%
    left_join(
      fc_data %>%
        transmute(
          Gene = SYMBOL,
          Time = as.character(Time),
          p.adjust = pval,
          log2FoldChange =as.numeric(log2FoldChange)
        ),
      by = c("Gene", "Time")
    ) %>%
    left_join(
      network_stats,
      by = "Gene"
    ) %>%
    mutate(
      Time = factor(Time, levels = times)
    )
  return(list(
    graph = ppi_graph,
    stats = network_stats,
    plot_data = plot_data,
    mapped = ppi_mapped
  ))
}

# plotting PPI network
plot_temporal_PPI <- function(
    ppi_results,
    times = c("0h", "2h", "6h", "24h", "48h"),
    fc_limit = NULL
) {
  graph <- ppi_results$graph
  set.seed(123)
  layout <- ggraph::create_layout(
    graph,
    layout = "fr"
  )
  node_positions <- as.data.frame(layout) %>%
    dplyr::select(name, x, y)
  edges <- igraph::as_data_frame(
    graph,
    what = "edges"
  ) %>%
    dplyr::left_join(
      node_positions %>%
        dplyr::rename(from = name),
      by = "from"
    ) %>%
    dplyr::left_join(
      node_positions %>%
        dplyr::rename(
          to = name,
          xend = x,
          yend = y
        ),
      by = "to"
    )
  nodes <- ppi_results$plot_data
  if (!"log2FoldChange" %in% names(nodes) &&
      "Log2Change" %in% names(nodes)) {
    nodes <- nodes %>%
      dplyr::rename(
        log2FoldChange = Log2Change
      )
  }
  if (!"log2FoldChange" %in% names(nodes)) {
    stop(
      "Column is missung in ppi_results$plot_data",
      "'log2FoldChange' bzw. 'Log2Change'."
    )
  }
  nodes <- nodes %>%
    dplyr::mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      Time = factor(as.character(Time), levels = times)
    ) %>%
    dplyr::left_join(
      node_positions,
      by = c("Gene" = "name")
    )
 edges_all <- tidyr::crossing(
    edges,
    Time = factor(times, levels = times)
  )
  if (is.null(fc_limit)) {

    fc_limit <- suppressWarnings(
      max(
        abs(nodes$log2FoldChange),
        na.rm = TRUE
      )
    )
  }
  if (!is.finite(fc_limit) || fc_limit == 0) {
    fc_limit <- 1
  }
 ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges_all,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      colour = "grey70",
      linewidth = 0.3,
      alpha = 0.6
    ) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(
        x = x,
        y = y,
        size = Degree,
        fill = log2FoldChange
      ),
      shape = 21,
      colour = "black",
      stroke = 0.3
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      limits = c(-fc_limit, fc_limit),
      oob = scales::squish,
      na.value = "grey85",
      name = expression(log[2]("Fold Change"))
    ) +
    ggplot2::scale_size_continuous(
      range = c(1.5, 5),
      name = "Degree"
    ) +
    ggplot2::facet_wrap(
      ~ Time,
      ncol = 3
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9)
    )
}
#-----------------------------------------------
setwd(dirs$code)
