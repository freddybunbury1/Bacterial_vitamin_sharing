input_path <- "data/intermediate/growth_dynamics_b12_methionine/od600_with_metadata_b12_methionine_2.csv"
output_dir <- "data/processed/growth_dynamics_b12_methionine"

per_well_output_path <- file.path(output_dir, "growth_curve_metrics_per_well_b12_methionine_2.csv")
blank_threshold_output_path <- file.path(output_dir, "blank_growth_thresholds_b12_methionine_2.csv")
strain_summary_output_path <- file.path(output_dir, "growth_metrics_by_strain_treatment_b12_methionine_2.csv")
treatment_contrast_output_path <- file.path(output_dir, "treatment_contrasts_by_strain_b12_methionine_2.csv")
auc_b0_vs_b10_source_path <- file.path(output_dir, "07A_mean_auc_B0_M0_vs_B10_M0_source.csv")
blank_model_coefficients_output_path <- file.path(output_dir, "blank_mixed_model_coefficients_b12_methionine_2.csv")
blank_model_decision_output_path <- file.path(output_dir, "blank_model_decision_b12_methionine_2.csv")

minimum_od_threshold <- 0.005
smoothing_span <- 0.2
blank_threshold_quantile <- 0.99

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required.")
}

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

trapz_auc <- function(x, y) {
  keep <- !(is.na(x) | is.na(y))
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 2) {
    return(NA_real_)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

compute_curve_metrics <- function(df, minimum_od_threshold, smoothing_span) {
  df <- df[order(df$growth_hours), , drop = FALSE]
  df <- df[!is.na(df$growth_hours) & !is.na(df$adjusted_OD600nm_model), , drop = FALSE]

  if (nrow(df) < 3) {
    return(data.frame(
      max_od_smoothed = NA_real_,
      time_of_max_od_smoothed = NA_real_,
      time_to_halfmax = NA_real_,
      auc_smoothed = NA_real_,
      mu_max = NA_real_,
      initial_od_smoothed = NA_real_,
      final_od_smoothed = NA_real_,
      n_timepoints = nrow(df)
    ))
  }

  adjusted_od_floor <- pmax(df$adjusted_OD600nm_model, minimum_od_threshold, na.rm = FALSE)
  log_od <- log(adjusted_od_floor)

  smoothed_log_od <- tryCatch(
    {
      if (nrow(df) > 5) {
        stats::predict(stats::loess(log_od ~ growth_hours, data = df, span = smoothing_span))
      } else {
        log_od
      }
    },
    error = function(e) log_od
  )

  smoothed_od <- exp(smoothed_log_od)
  max_idx <- which.max(smoothed_od)
  max_od_smoothed <- smoothed_od[max_idx]
  time_of_max_od_smoothed <- df$growth_hours[max_idx]
  halfmax_threshold <- 0.5 * max_od_smoothed
  halfmax_idx <- which(smoothed_od >= halfmax_threshold)[1]
  time_to_halfmax <- if (length(halfmax_idx) == 0) NA_real_ else df$growth_hours[halfmax_idx]
  auc_smoothed <- trapz_auc(df$growth_hours, smoothed_od)

  interval_mid_od <- (head(smoothed_od, -1) + tail(smoothed_od, -1)) / 2
  interval_in_mu_window <- interval_mid_od >= (0.10 * max_od_smoothed) &
    interval_mid_od <= (0.50 * max_od_smoothed)

  log_derivative <- diff(log(smoothed_od)) / diff(df$growth_hours)
  finite_derivative <- log_derivative[is.finite(log_derivative) & interval_in_mu_window]
  mu_max <- if (length(finite_derivative) == 0) NA_real_ else max(finite_derivative, na.rm = TRUE)

  data.frame(
    max_od_smoothed = max_od_smoothed,
    time_of_max_od_smoothed = time_of_max_od_smoothed,
    time_to_halfmax = time_to_halfmax,
    auc_smoothed = auc_smoothed,
    mu_max = mu_max,
    initial_od_smoothed = smoothed_od[1],
    final_od_smoothed = smoothed_od[length(smoothed_od)],
    n_timepoints = nrow(df)
  )
}

od_data <- utils::read.csv(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

od_data$date_time <- as.POSIXct(od_data$date_time, tz = "UTC")
od_data$condition_name <- paste0("B", od_data$b12_nM, "_M", od_data$methionine_uM)
od_data$curve_id <- paste(od_data$culture_plate, od_data$culture_well, sep = "__")
od_data$layout_id <- paste(od_data$array_plate, od_data$culture_well, sep = "__")

blank_data <- dplyr::filter(od_data, strain_id == "blank")

blank_model <- lme4::lmer(
  OD600nm ~ growth_hours + array_plate + factor(solution_number) + (1 | curve_id),
  data = blank_data,
  REML = TRUE
)

blank_model_coef <- as.data.frame(summary(blank_model)$coefficients)
blank_model_coef$term <- rownames(blank_model_coef)
rownames(blank_model_coef) <- NULL

blank_model_ci <- as.data.frame(stats::confint(blank_model, parm = "beta_", method = "Wald"))
blank_model_ci$term <- rownames(blank_model_ci)
rownames(blank_model_ci) <- NULL
names(blank_model_ci)[1:2] <- c("conf.low", "conf.high")

blank_model_results <- merge(blank_model_coef, blank_model_ci, by = "term", all.x = TRUE)
blank_model_results <- blank_model_results[, c(
  "term",
  "Estimate",
  "Std. Error",
  "t value",
  "conf.low",
  "conf.high"
)]
blank_model_results$ci_excludes_zero <- with(
  blank_model_results,
  !is.na(conf.low) & !is.na(conf.high) & ((conf.low > 0 & conf.high > 0) | (conf.low < 0 & conf.high < 0))
)

significant_terms <- blank_model_results$term[blank_model_results$ci_excludes_zero]
use_growth_hours <- "growth_hours" %in% significant_terms
use_array_plate <- any(grepl("^array_plate", significant_terms))
use_solution_number <- any(grepl("^factor\\(solution_number\\)", significant_terms))

blank_model_decision <- data.frame(
  full_model = "OD600nm ~ growth_hours + array_plate + factor(solution_number) + (1 | curve_id)",
  selected_model = "OD600nm ~ growth_hours + array_plate + (1 | curve_id)",
  growth_hours_selected = use_growth_hours,
  array_plate_selected = use_array_plate,
  solution_number_selected = use_solution_number,
  decision_reason = "Selected reduced model because growth_hours and array_plate showed confidence intervals excluding zero in the full model, while solution_number did not.",
  stringsAsFactors = FALSE
)

blank_model_selected <- lme4::lmer(
  OD600nm ~ growth_hours + array_plate + (1 | curve_id),
  data = blank_data,
  REML = TRUE
)

od_data$predicted_blank_od600nm <- stats::predict(
  blank_model_selected,
  newdata = od_data,
  re.form = NA,
  allow.new.levels = TRUE
)

od_data$adjusted_OD600nm_model <- od_data$OD600nm - od_data$predicted_blank_od600nm

split_indices <- split(seq_len(nrow(od_data)), od_data$curve_id, drop = TRUE)

curve_metrics_list <- lapply(split_indices, function(idx) {
  df_group <- od_data[idx, , drop = FALSE]

  group_key <- df_group[1, c(
    "curve_id",
    "layout_id",
    "culture_plate",
    "culture_plate_number",
    "culture_well",
    "solution_number",
    "array_plate",
    "technical_replicate",
    "divalent_salts_gl",
    "nacl_gl",
    "soytone_gl",
    "methionine_uM",
    "b12_nM",
    "condition_name",
    "strain_id",
    "broad_environment",
    "environmental_sample_id",
    "sample_location",
    "date_collected",
    "latitude_collected",
    "longitude_collected",
    "isolation_description",
    "growth_medium"
  ), drop = FALSE]

  metrics <- compute_curve_metrics(
    df_group,
    minimum_od_threshold = minimum_od_threshold,
    smoothing_span = smoothing_span
  )

  cbind(group_key, metrics)
})

curve_metrics <- do.call(rbind, curve_metrics_list)

blank_metrics <- dplyr::filter(curve_metrics, strain_id == "blank")

blank_thresholds <- data.frame(
  blank_threshold_quantile = blank_threshold_quantile,
  max_od_smoothed_threshold = stats::quantile(
    blank_metrics$max_od_smoothed,
    probs = blank_threshold_quantile,
    na.rm = TRUE
  ),
  auc_smoothed_threshold = stats::quantile(
    blank_metrics$auc_smoothed,
    probs = blank_threshold_quantile,
    na.rm = TRUE
  )
)

curve_metrics$growth_detected <- with(
  curve_metrics,
  !is.na(max_od_smoothed) &
    !is.na(auc_smoothed) &
    max_od_smoothed > blank_thresholds$max_od_smoothed_threshold[1] &
    auc_smoothed > blank_thresholds$auc_smoothed_threshold[1]
)

curve_metrics$growth_detected[curve_metrics$strain_id == "blank"] <- FALSE
curve_metrics$mu_max[!curve_metrics$growth_detected] <- NA_real_
curve_metrics$time_to_halfmax[!curve_metrics$growth_detected] <- NA_real_

strain_treatment_summary <- curve_metrics |>
  dplyr::filter(strain_id != "blank") |>
  dplyr::group_by(
    layout_id,
    array_plate,
    culture_well,
    strain_id,
    broad_environment,
    environmental_sample_id,
    sample_location,
    condition_name,
    b12_nM,
    methionine_uM
  ) |>
  dplyr::summarise(
    n_wells = dplyr::n(),
    n_growth_wells = sum(growth_detected, na.rm = TRUE),
    frac_growth_wells = mean(growth_detected, na.rm = TRUE),
    median_mu_max = stats::median(mu_max[growth_detected], na.rm = TRUE),
    median_auc_smoothed = stats::median(auc_smoothed, na.rm = TRUE),
    median_max_od_smoothed = stats::median(max_od_smoothed, na.rm = TRUE),
    median_time_of_max_od_smoothed = stats::median(time_of_max_od_smoothed[growth_detected], na.rm = TRUE),
    median_time_to_halfmax = stats::median(time_to_halfmax[growth_detected], na.rm = TRUE),
    .groups = "drop"
  )

summary_metric_columns <- c(
  "n_wells",
  "n_growth_wells",
  "frac_growth_wells",
  "median_mu_max",
  "median_auc_smoothed",
  "median_max_od_smoothed",
  "median_time_of_max_od_smoothed",
  "median_time_to_halfmax"
)

for (col in summary_metric_columns) {
  strain_treatment_summary[[col]][is.nan(strain_treatment_summary[[col]])] <- NA_real_
}

control_summary <- strain_treatment_summary |>
  dplyr::filter(condition_name == "B0_M0") |>
  dplyr::group_by(
    layout_id,
    array_plate
  ) |>
  dplyr::summarise(
    control_frac_growth_wells = stats::median(frac_growth_wells, na.rm = TRUE),
    control_median_mu_max = stats::median(median_mu_max, na.rm = TRUE),
    control_median_auc_smoothed = stats::median(median_auc_smoothed, na.rm = TRUE),
    control_median_max_od_smoothed = stats::median(median_max_od_smoothed, na.rm = TRUE),
    control_median_time_of_max_od_smoothed = stats::median(median_time_of_max_od_smoothed, na.rm = TRUE),
    control_median_time_to_halfmax = stats::median(median_time_to_halfmax, na.rm = TRUE),
    .groups = "drop"
  )

treatment_contrasts <- strain_treatment_summary |>
  dplyr::filter(condition_name != "B0_M0") |>
  dplyr::left_join(control_summary, by = c("layout_id", "array_plate")) |>
  dplyr::mutate(
    delta_frac_growth_wells = frac_growth_wells - control_frac_growth_wells,
    delta_mu_max = median_mu_max - control_median_mu_max,
    delta_auc_smoothed = median_auc_smoothed - control_median_auc_smoothed,
    delta_max_od_smoothed = median_max_od_smoothed - control_median_max_od_smoothed,
    delta_time_of_max_od_smoothed = median_time_of_max_od_smoothed - control_median_time_of_max_od_smoothed,
    delta_time_to_halfmax = median_time_to_halfmax - control_median_time_to_halfmax
  )

# Figure S3 uses only this B12/no-B12 AUC comparison. Keep the compact source
# export here instead of retaining the exploratory growth-figure script.
auc_b0_vs_b10 <- curve_metrics |>
  dplyr::filter(
    .data$strain_id != "blank",
    .data$condition_name %in% c("B0_M0", "B10_M0"),
    is.finite(.data$auc_smoothed)
  ) |>
  dplyr::group_by(.data$strain_id, .data$condition_name) |>
  dplyr::summarise(
    mean_auc = mean(.data$auc_smoothed),
    sd_auc = stats::sd(.data$auc_smoothed),
    n_replicates = dplyr::n(),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = "condition_name",
    values_from = c("mean_auc", "sd_auc", "n_replicates"),
    names_glue = "{.value}_{condition_name}"
  ) |>
  dplyr::filter(
    is.finite(.data$mean_auc_B0_M0),
    is.finite(.data$mean_auc_B10_M0)
  ) |>
  dplyr::mutate(
    delta_auc_B10_M0_vs_B0_M0 = .data$mean_auc_B10_M0 - .data$mean_auc_B0_M0,
    absolute_delta_auc = abs(.data$delta_auc_B10_M0_vs_B0_M0),
    auc_B0_M0_lower_raw = .data$mean_auc_B0_M0 - .data$sd_auc_B0_M0,
    auc_B0_M0_lower = pmax(.data$auc_B0_M0_lower_raw, 0.1),
    auc_B0_M0_upper = .data$mean_auc_B0_M0 + .data$sd_auc_B0_M0,
    auc_B10_M0_lower_raw = .data$mean_auc_B10_M0 - .data$sd_auc_B10_M0,
    auc_B10_M0_lower = pmax(.data$auc_B10_M0_lower_raw, 0.1),
    auc_B10_M0_upper = .data$mean_auc_B10_M0 + .data$sd_auc_B10_M0,
    label_top1pct = .data$absolute_delta_auc >= stats::quantile(
      .data$absolute_delta_auc, probs = 0.99, na.rm = TRUE
    )
  )

utils::write.csv(curve_metrics, per_well_output_path, row.names = FALSE, na = "")
utils::write.csv(blank_thresholds, blank_threshold_output_path, row.names = FALSE, na = "")
utils::write.csv(strain_treatment_summary, strain_summary_output_path, row.names = FALSE, na = "")
utils::write.csv(treatment_contrasts, treatment_contrast_output_path, row.names = FALSE, na = "")
utils::write.csv(auc_b0_vs_b10, auc_b0_vs_b10_source_path, row.names = FALSE, na = "")
utils::write.csv(blank_model_results, blank_model_coefficients_output_path, row.names = FALSE, na = "")
utils::write.csv(blank_model_decision, blank_model_decision_output_path, row.names = FALSE, na = "")

message("Wrote per-well metrics to: ", per_well_output_path)
message("Wrote blank thresholds to: ", blank_threshold_output_path)
message("Wrote strain-by-treatment summary to: ", strain_summary_output_path)
message("Wrote treatment contrasts to: ", treatment_contrast_output_path)
message("Wrote Figure S3 source data to: ", auc_b0_vs_b10_source_path)
message("Wrote blank mixed-model coefficients to: ", blank_model_coefficients_output_path)
message("Wrote blank model decision to: ", blank_model_decision_output_path)
message("Blank max OD threshold (99th percentile): ", round(blank_thresholds$max_od_smoothed_threshold[1], 4))
message("Blank AUC threshold (99th percentile): ", round(blank_thresholds$auc_smoothed_threshold[1], 4))
