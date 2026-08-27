#!/usr/bin/env Rscript

# Public end-to-end release/uptake/death analysis.
# Stage 1 estimates B12 and joins OD/viability; stage 2 fits the kinetic models.

section_name <- "b12_release_uptake_death"

intermediate_dir <- file.path("data/intermediate", section_name)
processed_dir <- file.path("data/processed", section_name)
estimator_processed_dir <- file.path("data/processed", "b12_concentration_estimator")
estimator_path <- file.path("R/utils", "b12_concentration_estimator.R")

culture_samples_path <- file.path(intermediate_dir, "culture_samples_b12_release_uptake_death.csv")
b12_assay_combined_path <- file.path(intermediate_dir, "b12_assay_combined_b12_release_uptake_death.csv")
od_combined_path <- file.path(intermediate_dir, "od600_combined_b12_release_uptake_death.csv")
viability_summary_input_path <- file.path(intermediate_dir, "viability_summary_b12_release_uptake_death.csv")

analysis_parameters_path <- file.path(processed_dir, "analysis_parameters_b12_release_uptake_death.csv")
b12_concentration_estimates_path <- file.path(processed_dir, "b12_concentration_estimates_b12_release_uptake_death.csv")
b12_concentration_estimates_uncorrected_path <- file.path(processed_dir, "b12_concentration_estimates_uncorrected_b12_release_uptake_death.csv")
b12_estimator_parameters_path <- file.path(estimator_processed_dir, "b12_estimator_parameters.csv")
b12_curve_library_path <- file.path(processed_dir, "b12_standard_curve_library_b12_release_uptake_death.csv")
b12_standards_used_path <- file.path(processed_dir, "b12_standards_used_b12_release_uptake_death.csv")
b12_corrected_assay_path <- file.path(processed_dir, "b12_assay_edge_corrected_b12_release_uptake_death.csv")
b12_adjacent_pair_summary_path <- file.path(processed_dir, "b12_adjacent_row_pair_summary_b12_release_uptake_death.csv")
b12_adjacent_pair_block_summary_path <- file.path(processed_dir, "b12_adjacent_row_pair_block_summary_b12_release_uptake_death.csv")
b12_adjacent_pair_model_tests_path <- file.path(processed_dir, "b12_adjacent_row_pair_quadratic_tests_b12_release_uptake_death.csv")
b12_edge_correction_pairs_path <- file.path(processed_dir, "b12_edge_correction_pairs_b12_release_uptake_death.csv")
b12_edge_correction_summary_path <- file.path(processed_dir, "b12_edge_correction_summary_b12_release_uptake_death.csv")
b12_results_path <- file.path(processed_dir, "b12_results_b12_release_uptake_death.csv")
b12_results_uncorrected_path <- file.path(processed_dir, "b12_results_uncorrected_b12_release_uptake_death.csv")
b12_od_viability_results_path <- file.path(processed_dir, "b12_od_viability_results_b12_release_uptake_death.csv")
b12_od_viability_results_filtered_path <- file.path(processed_dir, "b12_od_viability_results_filtered_b12_release_uptake_death.csv")
b12_summary_path <- file.path(processed_dir, "b12_summary_b12_release_uptake_death.csv")
b12_summary_uncorrected_path <- file.path(processed_dir, "b12_summary_uncorrected_b12_release_uptake_death.csv")
b12_sample_assay_edge_flags_path <- file.path(processed_dir, "b12_sample_assay_edge_flags_b12_release_uptake_death.csv")

standard_sample_types <- c("filtrate")
standard_cols_to_remove <- c(1)
add_offset_pm <- 0.1
edge_rows <- c("A", "P")
edge_correction <- "quadratic_edge_pairs"
upper_thresh_logb12_estimation_se <- 0.5
blank_od_filter_sd_multiplier <- 2
blank_dead_stain_filter_sd_multiplier <- 2
blank_live_stain_floor_sd_multiplier <- 1

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required.")
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  stop("Package 'stringr' is required.")
}

if (!requireNamespace("drc", quietly = TRUE)) {
  stop("Package 'drc' is required.")
}

source(estimator_path)

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(estimator_processed_dir, recursive = TRUE, showWarnings = FALSE)

analysis_parameters <- data.frame(
  parameter = c(
    "standard_sample_types",
    "standard_cols_to_remove",
    "add_offset_pm",
    "edge_rows",
    "edge_correction",
    "upper_thresh_logb12_estimation_se"
  ),
  value = c(
    paste(standard_sample_types, collapse = ";"),
    paste(standard_cols_to_remove, collapse = ";"),
    add_offset_pm,
    paste(edge_rows, collapse = ";"),
    edge_correction,
    upper_thresh_logb12_estimation_se
  ),
  stringsAsFactors = FALSE
)

culture_samples <- utils::read.csv(
  culture_samples_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

b12_assay_combined <- utils::read.csv(
  b12_assay_combined_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

od_combined <- utils::read.csv(
  od_combined_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

viability_summary_input <- utils::read.csv(
  viability_summary_input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

estimator_result <- b12_concentration_estimator(
  b12_assay_combined,
  standard_sample_types = standard_sample_types,
  standard_cols_to_remove = standard_cols_to_remove,
  add_offset_pm = add_offset_pm,
  edge_correction = edge_correction
)

estimator_result_uncorrected <- b12_concentration_estimator(
  b12_assay_combined,
  standard_sample_types = standard_sample_types,
  standard_cols_to_remove = standard_cols_to_remove,
  add_offset_pm = add_offset_pm,
  edge_correction = "none"
)

b12_concentration_estimates <- estimator_result$estimates
b12_concentration_estimates_uncorrected <- estimator_result_uncorrected$estimates

b12_concentration_estimates <- b12_concentration_estimates |>
  dplyr::mutate(
    passed_se_filter = is.finite(se_log_fm) &
      se_log_fm <= upper_thresh_logb12_estimation_se
  )

b12_concentration_estimates_uncorrected <- b12_concentration_estimates_uncorrected |>
  dplyr::mutate(
    passed_se_filter = is.finite(se_log_fm) &
      se_log_fm <= upper_thresh_logb12_estimation_se
  )

b12_sample_assay_edge_flags <- b12_assay_combined |>
  dplyr::filter(assay_type == "sample") |>
  dplyr::mutate(
    assay_row = stringr::str_extract(assay_well, "^[A-Z]+"),
    assay_col = as.integer(stringr::str_extract(assay_well, "\\d+")),
    is_edge_row = assay_row %in% edge_rows
  ) |>
  dplyr::group_by(sample_id, assay_batch, experiment_name, sample_time_minute) |>
  dplyr::summarise(
    assay_rows = paste(sort(unique(assay_row)), collapse = ";"),
    assay_wells = paste(sort(unique(assay_well)), collapse = ";"),
    n_assay_observations = dplyr::n(),
    n_edge_row_observations = sum(is_edge_row, na.rm = TRUE),
    any_edge_row = any(is_edge_row, na.rm = TRUE),
    all_edge_row = all(is_edge_row, na.rm = TRUE),
    .groups = "drop"
  )

b12_results <- culture_samples |>
  dplyr::mutate(sample_id = as.integer(sample_id)) |>
  dplyr::left_join(
    b12_concentration_estimates |>
      dplyr::mutate(sample_id = as.integer(sample_id)) |>
      dplyr::mutate(
        b12_pm = dplyr::if_else(passed_se_filter, b12_pm, NA_real_),
        b12_pm_ci_low = dplyr::if_else(passed_se_filter, b12_pm_ci_low, NA_real_),
        b12_pm_ci_high = dplyr::if_else(passed_se_filter, b12_pm_ci_high, NA_real_)
      ) |>
      dplyr::select(
        assay_batch,
        experiment_name,
        sample_id,
        sample_type,
        b12_pm,
        se_log_fm,
        b12_pm_ci_low,
        b12_pm_ci_high,
        n_obs,
        eff_info,
        passed_se_filter
      ),
    by = c("experiment_name", "sample_id")
  ) |>
  dplyr::left_join(
    b12_sample_assay_edge_flags |>
      dplyr::select(
        sample_id,
        any_edge_row,
        all_edge_row,
        n_edge_row_observations,
        n_assay_observations,
        assay_rows,
        assay_wells
      ),
    by = "sample_id"
  )

b12_results_uncorrected <- culture_samples |>
  dplyr::mutate(sample_id = as.integer(sample_id)) |>
  dplyr::left_join(
    b12_concentration_estimates_uncorrected |>
      dplyr::mutate(sample_id = as.integer(sample_id)) |>
      dplyr::mutate(
        b12_pm = dplyr::if_else(passed_se_filter, b12_pm, NA_real_),
        b12_pm_ci_low = dplyr::if_else(passed_se_filter, b12_pm_ci_low, NA_real_),
        b12_pm_ci_high = dplyr::if_else(passed_se_filter, b12_pm_ci_high, NA_real_)
      ) |>
      dplyr::select(
        assay_batch,
        experiment_name,
        sample_id,
        sample_type,
        b12_pm,
        se_log_fm,
        b12_pm_ci_low,
        b12_pm_ci_high,
        n_obs,
        eff_info,
        passed_se_filter
      ),
    by = c("experiment_name", "sample_id")
  ) |>
  dplyr::left_join(
    b12_sample_assay_edge_flags |>
      dplyr::select(
        sample_id,
        any_edge_row,
        all_edge_row,
        n_edge_row_observations,
        n_assay_observations,
        assay_rows,
        assay_wells
      ),
    by = "sample_id"
  )

b12_summary <- b12_results |>
  dplyr::group_by(
    experiment_type,
    strain_id,
    added_b12_pM,
    added_glucose_mM_C,
    sample_time_minute,
    killed_proportion
  ) |>
  dplyr::summarise(
    mean_b12 = mean(b12_pm, na.rm = TRUE),
    sd_b12 = stats::sd(b12_pm, na.rm = TRUE),
    n = dplyr::n(),
    n_estimated = sum(!is.na(b12_pm)),
    n_edge_row = sum(any_edge_row, na.rm = TRUE),
    .groups = "drop"
  )

b12_summary_uncorrected <- b12_results_uncorrected |>
  dplyr::group_by(
    experiment_type,
    strain_id,
    added_b12_pM,
    added_glucose_mM_C,
    sample_time_minute,
    killed_proportion
  ) |>
  dplyr::summarise(
    mean_b12 = mean(b12_pm, na.rm = TRUE),
    sd_b12 = stats::sd(b12_pm, na.rm = TRUE),
    n = dplyr::n(),
    n_estimated = sum(!is.na(b12_pm)),
    n_edge_row = sum(any_edge_row, na.rm = TRUE),
    .groups = "drop"
  )

join_keys <- c(
  "strain_id",
  "biological_replicate",
  "experiment_name",
  "experiment_type",
  "sample_time_minute",
  "sample_time_hour",
  "killed_proportion",
  "added_b12_pM",
  "added_glucose_mM_C"
)

b12_join_data <- b12_results |>
  dplyr::mutate(sample_time_hour = sample_time_minute / 60) |>
  dplyr::select(
    dplyr::all_of(join_keys),
    culture_well,
    filtrate_b12 = b12_pm,
    filtrate_b12_ci_low = b12_pm_ci_low,
    filtrate_b12_ci_high = b12_pm_ci_high,
    filtrate_b12_se_log_fm = se_log_fm,
    filtrate_b12_passed_se_filter = passed_se_filter,
    b12_n_obs = n_obs,
    b12_eff_info = eff_info,
    b12_any_edge_row = any_edge_row,
    b12_all_edge_row = all_edge_row
  )

od_join_data <- od_combined |>
  dplyr::select(
    dplyr::all_of(join_keys),
    culture_well,
    raw_od600 = od600
  )

blank_od_values <- od_join_data |>
  dplyr::filter(strain_id == "blank", is.finite(raw_od600))

blank_od_n <- nrow(blank_od_values)
blank_od_mean <- mean(blank_od_values$raw_od600, na.rm = TRUE)
blank_od_sd <- stats::sd(blank_od_values$raw_od600, na.rm = TRUE)
blank_od_filter_threshold <- blank_od_mean + blank_od_filter_sd_multiplier * blank_od_sd

viability_join_data <- viability_summary_input |>
  dplyr::mutate(
    stain_metric = dplyr::case_when(
      viability_status == "live" ~ "live_stain",
      viability_status == "dead" ~ "dead_stain",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(stain_metric)) |>
  dplyr::group_by(dplyr::across(dplyr::all_of(join_keys)), stain_metric) |>
  dplyr::summarise(
    stain_fluorescence = mean(fluorescence, na.rm = TRUE),
    stain_sd = stats::sd(fluorescence, na.rm = TRUE),
    stain_n = dplyr::n(),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = stain_metric,
    values_from = c(stain_fluorescence, stain_sd, stain_n),
    names_glue = "{stain_metric}_{.value}"
  ) |>
  dplyr::rename(
    live_stain = live_stain_stain_fluorescence,
    dead_stain = dead_stain_stain_fluorescence,
    live_stain_sd = live_stain_stain_sd,
    dead_stain_sd = dead_stain_stain_sd,
    live_stain_n = live_stain_stain_n,
    dead_stain_n = dead_stain_stain_n
  )

blank_dead_stain_values <- viability_join_data |>
  dplyr::filter(strain_id == "blank", is.finite(dead_stain))

blank_dead_stain_n <- nrow(blank_dead_stain_values)
blank_dead_stain_mean <- mean(blank_dead_stain_values$dead_stain, na.rm = TRUE)
blank_dead_stain_sd <- stats::sd(blank_dead_stain_values$dead_stain, na.rm = TRUE)
blank_dead_stain_filter_threshold <- blank_dead_stain_mean +
  blank_dead_stain_filter_sd_multiplier * blank_dead_stain_sd

blank_live_stain_values <- viability_join_data |>
  dplyr::filter(strain_id == "blank", is.finite(live_stain))

blank_live_stain_n <- nrow(blank_live_stain_values)
blank_live_stain_mean <- mean(blank_live_stain_values$live_stain, na.rm = TRUE)
blank_live_stain_sd <- stats::sd(blank_live_stain_values$live_stain, na.rm = TRUE)
blank_live_stain_floor <- blank_live_stain_mean +
  blank_live_stain_floor_sd_multiplier * blank_live_stain_sd

blank_measurement_parameters <- data.frame(
  parameter = c(
    "blank_od_filter_sd_multiplier",
    "blank_od_n",
    "blank_od_mean",
    "blank_od_sd",
    "blank_od_filter_threshold",
    "blank_od_filter_formula",
    "blank_dead_stain_filter_sd_multiplier",
    "blank_dead_stain_n",
    "blank_dead_stain_mean",
    "blank_dead_stain_sd",
    "blank_dead_stain_filter_threshold",
    "blank_dead_stain_filter_formula",
    "blank_live_stain_floor_sd_multiplier",
    "blank_live_stain_n",
    "blank_live_stain_mean",
    "blank_live_stain_sd",
    "blank_live_stain_floor",
    "stain_adjustment_formula",
    "dead_proportion_formula"
  ),
  value = c(
    blank_od_filter_sd_multiplier,
    blank_od_n,
    blank_od_mean,
    blank_od_sd,
    blank_od_filter_threshold,
    "strain_id == blank OR raw_od600 > blank_od_mean + blank_od_filter_sd_multiplier * blank_od_sd",
    blank_dead_stain_filter_sd_multiplier,
    blank_dead_stain_n,
    blank_dead_stain_mean,
    blank_dead_stain_sd,
    blank_dead_stain_filter_threshold,
    "strain_id == blank OR dead_stain > blank_dead_stain_mean + blank_dead_stain_filter_sd_multiplier * blank_dead_stain_sd",
    blank_live_stain_floor_sd_multiplier,
    blank_live_stain_n,
    blank_live_stain_mean,
    blank_live_stain_sd,
    blank_live_stain_floor,
    "adjusted_live_stain = floored_live_stain - pooled blank_live_stain_mean; adjusted_dead_stain = dead_stain - pooled blank_dead_stain_mean",
    "dead_proportion_floored_live = adjusted_live_stain / adjusted_dead_stain when adjusted_dead_stain > 0"
  ),
  stringsAsFactors = FALSE
)

b12_od_viability_results <- b12_join_data |>
  dplyr::full_join(
    od_join_data,
    by = c(join_keys, "culture_well")
  ) |>
  dplyr::full_join(
    viability_join_data,
    by = join_keys
  ) |>
  dplyr::mutate(
    od600 = raw_od600 - blank_od_mean,
    adjusted_od600 = od600,
    blank_od600_mean = blank_od_mean,
    blank_od600_filter_threshold = blank_od_filter_threshold,
    passed_blank_od_filter = dplyr::case_when(
      strain_id == "blank" ~ TRUE,
      is.na(raw_od600) ~ NA,
      TRUE ~ raw_od600 > blank_od_filter_threshold
    ),
    blank_live_stain_mean = blank_live_stain_mean,
    blank_live_stain_floor = blank_live_stain_floor,
    floored_live_stain = dplyr::if_else(
      is.na(live_stain),
      NA_real_,
      pmax(live_stain, blank_live_stain_floor, na.rm = TRUE)
    ),
    live_stain_below_floor = live_stain < floored_live_stain,
    adjusted_live_stain = floored_live_stain - blank_live_stain_mean,
    blank_dead_stain_mean = blank_dead_stain_mean,
    blank_dead_stain_filter_threshold = blank_dead_stain_filter_threshold,
    adjusted_dead_stain = dead_stain - blank_dead_stain_mean,
    passed_blank_dead_stain_filter = dplyr::case_when(
      strain_id == "blank" ~ TRUE,
      is.na(dead_stain) ~ NA,
      TRUE ~ dead_stain > blank_dead_stain_filter_threshold
    ),
    raw_dead_proportion = dplyr::if_else(dead_stain > 0, live_stain / dead_stain, NA_real_),
    dead_proportion_floored_live = dplyr::if_else(
      adjusted_dead_stain > 0,
      adjusted_live_stain / adjusted_dead_stain,
      NA_real_
    ),
    dead_proportion = dead_proportion_floored_live
  ) |>
  dplyr::arrange(
    experiment_name,
    experiment_type,
    sample_time_minute,
    strain_id,
    biological_replicate,
    killed_proportion,
    added_b12_pM,
    added_glucose_mM_C
  )

b12_od_viability_results_filtered <- b12_od_viability_results |>
  dplyr::filter(
    (is.na(passed_blank_dead_stain_filter) | passed_blank_dead_stain_filter)
  )

analysis_parameters <- dplyr::bind_rows(
  analysis_parameters,
  blank_measurement_parameters
)

utils::write.csv(analysis_parameters, analysis_parameters_path, row.names = FALSE, na = "")
utils::write.csv(b12_concentration_estimates, b12_concentration_estimates_path, row.names = FALSE, na = "")
utils::write.csv(b12_concentration_estimates_uncorrected, b12_concentration_estimates_uncorrected_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$estimator_parameters, b12_estimator_parameters_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$curve_library, b12_curve_library_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$standards_used, b12_standards_used_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$corrected_assay, b12_corrected_assay_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$adjacent_pair_diagnostics$pair_summary, b12_adjacent_pair_summary_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$adjacent_pair_diagnostics$block_summary, b12_adjacent_pair_block_summary_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$adjacent_pair_diagnostics$quadratic_model_tests, b12_adjacent_pair_model_tests_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$edge_correction_pairs, b12_edge_correction_pairs_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$edge_correction_summary, b12_edge_correction_summary_path, row.names = FALSE, na = "")
utils::write.csv(b12_results, b12_results_path, row.names = FALSE, na = "")
utils::write.csv(b12_results_uncorrected, b12_results_uncorrected_path, row.names = FALSE, na = "")
utils::write.csv(b12_od_viability_results, b12_od_viability_results_path, row.names = FALSE, na = "")
utils::write.csv(b12_od_viability_results_filtered, b12_od_viability_results_filtered_path, row.names = FALSE, na = "")
utils::write.csv(b12_summary, b12_summary_path, row.names = FALSE, na = "")
utils::write.csv(b12_summary_uncorrected, b12_summary_uncorrected_path, row.names = FALSE, na = "")
utils::write.csv(b12_sample_assay_edge_flags, b12_sample_assay_edge_flags_path, row.names = FALSE, na = "")

message("Wrote analysis parameters to: ", analysis_parameters_path)
message("Wrote B12 concentration estimates to: ", b12_concentration_estimates_path)
message("Wrote B12 results to: ", b12_results_path)
message("Wrote joined B12/OD/viability results to: ", b12_od_viability_results_path)
message("Wrote filtered B12/OD/viability results to: ", b12_od_viability_results_filtered_path)
message("Wrote B12 summary to: ", b12_summary_path)
message("Estimated sample groups: ", nrow(b12_concentration_estimates))
message("Culture sample rows: ", nrow(culture_samples))
message("Culture sample rows with B12 estimates: ", sum(!is.na(b12_results$b12_pm)))

# -----------------------------------------------------------------------------
# Kinetic modeling stage
# -----------------------------------------------------------------------------

section_name <- "b12_release_uptake_death"

processed_dir <- file.path("data/processed", section_name)

# Inputs
b12_results_path <- file.path(processed_dir, "b12_results_b12_release_uptake_death.csv")

# Canonical two-stage kinetic model used by the manuscript.
pipeline_path <- "R/utils/b12_two_stage_kinetic_pipeline.R"

# Modeling outputs
modeling_dir <- file.path(processed_dir, "modeling")

modeling_parameters_path <- file.path(modeling_dir, "analysis_parameters_b12_release_uptake_death_modeling.csv")
b12_model_summary_path <- file.path(modeling_dir, "b12_model_summary_b12_release_uptake_death.csv")
b12_delta_estimates_path <- file.path(modeling_dir, "b12_delta_estimates_b12_release_uptake_death.csv")
b12_b_estimates_path <- file.path(modeling_dir, "b12_b_estimates_b12_release_uptake_death.csv")
b12_fits_uptake_path <- file.path(modeling_dir, "b12_fits_uptake_b12_release_uptake_death.csv")
b12_fits_release_path <- file.path(modeling_dir, "b12_fits_release_b12_release_uptake_death.csv")
b12_mix_pred_uptake_path <- file.path(modeling_dir, "b12_mix_pred_uptake_b12_release_uptake_death.csv")
b12_mix_pred_release_path <- file.path(modeling_dir, "b12_mix_pred_release_b12_release_uptake_death.csv")
b12_with_prop_path <- file.path(modeling_dir, "b12_with_prop_b12_release_uptake_death.csv")

required_packages <- c("dplyr", "tidyr", "broom", "minpack.lm", "ggplot2")
missing <- required_packages[!vapply(
  required_packages,
  function(pkg) requireNamespace(pkg, quietly = TRUE),
  logical(1)
)]
if (length(missing) > 0) {
  stop("Required package(s) are not installed: ", paste(missing, collapse = ", "))
}

for (p in required_packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

dir.create(modeling_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(processed_dir)) {
  stop("Processed directory not found: ", processed_dir)
}

if (!file.exists(b12_results_path)) {
  stop("Input file not found: ", b12_results_path)
}

if (!file.exists(pipeline_path)) {
  stop("Modeling pipeline not found: ", pipeline_path)
}

# Load the single public modeling implementation.
source(pipeline_path)

# input data
b12_results <- utils::read.csv(
  b12_results_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# Fixed manuscript model settings.
group_by_glucose <- TRUE
group_by_replicate <- TRUE
fit_scale <- "log"

model_outputs <- build_b12_model_outputs(
  b12_results,
  group_by_glucose = group_by_glucose,
  group_by_replicate = group_by_replicate,
  fit_scale = fit_scale
)

write.csv(
  data.frame(
    parameter = c(
      "pipeline_file",
      "group_by_glucose",
      "group_by_replicate",
      "fit_scale",
      "pipeline_version"
    ),
    value = c(
      basename(pipeline_path),
      group_by_glucose,
      group_by_replicate,
      fit_scale,
      "ported_from_aps_march_modeling_section_2026_02_06"
    ),
    stringsAsFactors = FALSE
  ),
  modeling_parameters_path,
  row.names = FALSE
)

write.csv(
  data.frame(
    n_rows = nrow(b12_results),
    blank_offset = model_outputs$blank_offset,
    pipeline = basename(pipeline_path),
    group_by_glucose = group_by_glucose,
    group_by_replicate = group_by_replicate,
    fit_scale = fit_scale,
    stringsAsFactors = FALSE
  ),
  b12_model_summary_path,
  row.names = FALSE
)

write.csv(model_outputs$delta_estimates, b12_delta_estimates_path, row.names = FALSE)
write.csv(model_outputs$b_estimates, b12_b_estimates_path, row.names = FALSE)
write.csv(model_outputs$fits_uptake, b12_fits_uptake_path, row.names = FALSE)
write.csv(model_outputs$fits_release, b12_fits_release_path, row.names = FALSE)
write.csv(model_outputs$mix_pred_uptake, b12_mix_pred_uptake_path, row.names = FALSE)
write.csv(model_outputs$mix_pred_release, b12_mix_pred_release_path, row.names = FALSE)

# denominator and filtrate proportion for live/dead mix
full_kill_vals <- b12_results |>
  dplyr::filter(experiment_type == "live_dead_mix") |>
  dplyr::group_by(strain_id, biological_replicate, experiment_name, added_glucose_mM_C) |>
  dplyr::summarise(
    b12_full_kill = stats::median(
      dplyr::if_else(killed_proportion == 1, b12_pm, NA_real_),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

b12_with_prop <- b12_results |>
  dplyr::left_join(
    full_kill_vals,
    by = c("strain_id", "biological_replicate", "experiment_name", "added_glucose_mM_C")
  ) |>
  dplyr::mutate(
    filtrate_proportion = dplyr::if_else(
      experiment_type == "live_dead_mix" & !is.na(b12_full_kill) & b12_full_kill > 0,
      b12_pm / b12_full_kill,
      NA_real_
    )
  )

write.csv(b12_with_prop, b12_with_prop_path, row.names = FALSE)

message("Wrote canonical kinetic-model outputs to: ", modeling_dir)

# -----------------------------------------------------------------------------
# Compact compiled-figure source tables
# -----------------------------------------------------------------------------

figure_source_dir <- file.path(
  "data/processed/figure_source_data", "b12_release_uptake_death"
)
dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)

figure_b12_summary <- utils::read.csv(
  file.path(processed_dir, "b12_summary_b12_release_uptake_death.csv"),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
)
figure_joined_results <- utils::read.csv(
  file.path(processed_dir, "b12_od_viability_results_b12_release_uptake_death.csv"),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
)
figure_taxonomy <- utils::read.csv(
  file.path("data/intermediate/gtdbtk", "gtdbtk_bac120_taxonomy_with_strain_id.csv"),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
) |>
  dplyr::select("strain_id", "gtdb_order") |>
  dplyr::distinct(.data$strain_id, .keep_all = TRUE)
figure_trait_summary <- utils::read.csv(
  file.path("data/processed/b12_trait_survey", "b12_trait_survey_strain_summary.csv"),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
)
figure_traits <- figure_trait_summary |>
  dplyr::select("strain_id", "b12_trait_group2") |>
  dplyr::distinct(.data$strain_id, .keep_all = TRUE)

trait_group_levels <- c("retainer", "provider", "reclaimer")
facet_lookup <- figure_traits |>
  dplyr::mutate(
    trait_group2 = dplyr::case_when(
      .data$strain_id == "blank" ~ "blank",
      .data$b12_trait_group2 %in% trait_group_levels ~ .data$b12_trait_group2,
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::bind_rows(data.frame(
    strain_id = "blank", b12_trait_group2 = NA_character_,
    trait_group2 = "blank", stringsAsFactors = FALSE
  )) |>
  dplyr::filter(!is.na(.data$trait_group2)) |>
  dplyr::distinct(.data$strain_id, .keep_all = TRUE) |>
  dplyr::mutate(
    trait_group2 = factor(.data$trait_group2, levels = c("blank", trait_group_levels))
  )

add_trait_facets <- function(data) {
  out <- data |>
    dplyr::select(-dplyr::any_of(c(
      "b12_trait_group2", "trait_group2", "trait_group2_strain_id"
    ))) |>
    dplyr::left_join(facet_lookup, by = "strain_id") |>
    dplyr::mutate(
      trait_group2_strain_id = dplyr::if_else(
        .data$trait_group2 == "blank", "blank",
        paste(.data$trait_group2, .data$strain_id, sep = "_")
      )
    )
  facet_levels <- out |>
    dplyr::filter(!is.na(.data$trait_group2)) |>
    dplyr::distinct(.data$trait_group2, .data$strain_id, .data$trait_group2_strain_id) |>
    dplyr::arrange(.data$trait_group2, .data$strain_id) |>
    dplyr::pull(.data$trait_group2_strain_id)
  out |>
    dplyr::mutate(
      trait_group2_strain_id = factor(
        .data$trait_group2_strain_id, levels = unique(facet_levels)
      )
    ) |>
    dplyr::arrange(.data$trait_group2, .data$strain_id)
}

# Figure S9: dead-cell fluorescence against the killed-cell-density proxy.
live_dead_mix <- figure_joined_results |>
  dplyr::filter(.data$experiment_type == "live_dead_mix") |>
  dplyr::left_join(figure_taxonomy, by = "strain_id") |>
  dplyr::mutate(
    biological_replicate = factor(.data$biological_replicate),
    experiment_name = factor(.data$experiment_name),
    replicate_id = interaction(
      .data$experiment_name, .data$biological_replicate, sep = " / "
    )
  ) |>
  add_trait_facets() |>
  dplyr::arrange(
    .data$trait_group2, .data$strain_id, .data$experiment_name,
    .data$biological_replicate, .data$killed_proportion
  ) |>
  droplevels()

source_od_0min <- figure_joined_results |>
  dplyr::filter(
    .data$experiment_type == "live_dead_source",
    .data$sample_time_minute == 0,
    .data$killed_proportion == 0
  ) |>
  dplyr::mutate(biological_replicate = factor(.data$biological_replicate)) |>
  dplyr::select(
    "experiment_name", "strain_id", "biological_replicate",
    source_culture_od600_0min = "adjusted_od600"
  )

figure_s9_source <- live_dead_mix |>
  dplyr::left_join(
    source_od_0min,
    by = c("experiment_name", "strain_id", "biological_replicate")
  ) |>
  dplyr::mutate(
    killed_proportion_source_od600 =
      .data$killed_proportion * .data$source_culture_od600_0min
  )

# Figure S14E-F: paired glucose response normalized by live-stain signal.
release_endpoint_cultures <- figure_joined_results |>
  dplyr::filter(
    .data$experiment_type == "release_and_uptake",
    .data$added_b12_pM == 0,
    .data$sample_time_minute == 240,
    !is.na(.data$filtrate_b12),
    !is.na(.data$dead_stain)
  ) |>
  dplyr::left_join(figure_taxonomy, by = "strain_id") |>
  dplyr::mutate(added_glucose_mM_C = factor(.data$added_glucose_mM_C)) |>
  add_trait_facets() |>
  dplyr::filter(.data$filtrate_b12 > 0, .data$adjusted_live_stain > 0) |>
  dplyr::mutate(
    filtrate_b12_per_adjusted_live_stain =
      .data$filtrate_b12 / .data$adjusted_live_stain
  )

glucose_pair_source <- release_endpoint_cultures |>
  dplyr::filter(
    is.finite(.data$filtrate_b12_per_adjusted_live_stain),
    .data$filtrate_b12_per_adjusted_live_stain > 0,
    .data$added_glucose_mM_C %in% c(0, 10)
  ) |>
  dplyr::group_by(.data$strain_id, .data$b12_trait_group2, .data$added_glucose_mM_C) |>
  dplyr::summarise(
    mean_filtrate_b12_per_adjusted_live_stain = mean(
      .data$filtrate_b12_per_adjusted_live_stain, na.rm = TRUE
    ),
    n_cultures = dplyr::n(),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = "added_glucose_mM_C",
    values_from = c("mean_filtrate_b12_per_adjusted_live_stain", "n_cultures"),
    names_glue = "{.value}_glucose_{added_glucose_mM_C}"
  ) |>
  dplyr::filter(
    is.finite(.data$mean_filtrate_b12_per_adjusted_live_stain_glucose_0),
    .data$mean_filtrate_b12_per_adjusted_live_stain_glucose_0 > 0,
    is.finite(.data$mean_filtrate_b12_per_adjusted_live_stain_glucose_10),
    .data$mean_filtrate_b12_per_adjusted_live_stain_glucose_10 > 0
  )

glucose_fold_source <- glucose_pair_source |>
  dplyr::transmute(
    strain_id = .data$strain_id,
    b12_trait_group2 = factor(.data$b12_trait_group2, levels = trait_group_levels),
    glucose_0_value = .data$mean_filtrate_b12_per_adjusted_live_stain_glucose_0,
    glucose_10_value = .data$mean_filtrate_b12_per_adjusted_live_stain_glucose_10,
    glucose_10_over_0 = .data$glucose_10_value / .data$glucose_0_value,
    normalization = "adjusted live stain fluorescence"
  ) |>
  dplyr::filter(
    !is.na(.data$b12_trait_group2),
    is.finite(.data$glucose_10_over_0),
    .data$glucose_10_over_0 > 0
  )

# Figure S14A-D: strain-level four-hour endpoints.
endpoint_b12 <- figure_b12_summary |>
  dplyr::filter(
    .data$experiment_type == "release_and_uptake",
    .data$strain_id != "blank",
    .data$sample_time_minute == 240
  ) |>
  dplyr::left_join(figure_taxonomy, by = "strain_id") |>
  dplyr::left_join(figure_traits, by = "strain_id") |>
  dplyr::mutate(
    gtdb_order = dplyr::coalesce(.data$gtdb_order, "Unassigned"),
    b12_trait_group2 = dplyr::coalesce(.data$b12_trait_group2, "Unassigned"),
    b12_trait_group2 = factor(
      .data$b12_trait_group2, levels = c(trait_group_levels, "Unassigned")
    ),
    added_b12_pM = factor(
      .data$added_b12_pM, levels = sort(unique(.data$added_b12_pM))
    ),
    added_glucose_mM_C = factor(
      .data$added_glucose_mM_C, levels = sort(unique(.data$added_glucose_mM_C))
    )
  ) |>
  droplevels()

prepare_endpoint_measurement <- function(measurement) {
  measurement <- rlang::ensym(measurement)
  figure_joined_results |>
    dplyr::filter(
      .data$experiment_type == "release_and_uptake",
      .data$strain_id != "blank",
      .data$sample_time_hour == 4,
      .data$added_b12_pM == 0,
      is.finite(!!measurement)
    ) |>
    dplyr::select(
      "strain_id", "added_glucose_mM_C", measurement = !!measurement
    ) |>
    dplyr::left_join(figure_traits, by = "strain_id") |>
    dplyr::group_by(.data$strain_id, .data$b12_trait_group2, .data$added_glucose_mM_C) |>
    dplyr::summarise(
      mean_measurement = mean(.data$measurement, na.rm = TRUE),
      n_cultures = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      b12_trait_group2 = factor(.data$b12_trait_group2, levels = trait_group_levels),
      added_glucose_mM_C = factor(.data$added_glucose_mM_C, levels = c(0, 10))
    ) |>
    dplyr::filter(!is.na(.data$b12_trait_group2), is.finite(.data$mean_measurement)) |>
    droplevels()
}

write_figure_source <- function(data, filename) {
  utils::write.csv(
    data, file.path(figure_source_dir, filename), row.names = FALSE, na = ""
  )
}

write_figure_source(
  figure_s9_source,
  "06F_live_dead_mix_viability_stain_vs_killed_proportion_source_od600_trait_group2_strain_facets_source.csv"
)
write_figure_source(
  dplyr::filter(endpoint_b12, .data$added_b12_pM == 0),
  "12G_release_uptake_b12_240min_endpoint_all_strains_b12_0_trait_group2_source.csv"
)
write_figure_source(
  dplyr::filter(endpoint_b12, .data$added_b12_pM == 500),
  "12H_release_uptake_b12_240min_endpoint_all_strains_b12_500_trait_group2_source.csv"
)
write_figure_source(
  prepare_endpoint_measurement(adjusted_od600),
  "12I_release_uptake_adjusted_od600_240min_b12_0_trait_group2_source.csv"
)
write_figure_source(
  prepare_endpoint_measurement(adjusted_live_stain),
  "12K_release_uptake_adjusted_live_stain_240min_b12_0_trait_group2_log_y_source.csv"
)
write_figure_source(
  glucose_pair_source,
  "09E_filtrate_b12_per_adjusted_live_stain_glucose_10_vs_0_source.csv"
)
write_figure_source(
  glucose_fold_source,
  "09H_filtrate_b12_per_adjusted_live_stain_glucose_10_over_0_by_trait_group2_source.csv"
)

# The kinetic parameter comparison uses the same release/uptake-strain subset
# of the survey phenotype table. Generate it here because GTDB taxonomy is an
# explicit upstream product by this point in the pipeline.
selected_release_strains <- utils::read.csv(
  file.path(
    "data/raw/b12_release_uptake_death/metadata",
    "strains_selected_b12_release_uptake_20260520.csv"
  ),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
)
trait_taxonomy <- utils::read.csv(
  file.path("data/intermediate/gtdbtk", "gtdbtk_bac120_taxonomy_with_strain_id.csv"),
  check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA")
) |>
  dplyr::select("strain_id", "gtdb_phylum", "gtdb_class", "gtdb_order") |>
  dplyr::distinct(.data$strain_id, .keep_all = TRUE)
trait_release_source <- figure_trait_summary |>
  dplyr::filter(
    .data$n_samples >= 3L,
    .data$strain_id %in% selected_release_strains$strain_id,
    is.finite(.data$total_b12_gm), .data$total_b12_gm > 50,
    is.finite(.data$dead_proportion_gm), .data$dead_proportion_gm > 0,
    is.finite(.data$dead_proportion_gsd), .data$dead_proportion_gsd > 0,
    is.finite(.data$filtrate_proportion_gm), .data$filtrate_proportion_gm > 0,
    is.finite(.data$filtrate_proportion_gsd), .data$filtrate_proportion_gsd > 0,
    is.finite(.data$uptake_b12_mean)
  ) |>
  dplyr::left_join(trait_taxonomy, by = "strain_id") |>
  dplyr::mutate(
    gtdb_order = dplyr::coalesce(.data$gtdb_order, "Unassigned"),
    b12_trait_group2 = factor(.data$b12_trait_group2, levels = trait_group_levels)
  )
trait_figure_source_dir <- file.path(
  "data/processed/figure_source_data", "b12_trait_survey"
)
dir.create(trait_figure_source_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  trait_release_source,
  file.path(
    trait_figure_source_dir,
    "09E_fraction_dead_vs_extracellular_b12_fraction_release_uptake_strains_gtdb_order_source.csv"
  ),
  row.names = FALSE, na = ""
)

message("Wrote empirical compiled-figure source tables to: ", figure_source_dir)

# Compact source-export stage adapted from the approved modeling figure analysis.
section_name <- "b12_release_uptake_death"

figure_source_dir <- file.path("data/processed/figure_source_data", "b12_release_uptake_death_modeling")
dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)

b12_results_path <- file.path("data/processed", section_name, "b12_results_b12_release_uptake_death.csv")

b12_summary_path <- file.path("data/processed", section_name, "b12_summary_b12_release_uptake_death.csv")

b12_delta_estimates_path <- file.path("data/processed", section_name, "modeling", "b12_delta_estimates_b12_release_uptake_death.csv")

b12_fits_uptake_path <- file.path("data/processed", section_name, "modeling", "b12_fits_uptake_b12_release_uptake_death.csv")

b12_fits_release_path <- file.path("data/processed", section_name, "modeling", "b12_fits_release_b12_release_uptake_death.csv")

b12_mix_pred_uptake_path <- file.path("data/processed", section_name, "modeling", "b12_mix_pred_uptake_b12_release_uptake_death.csv")

b12_mix_pred_release_path <- file.path("data/processed", section_name, "modeling", "b12_mix_pred_release_b12_release_uptake_death.csv")

trait_survey_fraction_dead_release_path <- file.path("data/processed/figure_source_data", "b12_trait_survey", "09E_fraction_dead_vs_extracellular_b12_fraction_release_uptake_strains_gtdb_order_source.csv")

gtdb_metadata_path <- file.path("data/intermediate", "gtdbtk", "gtdbtk_bac120_taxonomy_with_strain_id.csv")

trait_metadata_path <- file.path("data/processed", "b12_trait_survey", "b12_trait_survey_strain_summary.csv")

read_csv_or_stop <- function(path) {
    if (!file.exists(path)) {
        stop("Missing input file: ", path)
    }
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

coerce_numeric <- function(df, cols) {
    for (col in cols) {
        if (col %in% names(df)) {
            df[[col]] <- as.numeric(df[[col]])
        }
    }
    df
}

plot_theme <- ggplot2::theme_bw(base_size = 11) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey95"), 
    plot.title = ggplot2::element_text(face = "bold", size = 12), axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, 
        vjust = 1))

b12_results <- read_csv_or_stop(b12_results_path)

b12_summary <- read_csv_or_stop(b12_summary_path)

b12_delta_estimates <- read_csv_or_stop(b12_delta_estimates_path)

b12_fits_uptake <- read_csv_or_stop(b12_fits_uptake_path)

b12_fits_release <- read_csv_or_stop(b12_fits_release_path)

b12_mix_pred_uptake <- read_csv_or_stop(b12_mix_pred_uptake_path)

b12_mix_pred_release <- read_csv_or_stop(b12_mix_pred_release_path)

trait_survey_fraction_dead_release <- read_csv_or_stop(trait_survey_fraction_dead_release_path)

gtdb_metadata <- if (file.exists(gtdb_metadata_path)) {
    dplyr::distinct(dplyr::select(read_csv_or_stop(gtdb_metadata_path), strain_id, gtdb_order))
} else {
    data.frame(strain_id = character(), gtdb_order = character())
}

trait_metadata <- if (file.exists(trait_metadata_path)) {
    dplyr::distinct(dplyr::select(read_csv_or_stop(trait_metadata_path), strain_id, b12_trait_group2))
} else {
    data.frame(strain_id = character(), b12_trait_group2 = character())
}

gtdb_order_levels <- sort(unique(stats::na.omit(gtdb_metadata$gtdb_order)))

gtdb_order_levels <- c(setdiff(gtdb_order_levels, c("Multiple", "Unassigned")), intersect(gtdb_order_levels, "Multiple"), 
    intersect(gtdb_order_levels, "Unassigned"))

gtdb_order_palette <- c("#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#b07aa1", "#76b7b2", "#edc948", "#af7aa1", "#ff9da7", 
    "#9c755f", "#bab0ac", "#86bc86", "#a0cbe8", "#ffbe7d", "#8cd17d")

gtdb_order_colors <- stats::setNames(rep(gtdb_order_palette, length.out = length(gtdb_order_levels)), gtdb_order_levels)

b12_trait_group2_levels <- c("retainer", "provider", "reclaimer")

b12_trait_group2_colors <- c(retainer = "#7570b3", provider = "#e7298a", reclaimer = "#d95f02", Unassigned = "grey55")

gtdb_order_color_scale <- function() {
    ggplot2::scale_color_manual(values = gtdb_order_colors, breaks = names(gtdb_order_colors), drop = TRUE, name = "GTDB order")
}

b12_trait_group2_color_scale <- function() {
    ggplot2::scale_color_manual(values = b12_trait_group2_colors, breaks = names(b12_trait_group2_colors), drop = TRUE, name = "Trait group")
}

order_by_trait_group2 <- function(dat) {
    strain_levels <- dplyr::pull(dplyr::arrange(dplyr::mutate(dplyr::distinct(dat, strain_id, b12_trait_group2), b12_trait_group2 = factor(as.character(b12_trait_group2), 
        levels = c(b12_trait_group2_levels, "Unassigned"))), b12_trait_group2, strain_id), strain_id)
    dplyr::mutate(dat, b12_trait_group2 = factor(as.character(b12_trait_group2), levels = c(b12_trait_group2_levels, "Unassigned")), 
        strain_facet = factor(strain_id, levels = unique(strain_levels)))
}

b12_results <- coerce_numeric(b12_results, c("added_b12_pM", "added_glucose_mM_C", "sample_time_minute", "killed_proportion"))

b12_delta_estimates <- coerce_numeric(b12_delta_estimates, c("added_glucose_mM_C", "biological_replicate", "delta"))

b12_fits_uptake <- coerce_numeric(b12_fits_uptake, c("k_per_cell", "k_per_cell_se"))

b12_fits_release <- coerce_numeric(b12_fits_release, c("k_per_cell", "k_per_cell_se", "E_per_cell", "E_per_cell_se", "r_per_cell", 
    "r_per_cell_se"))

b12_mix_pred_uptake <- coerce_numeric(b12_mix_pred_uptake, c("added_b12_pM", "added_glucose_mM_C", "killed_proportion", "sample_time_minute", 
    "predicted_b12", "b12_pm"))

b12_mix_pred_release <- coerce_numeric(b12_mix_pred_release, c("added_b12_pM", "added_glucose_mM_C", "killed_proportion", 
    "sample_time_minute", "predicted_b12", "b12_pm"))

trait_survey_fraction_dead_release <- coerce_numeric(trait_survey_fraction_dead_release, c("dead_proportion_gm", "dead_proportion_gsd", 
    "filtrate_proportion_gm", "filtrate_proportion_gsd"))

b12_results <- dplyr::mutate(dplyr::left_join(dplyr::left_join(b12_results, gtdb_metadata, by = "strain_id"), trait_metadata, 
    by = "strain_id"), strain_replicate = interaction(experiment_name, biological_replicate), gtdb_order = dplyr::if_else(is.na(gtdb_order) | 
    gtdb_order == "", "Unassigned", gtdb_order), b12_trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == 
    "", "Unassigned", b12_trait_group2))

blank_rows <- dplyr::filter(b12_results, strain_id == "blank", added_b12_pM == 0)

blank_offset <- mean(blank_rows$b12_pm, na.rm = TRUE)

b12_results <- dplyr::filter(b12_results, !is.na(b12_pm))

join_by_fit <- intersect(c("strain_id", "added_glucose_mM_C", "experiment_name", "biological_replicate"), names(b12_fits_uptake))

join_by_release <- intersect(c("strain_id", "added_glucose_mM_C", "experiment_name", "biological_replicate"), names(b12_fits_release))

join_by_delta <- intersect(c("strain_id", "added_glucose_mM_C", "experiment_name", "biological_replicate"), names(b12_delta_estimates))

geometric_mean <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0) {
        return(NA_real_)
    }
    exp(mean(log(x), na.rm = TRUE))
}

trapezoid_auc <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(unique(x)) < 2) {
        return(NA_real_)
    }
    dat <- dplyr::arrange(dplyr::summarise(dplyr::group_by(data.frame(x = x, y = y), x), y = mean(y, na.rm = TRUE), .groups = "drop"), 
        x)
    sum(diff(dat$x) * (head(dat$y, -1) + tail(dat$y, -1))/2)
}

glucose_shape_values <- c(`0` = 16, `10` = 17)

glucose_linetype_values <- c(`0` = "solid", `10` = "22")

add_glucose_factor <- function(dat) {
    dplyr::mutate(dat, added_glucose_mM_C_factor = factor(as.character(added_glucose_mM_C), levels = sort(unique(as.character(added_glucose_mM_C)))))
}

glucose_shape_scale <- function() {
    ggplot2::scale_shape_manual(values = glucose_shape_values, drop = TRUE, name = "Glucose (mM C)")
}

glucose_linetype_scale <- function() {
    ggplot2::scale_linetype_manual(values = glucose_linetype_values, drop = TRUE, name = "Glucose (mM C)")
}

summarise_timecourse <- function(dat) {
    add_glucose_factor(dplyr::mutate(dplyr::summarise(dplyr::group_by(dat, model_type, strain_id, gtdb_order, b12_trait_group2, 
        added_b12_pM, added_glucose_mM_C, sample_time_minute), b12_mean = mean(b12_pm, na.rm = TRUE), b12_sd = stats::sd(b12_pm, 
        na.rm = TRUE), n_obs = dplyr::n(), model_pred = mean(model_pred, na.rm = TRUE), .groups = "drop"), b12_sd = dplyr::if_else(is.na(b12_sd), 
        0, b12_sd), b12_ymin = pmax(b12_mean - b12_sd, 1e-06), b12_ymax = b12_mean + b12_sd))
}

summarise_mix_fraction <- function(dat) {
    full_kill_empirical <- dplyr::summarise(dplyr::group_by(dplyr::filter(dat, killed_proportion == 1, is.finite(b12_pm)), 
        model_type, strain_id, added_glucose_mM_C, experiment_name, biological_replicate), b12_full_kill_empirical = stats::median(b12_pm, 
        na.rm = TRUE), .groups = "drop")
    full_kill_predicted <- dplyr::summarise(dplyr::group_by(dplyr::filter(dat, killed_proportion == 1, is.finite(predicted_b12)), 
        model_type, strain_id, added_glucose_mM_C, experiment_name, biological_replicate), predicted_full_kill = stats::median(predicted_b12, 
        na.rm = TRUE), .groups = "drop")
    add_glucose_factor(dplyr::mutate(dplyr::summarise(dplyr::group_by(dplyr::mutate(dplyr::left_join(dplyr::left_join(dat, 
        full_kill_empirical, by = c("model_type", "strain_id", "added_glucose_mM_C", "experiment_name", "biological_replicate")), 
        full_kill_predicted, by = c("model_type", "strain_id", "added_glucose_mM_C", "experiment_name", "biological_replicate")), 
        b12_frac = dplyr::if_else(is.finite(b12_pm) & !is.na(b12_full_kill_empirical) & b12_full_kill_empirical > 0, b12_pm/b12_full_kill_empirical, 
            NA_real_), predicted_b12_frac = dplyr::if_else(is.finite(predicted_b12) & !is.na(predicted_full_kill) & predicted_full_kill > 
            0, predicted_b12/predicted_full_kill, NA_real_)), model_type, strain_id, gtdb_order, b12_trait_group2, added_glucose_mM_C, 
        killed_proportion), b12_mean = mean(b12_pm, na.rm = TRUE), b12_sd = stats::sd(b12_pm, na.rm = TRUE), n_obs = sum(is.finite(b12_frac)), 
        n_obs_raw = sum(is.finite(b12_pm)), b12_frac_mean = mean(b12_frac, na.rm = TRUE), b12_frac_sd = stats::sd(b12_frac, 
            na.rm = TRUE), predicted_b12 = mean(predicted_b12, na.rm = TRUE), predicted_b12_frac = mean(predicted_b12_frac, 
            na.rm = TRUE), n_pred = sum(is.finite(predicted_b12_frac)), .groups = "drop"), b12_mean = dplyr::if_else(is.finite(b12_mean), 
        b12_mean, NA_real_), b12_sd = dplyr::if_else(is.na(b12_sd), 0, b12_sd), b12_frac_mean = dplyr::if_else(is.finite(b12_frac_mean), 
        b12_frac_mean, NA_real_), b12_frac_sd = dplyr::if_else(is.na(b12_frac_sd), 0, b12_frac_sd), predicted_b12 = dplyr::if_else(is.finite(predicted_b12), 
        predicted_b12, NA_real_), predicted_b12_frac = dplyr::if_else(is.finite(predicted_b12_frac), predicted_b12_frac, 
        NA_real_), b12_frac_ymin = pmax(b12_frac_mean - b12_frac_sd, 1e-06), b12_frac_ymax = b12_frac_mean + b12_frac_sd))
}

make_timecourse_plot <- function(dat, title, subtitle, facet_trait = FALSE, show_glucose = FALSE, log_y = FALSE, show_b12_reference = TRUE) {
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = sample_time_minute, color = gtdb_order)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = b12_ymin, 
        ymax = b12_ymax), width = 0, alpha = 0.45, linewidth = 0.35)
    if (show_glucose) {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = b12_mean, shape = added_glucose_mM_C_factor), alpha = 0.75, size = 1.5) + 
            ggplot2::geom_line(ggplot2::aes(y = model_pred, linetype = added_glucose_mM_C_factor, group = interaction(strain_id, 
                added_glucose_mM_C_factor)), linewidth = 0.8) + glucose_shape_scale() + glucose_linetype_scale()
    }
    else {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = b12_mean), alpha = 0.75, size = 1.5) + ggplot2::geom_line(ggplot2::aes(y = model_pred, 
            group = strain_id), linewidth = 0.8)
    }
    if (show_b12_reference) {
        p <- p + ggplot2::geom_hline(yintercept = 500, linetype = "dashed", color = "grey70", linewidth = 0.4)
    }
    p <- p + gtdb_order_color_scale() + ggplot2::labs(title = title, subtitle = subtitle, x = "Time (minutes)", y = if (log_y) 
        "Filtrate B12 (pM, log scale)"
    else "Filtrate B12 (pM)", color = "GTDB order") + plot_theme
    if (facet_trait) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(b12_trait_group2, strain_facet), ncol = 6)
    }
    else {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(strain_id), ncol = 6)
    }
    if (log_y) {
        p <- p + ggplot2::scale_y_log10()
    }
    p
}

make_timecourse_trait_plot <- function(dat, title, subtitle, color_by = c("gtdb_order", "trait_group2"), show_glucose = FALSE, 
    show_b12_reference = TRUE, x_label = "Time (minutes)", y_label = "Filtrate B12 (pM)") {
    color_by <- match.arg(color_by)
    dat <- dplyr::mutate(dat, b12_trait_group2 = factor(as.character(b12_trait_group2), levels = c(b12_trait_group2_levels, 
        "Unassigned")))
    if (color_by == "trait_group2") {
        p <- ggplot2::ggplot(dat, ggplot2::aes(x = sample_time_minute, color = b12_trait_group2)) + b12_trait_group2_color_scale()
        color_label <- "Trait group"
    }
    else {
        p <- ggplot2::ggplot(dat, ggplot2::aes(x = sample_time_minute, color = gtdb_order)) + gtdb_order_color_scale()
        color_label <- "GTDB order"
    }
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = b12_ymin, ymax = b12_ymax, group = if (show_glucose) 
        interaction(strain_id, added_glucose_mM_C_factor)
    else strain_id), width = 0, alpha = 0.35, linewidth = 0.3)
    if (show_glucose) {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = b12_mean, shape = added_glucose_mM_C_factor, group = interaction(strain_id, 
            added_glucose_mM_C_factor)), alpha = 0.65, size = 1.5) + ggplot2::geom_line(ggplot2::aes(y = model_pred, linetype = added_glucose_mM_C_factor, 
            group = interaction(strain_id, added_glucose_mM_C_factor)), linewidth = 0.8, alpha = 0.85) + glucose_shape_scale() + 
            glucose_linetype_scale()
    }
    else {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = b12_mean, group = strain_id), alpha = 0.65, size = 1.5) + ggplot2::geom_line(ggplot2::aes(y = model_pred, 
            group = strain_id), linewidth = 0.8, alpha = 0.85)
    }
    if (show_b12_reference) {
        p <- p + ggplot2::geom_hline(yintercept = 500, linetype = "dashed", color = "grey70", linewidth = 0.4)
    }
    p + ggplot2::facet_wrap(ggplot2::vars(b12_trait_group2), nrow = 1) + ggplot2::labs(title = title, subtitle = subtitle, 
        x = x_label, y = y_label, color = color_label) + plot_theme
}

make_mix_plot <- function(dat, title, subtitle, facet_trait = FALSE, show_glucose = FALSE, log_y = FALSE) {
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = killed_proportion, color = gtdb_order)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = b12_frac_ymin, 
        ymax = b12_frac_ymax), width = 0.015, alpha = 0.45, linewidth = 0.35)
    if (show_glucose) {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = b12_frac_mean, shape = added_glucose_mM_C_factor), alpha = 0.75, size = 1.5) + 
            ggplot2::geom_line(ggplot2::aes(y = predicted_b12_frac, linetype = added_glucose_mM_C_factor, group = interaction(strain_id, 
                added_glucose_mM_C_factor)), linewidth = 0.8) + glucose_shape_scale() + glucose_linetype_scale()
    }
    else {
        p <- p + ggplot2::geom_point(ggplot2::aes(y = b12_frac_mean), alpha = 0.75, size = 1.5) + ggplot2::geom_line(ggplot2::aes(y = predicted_b12_frac, 
            group = strain_id), linewidth = 0.8)
    }
    p <- p + gtdb_order_color_scale() + ggplot2::labs(title = title, subtitle = subtitle, x = "Fraction dead", y = if (log_y) 
        "Extracellular B12 fraction (log scale)"
    else "Extracellular B12 fraction", color = "GTDB order") + plot_theme
    if (facet_trait) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(b12_trait_group2, strain_facet), ncol = 6)
    }
    else {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(strain_id), ncol = 6)
    }
    if (log_y) {
        p <- p + ggplot2::scale_y_log10()
    }
    p
}

make_mix_trait_plot <- function(dat, title, subtitle, color_by = c("gtdb_order", "trait_group2"), x_label = "Fraction dead", 
    y_label = "Extracellular B12 fraction") {
    color_by <- match.arg(color_by)
    dat <- dplyr::mutate(dat, b12_trait_group2 = factor(as.character(b12_trait_group2), levels = c(b12_trait_group2_levels, 
        "Unassigned")))
    if (color_by == "trait_group2") {
        p <- ggplot2::ggplot(dat, ggplot2::aes(x = killed_proportion, color = b12_trait_group2)) + b12_trait_group2_color_scale()
        color_label <- "Trait group"
    }
    else {
        p <- ggplot2::ggplot(dat, ggplot2::aes(x = killed_proportion, color = gtdb_order)) + gtdb_order_color_scale()
        color_label <- "GTDB order"
    }
    p + ggplot2::geom_errorbar(ggplot2::aes(ymin = b12_frac_ymin, ymax = b12_frac_ymax, group = strain_id), width = 0.015, 
        alpha = 0.35, linewidth = 0.3) + ggplot2::geom_point(ggplot2::aes(y = b12_frac_mean, group = strain_id), alpha = 0.65, 
        size = 1.5) + ggplot2::geom_line(ggplot2::aes(y = predicted_b12_frac, group = strain_id), linewidth = 0.8, alpha = 0.85) + 
        ggplot2::facet_wrap(ggplot2::vars(b12_trait_group2), nrow = 1) + ggplot2::labs(title = title, subtitle = subtitle, 
        x = x_label, y = y_label, color = color_label) + plot_theme
}

prepare_representative_plot_source <- function(dat, representatives) {
    representatives <- dplyr::arrange(dplyr::mutate(representatives, b12_trait_group2 = as.character(b12_trait_group2)), 
        factor(b12_trait_group2, levels = b12_trait_group2_levels), strain_id)
    representative_levels <- dplyr::pull(dplyr::mutate(representatives, representative_facet = paste(b12_trait_group2, strain_id, 
        sep = "_")), representative_facet)
    dplyr::mutate(dplyr::inner_join(dplyr::mutate(dat, b12_trait_group2 = as.character(b12_trait_group2)), dplyr::select(representatives, 
        b12_trait_group2, strain_id, selection_method, selection_description, selection_score, dplyr::any_of(c("distance_to_median_abs_auc_error", 
            "median_abs_auc_error", "empirical_auc", "uptake_auc", "abs_auc_error"))), by = c("b12_trait_group2", "strain_id")), 
        b12_trait_group2 = factor(b12_trait_group2, levels = b12_trait_group2_levels), representative_facet = factor(paste(as.character(b12_trait_group2), 
            strain_id, sep = "_"), levels = representative_levels))
}

uptake_raw <- dplyr::filter(dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::filter(b12_results, experiment_type == 
    "release_and_uptake", added_b12_pM > 0), b12_fits_uptake, by = join_by_fit), b12_delta_estimates, by = join_by_delta, 
    suffix = c("", "")), model_pred = blank_offset + (added_b12_pM + delta) * exp(-k_per_cell * sample_time_minute), model_type = "uptake_model"), 
    is.finite(sample_time_minute), is.finite(b12_pm), is.finite(model_pred))

uptake_summary <- summarise_timecourse(uptake_raw)

uptake_summary_glucose0 <- dplyr::filter(uptake_summary, added_glucose_mM_C == 0)

uptake_summary_glucose0_nonblank <- dplyr::filter(uptake_summary_glucose0, strain_id != "blank")

trait_group2_strain_facet_levels_01i <- c("NA", b12_trait_group2_levels)

uptake_summary_trait_group2_strain_facet_01i <- droplevels(dplyr::arrange(dplyr::mutate(dplyr::left_join(dplyr::filter(b12_summary, 
    experiment_type == "release_and_uptake", added_b12_pM == 500, added_glucose_mM_C == 0), trait_metadata, by = "strain_id"), 
    edge_correction = "quadratic_edge_pairs", trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == 
        "", "NA", b12_trait_group2), trait_group2 = factor(trait_group2, levels = trait_group2_strain_facet_levels_01i), 
    trait_group2_strain_id = paste(trait_group2, strain_id, sep = "_"), trait_group2_strain_id = factor(trait_group2_strain_id, 
        levels = dplyr::distinct(dplyr::arrange(data.frame(trait_group2 = trait_group2, strain_id = strain_id, trait_group2_strain_id = trait_group2_strain_id), 
            trait_group2, strain_id), trait_group2_strain_id)$trait_group2_strain_id)), trait_group2, strain_id, sample_time_minute))

uptake_model_trait_group2_strain_facet_01j <- droplevels(dplyr::arrange(dplyr::mutate(uptake_summary_glucose0, trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | 
    b12_trait_group2 == "" | b12_trait_group2 == "Unassigned", "NA", b12_trait_group2), trait_group2 = factor(trait_group2, 
    levels = trait_group2_strain_facet_levels_01i), trait_group2_strain_id = paste(trait_group2, strain_id, sep = "_"), trait_group2_strain_id = factor(trait_group2_strain_id, 
    levels = levels(uptake_summary_trait_group2_strain_facet_01i$trait_group2_strain_id))), trait_group2, strain_id, sample_time_minute))

uptake_overlay_trait_group2_strain_facet_01k <- dplyr::left_join(uptake_model_trait_group2_strain_facet_01j, dplyr::transmute(uptake_summary_trait_group2_strain_facet_01i, 
    strain_id, sample_time_minute, empirical_mean_b12 = mean_b12, empirical_sd_b12 = sd_b12), by = c("strain_id", "sample_time_minute"))

write.csv(uptake_overlay_trait_group2_strain_facet_01k, file.path(figure_source_dir, "uptake_fit_dynamics_source_01K.csv"), 
    row.names = FALSE)

release_uptake_raw <- dplyr::filter(dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::filter(b12_results, experiment_type == 
    "release_and_uptake", added_b12_pM == 0), b12_fits_release, by = join_by_release), b12_delta_estimates, by = join_by_delta, 
    suffix = c("", "")), model_pred = blank_offset + E_per_cell + (added_b12_pM + delta - E_per_cell) * exp(-k_per_cell * 
    sample_time_minute), model_type = "uptake_release_model"), is.finite(sample_time_minute), is.finite(b12_pm), is.finite(model_pred))

release_uptake_summary <- summarise_timecourse(release_uptake_raw)

release_uptake_summary_nonblank <- dplyr::filter(release_uptake_summary, strain_id != "blank")

release_uptake_summary_trait_group2_strain_facet_02o <- droplevels(dplyr::arrange(add_glucose_factor(dplyr::mutate(dplyr::left_join(dplyr::filter(b12_summary, 
    experiment_type == "release_and_uptake", added_b12_pM == 0, added_glucose_mM_C %in% c(0, 10)), trait_metadata, by = "strain_id"), 
    edge_correction = "quadratic_edge_pairs", trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == 
        "", "NA", b12_trait_group2), trait_group2 = factor(trait_group2, levels = trait_group2_strain_facet_levels_01i), 
    trait_group2_strain_id = paste(trait_group2, strain_id, sep = "_"), trait_group2_strain_id = factor(trait_group2_strain_id, 
        levels = levels(uptake_summary_trait_group2_strain_facet_01i$trait_group2_strain_id)))), trait_group2, strain_id, 
    added_glucose_mM_C, sample_time_minute))

write.csv(release_uptake_summary, file.path(figure_source_dir, "release_uptake_fit_dynamics_source_02A.csv"), row.names = FALSE)

write.csv(release_uptake_summary_trait_group2_strain_facet_02o, file.path(figure_source_dir, "release_uptake_fit_dynamics_source_02O.csv"), 
    row.names = FALSE)

release_fit_params_replicate <- add_glucose_factor(dplyr::filter(dplyr::mutate(dplyr::filter(tidyr::pivot_longer(dplyr::select(dplyr::mutate(dplyr::left_join(dplyr::left_join(b12_fits_release, 
    gtdb_metadata, by = "strain_id"), trait_metadata, by = "strain_id"), gtdb_order = dplyr::if_else(is.na(gtdb_order) | 
    gtdb_order == "", "Unassigned", gtdb_order), b12_trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == 
    "", "Unassigned", b12_trait_group2)), strain_id, added_glucose_mM_C, dplyr::any_of(c("experiment_name", "biological_replicate")), 
    k_per_cell, r_per_cell, E_per_cell, gtdb_order, b12_trait_group2), cols = dplyr::all_of(c("k_per_cell", "r_per_cell", 
    "E_per_cell")), names_to = "parameter", values_to = "estimate"), is.finite(estimate), estimate > 0), parameter = dplyr::recode(parameter, 
    k_per_cell = "k", r_per_cell = "r", E_per_cell = "E")), parameter %in% c("k", "r")))

release_fit_params <- add_glucose_factor(dplyr::summarise(dplyr::group_by(release_fit_params_replicate, strain_id, added_glucose_mM_C, 
    gtdb_order, b12_trait_group2, parameter), estimate = geometric_mean(estimate), n_parameter_fits = dplyr::n(), .groups = "drop"))

release_fit_params_k_all_glucose <- order_by_trait_group2(dplyr::filter(release_fit_params, strain_id != "blank", parameter == 
    "k", as.character(b12_trait_group2) %in% b12_trait_group2_levels))

release_fit_params_k_glucose0 <- order_by_trait_group2(dplyr::filter(release_fit_params, strain_id != "blank", added_glucose_mM_C == 
    0, parameter == "k"))

release_fit_params_r_all_glucose <- order_by_trait_group2(dplyr::filter(release_fit_params, strain_id != "blank", parameter == 
    "r", as.character(b12_trait_group2) %in% b12_trait_group2_levels))

release_fit_params_k_r_glucose_geomean <- order_by_trait_group2(dplyr::filter(dplyr::rename(tidyr::pivot_wider(dplyr::filter(dplyr::summarise(dplyr::group_by(dplyr::filter(release_fit_params, 
    strain_id != "blank", parameter %in% c("k", "r"), added_glucose_mM_C %in% c(0, 10), as.character(b12_trait_group2) %in% 
        b12_trait_group2_levels, is.finite(estimate), estimate > 0), strain_id, gtdb_order, b12_trait_group2, parameter), 
    estimate = geometric_mean(estimate), n_glucose_levels = dplyr::n_distinct(added_glucose_mM_C), glucose_values = paste(sort(unique(added_glucose_mM_C)), 
        collapse = ", "), n_parameter_fits = sum(n_parameter_fits, na.rm = TRUE), .groups = "drop"), n_glucose_levels == 
    2), names_from = parameter, values_from = c(estimate, n_glucose_levels, glucose_values, n_parameter_fits), names_glue = "{parameter}_{.value}"), 
    k = k_estimate, r = r_estimate), is.finite(k), is.finite(r)))

release_fit_params_r_by_glucose <- order_by_trait_group2(dplyr::filter(dplyr::rename(tidyr::pivot_wider(dplyr::select(dplyr::filter(release_fit_params_r_all_glucose, 
    added_glucose_mM_C %in% c(0, 10), is.finite(estimate), estimate > 0), strain_id, gtdb_order, b12_trait_group2, added_glucose_mM_C, 
    estimate, n_parameter_fits), names_from = added_glucose_mM_C, values_from = c(estimate, n_parameter_fits), names_glue = "glucose_{added_glucose_mM_C}_{.value}"), 
    r_glucose_0 = glucose_0_estimate, r_glucose_10 = glucose_10_estimate, n_parameter_fits_glucose_0 = glucose_0_n_parameter_fits, 
    n_parameter_fits_glucose_10 = glucose_10_n_parameter_fits), is.finite(r_glucose_0), is.finite(r_glucose_10)))

write.csv(release_fit_params_k_all_glucose, file.path(figure_source_dir, "fit_parameter_columns_source_03C.csv"), row.names = FALSE)

write.csv(release_fit_params_k_glucose0, file.path(figure_source_dir, "fit_parameter_columns_source_03E.csv"), row.names = FALSE)

write.csv(release_fit_params_r_all_glucose, file.path(figure_source_dir, "fit_parameter_columns_source_03H.csv"), row.names = FALSE)

write.csv(release_fit_params_k_r_glucose_geomean, file.path(figure_source_dir, "fit_parameter_columns_source_03N.csv"), row.names = FALSE)

write.csv(release_fit_params_r_by_glucose, file.path(figure_source_dir, "fit_parameter_columns_source_03O_03P.csv"), row.names = FALSE)

prepare_mix_source <- function(pred_df, model_type_label) {
    dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::filter(dplyr::mutate(pred_df, model_type = model_type_label), 
        experiment_type == "live_dead_mix", !is.na(predicted_b12), !is.na(killed_proportion)), gtdb_metadata, by = "strain_id"), 
        trait_metadata, by = "strain_id"), gtdb_order = dplyr::if_else(is.na(gtdb_order) | gtdb_order == "", "Unassigned", 
        gtdb_order), b12_trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == "", "Unassigned", b12_trait_group2), 
        strain_replicate = interaction(experiment_name, biological_replicate))
}

select_median_abs_auc_error_strains <- function(mix_summary) {
    auc_by_strain <- dplyr::filter(dplyr::mutate(dplyr::summarise(dplyr::group_by(dplyr::mutate(dplyr::filter(mix_summary, 
        strain_id != "blank", added_glucose_mM_C == 0, as.character(b12_trait_group2) %in% b12_trait_group2_levels, is.finite(killed_proportion), 
        is.finite(b12_frac_mean), is.finite(predicted_b12_frac)), b12_trait_group2 = as.character(b12_trait_group2)), b12_trait_group2, 
        strain_id, gtdb_order), empirical_auc = trapezoid_auc(killed_proportion, b12_frac_mean), uptake_auc = trapezoid_auc(killed_proportion, 
        predicted_b12_frac), n_fractions = dplyr::n_distinct(killed_proportion), .groups = "drop"), abs_auc_error = abs(uptake_auc - 
        empirical_auc)), is.finite(abs_auc_error))
    median_by_trait_group <- dplyr::summarise(dplyr::group_by(auc_by_strain, b12_trait_group2), median_abs_auc_error = stats::median(abs_auc_error, 
        na.rm = TRUE), .groups = "drop")
    dplyr::mutate(dplyr::ungroup(dplyr::slice(dplyr::arrange(dplyr::group_by(dplyr::mutate(dplyr::left_join(auc_by_strain, 
        median_by_trait_group, by = "b12_trait_group2"), distance_to_median_abs_auc_error = abs(abs_auc_error - median_abs_auc_error), 
        selection_score = abs_auc_error, b12_trait_group2 = factor(b12_trait_group2, levels = b12_trait_group2_levels)), 
        b12_trait_group2), distance_to_median_abs_auc_error, strain_id, .by_group = TRUE), 1)), selection_method = "median_abs_auc_error", 
        selection_description = "closest to trait-group median absolute uptake-model AUC error")
}

mix_uptake_raw <- prepare_mix_source(b12_mix_pred_uptake, "uptake_model")

mix_uptake_summary <- summarise_mix_fraction(mix_uptake_raw)

mix_uptake_summary_nonblank <- dplyr::filter(mix_uptake_summary, strain_id != "blank")

mix_release_raw <- prepare_mix_source(b12_mix_pred_release, "release_model")

mix_release_summary <- summarise_mix_fraction(mix_release_raw)

representative_median_abs_auc_error <- select_median_abs_auc_error_strains(mix_uptake_summary_nonblank)

uptake_representative_median_abs_auc_error <- prepare_representative_plot_source(uptake_summary_glucose0_nonblank, representative_median_abs_auc_error)

release_uptake_representative_median_abs_auc_error_all_glucose <- prepare_representative_plot_source(release_uptake_summary_nonblank, 
    representative_median_abs_auc_error)

mix_uptake_representative_median_abs_auc_error <- prepare_representative_plot_source(mix_uptake_summary_nonblank, representative_median_abs_auc_error)

mix_release_trait_group2_strain_facet_04m <- droplevels(dplyr::arrange(dplyr::mutate(mix_release_summary, trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | 
    b12_trait_group2 == "" | b12_trait_group2 == "Unassigned", "NA", b12_trait_group2), trait_group2 = factor(trait_group2, 
    levels = trait_group2_strain_facet_levels_01i), trait_group2_strain_id = paste(trait_group2, strain_id, sep = "_"), trait_group2_strain_id = factor(trait_group2_strain_id, 
    levels = levels(uptake_summary_trait_group2_strain_facet_01i$trait_group2_strain_id))), trait_group2, strain_id, killed_proportion))

write.csv(mix_release_trait_group2_strain_facet_04m, file.path(figure_source_dir, "dead_fraction_mix_source_04M_04N_04O.csv"), 
    row.names = FALSE)

write.csv(uptake_representative_median_abs_auc_error, file.path(figure_source_dir, "uptake_fit_dynamics_source_01H.csv"), 
    row.names = FALSE)

write.csv(release_uptake_representative_median_abs_auc_error_all_glucose, file.path(figure_source_dir, "release_uptake_fit_dynamics_source_02K.csv"), 
    row.names = FALSE)

write.csv(mix_uptake_representative_median_abs_auc_error, file.path(figure_source_dir, "dead_fraction_mix_source_04L.csv"), 
    row.names = FALSE)

prepare_prediction_error_source <- function(pred_df, killed_min = NULL, killed_max = NULL) {
    full_release <- dplyr::summarise(dplyr::group_by(dplyr::filter(pred_df, experiment_type == "live_dead_mix", strain_id != 
        "blank", killed_proportion == 1, is.finite(b12_pm)), strain_id, added_glucose_mM_C, experiment_name, biological_replicate), 
        scale_full_release = stats::median(b12_pm, na.rm = TRUE), .groups = "drop")
    error_by_fraction <- dplyr::filter(pred_df, experiment_type == "live_dead_mix", strain_id != "blank", killed_proportion < 
        1, is.finite(b12_pm), is.finite(predicted_b12), is.finite(null_predicted_b12))
    if (!is.null(killed_min)) {
        error_by_fraction <- dplyr::filter(error_by_fraction, killed_proportion > killed_min)
    }
    if (!is.null(killed_max)) {
        error_by_fraction <- dplyr::filter(error_by_fraction, killed_proportion < killed_max)
    }
    error_by_fraction <- dplyr::filter(dplyr::mutate(tidyr::pivot_longer(dplyr::filter(dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::left_join(error_by_fraction, 
        gtdb_metadata, by = "strain_id"), trait_metadata, by = "strain_id"), full_release, by = c("strain_id", "added_glucose_mM_C", 
        "experiment_name", "biological_replicate")), gtdb_order = dplyr::if_else(is.na(gtdb_order) | gtdb_order == "", "Unassigned", 
        gtdb_order), b12_trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == "", "Unassigned", b12_trait_group2), 
        b12_trait_group2 = factor(as.character(b12_trait_group2), levels = b12_trait_group2_levels)), !is.na(b12_trait_group2)), 
        cols = dplyr::all_of(c("null_predicted_b12", "predicted_b12")), names_to = "prediction_type", values_to = "predicted"), 
        prediction_type = dplyr::case_when(prediction_type == "null_predicted_b12" ~ "Death only", prediction_type == "predicted_b12" ~ 
            "Death + uptake", TRUE ~ prediction_type), prediction_type = factor(prediction_type, levels = c("Death only", 
            "Death + uptake")), scaled_abs_error = abs(predicted - b12_pm)/scale_full_release), is.finite(scaled_abs_error))
    add_glucose_factor(dplyr::mutate(dplyr::summarise(dplyr::group_by(error_by_fraction, strain_id, gtdb_order, b12_trait_group2, 
        added_glucose_mM_C, prediction_type), mean_abs_scaled_error = mean(scaled_abs_error, na.rm = TRUE), median_abs_scaled_error = stats::median(scaled_abs_error, 
        na.rm = TRUE), n_fractions = dplyr::n_distinct(killed_proportion), n_obs = dplyr::n(), .groups = "drop"), trait_group2 = b12_trait_group2))
}

prepare_auc_error_source <- function(pred_df) {
    full_release <- dplyr::summarise(dplyr::group_by(dplyr::filter(pred_df, experiment_type == "live_dead_mix", strain_id != 
        "blank", killed_proportion == 1, is.finite(b12_pm)), strain_id, added_glucose_mM_C, experiment_name, biological_replicate), 
        empirical_full_kill = stats::median(b12_pm, na.rm = TRUE), .groups = "drop")
    full_prediction <- dplyr::summarise(dplyr::group_by(dplyr::filter(pred_df, experiment_type == "live_dead_mix", strain_id != 
        "blank", killed_proportion == 1, is.finite(null_predicted_b12), is.finite(predicted_b12)), strain_id, added_glucose_mM_C, 
        experiment_name, biological_replicate), death_only_full_kill = stats::median(null_predicted_b12, na.rm = TRUE), uptake_full_kill = stats::median(predicted_b12, 
        na.rm = TRUE), .groups = "drop")
    auc_by_replicate <- dplyr::filter(dplyr::summarise(dplyr::group_by(dplyr::filter(dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::left_join(dplyr::left_join(dplyr::filter(pred_df, 
        experiment_type == "live_dead_mix", strain_id != "blank", is.finite(killed_proportion), is.finite(b12_pm), is.finite(null_predicted_b12), 
        is.finite(predicted_b12)), gtdb_metadata, by = "strain_id"), trait_metadata, by = "strain_id"), full_release, by = c("strain_id", 
        "added_glucose_mM_C", "experiment_name", "biological_replicate")), full_prediction, by = c("strain_id", "added_glucose_mM_C", 
        "experiment_name", "biological_replicate")), gtdb_order = dplyr::if_else(is.na(gtdb_order) | gtdb_order == "", "Unassigned", 
        gtdb_order), b12_trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == "", "Unassigned", b12_trait_group2), 
        b12_trait_group2 = factor(as.character(b12_trait_group2), levels = b12_trait_group2_levels), empirical_fraction = dplyr::if_else(!is.na(empirical_full_kill) & 
            empirical_full_kill > 0, b12_pm/empirical_full_kill, NA_real_), death_only_fraction = dplyr::if_else(!is.na(death_only_full_kill) & 
            death_only_full_kill > 0, null_predicted_b12/death_only_full_kill, NA_real_), uptake_fraction = dplyr::if_else(!is.na(uptake_full_kill) & 
            uptake_full_kill > 0, predicted_b12/uptake_full_kill, NA_real_)), !is.na(b12_trait_group2), is.finite(empirical_fraction), 
        is.finite(death_only_fraction), is.finite(uptake_fraction)), strain_id, gtdb_order, b12_trait_group2, added_glucose_mM_C, 
        experiment_name, biological_replicate), empirical_auc = trapezoid_auc(killed_proportion, empirical_fraction), death_only_auc = trapezoid_auc(killed_proportion, 
        death_only_fraction), uptake_auc = trapezoid_auc(killed_proportion, uptake_fraction), n_fractions = dplyr::n_distinct(killed_proportion), 
        min_fraction = min(killed_proportion, na.rm = TRUE), max_fraction = max(killed_proportion, na.rm = TRUE), .groups = "drop"), 
        is.finite(empirical_auc), is.finite(death_only_auc), is.finite(uptake_auc))
    add_glucose_factor(dplyr::mutate(dplyr::summarise(dplyr::group_by(dplyr::mutate(tidyr::pivot_longer(auc_by_replicate, 
        cols = dplyr::all_of(c("death_only_auc", "uptake_auc")), names_to = "prediction_type", values_to = "model_auc"), 
        prediction_type = dplyr::case_when(prediction_type == "death_only_auc" ~ "Death only", prediction_type == "uptake_auc" ~ 
            "Death + uptake", TRUE ~ prediction_type), prediction_type = factor(prediction_type, levels = c("Death only", 
            "Death + uptake")), signed_auc_error = model_auc - empirical_auc, abs_auc_error = abs(signed_auc_error)), strain_id, 
        gtdb_order, b12_trait_group2, added_glucose_mM_C, prediction_type), empirical_auc = mean(empirical_auc, na.rm = TRUE), 
        model_auc = mean(model_auc, na.rm = TRUE), mean_signed_auc_error = mean(signed_auc_error, na.rm = TRUE), mean_abs_scaled_error = mean(abs_auc_error, 
            na.rm = TRUE), median_abs_scaled_error = stats::median(abs_auc_error, na.rm = TRUE), n_replicates = dplyr::n(), 
        n_fractions_min = min(n_fractions, na.rm = TRUE), n_fractions_max = max(n_fractions, na.rm = TRUE), min_fraction = min(min_fraction, 
            na.rm = TRUE), max_fraction = max(max_fraction, na.rm = TRUE), .groups = "drop"), trait_group2 = b12_trait_group2))
}

auc_error_06e <- prepare_auc_error_source(b12_mix_pred_uptake)

auc_error_06j <- dplyr::filter(auc_error_06e, as.character(prediction_type) == "Death + uptake")

auc_mae_observed_10a <- mean(abs(auc_error_06j$model_auc - auc_error_06j$empirical_auc))

auc_mae_permutation_10a <- data.frame(permutation = seq_len(1000), mae = replicate(1000, mean(abs(sample(auc_error_06j$model_auc) - 
    auc_error_06j$empirical_auc))))

auc_mae_permutation_p_10a <- (1 + sum(auc_mae_permutation_10a$mae <= auc_mae_observed_10a))/(nrow(auc_mae_permutation_10a) + 
    1)

auc_empirical_mean_10a <- mean(auc_error_06j$empirical_auc)

auc_model_mean_10a <- mean(auc_error_06j$model_auc)

auc_ccc_10a <- 2 * mean((auc_error_06j$empirical_auc - auc_empirical_mean_10a) * (auc_error_06j$model_auc - auc_model_mean_10a))/(mean((auc_error_06j$empirical_auc - 
    auc_empirical_mean_10a)^2) + mean((auc_error_06j$model_auc - auc_model_mean_10a)^2) + (auc_empirical_mean_10a - auc_model_mean_10a)^2)

auc_mae_permutation_summary_10a <- data.frame(n_strains = nrow(auc_error_06j), n_permutations = nrow(auc_mae_permutation_10a), 
    concordance_correlation_coefficient = auc_ccc_10a, observed_mae = auc_mae_observed_10a, permutation_p_value = auc_mae_permutation_p_10a, 
    seed = 20260804)

write.csv(auc_error_06j, file.path(figure_source_dir, "empirical_vs_model_auc_source_07J.csv"), row.names = FALSE)

write.csv(auc_mae_permutation_summary_10a, file.path(figure_source_dir, "auc_mae_strain_label_permutation_summary_10A.csv"), 
    row.names = FALSE)

e_by_strain_glucose <- dplyr::summarise(dplyr::group_by(dplyr::filter(b12_fits_release, strain_id != "blank", added_glucose_mM_C %in% 
    c(0, 10), is.finite(E_per_cell), E_per_cell > 0), strain_id, added_glucose_mM_C), E_geomean = geometric_mean(E_per_cell), 
    n_E_fits = dplyr::n(), .groups = "drop")

e_by_strain_all_glucose <- dplyr::filter(dplyr::summarise(dplyr::group_by(e_by_strain_glucose, strain_id), E_geomean = geometric_mean(E_geomean), 
    n_glucose_levels = dplyr::n_distinct(added_glucose_mM_C), glucose_values = paste(sort(unique(added_glucose_mM_C)), collapse = ", "), 
    n_E_fits = sum(n_E_fits, na.rm = TRUE), .groups = "drop"), n_glucose_levels == 2)

prepare_e_fraction_dead_source <- function(e_summary, e_source_label) {
    droplevels(dplyr::filter(dplyr::mutate(dplyr::left_join(trait_survey_fraction_dead_release, e_summary, by = "strain_id"), 
        E_source = e_source_label, b12_trait_group2 = dplyr::if_else(is.na(b12_trait_group2) | b12_trait_group2 == "", "Unassigned", 
            b12_trait_group2), b12_trait_group2 = factor(as.character(b12_trait_group2), levels = c(b12_trait_group2_levels, 
            "Unassigned"))), !is.na(b12_trait_group2), is.finite(dead_proportion_gm), dead_proportion_gm > 0, is.finite(dead_proportion_gsd), 
        dead_proportion_gsd > 0, is.finite(filtrate_proportion_gm), filtrate_proportion_gm > 0, is.finite(filtrate_proportion_gsd), 
        filtrate_proportion_gsd > 0, is.finite(E_geomean), E_geomean > 0))
}

fraction_dead_e_08a <- prepare_e_fraction_dead_source(e_by_strain_all_glucose, "Geometric mean of glucose-specific E values from added_glucose_mM_C = 0 and 10")

fraction_dead_k_over_r_09b <- dplyr::filter(dplyr::mutate(fraction_dead_e_08a, k_over_r = 1/E_geomean, ec_b12_fraction_per_dead_cell_fraction = filtrate_proportion_gm/dead_proportion_gm), 
    is.finite(k_over_r), k_over_r > 0, is.finite(ec_b12_fraction_per_dead_cell_fraction), ec_b12_fraction_per_dead_cell_fraction > 
        0)

regression_09b <- stats::lm(log10(ec_b12_fraction_per_dead_cell_fraction) ~ log10(k_over_r), data = fraction_dead_k_over_r_09b)

regression_summary_09b <- summary(regression_09b)

regression_stats_09b <- data.frame(n_strains = nrow(fraction_dead_k_over_r_09b), intercept = unname(stats::coef(regression_09b)[1]), 
    slope = unname(stats::coef(regression_09b)[2]), r_squared = regression_summary_09b$r.squared, slope_p_value = regression_summary_09b$coefficients[2, 
        "Pr(>|t|)"])

write.csv(fraction_dead_k_over_r_09b, file.path(figure_source_dir, "ec_b12_fraction_per_dead_cell_fraction_vs_k_over_r_source_09B.csv"), 
    row.names = FALSE)

write.csv(regression_stats_09b, file.path(figure_source_dir, "ec_b12_fraction_per_dead_cell_fraction_vs_k_over_r_regression_09B.csv"), 
    row.names = FALSE)
