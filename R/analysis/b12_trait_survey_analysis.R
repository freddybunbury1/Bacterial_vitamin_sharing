section_name <- "b12_trait_survey"

intermediate_dir <- file.path("data/intermediate", section_name)
processed_dir <- file.path("data/processed", section_name)
figure_source_dir <- file.path("data/processed/figure_source_data", section_name)
estimator_processed_dir <- file.path("data/processed", "b12_concentration_estimator")
helper_path <- file.path("R/utils", "b12_trait_survey_helpers.R")
estimator_path <- file.path("R/utils", "b12_concentration_estimator.R")

culture_samples_path <- file.path(intermediate_dir, "culture_samples_b12_trait_survey.csv")
od_viability_summary_path <- file.path(intermediate_dir, "od600_viability_summary_b12_trait_survey.csv")
b12_assay_combined_path <- file.path(intermediate_dir, "b12_assay_combined_b12_trait_survey.csv")

analysis_parameters_path <- file.path(processed_dir, "analysis_parameters_b12_trait_survey.csv")
b12_estimator_output_path <- file.path(processed_dir, "b12_estimator_output_b12_trait_survey.csv")
b12_trait_survey_results_path <- file.path(processed_dir, "b12_trait_survey_results.csv")
b12_trait_survey_results_filtered_path <- file.path(processed_dir, "b12_trait_survey_results_filtered.csv")
b12_trait_survey_time_summary_path <- file.path(processed_dir, "b12_trait_survey_time_summary.csv")
b12_trait_survey_strain_summary_path <- file.path(processed_dir, "b12_trait_survey_strain_summary.csv")
b12_trait_survey_gwas_phenotypes_path <- file.path(processed_dir, "b12_trait_survey_gwas_phenotypes.csv")
filtrate_fraction_dead_uptake_log10_model_results_path <- file.path(
  processed_dir,
  "filtrate_fraction_dead_uptake_log10_lm_results.csv"
)
b12_trait_survey_estimator_diagnostics_path <- file.path(processed_dir, "b12_trait_survey_estimator_diagnostics.csv")
b12_estimator_parameters_path <- file.path(estimator_processed_dir, "b12_estimator_parameters.csv")

standard_sample_types <- c("filtrate")
standard_cols_to_remove <- c(1)
add_offset_pm <- 0.1
edge_correction <- "quadratic_edge_pairs"
upper_thresh_logb12_estimation_se <- 0.5
blank_od_filter_sd_multiplier <- 2
blank_live_stain_floor_sd_multiplier <- 1
added_b12 <- 500
uptake_0hr_error_thresh <- 300
b12_uptake_baseline <- "mean_of_both"
b12_trait_group_total_b12_threshold_pm <- 50
b12_trait_group_uptake_b12_threshold_pm <- 250
b12_trait_group_filtrate_proportion_threshold <- 0.05

source(helper_path)
source(estimator_path)
check_required_packages(c("dplyr", "tidyr", "purrr", "stringr", "drc"))

allowed_b12_uptake_baselines <- c("uptake_0hr", "filtrate_plus_added", "mean_of_both")
if (!b12_uptake_baseline %in% allowed_b12_uptake_baselines) {
  stop(
    "b12_uptake_baseline must be one of: ",
    paste(allowed_b12_uptake_baselines, collapse = ", ")
  )
}

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(estimator_processed_dir, recursive = TRUE, showWarnings = FALSE)

culture_samples <- utils::read.csv(
  culture_samples_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

od_viab_summary <- utils::read.csv(
  od_viability_summary_path,
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

estimator_result <- b12_concentration_estimator(
  b12_assay_combined,
  standard_sample_types = standard_sample_types,
  standard_cols_to_remove = standard_cols_to_remove,
  add_offset_pm = add_offset_pm,
  edge_correction = edge_correction
)

b12_concentration_estimates <- estimator_result$estimates

b12_assay_flags <- estimator_result$corrected_assay |>
  dplyr::filter(assay_type == "sample") |>
  dplyr::group_by(assay_batch, experiment_name, sample_id, sample_type) |>
  dplyr::summarise(
    assay_rows = paste(sort(unique(assay_row)), collapse = ";"),
    assay_wells = paste(sort(unique(assay_well)), collapse = ";"),
    n_assay_observations = dplyr::n(),
    n_edge_corrected_observations = sum(edge_correction_applied, na.rm = TRUE),
    any_edge_correction_applied = any(edge_correction_applied, na.rm = TRUE),
    edge_correction_pairs = paste(sort(unique(stats::na.omit(edge_correction_pair))), collapse = ";"),
    .groups = "drop"
  )

b12_estimator_output <- b12_concentration_estimates |>
  dplyr::left_join(
    b12_assay_flags,
    by = c("assay_batch", "experiment_name", "sample_id", "sample_type")
  ) |>
  dplyr::mutate(
    passed_se_filter = is.finite(se_log_fm) & se_log_fm <= upper_thresh_logb12_estimation_se
  )

b12_summary_wide <- b12_estimator_output |>
  dplyr::filter(passed_se_filter) |>
  dplyr::select(
    sample_id,
    sample_type,
    b12_pm,
    hit_lower,
    hit_upper,
    response_lower_censored,
    response_upper_censored,
    response_any_censored
  ) |>
  tidyr::pivot_wider(
    names_from = sample_type,
    values_from = c(
      b12_pm,
      hit_lower,
      hit_upper,
      response_lower_censored,
      response_upper_censored,
      response_any_censored
    ),
    names_glue = "{sample_type}_{.value}"
  ) |>
  dplyr::rename(
    total_b12 = total_b12_pm,
    filtrate_b12 = filtrate_b12_pm,
    uptake_1hr_b12 = uptake_1hr_b12_pm,
    uptake_0hr_b12 = uptake_0hr_b12_pm,
    total_b12_hit_lower = total_hit_lower,
    filtrate_b12_hit_lower = filtrate_hit_lower,
    uptake_1hr_b12_hit_lower = uptake_1hr_hit_lower,
    uptake_0hr_b12_hit_lower = uptake_0hr_hit_lower,
    total_b12_hit_upper = total_hit_upper,
    filtrate_b12_hit_upper = filtrate_hit_upper,
    uptake_1hr_b12_hit_upper = uptake_1hr_hit_upper,
    uptake_0hr_b12_hit_upper = uptake_0hr_hit_upper,
    total_b12_response_lower_censored = total_response_lower_censored,
    filtrate_b12_response_lower_censored = filtrate_response_lower_censored,
    uptake_1hr_b12_response_lower_censored = uptake_1hr_response_lower_censored,
    uptake_0hr_b12_response_lower_censored = uptake_0hr_response_lower_censored,
    total_b12_response_upper_censored = total_response_upper_censored,
    filtrate_b12_response_upper_censored = filtrate_response_upper_censored,
    uptake_1hr_b12_response_upper_censored = uptake_1hr_response_upper_censored,
    uptake_0hr_b12_response_upper_censored = uptake_0hr_response_upper_censored,
    total_b12_response_any_censored = total_response_any_censored,
    filtrate_b12_response_any_censored = filtrate_response_any_censored,
    uptake_1hr_b12_response_any_censored = uptake_1hr_response_any_censored,
    uptake_0hr_b12_response_any_censored = uptake_0hr_response_any_censored
  )

blank_od_values <- od_viab_summary |>
  dplyr::left_join(
    culture_samples |>
      dplyr::select(sample_id, strain_id),
    by = "sample_id"
  ) |>
  dplyr::filter(strain_id == "blank", is.finite(raw_od600))

blank_od_n <- nrow(blank_od_values)
blank_od_mean <- mean(blank_od_values$raw_od600, na.rm = TRUE)
blank_od_sd <- stats::sd(blank_od_values$raw_od600, na.rm = TRUE)
blank_od_filter_threshold <- blank_od_mean + blank_od_filter_sd_multiplier * blank_od_sd

blank_dead_stain_summary <- od_viab_summary |>
  dplyr::left_join(
    culture_samples |>
      dplyr::select(sample_id, strain_id, experiment_name, sample_time_hour),
    by = "sample_id"
  ) |>
  dplyr::filter(strain_id == "blank", is.finite(dead_stain)) |>
  dplyr::summarise(
    blank_dead_stain_n = dplyr::n(),
    blank_dead_stain_mean = mean(dead_stain, na.rm = TRUE),
    blank_dead_stain_sd = stats::sd(dead_stain, na.rm = TRUE),
    .groups = "drop"
  )

blank_live_stain_floor <- od_viab_summary |>
  dplyr::left_join(
    culture_samples |>
      dplyr::select(sample_id, strain_id, experiment_name, sample_time_hour),
    by = "sample_id"
  ) |>
  dplyr::filter(strain_id == "blank", is.finite(live_stain)) |>
  dplyr::summarise(
    blank_live_stain_n = dplyr::n(),
    blank_live_stain_mean = mean(live_stain, na.rm = TRUE),
    blank_live_stain_sd = stats::sd(live_stain, na.rm = TRUE),
    blank_live_stain_floor =
      blank_live_stain_mean + blank_live_stain_floor_sd_multiplier * blank_live_stain_sd
  )

make_biological_replicate <- function(experiment_name) {
  dplyr::dense_rank(experiment_name)
}

b12_trait_survey_results <- culture_samples |>
  dplyr::left_join(od_viab_summary, by = "sample_id") |>
  dplyr::left_join(b12_summary_wide, by = "sample_id") |>
  dplyr::cross_join(
    blank_dead_stain_summary |>
      dplyr::select(blank_dead_stain_mean)
  ) |>
  dplyr::cross_join(
    blank_live_stain_floor |>
      dplyr::select(blank_live_stain_mean, blank_live_stain_floor)
  ) |>
  dplyr::mutate(
    biological_replicate = make_biological_replicate(experiment_name),
    od600 = raw_od600 - blank_od_mean,
    adjusted_od600 = od600,
    floored_live_stain = dplyr::if_else(
      is.na(live_stain),
      NA_real_,
      pmax(live_stain, blank_live_stain_floor, na.rm = TRUE)
    ),
    live_stain_below_floor = live_stain < floored_live_stain,
    adjusted_live_stain = floored_live_stain - blank_live_stain_mean,
    adjusted_dead_stain = dead_stain - blank_dead_stain_mean,
    raw_dead_proportion = dplyr::if_else(dead_stain > 0, live_stain / dead_stain, NA_real_),
    dead_proportion_floored_live = dplyr::if_else(
      adjusted_dead_stain > 0,
      adjusted_live_stain / adjusted_dead_stain,
      NA_real_
    ),
    dead_proportion = dead_proportion_floored_live,
    filtrate_proportion = dplyr::if_else(total_b12 > 0, filtrate_b12 / total_b12, NA_real_),
    uptake_0hr_expected_b12 = filtrate_b12 + added_b12,
    uptake_ok = is.finite(uptake_0hr_b12) &
      is.finite(uptake_1hr_b12) &
      is.finite(filtrate_b12) &
      abs(uptake_0hr_b12 - uptake_0hr_expected_b12) < uptake_0hr_error_thresh,
    uptake_baseline_b12 = dplyr::case_when(
      b12_uptake_baseline == "uptake_0hr" ~ uptake_0hr_b12,
      b12_uptake_baseline == "filtrate_plus_added" ~ uptake_0hr_expected_b12,
      b12_uptake_baseline == "mean_of_both" ~ rowMeans(
        dplyr::pick(uptake_0hr_b12, uptake_0hr_expected_b12),
        na.rm = FALSE
      ),
      TRUE ~ NA_real_
    ),
    uptake_b12 = dplyr::if_else(uptake_ok, uptake_baseline_b12 - uptake_1hr_b12, NA_real_),
    total_b12_per_od = dplyr::if_else(adjusted_od600 > 0, total_b12 / adjusted_od600, NA_real_),
    filtrate_b12_per_od = dplyr::if_else(adjusted_od600 > 0, filtrate_b12 / adjusted_od600, NA_real_),
    uptake_per_od = dplyr::if_else(uptake_ok & adjusted_od600 > 0, uptake_b12 / adjusted_od600, NA_real_),
    passed_blank_od_filter = strain_id == "blank" | raw_od600 > blank_od_filter_threshold
  ) |>
  dplyr::select(
    isolation_source,
    strain_id,
    array_plate,
    culture_plate_number,
    culture_well,
    layout_row_id,
    experiment_name,
    biological_replicate,
    sample_time_hour,
    sample_id,
    raw_od600,
    od600,
    adjusted_od600,
    dead_stain,
    live_stain,
    blank_live_stain_mean,
    blank_live_stain_floor,
    floored_live_stain,
    live_stain_below_floor,
    adjusted_live_stain,
    blank_dead_stain_mean,
    adjusted_dead_stain,
    total_b12,
    filtrate_b12,
    uptake_0hr_b12,
    uptake_1hr_b12,
    uptake_0hr_expected_b12,
    uptake_baseline_b12,
    uptake_ok,
    uptake_b12,
    total_b12_per_od,
    filtrate_b12_per_od,
    uptake_per_od,
    raw_dead_proportion,
    dead_proportion,
    dead_proportion_floored_live,
    filtrate_proportion,
    total_b12_hit_lower,
    filtrate_b12_hit_lower,
    uptake_0hr_b12_hit_lower,
    uptake_1hr_b12_hit_lower,
    total_b12_hit_upper,
    filtrate_b12_hit_upper,
    uptake_0hr_b12_hit_upper,
    uptake_1hr_b12_hit_upper,
    total_b12_response_lower_censored,
    filtrate_b12_response_lower_censored,
    uptake_0hr_b12_response_lower_censored,
    uptake_1hr_b12_response_lower_censored,
    total_b12_response_upper_censored,
    filtrate_b12_response_upper_censored,
    uptake_0hr_b12_response_upper_censored,
    uptake_1hr_b12_response_upper_censored,
    total_b12_response_any_censored,
    filtrate_b12_response_any_censored,
    uptake_0hr_b12_response_any_censored,
    uptake_1hr_b12_response_any_censored,
    passed_blank_od_filter
  ) |>
  dplyr::arrange(experiment_name, sample_time_hour, culture_plate_number, culture_well)

b12_trait_survey_results_filtered <- b12_trait_survey_results |>
  dplyr::filter(passed_blank_od_filter)

geo_mean_positive <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) {
    return(NA_real_)
  }
  exp(mean(log(x)))
}

geo_sd_positive <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 2) {
    return(NA_real_)
  }
  exp(stats::sd(log(x)))
}

prop_true <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x)
}

classify_b12_trait_group1 <- function(total_b12, uptake_b12) {
  dplyr::case_when(
    !is.finite(total_b12) | !is.finite(uptake_b12) ~ NA_character_,
    total_b12 < b12_trait_group_total_b12_threshold_pm &
      uptake_b12 < b12_trait_group_uptake_b12_threshold_pm ~ "bystander",
    total_b12 < b12_trait_group_total_b12_threshold_pm &
      uptake_b12 >= b12_trait_group_uptake_b12_threshold_pm ~ "scavenger",
    total_b12 >= b12_trait_group_total_b12_threshold_pm &
      uptake_b12 >= b12_trait_group_uptake_b12_threshold_pm ~ "reclaimer",
    total_b12 >= b12_trait_group_total_b12_threshold_pm &
      uptake_b12 < b12_trait_group_uptake_b12_threshold_pm ~ "sharer",
    TRUE ~ NA_character_
  )
}

classify_b12_trait_group2 <- function(total_b12, filtrate_proportion, uptake_b12) {
  dplyr::case_when(
    !is.finite(total_b12) | total_b12 < b12_trait_group_total_b12_threshold_pm ~ NA_character_,
    !is.finite(filtrate_proportion) | !is.finite(uptake_b12) ~ NA_character_,
    uptake_b12 >= b12_trait_group_uptake_b12_threshold_pm ~ "reclaimer",
    filtrate_proportion < b12_trait_group_filtrate_proportion_threshold &
      uptake_b12 < b12_trait_group_uptake_b12_threshold_pm ~ "retainer",
    filtrate_proportion >= b12_trait_group_filtrate_proportion_threshold &
      uptake_b12 < b12_trait_group_uptake_b12_threshold_pm ~ "provider",
    TRUE ~ NA_character_
  )
}

raw_measurements <- c(
  "od600",
  "adjusted_od600",
  "live_stain",
  "floored_live_stain",
  "adjusted_live_stain",
  "dead_stain",
  "adjusted_dead_stain",
  "total_b12",
  "filtrate_b12",
  "uptake_0hr_b12",
  "uptake_1hr_b12"
)

geometric_summary_vars <- c(
  raw_measurements,
  "total_b12_per_od",
  "filtrate_b12_per_od",
  "raw_dead_proportion",
  "dead_proportion",
  "dead_proportion_floored_live",
  "filtrate_proportion"
)

arithmetic_summary_vars <- c(
  "uptake_b12",
  "uptake_per_od",
  "raw_dead_proportion",
  "dead_proportion",
  "dead_proportion_floored_live",
  "filtrate_proportion"
)

time_group_cols <- c(
  "isolation_source",
  "strain_id",
  "sample_time_hour"
)

strain_group_cols <- c(
  "isolation_source",
  "strain_id"
)

b12_trait_survey_time_summary <- b12_trait_survey_results_filtered |>
  dplyr::filter(strain_id != "blank") |>
  dplyr::group_by(dplyr::across(dplyr::all_of(time_group_cols))) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(geometric_summary_vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
    dplyr::across(dplyr::all_of(geometric_summary_vars), geo_mean_positive, .names = "{.col}_gm"),
    dplyr::across(dplyr::all_of(geometric_summary_vars), geo_sd_positive, .names = "{.col}_gsd"),
    dplyr::across(dplyr::all_of(arithmetic_summary_vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
    dplyr::across(dplyr::all_of(arithmetic_summary_vars), ~ stats::sd(.x, na.rm = TRUE), .names = "{.col}_sd"),
    total_b12_prop_hit_lower = prop_true(total_b12_hit_lower),
    filtrate_b12_prop_hit_lower = prop_true(filtrate_b12_hit_lower),
    uptake_0hr_b12_prop_hit_lower = prop_true(uptake_0hr_b12_hit_lower),
    uptake_1hr_b12_prop_hit_lower = prop_true(uptake_1hr_b12_hit_lower),
    total_b12_prop_hit_upper = prop_true(total_b12_hit_upper),
    filtrate_b12_prop_hit_upper = prop_true(filtrate_b12_hit_upper),
    uptake_0hr_b12_prop_hit_upper = prop_true(uptake_0hr_b12_hit_upper),
    uptake_1hr_b12_prop_hit_upper = prop_true(uptake_1hr_b12_hit_upper),
    total_b12_prop_response_censored = prop_true(total_b12_response_any_censored),
    filtrate_b12_prop_response_censored = prop_true(filtrate_b12_response_any_censored),
    uptake_0hr_b12_prop_response_censored = prop_true(uptake_0hr_b12_response_any_censored),
    uptake_1hr_b12_prop_response_censored = prop_true(uptake_1hr_b12_response_any_censored),
    n_samples = dplyr::n(),
    n_biological_replicates = dplyr::n_distinct(biological_replicate),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    b12_trait_group1 = classify_b12_trait_group1(total_b12_gm, uptake_b12_mean),
    b12_trait_group2 = classify_b12_trait_group2(total_b12_gm, filtrate_proportion_gm, uptake_b12_mean)
  )

b12_trait_survey_strain_summary <- b12_trait_survey_results_filtered |>
  dplyr::filter(strain_id != "blank") |>
  dplyr::group_by(dplyr::across(dplyr::all_of(strain_group_cols))) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(geometric_summary_vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
    dplyr::across(dplyr::all_of(geometric_summary_vars), geo_mean_positive, .names = "{.col}_gm"),
    dplyr::across(dplyr::all_of(geometric_summary_vars), geo_sd_positive, .names = "{.col}_gsd"),
    dplyr::across(dplyr::all_of(arithmetic_summary_vars), ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
    dplyr::across(dplyr::all_of(arithmetic_summary_vars), ~ stats::sd(.x, na.rm = TRUE), .names = "{.col}_sd"),
    total_b12_prop_hit_lower = prop_true(total_b12_hit_lower),
    filtrate_b12_prop_hit_lower = prop_true(filtrate_b12_hit_lower),
    uptake_0hr_b12_prop_hit_lower = prop_true(uptake_0hr_b12_hit_lower),
    uptake_1hr_b12_prop_hit_lower = prop_true(uptake_1hr_b12_hit_lower),
    total_b12_prop_hit_upper = prop_true(total_b12_hit_upper),
    filtrate_b12_prop_hit_upper = prop_true(filtrate_b12_hit_upper),
    uptake_0hr_b12_prop_hit_upper = prop_true(uptake_0hr_b12_hit_upper),
    uptake_1hr_b12_prop_hit_upper = prop_true(uptake_1hr_b12_hit_upper),
    total_b12_prop_response_censored = prop_true(total_b12_response_any_censored),
    filtrate_b12_prop_response_censored = prop_true(filtrate_b12_response_any_censored),
    uptake_0hr_b12_prop_response_censored = prop_true(uptake_0hr_b12_response_any_censored),
    uptake_1hr_b12_prop_response_censored = prop_true(uptake_1hr_b12_response_any_censored),
    n_samples = dplyr::n(),
    n_biological_replicates = dplyr::n_distinct(biological_replicate),
    n_sample_times = dplyr::n_distinct(sample_time_hour),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    b12_trait_group1 = classify_b12_trait_group1(total_b12_gm, uptake_b12_mean),
    b12_trait_group2 = classify_b12_trait_group2(total_b12_gm, filtrate_proportion_gm, uptake_b12_mean)
  )

b12_trait_survey_gwas_phenotypes <- b12_trait_survey_strain_summary |>
  dplyr::mutate(
    log10_total_b12_gm = dplyr::if_else(
      total_b12_gm > 0,
      log10(total_b12_gm),
      NA_real_
    ),
    log10_total_b12_per_od_gm = dplyr::if_else(
      total_b12_per_od_gm > 0,
      log10(total_b12_per_od_gm),
      NA_real_
    )
  ) |>
  dplyr::select(
    strain_id,
    total_b12_gm,
    log10_total_b12_gm,
    total_b12_gsd,
    total_b12_per_od_gm,
    log10_total_b12_per_od_gm,
    total_b12_per_od_gsd,
    total_b12_prop_hit_lower,
    uptake_b12_mean,
    uptake_b12_sd,
    filtrate_proportion_gm,
    dead_proportion_gm,
    n_samples,
    n_biological_replicates,
    n_sample_times
  )

diagnostic_to_long <- function(dat, diagnostic_type) {
  dat |>
    dplyr::mutate(diagnostic_type = diagnostic_type) |>
    dplyr::relocate(diagnostic_type)
}

b12_trait_survey_estimator_diagnostics <- dplyr::bind_rows(
  diagnostic_to_long(estimator_result$edge_correction_summary, "edge_correction_summary"),
  diagnostic_to_long(estimator_result$adjacent_pair_diagnostics$pair_summary, "adjacent_row_pair_summary"),
  diagnostic_to_long(estimator_result$adjacent_pair_diagnostics$block_summary, "adjacent_row_pair_block_summary"),
  diagnostic_to_long(estimator_result$adjacent_pair_diagnostics$quadratic_model_tests, "adjacent_row_pair_quadratic_model_tests")
)

make_filtrate_fraction_dead_uptake_log10_lm_table <- function(
  dat,
  data_level,
  input_file,
  filtrate_fraction_col,
  fraction_dead_col,
  uptake_col
) {
  model_data <- dat |>
    dplyr::filter(
      is.finite(.data[[filtrate_fraction_col]]),
      .data[[filtrate_fraction_col]] > 0,
      is.finite(.data[[fraction_dead_col]]),
      .data[[fraction_dead_col]] > 0,
      is.finite(.data[[uptake_col]])
    ) |>
    dplyr::mutate(
      log10_filtrate_fraction = log10(.data[[filtrate_fraction_col]]),
      log10_fraction_dead = log10(.data[[fraction_dead_col]]),
      uptake_b12_linear = .data[[uptake_col]]
    )

  model <- stats::lm(
    log10_filtrate_fraction ~ log10_fraction_dead * uptake_b12_linear,
    data = model_data
  )
  model_summary <- summary(model)
  coefficient_table <- as.data.frame(coef(model_summary), check.names = FALSE)
  coefficient_table$term <- row.names(coefficient_table)
  row.names(coefficient_table) <- NULL

  coefficient_results <- coefficient_table |>
    dplyr::transmute(
      data_level = data_level,
      input_file = input_file,
      model_formula = "log10(filtrate_fraction) ~ log10(fraction_dead) * uptake_b12",
      n = stats::nobs(model),
      adjusted_r_squared = unname(model_summary$adj.r.squared),
      result_type = "coefficient",
      backward_step = NA_character_,
      term = term,
      estimate = .data[["Estimate"]],
      std_error = .data[["Std. Error"]],
      statistic = .data[["t value"]],
      p_value = .data[["Pr(>|t|)"]],
      decision = NA_character_
    )

  interaction_test <- as.data.frame(stats::drop1(model, test = "F"), check.names = FALSE)
  interaction_test$term <- row.names(interaction_test)
  row.names(interaction_test) <- NULL
  interaction_term <- "log10_fraction_dead:uptake_b12_linear"
  interaction_test <- interaction_test |>
    dplyr::filter(.data$term == interaction_term) |>
    dplyr::transmute(
      data_level = data_level,
      input_file = input_file,
      model_formula = "log10(filtrate_fraction) ~ log10(fraction_dead) * uptake_b12",
      n = stats::nobs(model),
      adjusted_r_squared = unname(model_summary$adj.r.squared),
      result_type = "backward_elimination",
      backward_step = "test_interaction",
      term = term,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = .data[["F value"]],
      p_value = .data[["Pr(>F)"]],
      decision = dplyr::if_else(
        .data[["Pr(>F)"]] < 0.05,
        "retain_interaction_stop_hierarchical_elimination",
        "remove_interaction_continue_to_main_effects"
      )
    )

  dplyr::bind_rows(coefficient_results, interaction_test)
}

filtrate_fraction_dead_uptake_log10_model_results <- dplyr::bind_rows(
  make_filtrate_fraction_dead_uptake_log10_lm_table(
    b12_trait_survey_results_filtered,
    data_level = "original_filtered_samples",
    input_file = b12_trait_survey_results_filtered_path,
    filtrate_fraction_col = "filtrate_proportion",
    fraction_dead_col = "dead_proportion",
    uptake_col = "uptake_b12"
  ),
  make_filtrate_fraction_dead_uptake_log10_lm_table(
    b12_trait_survey_time_summary,
    data_level = "strain_id_sample_time_hour_summary",
    input_file = b12_trait_survey_time_summary_path,
    filtrate_fraction_col = "filtrate_proportion_gm",
    fraction_dead_col = "dead_proportion_gm",
    uptake_col = "uptake_b12_mean"
  ),
  make_filtrate_fraction_dead_uptake_log10_lm_table(
    b12_trait_survey_strain_summary,
    data_level = "strain_id_summary",
    input_file = b12_trait_survey_strain_summary_path,
    filtrate_fraction_col = "filtrate_proportion_gm",
    fraction_dead_col = "dead_proportion_gm",
    uptake_col = "uptake_b12_mean"
  )
)

analysis_parameters <- dplyr::bind_rows(
  tibble::tibble(
    parameter = c(
      "standard_sample_types",
      "standard_cols_to_remove",
      "add_offset_pm",
      "edge_correction",
      "upper_thresh_logb12_estimation_se",
      "added_b12",
      "b12_uptake_baseline",
      "b12_uptake_baseline_formula_uptake_0hr",
      "b12_uptake_baseline_formula_filtrate_plus_added",
      "b12_uptake_baseline_formula_mean_of_both",
      "uptake_0hr_error_thresh",
      "uptake_0hr_error_thresh_definition",
      "uptake_ok_formula",
      "b12_trait_group_total_b12_threshold_pm",
      "b12_trait_group_uptake_b12_threshold_pm",
      "b12_trait_group_filtrate_proportion_threshold",
      "b12_trait_group1_formula",
      "b12_trait_group2_formula",
      "blank_od_n",
      "blank_od_mean",
      "blank_od_sd",
      "blank_od_filter_sd_multiplier",
      "blank_od_filter_threshold",
      "blank_od_filter_formula",
      "blank_live_stain_floor_sd_multiplier",
      "blank_live_stain_floor_grouping",
      "blank_live_stain_floor_n",
      "blank_live_stain_mean",
      "blank_live_stain_floor",
      "blank_dead_stain_n",
      "blank_dead_stain_mean",
      "blank_dead_stain_sd",
      "stain_adjustment_formula",
      "dead_proportion_formula"
    ),
    value = c(
      paste(standard_sample_types, collapse = ";"),
      paste(standard_cols_to_remove, collapse = ";"),
      add_offset_pm,
      edge_correction,
      upper_thresh_logb12_estimation_se,
      added_b12,
      b12_uptake_baseline,
      "uptake_b12 = uptake_0hr_b12 - uptake_1hr_b12",
      "uptake_b12 = filtrate_b12 + added_b12 - uptake_1hr_b12",
      "uptake_b12 = mean(uptake_0hr_b12, filtrate_b12 + added_b12) - uptake_1hr_b12",
      uptake_0hr_error_thresh,
      "Maximum allowed absolute difference, in pM, between direct uptake_0hr_b12 and physically expected filtrate_b12 + added_b12",
      "finite uptake_0hr_b12, uptake_1hr_b12, and filtrate_b12; abs(uptake_0hr_b12 - (filtrate_b12 + added_b12)) < uptake_0hr_error_thresh",
      b12_trait_group_total_b12_threshold_pm,
      b12_trait_group_uptake_b12_threshold_pm,
      b12_trait_group_filtrate_proportion_threshold,
      "Uses total_b12_gm and uptake_b12_mean: bystander total < threshold and uptake < threshold; scavenger total < threshold and uptake >= threshold; reclaimer total >= threshold and uptake >= threshold; sharer total >= threshold and uptake < threshold",
      "Uses total_b12_gm, filtrate_proportion_gm, and uptake_b12_mean: NA if total < threshold; reclaimer uptake >= threshold; retainer filtrate proportion < threshold and uptake < threshold; provider filtrate proportion >= threshold and uptake < threshold",
      blank_od_n,
      blank_od_mean,
      blank_od_sd,
      blank_od_filter_sd_multiplier,
      blank_od_filter_threshold,
      "strain_id == blank OR raw_od600 > blank_od_mean + blank_od_filter_sd_multiplier * blank_od_sd",
      blank_live_stain_floor_sd_multiplier,
      "pooled across experiment_name and sample_time_hour",
      blank_live_stain_floor$blank_live_stain_n,
      blank_live_stain_floor$blank_live_stain_mean,
      blank_live_stain_floor$blank_live_stain_floor,
      blank_dead_stain_summary$blank_dead_stain_n,
      blank_dead_stain_summary$blank_dead_stain_mean,
      blank_dead_stain_summary$blank_dead_stain_sd,
      "adjusted_live_stain = floored_live_stain - pooled blank_live_stain_mean; adjusted_dead_stain = dead_stain - pooled blank_dead_stain_mean",
      "dead_proportion and dead_proportion_floored_live = adjusted_live_stain / adjusted_dead_stain when adjusted_dead_stain > 0"
    )
  ),
  estimator_result$estimator_parameters |>
    dplyr::mutate(parameter = paste0("estimator_", parameter))
)

utils::write.csv(analysis_parameters, analysis_parameters_path, row.names = FALSE, na = "")
utils::write.csv(estimator_result$estimator_parameters, b12_estimator_parameters_path, row.names = FALSE, na = "")
utils::write.csv(b12_estimator_output, b12_estimator_output_path, row.names = FALSE, na = "")
utils::write.csv(b12_trait_survey_results, b12_trait_survey_results_path, row.names = FALSE, na = "")
utils::write.csv(b12_trait_survey_results_filtered, b12_trait_survey_results_filtered_path, row.names = FALSE, na = "")
utils::write.csv(b12_trait_survey_time_summary, b12_trait_survey_time_summary_path, row.names = FALSE, na = "")
utils::write.csv(b12_trait_survey_strain_summary, b12_trait_survey_strain_summary_path, row.names = FALSE, na = "")
utils::write.csv(b12_trait_survey_gwas_phenotypes, b12_trait_survey_gwas_phenotypes_path, row.names = FALSE, na = "")
utils::write.csv(
  filtrate_fraction_dead_uptake_log10_model_results,
  filtrate_fraction_dead_uptake_log10_model_results_path,
  row.names = FALSE,
  na = ""
)
utils::write.csv(b12_trait_survey_estimator_diagnostics, b12_trait_survey_estimator_diagnostics_path, row.names = FALSE, na = "")

# Compact source tables consumed by Figures 3 and S10. These are deterministic
# views of the strain summary; the former exploratory figure script is not part
# of the public pipeline.
trait_group_levels <- c("retainer", "provider", "reclaimer")
trait_figure_strains <- b12_trait_survey_strain_summary |>
  dplyr::filter(.data$n_samples >= 3L) |>
  dplyr::mutate(
    b12_trait_group2 = factor(.data$b12_trait_group2, levels = trait_group_levels)
  )

fraction_dead_extracellular <- trait_figure_strains |>
  dplyr::filter(
    is.finite(.data$total_b12_gm), .data$total_b12_gm > 50,
    is.finite(.data$dead_proportion_gm), .data$dead_proportion_gm > 0,
    is.finite(.data$dead_proportion_gsd), .data$dead_proportion_gsd > 0,
    is.finite(.data$filtrate_proportion_gm), .data$filtrate_proportion_gm > 0,
    is.finite(.data$filtrate_proportion_gsd), .data$filtrate_proportion_gsd > 0,
    is.finite(.data$uptake_b12_mean)
  )

extracellular_over_dead_uptake <- fraction_dead_extracellular |>
  dplyr::mutate(
    extracellular_over_dead = .data$filtrate_proportion_gm / .data$dead_proportion_gm,
    log10_extracellular_over_dead = log10(.data$extracellular_over_dead),
    extracellular_over_dead_gsd = exp(sqrt(
      log(.data$filtrate_proportion_gsd)^2 + log(.data$dead_proportion_gsd)^2
    ))
  ) |>
  dplyr::filter(
    is.finite(.data$extracellular_over_dead),
    .data$extracellular_over_dead > 0,
    is.finite(.data$log10_extracellular_over_dead)
  )

trait_group2_boxplot_source <- trait_figure_strains |>
  dplyr::filter(!is.na(.data$b12_trait_group2)) |>
  dplyr::mutate(
    extracellular_b12_over_dead_fraction =
      .data$filtrate_proportion_gm / .data$dead_proportion_gm
  ) |>
  dplyr::filter(
    is.finite(.data$extracellular_b12_over_dead_fraction),
    .data$extracellular_b12_over_dead_fraction > 0
  )

utils::write.csv(
  fraction_dead_extracellular,
  file.path(
    figure_source_dir,
    "15B_dead_cell_fraction_vs_ec_b12_fraction_trait_group2_source.csv"
  ),
  row.names = FALSE, na = ""
)
utils::write.csv(
  extracellular_over_dead_uptake,
  file.path(
    figure_source_dir,
    "15C_ec_b12_fraction_over_dead_cell_fraction_vs_uptake_trait_group2_source.csv"
  ),
  row.names = FALSE, na = ""
)
utils::write.csv(
  trait_group2_boxplot_source,
  file.path(figure_source_dir, "12D_trait_group2_boxplot_source.csv"),
  row.names = FALSE, na = ""
)

message("Wrote analysis parameters to: ", analysis_parameters_path)
message("Wrote B12 estimator output to: ", b12_estimator_output_path)
message("Wrote unfiltered trait survey results to: ", b12_trait_survey_results_path)
message("Wrote OD-filtered trait survey results to: ", b12_trait_survey_results_filtered_path)
message("Wrote time summary to: ", b12_trait_survey_time_summary_path)
message("Wrote strain summary to: ", b12_trait_survey_strain_summary_path)
message("Wrote GWAS phenotype table to: ", b12_trait_survey_gwas_phenotypes_path)
message("Wrote filtrate-fraction/death/uptake model results to: ", filtrate_fraction_dead_uptake_log10_model_results_path)
message("Wrote estimator diagnostics to: ", b12_trait_survey_estimator_diagnostics_path)
message("Wrote compiled-figure source tables to: ", figure_source_dir)
message("Unfiltered rows: ", nrow(b12_trait_survey_results))
message("Filtered rows: ", nrow(b12_trait_survey_results_filtered))
