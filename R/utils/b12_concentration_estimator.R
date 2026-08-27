# R/b12_concentration_estimator.R

b12_concentration_estimator <- function(
    b12_assay_combined,
    standard_sample_types = c("filtrate"),
    standard_cols_to_remove = c(1),
    add_offset_pm = 0.1,
    edge_correction = c("none", "quadratic_edge_pairs"),
    edge_outer_inner_rows = c(A = "B", P = "O"),
    min_edge_correction_pairs = 8,
    confidence_level = 0.95,
    hit_upper_tolerance = 1e-6,
    response_bound_fraction = 0.05,
    verbose = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(stringr)
    library(drc)
  })

  edge_correction <- match.arg(edge_correction)
  ci_multiplier <- stats::qnorm(1 - (1 - confidence_level) / 2)

  if (!is.numeric(response_bound_fraction) ||
      length(response_bound_fraction) != 1 ||
      !is.finite(response_bound_fraction) ||
      response_bound_fraction < 0 ||
      response_bound_fraction >= 0.5) {
    stop("response_bound_fraction must be a finite numeric value >= 0 and < 0.5.")
  }

  f4pl <- function(x, a, d, c, b) d + (a - d) / (1 + exp(b * (x - c)))
  f4pl_prime <- function(x, a, d, c, b) {
    z <- exp(b * (x - c))
    -(a - d) * b * z / (1 + z)^2
  }

  req_cols <- c(
    "assay_batch", "experiment_name", "assay_type", "assay_well", "assay_plate_number",
    "measure_time_hours", "fluorescence",
    "sample_volume_added", "assay_volume_added", "storage_dilution_factor",
    "sample_type", "sample_id", "sample_time_hour", "strain_id",
    "sample_b12_conc_pm", "standard_bacteria_strain_id", "assay_dilution"
  )
  missing <- setdiff(req_cols, names(b12_assay_combined))
  if (length(missing) > 0) {
    stop("b12_assay_combined is missing required columns: ", paste(missing, collapse = ", "))
  }

  df <- dplyr::as_tibble(b12_assay_combined)

  if (!("standard_processing" %in% names(df))) {
    df$standard_processing <- NA_character_
  }

  if (!("log_fluorescence" %in% names(df))) {
    df$log_fluorescence <- log(df$fluorescence)
  } else {
    df$log_fluorescence <- ifelse(
      is.finite(df$log_fluorescence),
      df$log_fluorescence,
      log(df$fluorescence)
    )
  }

  df$assay_dilution <- ifelse(
    is.finite(df$assay_dilution),
    df$assay_dilution,
    ((df$sample_volume_added + df$assay_volume_added) / df$sample_volume_added) *
      df$storage_dilution_factor
  )

  if (sum(df$assay_type == "standard", na.rm = TRUE) == 0) {
    stop("No standards found (assay_type=='standard') in b12_assay_combined.")
  }

  df <- df %>%
    dplyr::mutate(
      assay_row = stringr::str_extract(assay_well, "^[A-Z]+"),
      assay_col = as.integer(stringr::str_extract(assay_well, "\\d+")),
      fluorescence_raw = fluorescence,
      log_fluorescence_raw = log_fluorescence,
      edge_correction_applied = FALSE,
      edge_correction_pair = NA_character_,
      adjusted_b12_conc_pm = dplyr::if_else(
        assay_type == "standard",
        sample_b12_conc_pm / assay_dilution + add_offset_pm,
        NA_real_
      ),
      log_b12_fm = dplyr::if_else(
        assay_type == "standard" &
          is.finite(adjusted_b12_conc_pm) &
          adjusted_b12_conc_pm > 0,
        log(adjusted_b12_conc_pm * 1000),
        NA_real_
      )
    )

  build_adjacent_pair_diagnostics <- function(dat) {
    row_pairs <- tibble::tibble(
      row1 = LETTERS[seq(1, 15, by = 2)],
      row2 = LETTERS[seq(2, 16, by = 2)],
      row_pair = paste0(row1, "_", row2)
    )

    pair_side_1 <- dat %>%
      dplyr::filter(
        assay_type == "standard",
        is.finite(log_fluorescence_raw),
        is.finite(log_b12_fm),
        !is.na(adjusted_b12_conc_pm)
      ) %>%
      dplyr::inner_join(row_pairs, by = c("assay_row" = "row1")) %>%
      dplyr::mutate(side = "row1")

    pair_side_2 <- dat %>%
      dplyr::filter(
        assay_type == "standard",
        is.finite(log_fluorescence_raw),
        is.finite(log_b12_fm),
        !is.na(adjusted_b12_conc_pm)
      ) %>%
      dplyr::inner_join(row_pairs, by = c("assay_row" = "row2")) %>%
      dplyr::mutate(side = "row2")

    pair_long <- dplyr::bind_rows(pair_side_1, pair_side_2) %>%
      dplyr::mutate(
        standard_processing = dplyr::coalesce(standard_processing, "not_applicable"),
        group_id = paste(assay_batch, measure_time_hours, sep = "__"),
        log10_conc = log10(adjusted_b12_conc_pm)
      )

    paired <- pair_long %>%
      dplyr::group_by(
        assay_batch,
        measure_time_hours,
        standard_processing,
        group_id,
        row_pair,
        standard_bacteria_strain_id,
        sample_type,
        adjusted_b12_conc_pm,
        log10_conc,
        side
      ) %>%
      dplyr::summarise(log_fluorescence = mean(log_fluorescence_raw), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = side, values_from = log_fluorescence) %>%
      dplyr::filter(is.finite(row1), is.finite(row2)) %>%
      dplyr::mutate(
        log_diff_row2_minus_row1 = row2 - row1,
        fold_row2_over_row1 = exp(log_diff_row2_minus_row1)
      )

    pair_summary <- paired %>%
      dplyr::group_by(row_pair) %>%
      dplyr::summarise(
        n_pairs = dplyr::n(),
        mean_log_diff = mean(log_diff_row2_minus_row1, na.rm = TRUE),
        mean_fold_row2_over_row1 = exp(mean_log_diff),
        pct_row2_vs_row1 = (mean_fold_row2_over_row1 - 1) * 100,
        .groups = "drop"
      )

    block_summary <- paired %>%
      dplyr::group_by(group_id, assay_batch, measure_time_hours, row_pair) %>%
      dplyr::summarise(log_diff = mean(log_diff_row2_minus_row1, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(row_pair) %>%
      dplyr::summarise(
        n_blocks = dplyr::n(),
        mean_log_diff = mean(log_diff, na.rm = TRUE),
        mean_fold_row2_over_row1 = exp(mean_log_diff),
        pct_row2_vs_row1 = (mean_fold_row2_over_row1 - 1) * 100,
        t_statistic = if (dplyr::n() > 1) unname(stats::t.test(log_diff)$statistic) else NA_real_,
        p_value = if (dplyr::n() > 1) stats::t.test(log_diff)$p.value else NA_real_,
        ci_low_pct = if (dplyr::n() > 1) (exp(stats::t.test(log_diff)$conf.int[1]) - 1) * 100 else NA_real_,
        ci_high_pct = if (dplyr::n() > 1) (exp(stats::t.test(log_diff)$conf.int[2]) - 1) * 100 else NA_real_,
        .groups = "drop"
      )

    model_tests <- tibble::tibble(
      comparison = character(),
      df = numeric(),
      sum_of_squares = numeric(),
      f_statistic = numeric(),
      p_value = numeric()
    )

    if (nrow(paired) > 0 && dplyr::n_distinct(paired$row_pair) > 1) {
      model_data <- paired %>%
        dplyr::mutate(
          row_pair = factor(row_pair),
          group_id = factor(group_id),
          z_log10_conc = as.numeric(scale(log10_conc, center = TRUE, scale = FALSE))
        )

      m_group_only <- stats::lm(log_diff_row2_minus_row1 ~ group_id, data = model_data)
      m_const <- stats::lm(log_diff_row2_minus_row1 ~ row_pair + group_id, data = model_data)
      m_quad_common <- stats::lm(
        log_diff_row2_minus_row1 ~ row_pair + z_log10_conc + I(z_log10_conc^2) + group_id,
        data = model_data
      )
      m_quad_pair <- stats::lm(
        log_diff_row2_minus_row1 ~ row_pair * (z_log10_conc + I(z_log10_conc^2)) + group_id,
        data = model_data
      )
      a <- as.data.frame(stats::anova(m_group_only, m_const, m_quad_common, m_quad_pair))
      model_tests <- tibble::tibble(
        comparison = c(
          "group_id_only_baseline",
          "row_pair_offsets",
          "shared_quadratic_concentration_effect",
          "row_pair_specific_quadratic_concentration_effect"
        ),
        df = a$Df,
        sum_of_squares = a$`Sum of Sq`,
        f_statistic = a$F,
        p_value = a$`Pr(>F)`
      )
    }

    list(
      paired = paired,
      pair_summary = pair_summary,
      block_summary = block_summary,
      quadratic_model_tests = model_tests
    )
  }

  fit_edge_correction <- function(dat) {
    edge_pairs <- tibble::tibble(
      outer_row = names(edge_outer_inner_rows),
      inner_row = unname(edge_outer_inner_rows),
      edge_pair = paste0(outer_row, "_to_", inner_row)
    )

    outer <- dat %>%
      dplyr::filter(
        assay_type == "standard",
        assay_row %in% edge_pairs$outer_row,
        is.finite(log_fluorescence_raw),
        is.finite(adjusted_b12_conc_pm)
      ) %>%
      dplyr::inner_join(edge_pairs, by = c("assay_row" = "outer_row")) %>%
      dplyr::mutate(side = "outer")

    inner <- dat %>%
      dplyr::filter(
        assay_type == "standard",
        assay_row %in% edge_pairs$inner_row,
        is.finite(log_fluorescence_raw),
        is.finite(adjusted_b12_conc_pm)
      ) %>%
      dplyr::inner_join(edge_pairs, by = c("assay_row" = "inner_row")) %>%
      dplyr::mutate(side = "inner")

    paired <- dplyr::bind_rows(outer, inner) %>%
      dplyr::mutate(standard_processing = dplyr::coalesce(standard_processing, "not_applicable")) %>%
      dplyr::group_by(
        assay_batch,
        measure_time_hours,
        standard_processing,
        edge_pair,
        standard_bacteria_strain_id,
        sample_type,
        adjusted_b12_conc_pm,
        side
      ) %>%
      dplyr::summarise(log_fluorescence = mean(log_fluorescence_raw), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = side, values_from = log_fluorescence) %>%
      dplyr::filter(is.finite(outer), is.finite(inner)) %>%
      dplyr::mutate(
        z_log_outer = as.numeric(scale(outer, center = TRUE, scale = FALSE)),
        log_diff_inner_minus_outer = inner - outer,
        fold_inner_over_outer = exp(log_diff_inner_minus_outer)
      )

    if (nrow(paired) < min_edge_correction_pairs || dplyr::n_distinct(paired$edge_pair) < 1) {
      stop("Not enough paired edge standards to fit quadratic edge correction.")
    }

    outer_center <- attr(scale(paired$outer, center = TRUE, scale = FALSE), "scaled:center")
    paired$z_log_outer <- paired$outer - outer_center
    paired$edge_pair <- factor(paired$edge_pair)

    model <- stats::lm(inner ~ edge_pair * (z_log_outer + I(z_log_outer^2)), data = paired)

    edge_summary <- paired %>%
      dplyr::group_by(edge_pair) %>%
      dplyr::summarise(
        n_pairs = dplyr::n(),
        mean_fold_inner_over_outer = exp(mean(log_diff_inner_minus_outer, na.rm = TRUE)),
        pct_inner_vs_outer = (mean_fold_inner_over_outer - 1) * 100,
        .groups = "drop"
      )

    list(
      model = model,
      paired = paired,
      summary = edge_summary,
      outer_center = outer_center
    )
  }

  adjacent_pair_diagnostics <- build_adjacent_pair_diagnostics(df)
  edge_correction_result <- NULL

  if (edge_correction == "quadratic_edge_pairs") {
    edge_correction_result <- fit_edge_correction(df)

    edge_pairs <- tibble::tibble(
      assay_row = names(edge_outer_inner_rows),
      edge_pair = paste0(names(edge_outer_inner_rows), "_to_", unname(edge_outer_inner_rows))
    )

    correction_index <- df %>%
      dplyr::select(assay_row) %>%
      dplyr::mutate(row_index = dplyr::row_number()) %>%
      dplyr::inner_join(edge_pairs, by = "assay_row") %>%
      dplyr::filter(is.finite(df$log_fluorescence_raw[row_index]))

    if (nrow(correction_index) > 0) {
      new_data <- tibble::tibble(
        edge_pair = factor(
          correction_index$edge_pair,
          levels = levels(edge_correction_result$paired$edge_pair)
        ),
        z_log_outer = df$log_fluorescence_raw[correction_index$row_index] -
          edge_correction_result$outer_center
      )
      corrected_log <- as.numeric(stats::predict(edge_correction_result$model, newdata = new_data))

      df$log_fluorescence[correction_index$row_index] <- corrected_log
      df$fluorescence[correction_index$row_index] <- exp(corrected_log)
      df$edge_correction_applied[correction_index$row_index] <- TRUE
      df$edge_correction_pair[correction_index$row_index] <- correction_index$edge_pair
    }
  }

  standards_used <- df %>%
    dplyr::filter(
      assay_type == "standard",
      sample_type %in% standard_sample_types,
      is.finite(log_b12_fm),
      is.finite(log_fluorescence),
      !is.na(standard_bacteria_strain_id)
    ) %>%
    dplyr::filter(!assay_col %in% standard_cols_to_remove)

  if (nrow(standards_used) == 0) {
    stop("No standards after filtering by standard_sample_types, standard_cols_to_remove, and finite logs.")
  }

  curve_library <- standards_used %>%
    dplyr::group_by(assay_batch, measure_time_hours, standard_bacteria_strain_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model = purrr::map(data, ~tryCatch(
        drc::drm(log_fluorescence ~ log_b12_fm, data = .x, fct = drc::L.4()),
        error = function(e) NULL
      ))
    ) %>%
    dplyr::filter(!purrr::map_lgl(model, is.null)) %>%
    dplyr::mutate(
      coefs = purrr::map(model, ~as.list(coef(.x))),
      b_slope = purrr::map_dbl(coefs, ~.x[["b:(Intercept)"]]),
      d_min = purrr::map_dbl(coefs, ~.x[["c:(Intercept)"]]),
      a_max = purrr::map_dbl(coefs, ~.x[["d:(Intercept)"]]),
      c_ec50 = purrr::map_dbl(coefs, ~.x[["e:(Intercept)"]]),
      sigma2 = purrr::map2_dbl(model, data, ~mean(residuals(.x)^2, na.rm = TRUE)),
      x_min_std = purrr::map_dbl(data, ~min(.x$log_b12_fm, na.rm = TRUE)),
      x_max_std = purrr::map_dbl(data, ~max(.x$log_b12_fm, na.rm = TRUE)),
      response_range = a_max - d_min,
      log_fluorescence_response_lower = d_min + response_bound_fraction * response_range,
      log_fluorescence_response_upper = d_min + (1 - response_bound_fraction) * response_range,
      valid_response_bounds = is.finite(log_fluorescence_response_lower) &
        is.finite(log_fluorescence_response_upper) &
        log_fluorescence_response_lower <= log_fluorescence_response_upper
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      assay_batch, measure_time_hours, standard_bacteria_strain_id,
      a_max, d_min, c_ec50, b_slope, sigma2, x_min_std, x_max_std,
      response_range, log_fluorescence_response_lower, log_fluorescence_response_upper,
      valid_response_bounds
    )

  if (nrow(curve_library) == 0) {
    stop("No curves successfully fit. Inspect standards or loosen filters.")
  }

  samples_base <- df %>%
    dplyr::filter(assay_type == "sample", is.finite(log_fluorescence)) %>%
    dplyr::select(-dplyr::any_of("standard_bacteria_strain_id")) %>%
    dplyr::mutate(
      obs_dil = ((sample_volume_added + assay_volume_added) / sample_volume_added) *
        storage_dilution_factor,
      log_dil = log(obs_dil)
    ) %>%
    dplyr::select(
      assay_batch, experiment_name, sample_id, sample_time_hour, sample_type, strain_id,
      measure_time_hours, log_fluorescence, log_dil,
      fluorescence_raw, log_fluorescence_raw, edge_correction_applied, edge_correction_pair
    )

  samples_joined <- samples_base %>%
    dplyr::inner_join(
      curve_library,
      by = c("assay_batch", "measure_time_hours"),
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(curve_alpha_raw = 1 / sigma2) %>%
    dplyr::group_by(
      assay_batch, experiment_name, sample_id, sample_time_hour, sample_type, strain_id,
      measure_time_hours, log_fluorescence, log_dil
    ) %>%
    dplyr::mutate(curve_alpha = curve_alpha_raw / sum(curve_alpha_raw, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-curve_alpha_raw) %>%
    dplyr::mutate(
      response_lower_censored_curve = valid_response_bounds &
        log_fluorescence < log_fluorescence_response_lower,
      response_upper_censored_curve = valid_response_bounds &
        log_fluorescence > log_fluorescence_response_upper,
      response_any_censored_curve = response_lower_censored_curve |
        response_upper_censored_curve,
      log_fluorescence_bounded = dplyr::case_when(
        response_lower_censored_curve ~ log_fluorescence_response_lower,
        response_upper_censored_curve ~ log_fluorescence_response_upper,
        TRUE ~ log_fluorescence
      )
    )

  estimate_xsource_joint <- function(d) {
    lower <- max(d$x_min_std + d$log_dil, na.rm = TRUE)
    upper <- min(d$x_max_std + d$log_dil, na.rm = TRUE)

    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      return(tibble::tibble(
        x_source_log_fm = NA_real_, se_log_fm = NA_real_,
        ci_low_log_fm = NA_real_, ci_high_log_fm = NA_real_,
        b12_pm = NA_real_, b12_pm_ci_low = NA_real_, b12_pm_ci_high = NA_real_,
        n_obs = NA_integer_, eff_info = NA_real_, hit_lower = NA, hit_upper = NA,
        response_lower_censored = NA, response_upper_censored = NA,
        response_any_censored = NA,
        response_lower_all_censored = NA, response_upper_all_censored = NA,
        response_any_all_censored = NA,
        response_lower_censored_weight = NA_real_,
        response_upper_censored_weight = NA_real_,
        response_any_censored_weight = NA_real_
      ))
    }

    start <- median(d$c_ec50 + d$log_dil, na.rm = TRUE)

    obj <- function(x_source) {
      x_well <- x_source - d$log_dil
      y_hat <- f4pl(x_well, d$a_max, d$d_min, d$c_ec50, d$b_slope)
      w <- d$curve_alpha
      sum(w * (d$log_fluorescence_bounded - y_hat)^2, na.rm = TRUE)
    }

    res <- tryCatch(
      stats::optim(par = start, fn = obj, method = "Brent", lower = lower, upper = upper),
      error = function(e) list(par = NA_real_)
    )
    x_hat <- res$par

    if (!is.finite(x_hat)) {
      return(tibble::tibble(
        x_source_log_fm = NA_real_, se_log_fm = NA_real_,
        ci_low_log_fm = NA_real_, ci_high_log_fm = NA_real_,
        b12_pm = NA_real_, b12_pm_ci_low = NA_real_, b12_pm_ci_high = NA_real_,
        n_obs = dplyr::n_distinct(paste(d$measure_time_hours, d$log_dil, d$log_fluorescence)),
        eff_info = NA_real_, hit_lower = NA, hit_upper = NA,
        response_lower_censored = any(d$response_lower_censored_curve, na.rm = TRUE),
        response_upper_censored = any(d$response_upper_censored_curve, na.rm = TRUE),
        response_any_censored = any(d$response_any_censored_curve, na.rm = TRUE),
        response_lower_all_censored = all(d$response_lower_censored_curve),
        response_upper_all_censored = all(d$response_upper_censored_curve),
        response_any_all_censored = all(d$response_any_censored_curve),
        response_lower_censored_weight = sum(d$curve_alpha * d$response_lower_censored_curve, na.rm = TRUE) /
          sum(d$curve_alpha, na.rm = TRUE),
        response_upper_censored_weight = sum(d$curve_alpha * d$response_upper_censored_curve, na.rm = TRUE) /
          sum(d$curve_alpha, na.rm = TRUE),
        response_any_censored_weight = sum(d$curve_alpha * d$response_any_censored_curve, na.rm = TRUE) /
          sum(d$curve_alpha, na.rm = TRUE)
      ))
    }

    x_well_hat <- x_hat - d$log_dil
    fp_hat <- f4pl_prime(x_well_hat, d$a_max, d$d_min, d$c_ec50, d$b_slope)
    info <- sum(d$curve_alpha * (fp_hat^2) / d$sigma2, na.rm = TRUE)
    se <- if (is.finite(info) && info > 0) 1 / sqrt(info) else NA_real_

    tibble::tibble(
      x_source_log_fm = x_hat,
      se_log_fm = se,
      ci_low_log_fm = x_hat - ci_multiplier * se,
      ci_high_log_fm = x_hat + ci_multiplier * se,
      b12_pm = exp(x_hat) / 1000,
      b12_pm_ci_low = exp(x_hat - ci_multiplier * se) / 1000,
      b12_pm_ci_high = exp(x_hat + ci_multiplier * se) / 1000,
      n_obs = dplyr::n_distinct(paste(d$measure_time_hours, d$log_dil, d$log_fluorescence)),
      eff_info = info,
      hit_lower = abs(x_hat - lower) < hit_upper_tolerance,
      hit_upper = abs(x_hat - upper) < hit_upper_tolerance,
      response_lower_censored = any(d$response_lower_censored_curve, na.rm = TRUE),
      response_upper_censored = any(d$response_upper_censored_curve, na.rm = TRUE),
      response_any_censored = any(d$response_any_censored_curve, na.rm = TRUE),
      response_lower_all_censored = all(d$response_lower_censored_curve),
      response_upper_all_censored = all(d$response_upper_censored_curve),
      response_any_all_censored = all(d$response_any_censored_curve),
      response_lower_censored_weight = sum(d$curve_alpha * d$response_lower_censored_curve, na.rm = TRUE) /
        sum(d$curve_alpha, na.rm = TRUE),
      response_upper_censored_weight = sum(d$curve_alpha * d$response_upper_censored_curve, na.rm = TRUE) /
        sum(d$curve_alpha, na.rm = TRUE),
      response_any_censored_weight = sum(d$curve_alpha * d$response_any_censored_curve, na.rm = TRUE) /
        sum(d$curve_alpha, na.rm = TRUE)
    )
  }

  b12_concentration_estimates <- samples_joined %>%
    dplyr::group_by(assay_batch, experiment_name, sample_id, sample_time_hour, sample_type, strain_id) %>%
    tidyr::nest() %>%
    dplyr::mutate(est = purrr::map(data, estimate_xsource_joint)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(est) %>%
    dplyr::ungroup()

  if (verbose) {
    message(
      "b12_concentration_estimator: fitted ", nrow(curve_library),
      " standard curves; estimated ", nrow(b12_concentration_estimates),
      " sample concentrations; edge_correction=", edge_correction, "."
    )
  }

  list(
    estimates = b12_concentration_estimates,
    estimator_parameters = tibble::tibble(
      parameter = c(
        "standard_sample_types",
        "standard_cols_to_remove",
        "add_offset_pm",
        "edge_correction",
        "edge_outer_inner_rows",
        "min_edge_correction_pairs",
        "confidence_level",
        "ci_multiplier",
        "hit_upper_tolerance",
        "response_bound_fraction",
        "standard_curve_model",
        "standard_curve_grouping",
        "sample_curve_join_keys",
        "standard_exclusion_rule",
        "sample_exclusion_rule",
        "response_bound_rule",
        "response_censoring_rule",
        "hit_lower_rule",
        "hit_upper_rule",
        "na_estimate_rule"
      ),
      value = c(
        paste(standard_sample_types, collapse = ";"),
        paste(standard_cols_to_remove, collapse = ";"),
        add_offset_pm,
        edge_correction,
        paste(paste0(names(edge_outer_inner_rows), "_to_", unname(edge_outer_inner_rows)), collapse = ";"),
        min_edge_correction_pairs,
        confidence_level,
        ci_multiplier,
        hit_upper_tolerance,
        response_bound_fraction,
        "drc::L.4 four-parameter logistic model on log_fluorescence ~ log_b12_fm",
        "assay_batch + measure_time_hours + standard_bacteria_strain_id",
        "assay_batch + measure_time_hours",
        "standards must match standard_sample_types, have finite log_b12_fm and log_fluorescence, non-missing standard_bacteria_strain_id, and not be in standard_cols_to_remove",
        "samples must have assay_type == sample and finite corrected log_fluorescence",
        "per-curve response bounds are d_min + response_bound_fraction * (a_max - d_min) and d_min + (1 - response_bound_fraction) * (a_max - d_min)",
        "sample log_fluorescence is clamped to the per-curve response bounds before joint concentration optimization; response_*_censored flags indicate whether raw corrected fluorescence was clamped",
        "hit_lower is TRUE when the bounded optimizer returns the lower concentration bound within hit_upper_tolerance",
        "hit_upper is TRUE when the bounded optimizer returns the upper concentration bound within hit_upper_tolerance",
        "sample estimate is NA if no valid optimization range exists, optimization fails, or no fitted curve joins by assay_batch + measure_time_hours"
      )
    ),
    curve_library = curve_library,
    standards_used = standards_used,
    samples_joined = samples_joined,
    corrected_assay = df,
    edge_correction = edge_correction,
    edge_correction_model = if (!is.null(edge_correction_result)) edge_correction_result$model else NULL,
    edge_correction_pairs = if (!is.null(edge_correction_result)) edge_correction_result$paired else NULL,
    edge_correction_summary = if (!is.null(edge_correction_result)) edge_correction_result$summary else NULL,
    adjacent_pair_diagnostics = adjacent_pair_diagnostics
  )
}
