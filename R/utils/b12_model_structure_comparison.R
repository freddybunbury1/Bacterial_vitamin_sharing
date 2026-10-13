suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
if (!exists("build_kinetic_od_reference", mode = "function")) {
  source("R/utils/b12_two_stage_kinetic_pipeline.R")
}

# Keep this utility stable when sourced after packages that mask common dplyr
# verbs (notably MASS::select in the compiled analysis session).
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
group_by <- dplyr::group_by
distinct <- dplyr::distinct
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
pull <- dplyr::pull
across <- dplyr::across
n <- dplyr::n
if_else <- dplyr::if_else

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

trapezoid_auc <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

read_model_inputs <- function(repo_root = ".") {
  b12_path <- file.path(repo_root, "data/processed/b12_release_uptake_death/b12_results_b12_release_uptake_death.csv")
  od_path <- file.path(repo_root, "data/processed/b12_release_uptake_death/b12_od_viability_results_b12_release_uptake_death.csv")
  if (!file.exists(b12_path) || !file.exists(od_path)) stop("Required processed inputs are missing. Run the manuscript preprocessing/analysis first.")

  b12 <- read.csv(b12_path, check.names = FALSE, stringsAsFactors = FALSE)
  od <- read.csv(od_path, check.names = FALSE, stringsAsFactors = FALSE)

  # The same mean initial/final culture density used for the canonical rates.
  # Also use this density when transferring uptake to the live fraction, so
  # both model structures implement the main model's x * (1-f) assumption.
  reference_path <- file.path(repo_root, "data/processed/b12_release_uptake_death/modeling",
                              "b12_rate_od_reference_b12_release_uptake_death.csv")
  if (!file.exists(reference_path)) stop("Run the canonical first-order stage to export the shared OD reference first.")
  reference <- read.csv(reference_path, check.names = FALSE)
  ru_od <- reference |>
    select(strain_id, experiment_name, biological_replicate, added_glucose_mM_C,
           adjusted_od600 = od600_reference) |>
    mutate(experiment_type = "release_and_uptake")

  live_dead_od <- ru_od |>
    filter(added_glucose_mM_C == 0) |>
    mutate(experiment_type = "live_dead_mix")

  b12 <- b12 |>
    left_join(
      bind_rows(ru_od, live_dead_od),
      by = c("strain_id", "experiment_name", "biological_replicate", "added_glucose_mM_C", "experiment_type")
    ) |>
    mutate(
      od_valid = is.finite(adjusted_od600) & adjusted_od600 > 0,
      row_id = row_number()
    )

  trait_path <- file.path(repo_root, "data/processed/b12_trait_survey/b12_trait_survey_strain_summary.csv")
  if (file.exists(trait_path)) {
    trait <- read.csv(trait_path, check.names = FALSE, stringsAsFactors = FALSE) |>
      select(strain_id, b12_trait_group2) |>
      distinct()
    b12 <- b12 |> left_join(trait, by = "strain_id")
  } else {
    b12$b12_trait_group2 <- NA_character_
  }

  blank_offset <- b12 |>
    filter(strain_id == "blank", added_b12_pM == 0) |>
    summarise(value = safe_mean(b12_pm)) |>
    pull(value)

  quota_blank_sd_pM <- b12 |>
    filter(strain_id == "blank", added_b12_pM == 500, is.finite(b12_pm)) |>
    summarise(value = sd(b12_pm)) |>
    pull(value)

  list(data = b12, blank_offset = blank_offset, quota_blank_sd_pM = quota_blank_sd_pM)
}

group_variables <- function(spec) {
  vars <- "strain_id"
  if (spec$glucose_parameterization == "separate") vars <- c(vars, "added_glucose_mM_C")
  if (spec$replicate_parameterization == "separate") vars <- c(vars, "experiment_name", "biological_replicate")
  vars
}

release_initial_table <- function(dat, group_vars, treatment, blank_offset) {
  base <- dat |>
    filter(experiment_type == "release_and_uptake", added_b12_pM == 0, is.finite(b12_pm), is.finite(sample_time_minute))
  if (treatment == "zero" || treatment == "estimated_in_model") return(NULL)
  if (treatment == "global_fixed") {
    earliest <- min(base$sample_time_minute, na.rm = TRUE)
    return(data.frame(release_initial_b12_pM = max(safe_mean(base$b12_pm[base$sample_time_minute == earliest] - blank_offset), 0)))
  }
  base |>
    group_by(across(all_of(group_vars))) |>
    filter(sample_time_minute == min(sample_time_minute, na.rm = TRUE)) |>
    summarise(release_initial_b12_pM = max(safe_mean(b12_pm - blank_offset), 0), .groups = "drop")
}

# B is concentration above the additive assay offset. Q is cumulative uptake
# per unit OD. Uptake conserves B + od * Q; release adds od * r * time.
quota_predict_state <- function(times, b0, od, pars, release_per_od = 0,
                                q0 = 0, rtol = 1e-8, atol = 1e-8,
                                force_adaptive = FALSE) {
  k <- pars[["k_per_od"]]
  capacity <- pars[["Qmax_per_od"]]
  if (any(!is.finite(c(times, b0, od, k, capacity, release_per_od, q0))) ||
      any(times < 0) || min(b0, od, k, release_per_od, q0) < 0 ||
      capacity <= 0 || q0 > capacity || rtol <= 0 || atol <= 0) {
    stop("Invalid quota time, initial condition, density or parameter.")
  }
  if (!length(times)) return(data.frame(time = numeric(), B = numeric(), Q = numeric()))
  if (all(times == 0)) {
    B <- rep(b0, length(times)); Q <- rep(q0, length(times))
  } else if (od == 0 || k == 0 || q0 == capacity) {
    B <- b0 + od * release_per_od * times
    Q <- if (od == 0 && k > 0) {
      q0 + (capacity - q0) * (-expm1(-k * b0 / capacity * times))
    } else rep(q0, length(times))
  } else if (release_per_od == 0 && !force_adaptive) {
    available <- od * (capacity - q0)
    delta <- available - b0
    rate <- k / capacity
    if (b0 == 0) {
      B <- rep(0, length(times))
    } else if (delta >= 0) {
      z <- rate * delta * times
      # phi(z)=(1-exp(-z))/z has limit 1 at zero. This includes the
      # delta=0 solution without subtracting nearly equal concentrations.
      phi <- rep(1, length(z))
      nz <- z != 0
      phi[nz] <- -expm1(-z[nz]) / z[nz]
      B <- b0 * exp(-z) / (1 + b0 * rate * times * phi)
    } else {
      floor_B <- -delta
      z <- rate * floor_B * times
      B <- floor_B / (floor_B / b0 + available / b0 * (-expm1(-z)))
    }
    Q <- q0 + (b0 - B) / od
  } else {
    if (!requireNamespace("deSolve", quietly = TRUE)) stop("Install deSolve for quota release predictions.")
    target <- sort(unique(c(0, times)))
    deriv <- function(time, state, parameters) {
      flux <- k * max(state[["B"]], 0) * max(1 - state[["Q"]] / capacity, 0)
      list(c(B = od * (release_per_od - flux), Q = flux))
    }
    solution <- deSolve::ode(
      y = c(B = b0, Q = q0), times = target, func = deriv, parms = NULL,
      method = "lsoda", rtol = rtol, atol = atol, maxsteps = 100000L
    )
    if (attr(solution, "istate")[[1]] < 0 || nrow(solution) != length(target)) {
      stop("Quota adaptive solver did not complete the requested trajectory.")
    }
    B <- solution[match(times, solution[, "time"]), "B"]
    Q <- solution[match(times, solution[, "time"]), "Q"]
  }
  expected_mass <- b0 + od * q0 + od * release_per_od * times
  b_tol <- 100 * (atol + rtol * pmax(1, expected_mass))
  q_tol <- 100 * (atol + rtol * pmax(1, capacity))
  if (any(!is.finite(c(B, Q))) || any(B < -b_tol) ||
      any(Q < q0 - q_tol | Q > capacity + q_tol) ||
      any(abs(B + od * Q - expected_mass) > b_tol + od * q_tol)) {
    stop("Quota trajectory failed nonnegativity, capacity or mass-balance checks.")
  }
  if (release_per_od == 0 &&
      (any(B > b0 + b_tol) || any(B < max(b0 - od * (capacity - q0), 0) - b_tol))) {
    stop("Quota uptake trajectory exceeds available uptake capacity.")
  }
  data.frame(time = times, B = B, Q = Q)
}

predict_model_b12 <- function(times, b0, od, structure, pars, release_per_od = 0, q0 = 0) {
  if (!length(times)) return(numeric())
  if (structure == "first_order") {
    rate <- pars[["k_per_od"]] * od
    if (release_per_od == 0) return(b0 * exp(-rate * times))
    equilibrium <- release_per_od / pars[["k_per_od"]]
    return(equilibrium + (b0 - equilibrium) * exp(-rate * times))
  }
  if (structure != "quota_saturation") stop("Unknown uptake structure: ", structure)
  quota_predict_state(times, b0, od, pars, release_per_od, q0)$B
}

parameter_names <- function(structure, include_release = FALSE, include_release_initial = FALSE) {
  uptake <- switch(structure,
    first_order = "k_per_od",
    quota_saturation = c("k_per_od", "Qmax_per_od")
  )
  c(uptake, if (include_release) "r_per_od", if (include_release_initial) "release_initial_b12_pM")
}

initial_parameters <- function(names, dat) {
  vals <- c(k_per_od = 0.05, Qmax_per_od = 500,
            r_per_od = 0.1, release_initial_b12_pM = max(safe_mean(dat$b12_pm[dat$added_b12_pM == 0]) - min(dat$b12_pm, na.rm = TRUE), 0.01))
  log(pmax(vals[names], 1e-8))
}

predict_kinetic_group <- function(dat, structure, pars, blank_offset, release_initial_b12_pM, include_release) {
  pred <- numeric(nrow(dat))
  split_rows <- split(seq_len(nrow(dat)), interaction(dat$added_b12_pM, dat$adjusted_od600, drop = TRUE))
  for (idx in split_rows) {
    b0 <- dat$added_b12_pM[idx[[1]]] + release_initial_b12_pM
    pred[idx] <- blank_offset + predict_model_b12(
      dat$sample_time_minute[idx], b0, dat$adjusted_od600[idx[[1]]], structure, pars,
      release_per_od = if (include_release) pars[["r_per_od"]] else 0
    )
  }
  pred
}

fit_parameters <- function(dat, structure, fit_scale, blank_offset, release_initial_value,
                           include_release, estimate_release_initial, fixed_uptake = NULL,
                           quota_capacity_floor_pM = 0) {
  names <- parameter_names(structure, include_release, estimate_release_initial)
  if (!is.null(fixed_uptake)) names <- setdiff(names, names(fixed_uptake))
  if (!length(names)) return(list(ok = TRUE, parameters = fixed_uptake, objective = NA_real_, convergence = 0L))
  lower_values <- setNames(rep(1e-8, length(names)), names)
  if (structure == "quota_saturation" && "Qmax_per_od" %in% names &&
      is.finite(quota_capacity_floor_pM) && quota_capacity_floor_pM > 0) {
    group_od <- safe_mean(dat$adjusted_od600)
    if (is.finite(group_od) && group_od > 0) {
      lower_values[["Qmax_per_od"]] <- quota_capacity_floor_pM / group_od
    }
  }
  start <- pmax(initial_parameters(names, dat), log(lower_values))
  objective <- function(theta) {
    pars <- c(fixed_uptake, setNames(exp(theta), names))
    release_initial_b12_pM <- if (estimate_release_initial) pars[["release_initial_b12_pM"]] else release_initial_value
    pred <- tryCatch(
      predict_kinetic_group(dat, structure, pars, blank_offset, release_initial_b12_pM, include_release),
      error = function(e) rep(NA_real_, nrow(dat))
    )
    if (any(!is.finite(pred)) || any(pred <= 0)) return(.Machine$double.xmax / 1e100)
    residual <- if (fit_scale == "log") log(dat$b12_pm) - log(pred) else dat$b12_pm - pred
    sum(residual^2)
  }
  fit <- tryCatch({
    if (length(start) == 1L) {
      # A broad single optimize() call assumes unimodality. Screen the whole
      # log interval and refine every sampled local minimum instead.
      grid <- seq(-20, 20, length.out = 161)
      values <- vapply(grid, objective, numeric(1))
      local <- which(values[2:160] <= values[1:159] & values[2:160] <= values[3:161]) + 1L
      candidates <- lapply(local, function(i) optimize(objective,
        interval = grid[c(i - 1L, i + 1L)], tol = 1e-8))
      candidates <- c(candidates, list(list(minimum=grid[1], objective=values[1]),
                                       list(minimum=grid[161], objective=values[161])))
      one <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "objective"))]]
      list(par = one$minimum, value = one$objective, convergence = 0L)
    } else {
      upper_values <- c(
        k_per_od = 1e3,
        Qmax_per_od = 1e6, r_per_od = 1e4, release_initial_b12_pM = 1e4
      )
      initial_fit <- optim(
        start, objective, method = "L-BFGS-B",
        lower = log(unname(lower_values[names])),
        upper = log(unname(upper_values[names])),
        control = list(maxit = 1500)
      )
      # Near-flat uptake curves can trigger L-BFGS-B line-search failures.
      # Retry only failed quota fits with the same objective and log bounds.
      if (structure == "quota_saturation" && initial_fit$convergence != 0) {
        retry <- nlminb(
          initial_fit$par, objective,
          lower = log(unname(lower_values[names])),
          upper = log(unname(upper_values[names])),
          control = list(iter.max = 1500, eval.max = 3000, rel.tol = 1e-10)
        )
        if (retry$convergence == 0 && is.finite(retry$objective) &&
            retry$objective <= initial_fit$value + 1e-8) {
          initial_fit <- list(par = retry$par, value = retry$objective, convergence = 0L)
        }
      }
      initial_fit
    }
  }, error = identity)
  if (inherits(fit, "error")) return(list(ok = FALSE, message = conditionMessage(fit)))
  if (!is.finite(fit$value) || fit$value >= .Machine$double.xmax / 1e100) {
    return(list(ok = FALSE, message = "invalid_model_predictions"))
  }
  list(ok = fit$convergence == 0, parameters = c(fixed_uptake, setNames(exp(fit$par), names)),
       objective = fit$value, convergence = fit$convergence,
       message = if (fit$convergence == 0) "ok" else "optimizer_nonconvergence")
}

fit_one_group <- function(dat, spec, blank_offset, release_initial_value, quota_capacity_floor_pM = 0) {
  dat <- dat |>
    filter(experiment_type == "release_and_uptake", od_valid, is.finite(b12_pm), b12_pm > 0, is.finite(sample_time_minute))
  min_n <- 5L
  uptake_dat <- filter(dat, added_b12_pM > 0)
  release_dat <- filter(dat, added_b12_pM == 0)
  if (nrow(uptake_dat) < min_n || nrow(release_dat) < 3L) return(list(ok = FALSE, message = "insufficient_valid_kinetic_rows"))
  estimate_release_initial <- spec$release_initial_treatment == "estimated_in_model"
  if (spec$pipeline == "two_stage") {
    up <- fit_parameters(uptake_dat, spec$uptake_structure, spec$fit_scale, blank_offset,
                         0, FALSE, FALSE,
                         quota_capacity_floor_pM = quota_capacity_floor_pM)
    if (!up$ok) return(up)
    uptake_names <- parameter_names(spec$uptake_structure)
    release_initial_second <- release_initial_value
    rel <- fit_parameters(release_dat, spec$uptake_structure, spec$fit_scale, blank_offset,
                          release_initial_second, TRUE, FALSE, up$parameters[uptake_names],
                          quota_capacity_floor_pM = quota_capacity_floor_pM)
    if (!rel$ok) return(rel)
    pars <- c(rel$parameters, release_initial_b12_pM = release_initial_second)
    objective <- up$objective + rel$objective
  } else {
    joint <- fit_parameters(dat, spec$uptake_structure, spec$fit_scale, blank_offset,
                            release_initial_value, TRUE, estimate_release_initial,
                            quota_capacity_floor_pM = quota_capacity_floor_pM)
    if (!joint$ok) return(joint)
    pars <- joint$parameters
    if (!"release_initial_b12_pM" %in% names(pars)) pars <- c(pars, release_initial_b12_pM = release_initial_value)
    objective <- joint$objective
  }
  if (spec$uptake_structure == "first_order") pars <- c(pars, E_pM = pars[["r_per_od"]] / pars[["k_per_od"]])
  list(ok = TRUE, parameters = pars, objective = objective, message = "ok", n_fit = nrow(dat))
}

build_group_key <- function(dat, group_vars) do.call(paste, c(dat[group_vars], sep = "|"))

use_canonical_first_order_summary <- function(source, canonical, keys, target, value) {
  if (anyDuplicated(canonical[keys])) stop("Duplicate canonical summary coordinates.")
  values <- canonical[c(keys, value)]
  names(values)[length(names(values))] <- ".canonical_value"
  source <- left_join(source, values, by = keys)
  first <- source$model_structure == "first_order"
  source[[target]][first] <- source$.canonical_value[first]
  # Blank controls are empirical displays, not OD-normalized model fits.
  source[[target]][source$strain_id == "blank"] <- NA_real_
  source$.canonical_value <- NULL
  source
}

# Adapter only: canonical first-order fits and trajectories are passed in,
# never fitted or predicted again for the supplementary figures.
reuse_first_order_result <- function(spec, fits, release_initial, uptake_raw,
                                    release_raw, mix_raw, od_reference, blank_offset) {
  keys <- group_variables(spec)
  fitted <- fits |>
    filter(strain_id != "blank") |>
    left_join(release_initial, by = keys) |>
    mutate(group_key = build_group_key(pick(all_of(keys)), keys))
  parameters <- fitted |>
    transmute(across(all_of(keys)), model_spec_id = spec$model_spec_id, group_key,
              k_per_od = k_per_od600, r_per_od = r_per_od600,
              release_initial_b12_pM, E_pM = E_per_cell) |>
    pivot_longer(c(k_per_od, r_per_od, release_initial_b12_pM, E_pM),
                 names_to = "parameter", values_to = "estimate")
  kinetic <- bind_rows(uptake_raw, release_raw) |>
    filter(strain_id != "blank") |>
    mutate(.group_key = build_group_key(pick(all_of(keys)), keys),
           adjusted_od600 = od600_reference, od_valid = is.finite(adjusted_od600) & adjusted_od600 > 0,
           fitted_b12 = model_pred, model_spec_id = spec$model_spec_id)
  status <- fitted |>
    select(all_of(keys), group_key) |>
    left_join(kinetic |> count(.group_key, name = "n_fit"), by = c("group_key" = ".group_key")) |>
    mutate(model_spec_id = spec$model_spec_id, fit_ok = TRUE,
           fit_message = "reused_canonical_first_order", objective = NA_real_)
  predictions <- mix_raw |>
    filter(strain_id != "blank") |>
    left_join(od_reference |> select(all_of(keys), od600_reference, measured_source_od600 = od600_t0), by = keys) |>
    mutate(.group_key = build_group_key(pick(all_of(keys)), keys),
           adjusted_od600 = od600_reference, source_od600 = od600_reference,
           od_valid = is.finite(adjusted_od600) & adjusted_od600 > 0,
           full_kill_b12 = b_full_corrected + blank_offset,
           prediction_ok = is.finite(predicted_b12), model_spec_id = spec$model_spec_id,
           error = b12_pm - predicted_b12,
           normalized_error = error / full_kill_b12,
           absolute_normalized_error = abs(normalized_error),
           squared_normalized_error = normalized_error^2)
  list(status = status, parameters = parameters, kinetic_predictions = kinetic,
       predictions = predictions, metrics = summarize_prediction_metrics(predictions, spec))
}

summarize_prediction_metrics <- function(pred_rows, spec) {
  auc_by_curve <- pred_rows |>
    filter(prediction_ok, is.finite(full_kill_b12), full_kill_b12 > 0) |>
    mutate(observed_b12_fraction = b12_pm / full_kill_b12,
           predicted_b12_fraction = predicted_b12 / full_kill_b12) |>
    group_by(strain_id, experiment_name, biological_replicate) |>
    summarise(observed_auc = trapezoid_auc(killed_proportion, observed_b12_fraction),
              predicted_auc = trapezoid_auc(killed_proportion, predicted_b12_fraction),
              auc_difference = predicted_auc - observed_auc,
              absolute_auc_difference = abs(auc_difference), .groups = "drop")
  pred_rows |>
    filter(prediction_ok, is.finite(normalized_error)) |>
    summarise(n_prediction = n(), normalized_mae = mean(absolute_normalized_error),
              normalized_rmse = sqrt(mean(squared_normalized_error)), .groups = "drop") |>
    mutate(mean_absolute_auc_difference = safe_mean(auc_by_curve$absolute_auc_difference),
           mean_signed_auc_difference = safe_mean(auc_by_curve$auc_difference),
           model_spec_id = spec$model_spec_id, .before = 1) |>
    mutate(normalized_mae = if_else(n_prediction > 0, normalized_mae, NA_real_),
           normalized_rmse = if_else(n_prediction > 0, normalized_rmse, NA_real_))
}

run_model_spec <- function(spec, data, blank_offset, quota_blank_sd_pM = 0, max_groups = Inf) {
  if (spec$uptake_structure == "quota_saturation" && !requireNamespace("deSolve", quietly = TRUE)) {
    stop("Install deSolve before fitting the supplementary quota model.")
  }
  group_vars <- group_variables(spec)
  kinetic <- data |>
    filter(experiment_type == "release_and_uptake", strain_id != "blank")
  dt <- release_initial_table(data, group_vars, spec$release_initial_treatment, blank_offset)
  if (!is.null(dt) && length(group_vars) && spec$release_initial_treatment == "group_specific_fixed") kinetic <- left_join(kinetic, dt, by = group_vars)
  if (spec$release_initial_treatment == "global_fixed") kinetic$release_initial_b12_pM <- dt$release_initial_b12_pM[[1]]
  if (spec$release_initial_treatment %in% c("zero", "estimated_in_model")) kinetic$release_initial_b12_pM <- 0
  kinetic$.group_key <- build_group_key(kinetic, group_vars)
  groups <- split(kinetic, kinetic$.group_key)
  if (is.finite(max_groups)) groups <- head(groups, max_groups)

  fits <- lapply(names(groups), function(key) {
    g <- groups[[key]]
    quota_floor_multiplier <- spec$quota_floor_sd_multiplier %||% 0
    quota_capacity_floor_pM <- quota_floor_multiplier * quota_blank_sd_pM
    result <- fit_one_group(
      g, spec, blank_offset, unique(g$release_initial_b12_pM)[[1]], quota_capacity_floor_pM
    )
    identifiers <- g[1, group_vars, drop = FALSE]
    rownames(identifiers) <- NULL
    status <- cbind(identifiers, data.frame(model_spec_id = spec$model_spec_id, group_key = key,
      fit_ok = unname(result$ok), fit_message = unname(result$message %||% "unknown"), n_fit = unname(result$n_fit %||% nrow(g)),
      objective = result$objective %||% NA_real_))
    parameters <- if (result$ok) cbind(identifiers, data.frame(model_spec_id = spec$model_spec_id,
      group_key = key, parameter = unname(names(result$parameters)), estimate = unname(as.numeric(result$parameters)))) else data.frame()
    list(status = status, parameters = parameters, pars = if (result$ok) result$parameters else NULL)
  })
  status <- bind_rows(lapply(fits, `[[`, "status"))
  parameters <- bind_rows(lapply(fits, `[[`, "parameters"))
  pars_by_key <- setNames(lapply(fits, `[[`, "pars"), names(groups))

  kinetic_predictions <- bind_rows(lapply(names(groups), function(key) {
    g <- groups[[key]]
    pars <- pars_by_key[[key]]
    if (is.null(pars)) return(data.frame())
    if (spec$prediction_model == "null_no_uptake") {
      g$fitted_b12 <- blank_offset + g$added_b12_pM
    } else if (spec$prediction_model == "uptake_only") {
      g$fitted_b12 <- predict_kinetic_group(
        g, spec$uptake_structure, pars, blank_offset,
        0, include_release = FALSE
      )
    } else if (spec$pipeline == "two_stage") {
      g$fitted_b12 <- NA_real_
      uptake_idx <- g$added_b12_pM > 0
      release_idx <- !uptake_idx
      if (any(uptake_idx)) {
        g$fitted_b12[uptake_idx] <- predict_kinetic_group(
          g[uptake_idx, ], spec$uptake_structure, pars, blank_offset,
          0, include_release = FALSE
        )
      }
      if (any(release_idx)) {
        g$fitted_b12[release_idx] <- predict_kinetic_group(
          g[release_idx, ], spec$uptake_structure, pars, blank_offset,
          pars[["release_initial_b12_pM"]], include_release = TRUE
        )
      }
    } else {
      g$fitted_b12 <- predict_kinetic_group(
        g, spec$uptake_structure, pars, blank_offset,
        pars[["release_initial_b12_pM"]], include_release = TRUE
      )
    }
    g$model_spec_id <- spec$model_spec_id
    g
  }))

  pred_data <- data |>
    filter(experiment_type == "live_dead_mix", strain_id != "blank")
  pred_data$.group_key <- build_group_key(pred_data, group_vars)

  full_kill_keys <- c("strain_id", "experiment_name", "biological_replicate", "added_glucose_mM_C")
  b_table <- pred_data |>
    filter(killed_proportion == 1, od_valid, is.finite(b12_pm)) |>
    group_by(across(all_of(full_kill_keys))) |>
    summarise(full_kill_b12 = safe_mean(b12_pm), source_od600 = safe_mean(adjusted_od600), .groups = "drop")
  pred_data <- left_join(pred_data, b_table, by = full_kill_keys)
  # S13 can omit release without changing the S12 uptake and S14 release fits.
  mix_prediction_model <- spec$live_dead_prediction_model %||% spec$prediction_model

  pred_rows <- lapply(seq_len(nrow(pred_data)), function(i) {
    row <- pred_data[i, ]
    pars <- pars_by_key[[row$.group_key]]
    if (is.null(pars) || !row$od_valid || !is.finite(row$full_kill_b12) || row$full_kill_b12 <= 0) {
      row$predicted_b12 <- NA_real_; row$prediction_ok <- FALSE
    } else {
      release_initial_value <- pars[["release_initial_b12_pM"]]
      b_per_od <- max((row$full_kill_b12 - blank_offset) / row$source_od600, 0)
      release_initial_b12_pM <- if (mix_prediction_model == "release_plus_uptake") release_initial_value else 0
      b0 <- release_initial_b12_pM + b_per_od * row$adjusted_od600 * row$killed_proportion
      live_od <- row$adjusted_od600 * (1 - row$killed_proportion)
      if (mix_prediction_model == "null_no_uptake") value <- b0
      else value <- predict_model_b12(row$sample_time_minute, b0, live_od, spec$uptake_structure, pars,
        release_per_od = if (mix_prediction_model == "release_plus_uptake") pars[["r_per_od"]] else 0)
      row$predicted_b12 <- blank_offset + value
      row$prediction_ok <- is.finite(row$predicted_b12)
    }
    row$model_spec_id <- spec$model_spec_id
    row
  }) |>
    bind_rows() |>
    mutate(
      error = b12_pm - predicted_b12,
      normalized_error = if_else(is.finite(full_kill_b12) & full_kill_b12 > 0, error / full_kill_b12, NA_real_),
      absolute_normalized_error = abs(normalized_error),
      squared_normalized_error = normalized_error^2
    )

  metrics <- summarize_prediction_metrics(pred_rows, spec)
  list(
    status = status,
    parameters = parameters,
    kinetic_predictions = kinetic_predictions,
    predictions = pred_rows,
    metrics = metrics
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x


# Recheck the fitted trajectories against a tighter independent integration.
# These diagnostics are exported during the normal supplementary fit workflow.
validate_quota_result <- function(result, blank_offset, live_dead_prediction_model = "release_plus_uptake") {
  records <- list()
  check_curve <- function(key, setting, times, b0, od, pars, release) {
    current <- quota_predict_state(times, b0, od, pars, release)
    reference <- quota_predict_state(times, b0, od, pars, release,
                                     rtol = 1e-10, atol = 1e-10,
                                     force_adaptive = TRUE)
    b_error <- max(abs(current$B - reference$B))
    q_error <- max(abs(current$Q - reference$Q))
    mass_error <- max(abs(current$B + od * current$Q -
                           (b0 + od * release * times)))
    ok <- b_error <= 1e-4 + 1e-7 * max(1, abs(reference$B)) &&
      q_error <= 1e-4 + 1e-7 * pars[['Qmax_per_od']]
    records[[length(records) + 1L]] <<- data.frame(
      group_key = key, setting = setting, max_B_difference_pM = b_error,
      max_Q_difference_per_od = q_error, mass_balance_error_pM = mass_error,
      minimum_B_pM = min(current$B),
      maximum_capacity_excess_per_od = max(current$Q - pars[['Qmax_per_od']]),
      numerical_check_ok = ok
    )
  }
  for (key in unique(result$parameters$group_key)) {
    p <- result$parameters[result$parameters$group_key == key, ]
    pars <- setNames(p$estimate, p$parameter)
    g <- result$kinetic_predictions[result$kinetic_predictions$.group_key == key, ]
    for (added in unique(g$added_b12_pM)) {
      d <- g[g$added_b12_pM == added & g$od_valid, ]
      if (!nrow(d)) next
      check_curve(key, if (added > 0) 'uptake' else 'release',
                  sort(unique(c(0, d$sample_time_minute))),
                  if (added > 0) added else pars[['release_initial_b12_pM']],
                  d$adjusted_od600[1], pars,
                  if (added > 0) 0 else pars[['r_per_od']])
    }
    mix <- result$predictions[result$predictions$.group_key == key &
                                result$predictions$prediction_ok, ]
    for (i in seq_len(nrow(mix))) {
      d <- mix[i, ]
      b <- max(d$full_kill_b12 - blank_offset, 0) / d$source_od600
      include_release <- live_dead_prediction_model == 'release_plus_uptake'
      b0 <- (if (include_release) pars[['release_initial_b12_pM']] else 0) +
        b * d$adjusted_od600 * d$killed_proportion
      check_curve(key, paste0('live_dead_', d$killed_proportion),
                  c(0, d$sample_time_minute), b0,
                  d$adjusted_od600 * (1 - d$killed_proportion), pars,
                  if (include_release) pars[['r_per_od']] else 0)
    }
  }
  diagnostics <- bind_rows(records)
  if (nrow(diagnostics) && any(!diagnostics$numerical_check_ok)) {
    stop('Fitted quota trajectories failed tolerance-convergence validation.')
  }
  diagnostics
}

write_model_spec_result <- function(result, spec, output_dir, blank_offset) {
  if (spec$uptake_structure == 'quota_saturation') {
    result$numerical_checks <- validate_quota_result(
      result, blank_offset, spec$live_dead_prediction_model %||% spec$prediction_model
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(result)) {
    write.csv(result[[name]], file.path(output_dir, paste0(spec$model_spec_id, '_', name, '.csv')),
              row.names = FALSE, na = '')
  }
  invisible(result)
}
