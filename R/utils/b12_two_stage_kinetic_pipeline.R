############################################################
############ B12 Two-Stage Kinetic Model Pipeline
############################################################

library(minpack.lm)
library(dplyr)
library(broom)
library(tidyr)

# Locate candidate basins before local least-squares refinement. A single
# starting value or a single broad optimize() interval can miss a better fit.
screen_scalar_fit <- function(objective, lower, upper) {
  grid <- seq(log(lower), log(upper), length.out = 161)
  values <- vapply(grid, objective, numeric(1))
  local <- which(values[2:160] <= values[1:159] & values[2:160] <= values[3:161]) + 1L
  candidates <- lapply(local, function(i) optimize(objective,
    interval = grid[c(i-1L, i+1L)], tol = 1e-10))
  candidates <- c(candidates, list(list(minimum=grid[1], objective=values[1]),
                                   list(minimum=grid[161], objective=values[161])))
  best <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "objective"))]]
  list(start = exp(best$minimum), objective = best$objective)
}

refine_screened_fit <- function(fit_one, starts, screen) {
  fits <- lapply(unique(starts), function(start) tryCatch(fit_one(start), error=function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("No scalar least-squares refinement converged.")
  best <- fits[[which.min(vapply(fits, deviance, numeric(1)))]]
  if (deviance(best) > screen$objective + 1e-6 * max(1, screen$objective)) {
    stop("Scalar refinement is worse than the screened solution.")
  }
  best
}

# The legacy k_per_cell and r_per_cell fields are culture-level coefficients.
# Retain them for predictions; expose explicitly named OD-specific summaries.
build_kinetic_od_reference <- function(od_results) {
  source_keys <- c("strain_id", "experiment_name", "biological_replicate")
  fit_keys <- c(source_keys, "added_glucose_mM_C")
  required <- c(fit_keys, "experiment_type", "sample_time_minute",
                "killed_proportion", "added_b12_pM", "adjusted_od600")
  if (!all(required %in% names(od_results))) {
    stop("OD reference input lacks required columns.")
  }

  # The untreated live/dead source is the common starting culture, before
  # splitting into the glucose and added-B12 conditions.
  initial <- od_results |>
    dplyr::filter(experiment_type == "live_dead_source",
                  sample_time_minute == 0, killed_proportion == 0,
                  added_glucose_mM_C == 0, added_b12_pM == 0) |>
    dplyr::select(dplyr::all_of(source_keys), od600_t0 = adjusted_od600)
  endpoint <- od_results |>
    dplyr::filter(experiment_type == "release_and_uptake",
                  sample_time_minute == 240, added_b12_pM %in% c(0, 500))
  if (anyDuplicated(initial[source_keys]) ||
      anyDuplicated(endpoint[c(fit_keys, "added_b12_pM")])) {
    stop("Duplicate culture keys in initial/final OD measurements.")
  }

  reference <- endpoint |>
    dplyr::group_by(dplyr::across(dplyr::all_of(fit_keys))) |>
    dplyr::summarise(
      # Average matched release and uptake wells to use one density for k and r.
      od600_t240 = mean(adjusted_od600),
      n_od600_t240 = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::left_join(initial, by = source_keys) |>
    dplyr::mutate(
      od600_reference = (od600_t0 + od600_t240) / 2,
      rate_od_status = dplyr::case_when(
        strain_id == "blank" ~ "blank_control",
        n_od600_t240 != 2L ~ "missing_final_well",
        !is.finite(od600_t0) | !is.finite(od600_t240) ~ "missing_od",
        od600_t0 <= 0 | od600_t240 <= 0 ~ "nonpositive_od",
        TRUE ~ "valid"
      )
    )
  invalid <- dplyr::filter(reference, !rate_od_status %in% c("valid", "blank_control"))
  if (nrow(invalid)) {
    stop("Cannot normalize rates: invalid initial/final OD for ",
         paste(unique(invalid$strain_id), collapse = ", "))
  }
  reference
}

add_od_specific_rates <- function(fits, od_reference) {
  fit_keys <- c("strain_id", "experiment_name", "biological_replicate",
                "added_glucose_mM_C")
  if (anyDuplicated(od_reference[fit_keys])) stop("Duplicate OD reference keys.")
  out <- dplyr::left_join(fits, od_reference, by = fit_keys)
  if (any(out$strain_id != "blank" &
          (is.na(out$rate_od_status) | out$rate_od_status != "valid"))) {
    stop("Every non-blank fitted culture must have a valid matched OD reference.")
  }
  denominator <- ifelse(out$rate_od_status == "valid", out$od600_reference, NA_real_)
  for (parameter in c("k", "r")) {
    legacy <- paste0(parameter, "_per_cell")
    if (legacy %in% names(out)) {
      out[[paste0(parameter, "_per_od600")]] <- out[[legacy]] / denominator
      # Conditional fit SE: uncertainty in measured OD is not propagated.
      out[[paste0(parameter, "_per_od600_se")]] <- out[[paste0(legacy, "_se")]] / denominator
    }
  }
  out
}

############################################################
############ Uptake-only model (estimate k on B0 > 0)
############################################################

fit_uptake_model <- function(b12_results, blank_offset, group_vars,
                             fit_scale = c("log", "linear")){
  
  fit_scale <- match.arg(fit_scale)
  
  fits_uptake <- b12_results %>%
    dplyr::filter(
      experiment_type == "release_and_uptake",
      added_b12_pM > 0
    ) %>%
    dplyr::mutate(
      y = as.numeric(b12_pm),
      x = as.numeric(sample_time_minute),
      B0_nominal = as.numeric(added_b12_pM)
    ) %>%
    dplyr::filter(is.finite(y), y > 0, is.finite(x), x >= 0) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::group_modify(~{
      
      dat <- .x
      if(nrow(dat) < 4) return(data.frame())
      
      x_max <- max(dat$x)
      dat$x_scaled <- dat$x / x_max
      
      fit <- tryCatch({
        k_scaled_lower <- 1e-4 * x_max
        objective <- function(log_rate) {
          predicted <- blank_offset + dat$B0_nominal * exp(-exp(log_rate) * dat$x)
          residual <- if (fit_scale == "log") log(dat$y) - log(predicted) else dat$y - predicted
          sum(residual^2)
        }
        screen <- screen_scalar_fit(objective, 1e-4, 1e3)
        fit_one <- function(initial) {
        
        if (fit_scale == "log") {
          nlsLM(
            log(y) ~ log(blank_offset +
                           B0_nominal *
                           exp(-k_scaled * x_scaled)),
            data = dat,
            start = list(k_scaled = initial),
            lower = c(k_scaled = k_scaled_lower),
            control = nls.lm.control(maxiter = 500)
          )
        } else {
          nlsLM(
            y ~ blank_offset +
              B0_nominal *
              exp(-k_scaled * x_scaled),
            data = dat,
            start = list(k_scaled = initial),
            lower = c(k_scaled = k_scaled_lower),
            control = nls.lm.control(maxiter = 500)
          )
        }
        }
        refine_screened_fit(fit_one,
          c(0.1, max(screen$start * x_max, k_scaled_lower * (1 + 1e-8))), screen)
        
      }, error = function(e){
        message("Uptake fit failed: ", unique(dat$strain_id))
        return(NULL)
      })
      
      if(is.null(fit)) return(data.frame())
      
      out <- broom::tidy(fit)
      
      out$term[out$term == "k_scaled"] <- "k_per_cell"
      out$estimate[out$term == "k_per_cell"] <- out$estimate[out$term == "k_per_cell"] / x_max
      out$std.error[out$term == "k_per_cell"] <- out$std.error[out$term == "k_per_cell"] / x_max

      out
      
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(group_vars), term, estimate, std.error) %>%
    tidyr::pivot_wider(
      names_from = term,
      values_from = c(estimate, std.error),
      names_sep = "_"
    ) %>%
    dplyr::rename(
      k_per_cell = estimate_k_per_cell,
      k_per_cell_se = std.error_k_per_cell
    )
  
  return(fits_uptake)
  
}

############################################################
############ Release model (estimate r with k fixed, B0 = 0)
############################################################
fit_release_model <- function(b12_results, release_initial_b12_estimates, fits_uptake,
                              blank_offset, group_vars,
                              fit_scale = c("log", "linear")){
  
  fit_scale <- match.arg(fit_scale)
  
  fits_release <- b12_results %>%
    dplyr::filter(
      experiment_type == "release_and_uptake",
      added_b12_pM == 0
    ) %>%
    dplyr::mutate(
      y = as.numeric(b12_pm),
      x = as.numeric(sample_time_minute)
    ) %>%
    dplyr::filter(is.finite(y), y > 0, is.finite(x), x >= 0) %>%
    dplyr::left_join(release_initial_b12_estimates, by = group_vars) %>%
    dplyr::left_join(fits_uptake, by = group_vars) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::group_modify(~{
      
      dat <- .x
      if(nrow(dat) < 3) return(data.frame())
      
      k_val <- unique(dat$k_per_cell)
      if(length(k_val) != 1 || is.na(k_val)) return(data.frame())
      
      fit <- tryCatch({
        objective <- function(log_rate) {
          equilibrium <- exp(log_rate) / k_val
          predicted <- blank_offset + equilibrium +
            (dat$release_initial_b12_pM - equilibrium) * exp(-k_val * dat$x)
          residual <- if (fit_scale == "log") log(dat$y) - log(predicted) else dat$y - predicted
          sum(residual^2)
        }
        screen <- screen_scalar_fit(objective, 1e-4, 1e4)
        fit_one <- function(initial) {
        
        if (fit_scale == "log") {
          nlsLM(
            log(y) ~ log(
              blank_offset +
                (r_per_cell / k_val) +
                (release_initial_b12_pM - r_per_cell / k_val) *
                exp(-k_val * x)
            ),
            data = dat,
            start = list(r_per_cell = initial),
            lower = c(r_per_cell = 1e-4),
            control = nls.lm.control(maxiter = 500)
          )
        } else {
          nlsLM(
            y ~ blank_offset +
              (r_per_cell / k_val) +
              (release_initial_b12_pM - r_per_cell / k_val) *
              exp(-k_val * x),
            data = dat,
            start = list(r_per_cell = initial),
            lower = c(r_per_cell = 1e-4),
            control = nls.lm.control(maxiter = 500)
          )
        }
        }
        refine_screened_fit(fit_one, c(0.1, max(screen$start, 1e-4 * (1 + 1e-8))), screen)
        
      }, error = function(e){
        message("Release fit failed")
        return(NULL)
      })
      
      if(is.null(fit)) return(data.frame())
      
      out <- broom::tidy(fit)
      out
      
    }) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(group_vars), term, estimate, std.error) %>%
    tidyr::pivot_wider(
      names_from = term,
      values_from = c(estimate, std.error),
      names_sep = "_"
    ) %>%
    dplyr::rename(
      r_per_cell = estimate_r_per_cell,
      r_per_cell_se = std.error_r_per_cell
    ) %>%
    dplyr::left_join(fits_uptake, by = group_vars) %>%
    dplyr::mutate(
      E_per_cell = r_per_cell / k_per_cell
    )
  
  return(fits_release)
  
}
############################################################
############ Mix predictions (uptake-only)
############################################################

predict_mix_uptake <- function(b12_results, fits_uptake,
                               b_estimates, blank_offset, group_vars){
  
  mix_pred_uptake <- b12_results %>%
    dplyr::filter(experiment_type == "live_dead_mix") %>%
    dplyr::left_join(fits_uptake, by = group_vars) %>%
    dplyr::left_join(b_estimates, by = group_vars) %>%
    dplyr::mutate(
      f_dead = killed_proportion,
      f_live = 1 - killed_proportion,
      B0 = b_per_cell * f_dead,
      predicted_b12 = blank_offset +
        B0 * exp(-k_per_cell *
                   f_live *
                   sample_time_minute),
      null_predicted_b12 = blank_offset + B0
    )
  
  return(mix_pred_uptake)
  
}

############################################################
############ Master pipeline
############################################################

build_b12_model_outputs <- function(
  b12_results,
  group_by_glucose = TRUE,
  group_by_replicate = FALSE,
  fit_scale = c("log", "linear")
){
  
  fit_scale <- match.arg(fit_scale)
  
  if(group_by_glucose){
    group_vars <- c("strain_id","added_glucose_mM_C")
  } else {
    group_vars <- c("strain_id")
  }
  if(group_by_replicate){
    group_vars <- c(group_vars, "experiment_name", "biological_replicate")
  }
  
  join_vars <- group_vars
  
  ############################################################
  ############ Assay floor
  ############################################################
  
  blank_offset <- b12_results %>%
    dplyr::filter(strain_id == "blank",
                  added_b12_pM == 0) %>%
    dplyr::summarise(mean_blank = mean(b12_pm, na.rm = TRUE)) %>%
    dplyr::pull(mean_blank)
  
  ############################################################
  ############ Release-specific initial B12 estimate
  ############################################################
  
  release_initial_b12_estimates <- b12_results %>%
    dplyr::filter(experiment_type == "release_and_uptake",
                  added_b12_pM == 0) %>%
    dplyr::mutate(y = b12_pm,
                  x = sample_time_minute,
                  B0_nominal = added_b12_pM) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::slice_min(x, n = 1, with_ties = TRUE) %>%
    dplyr::summarise(
      release_initial_b12_pM = mean(y - blank_offset - B0_nominal,
                   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(release_initial_b12_pM = pmax(release_initial_b12_pM, 0))
  
  ############################################################
  ############ Dead-cell release estimate
  ############################################################
  
  b_estimates <- b12_results %>%
    dplyr::filter(experiment_type == "live_dead_mix",
                  killed_proportion == 1) %>%
    dplyr::mutate(b12_full_corrected = b12_pm - blank_offset) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarise(
      b_full_corrected = mean(b12_full_corrected, na.rm = TRUE),
      b_per_cell = mean(b12_full_corrected, na.rm = TRUE),
      .groups = "drop"
    )
  
  ############################################################
  ############ Fit models
  ############################################################
  
  fits_uptake <- fit_uptake_model(
    b12_results,
    blank_offset,
    group_vars,
    fit_scale = fit_scale
  )
  
  fits_release <- fit_release_model(
    b12_results,
    release_initial_b12_estimates,
    fits_uptake,
    blank_offset,
    group_vars,
    fit_scale = fit_scale
  )
  
  ############################################################
  ############ Predict mixes
  ############################################################
  
  mix_pred_uptake <- predict_mix_uptake(
    b12_results,
    fits_uptake,
    b_estimates,
    blank_offset,
    group_vars
  )
  
  ############################################################
  ############ Output
  ############################################################
  
  list(
    blank_offset = blank_offset,
    release_initial_b12_estimates = release_initial_b12_estimates,
    b_estimates = b_estimates,
    fits_uptake = fits_uptake,
    fits_release = fits_release,
    mix_pred_uptake = mix_pred_uptake,
    group_vars = group_vars,
    join_vars = join_vars,
    group_by_glucose = group_by_glucose,
    group_by_replicate = group_by_replicate,
    fit_scale = fit_scale
  )
  
}
