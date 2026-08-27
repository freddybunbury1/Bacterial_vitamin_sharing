############################################################
############ B12 Two-Stage Kinetic Model Pipeline
############################################################

library(minpack.lm)
library(dplyr)
library(broom)
library(tidyr)

############################################################
############ Uptake-only model (estimate k on B0 > 0)
############################################################

fit_uptake_model <- function(b12_results, delta_estimates, blank_offset, group_vars,
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
    dplyr::filter(!is.na(y), !is.na(x)) %>%
    dplyr::left_join(delta_estimates, by = group_vars) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::group_modify(~{
      
      dat <- .x
      if(nrow(dat) < 4) return(data.frame())
      
      x_max <- max(dat$x)
      dat$x_scaled <- dat$x / x_max
      
      fit <- tryCatch({
        k_scaled_lower <- 1e-4 * x_max
        
        if (fit_scale == "log") {
          nlsLM(
            log(y) ~ log(blank_offset +
                           (B0_nominal + delta) *
                           exp(-k_scaled * x_scaled)),
            data = dat,
            start = list(k_scaled = 0.1),
            lower = c(k_scaled = k_scaled_lower),
            control = nls.lm.control(maxiter = 500)
          )
        } else {
          nlsLM(
            y ~ blank_offset +
              (B0_nominal + delta) *
              exp(-k_scaled * x_scaled),
            data = dat,
            start = list(k_scaled = 0.1),
            lower = c(k_scaled = k_scaled_lower),
            control = nls.lm.control(maxiter = 500)
          )
        }
        
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
fit_release_model <- function(b12_results, delta_estimates, fits_uptake,
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
    dplyr::filter(!is.na(y), !is.na(x)) %>%
    dplyr::left_join(delta_estimates, by = group_vars) %>%
    dplyr::left_join(fits_uptake, by = group_vars) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::group_modify(~{
      
      dat <- .x
      if(nrow(dat) < 3) return(data.frame())
      
      k_val <- unique(dat$k_per_cell)
      if(length(k_val) != 1 || is.na(k_val)) return(data.frame())
      
      fit <- tryCatch({
        
        if (fit_scale == "log") {
          nlsLM(
            log(y) ~ log(
              blank_offset +
                (r_per_cell / k_val) +
                (delta - r_per_cell / k_val) *
                exp(-k_val * x)
            ),
            data = dat,
            start = list(r_per_cell = 0.1),
            lower = c(r_per_cell = 1e-4),
            control = nls.lm.control(maxiter = 500)
          )
        } else {
          nlsLM(
            y ~ blank_offset +
              (r_per_cell / k_val) +
              (delta - r_per_cell / k_val) *
              exp(-k_val * x),
            data = dat,
            start = list(r_per_cell = 0.1),
            lower = c(r_per_cell = 1e-4),
            control = nls.lm.control(maxiter = 500)
          )
        }
        
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

predict_mix_uptake <- function(b12_results, fits_uptake, delta_estimates,
                               b_estimates, blank_offset, group_vars){
  
  mix_pred_uptake <- b12_results %>%
    dplyr::filter(experiment_type == "live_dead_mix") %>%
    dplyr::left_join(fits_uptake, by = group_vars) %>%
    dplyr::left_join(delta_estimates, by = group_vars) %>%
    dplyr::left_join(b_estimates, by = group_vars) %>%
    dplyr::mutate(
      f_dead = killed_proportion,
      f_live = 1 - killed_proportion,
      B0 = delta + b_per_cell * f_dead,
      predicted_b12 = blank_offset +
        B0 * exp(-k_per_cell *
                   f_live *
                   sample_time_minute),
      null_predicted_b12 = blank_offset + B0
    )
  
  return(mix_pred_uptake)
  
}

############################################################
############ Mix predictions (release+uptake)
############################################################

predict_mix_release <- function(b12_results, fits_release, delta_estimates,
                                b_estimates, blank_offset, group_vars){
  
  mix_pred_release <- b12_results %>%
    dplyr::filter(experiment_type == "live_dead_mix") %>%
    dplyr::left_join(fits_release, by = group_vars) %>%
    dplyr::left_join(delta_estimates, by = group_vars) %>%
    dplyr::left_join(b_estimates, by = group_vars) %>%
    dplyr::mutate(
      f_dead = killed_proportion,
      f_live = 1 - killed_proportion,
      B0 = delta + b_per_cell * f_dead,
      predicted_b12 = blank_offset +
        E_per_cell +
        (B0 - E_per_cell) *
        exp(-k_per_cell *
              f_live *
              sample_time_minute),
      null_predicted_b12 = blank_offset + B0
    )
  
  return(mix_pred_release)
  
}

############################################################
############ Master pipeline (same interface)
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
  ############ Delta estimate
  ############################################################
  
  delta_estimates <- b12_results %>%
    dplyr::filter(experiment_type == "release_and_uptake",
                  added_b12_pM == 0) %>%
    dplyr::mutate(y = b12_pm,
                  x = sample_time_minute,
                  B0_nominal = added_b12_pM) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::slice_min(x, n = 1, with_ties = TRUE) %>%
    dplyr::summarise(
      delta = mean(y - blank_offset - B0_nominal,
                   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(delta = pmax(delta, 0))
  
  ############################################################
  ############ Dead-cell release estimate
  ############################################################
  
  b_estimates <- b12_results %>%
    dplyr::filter(experiment_type == "live_dead_mix",
                  killed_proportion == 1) %>%
    dplyr::mutate(b12_full_corrected = b12_pm - blank_offset) %>%
    dplyr::left_join(delta_estimates, by = group_vars) %>%
    dplyr::group_by(across(all_of(group_vars))) %>%
    dplyr::summarise(
      b_full_corrected = mean(b12_full_corrected, na.rm = TRUE),
      b_per_cell = mean(b12_full_corrected - delta, na.rm = TRUE),
      .groups = "drop"
    )
  
  ############################################################
  ############ Fit models
  ############################################################
  
  fits_uptake <- fit_uptake_model(
    b12_results,
    delta_estimates,
    blank_offset,
    group_vars,
    fit_scale = fit_scale
  )
  
  fits_release <- fit_release_model(
    b12_results,
    delta_estimates,
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
    delta_estimates,
    b_estimates,
    blank_offset,
    group_vars
  )
  
  mix_pred_release <- predict_mix_release(
    b12_results,
    fits_release,
    delta_estimates,
    b_estimates,
    blank_offset,
    group_vars
  )
  
  ############################################################
  ############ Output (identical structure)
  ############################################################
  
  list(
    blank_offset = blank_offset,
    delta_estimates = delta_estimates,
    b_estimates = b_estimates,
    fits_uptake = fits_uptake,
    fits_release = fits_release,
    mix_pred_uptake = mix_pred_uptake,
    mix_pred_release = mix_pred_release,
    group_vars = group_vars,
    join_vars = join_vars,
    group_by_glucose = group_by_glucose,
    group_by_replicate = group_by_replicate,
    fit_scale = fit_scale
  )
  
}
