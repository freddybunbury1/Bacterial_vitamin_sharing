# Run from the repository root after rebuilding the analysis and figures:
# Rscript --vanilla tests/test_od_specific_rates.R
suppressPackageStartupMessages(source("R/utils/b12_two_stage_kinetic_pipeline.R"))
section <- "data/processed/b12_release_uptake_death"
modeling <- file.path(section, "modeling")
sources <- "data/processed/figure_source_data/b12_release_uptake_death_modeling"
read_data <- function(path) read.csv(path, check.names = FALSE)
equal <- function(a, b) stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE)))
expect_error <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
keys <- c("strain_id", "experiment_name", "biological_replicate", "added_glucose_mM_C")
od <- read_data(file.path(section, "b12_od_viability_results_b12_release_uptake_death.csv"))
reference <- build_kinetic_od_reference(od)
stopifnot(nrow(reference) == 192L, sum(reference$rate_od_status == "valid") == 184L)
equal(reference$od600_reference, (reference$od600_t0 + reference$od600_t240) / 2)

# Independently match every denominator to its original measurement rows.
for (i in seq_len(nrow(reference))) {
  row <- reference[i, ]
  matched <- od$strain_id == row$strain_id & od$experiment_name == row$experiment_name &
    od$biological_replicate == row$biological_replicate
  initial <- subset(od, matched & experiment_type == "live_dead_source" &
                      sample_time_minute == 0 & killed_proportion == 0 &
                      added_glucose_mM_C == 0 & added_b12_pM == 0)$adjusted_od600
  final <- subset(od, matched & experiment_type == "release_and_uptake" &
                    sample_time_minute == 240 & added_glucose_mM_C == row$added_glucose_mM_C &
                    added_b12_pM %in% c(0, 500))$adjusted_od600
  stopifnot(length(initial) == 1L, length(final) == 2L)
  equal(row$od600_reference, (initial + mean(final)) / 2)
}

uptake <- read_data(file.path(modeling, "b12_fits_uptake_b12_release_uptake_death.csv"))
release <- read_data(file.path(modeling, "b12_fits_release_b12_release_uptake_death.csv"))
for (fits in list(uptake, release)) {
  valid <- subset(fits, strain_id != "blank")
  stopifnot(nrow(valid) == 184L, all(is.finite(valid$od600_reference)),
            all(valid$od600_reference > 0), !anyDuplicated(fits[keys]))
  equal(valid$k_per_od600 * valid$od600_reference, valid$k_per_cell)
  equal(valid$k_per_od600_se * valid$od600_reference, valid$k_per_cell_se)
  stopifnot(all(is.na(subset(fits, strain_id == "blank")$k_per_od600)))
}
valid_release <- subset(release, strain_id != "blank")
equal(valid_release$r_per_od600 * valid_release$od600_reference, valid_release$r_per_cell)
equal(valid_release$r_per_od600_se * valid_release$od600_reference, valid_release$r_per_cell_se)
equal(valid_release$k_per_od600 / valid_release$r_per_od600,
      valid_release$k_per_cell / valid_release$r_per_cell)

# Bad data must not silently switch to a one-timepoint denominator or duplicate fits.
initial_row <- which(od$strain_id != "blank" & od$experiment_type == "live_dead_source" &
                       od$killed_proportion == 0)[1]
final_row <- which(od$strain_id == od$strain_id[initial_row] &
                    od$experiment_type == "release_and_uptake" &
                    od$sample_time_minute == 240)[1]
expect_error(build_kinetic_od_reference(rbind(od, od[initial_row, ])))
expect_error(build_kinetic_od_reference(rbind(od, od[final_row, ])))
expect_error(build_kinetic_od_reference(od[-initial_row, ]))
expect_error(build_kinetic_od_reference(od[-final_row, ]))
bad <- od
bad$adjusted_od600[initial_row] <- 0
expect_error(build_kinetic_od_reference(bad))
bad$adjusted_od600[initial_row] <- NA_real_
expect_error(build_kinetic_od_reference(bad))
legacy_uptake <- uptake[c(keys, "k_per_cell", "k_per_cell_se")]
expect_error(add_od_specific_rates(legacy_uptake, subset(reference, strain_id == "blank")))
equal(add_od_specific_rates(legacy_uptake, reference)$k_per_od600, uptake$k_per_od600)

# Figure summaries must average replicate-normalized rates, not divide strain means.
geomean <- function(x) exp(mean(log(x)))
summary <- aggregate(cbind(k_per_od600, r_per_od600) ~ strain_id + added_glucose_mM_C,
                     valid_release, geomean)
k_source <- read_data(file.path(sources, "fit_parameter_columns_source_03C.csv"))
r_source <- read_data(file.path(sources, "fit_parameter_columns_source_03H.csv"))
k_check <- merge(summary, k_source, by = c("strain_id", "added_glucose_mM_C"))
r_check <- merge(summary, r_source, by = c("strain_id", "added_glucose_mM_C"))
stopifnot(nrow(k_check) == 46L, nrow(r_check) == 46L,
          all(k_check$n_parameter_fits == 4), all(r_check$n_parameter_fits == 4))
equal(k_check$k_per_od600, k_check$estimate)
equal(r_check$r_per_od600, r_check$estimate)

s16 <- read_data(file.path(sources, "FigureS16_first_order_rates_source.csv"))
fig4b <- read_data(file.path(sources, "fit_parameter_columns_source_03E.csv"))
fig5c <- read_data(file.path(sources, "fit_parameter_columns_source_03N.csv"))
shared <- Reduce(function(a, b) merge(a, b, by = "strain_id"), list(s16, fig4b, fig5c))
stopifnot(nrow(shared) == 23L)
equal(shared$uptake_rate_k, shared$estimate)
equal(shared$uptake_rate_k, shared$k)
equal(shared$release_rate_r, shared$r)
ratio <- read_data(file.path(sources, "ec_b12_fraction_per_dead_cell_fraction_vs_k_over_r_source_09B.csv"))
shared_ratio <- merge(fig5c, ratio[c("strain_id", "k_over_r")], by = "strain_id")
stopifnot(nrow(shared_ratio) == 23L)
equal(shared_ratio$k / shared_ratio$r, shared_ratio$k_over_r)
cat("PASS: OD matching, 184 normalized culture fits, invalid-input guards, and agreement across Figures 4B, 5B-D and S16.\n")
