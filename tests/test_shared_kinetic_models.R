# Run after rebuilding from the release root.
suppressPackageStartupMessages(source('R/utils/b12_model_structure_comparison.R'))
src <- 'data/processed/figure_source_data/b12_release_uptake_death_modeling'
model <- 'data/processed/b12_release_uptake_death/modeling'
stopifnot(!exists('predict_mix_release', mode='function'),
          !file.exists(file.path(model,'b12_mix_pred_release_b12_release_uptake_death.csv')),
          !file.exists(file.path(model,'b12_with_prop_b12_release_uptake_death.csv')),
          !file.exists(file.path(src,'dead_fraction_mix_source_04M_04N_04O.csv')))
read_data <- function(path) read.csv(path, check.names=FALSE)
equal <- function(x,y) stopifnot(identical(as.numeric(x),as.numeric(y)))
check_pair <- function(main_file,supp_file,keys,fields) {
  a <- read_data(file.path(src,main_file))
  b <- subset(read_data(file.path(src,supp_file)),model_structure=='first_order')
  stopifnot(!anyDuplicated(a[keys]),!anyDuplicated(b[keys]))
  z <- merge(a,b,by=keys,suffixes=c('.main','.supp'))
  stopifnot(nrow(z)==nrow(a))
  for (f in fields) {
    x <- if(f[1] %in% names(z)) f[1] else paste0(f[1],'.main')
    y <- if(f[2] %in% names(z)) f[2] else paste0(f[2],'.supp')
    equal(z[[x]],z[[y]])
  }
  cat('PASS:',main_file,'is an exact subset of',supp_file,'at',nrow(z),'coordinates.\n')
}
check_pair('uptake_fit_dynamics_source_01H.csv','supp_model_structures_uptake_dynamics.csv',
           c('strain_id','sample_time_minute'),
           list(c('b12_mean','empirical_mean_b12'),c('b12_sd','empirical_sd_b12'),c('model_pred','model_pred')))
check_pair('dead_fraction_mix_source_04L.csv','supp_model_structures_live_dead_mix.csv',
           c('strain_id','killed_proportion'),
           list(c('b12_frac_mean','empirical_fraction_mean'),c('b12_frac_sd','empirical_fraction_sd'),c('predicted_b12_frac','predicted_fraction')))
check_pair('release_uptake_fit_dynamics_source_02K.csv','supp_model_structures_release_dynamics.csv',
           c('strain_id','added_glucose_mM_C','sample_time_minute'),
           list(c('b12_mean','empirical_mean_b12'),c('b12_sd','empirical_sd_b12'),c('model_pred','model_pred')))
# Full first-order time-course summaries, beyond the three selected strains.
for (pair in list(c('uptake_fit_dynamics_source_01K.csv','supp_model_structures_uptake_dynamics.csv'),
                  c('release_uptake_fit_dynamics_source_02A.csv','supp_model_structures_release_dynamics.csv'))) {
  a <- subset(read_data(file.path(src,pair[1])),strain_id!='blank' & added_glucose_mM_C %in% c(0,10))
  if(grepl('uptake_fit_dynamics_source',pair[1])) a <- subset(a,added_glucose_mM_C==0)
  b <- subset(read_data(file.path(src,pair[2])),model_structure=='first_order' & strain_id!='blank')
  keys <- c('strain_id','sample_time_minute')
  if('added_glucose_mM_C' %in% names(b)) keys <- c(keys,'added_glucose_mM_C')
  z <- merge(a,b,by=keys,suffixes=c('.main','.supp'))
  stopifnot(nrow(z)==nrow(a));equal(z$model_pred.main,z$model_pred.supp)
}
keys <- c('strain_id','added_glucose_mM_C','experiment_name','biological_replicate')
fits <- subset(read_data(file.path(model,'b12_fits_release_b12_release_uptake_death.csv')),strain_id!='blank')
parameters <- read_data(file.path(model,'supp_first_order_parameters.csv'))
for (f in list(c('k_per_od','k_per_od600'),c('r_per_od','r_per_od600'),c('E_pM','E_per_cell'))) {
  z <- merge(fits,subset(parameters,parameter==f[1]),by=keys)
  stopifnot(nrow(z)==nrow(fits));equal(z[[f[2]]],z$estimate)
}
status <- read_data(file.path(model,'supp_first_order_status.csv'))
stopifnot(all(status$fit_message=='reused_canonical_first_order'))
mix <- subset(read_data(file.path(model,'b12_mix_pred_uptake_b12_release_uptake_death.csv')),
              strain_id!='blank' & is.finite(predicted_b12))
supp <- read_data(file.path(model,'supp_first_order_predictions.csv'))
z <- merge(mix,supp,by=c('experiment_name','sample_id'),suffixes=c('.main','.supp'))
stopifnot(nrow(z)==nrow(mix));equal(z$predicted_b12.main,z$predicted_b12.supp)
stopifnot(all(fits$k_per_cell>=1e-4*(1-1e-10)),all(fits$r_per_cell>=1e-4*(1-1e-10)))
inputs <- read_model_inputs('.')
quota <- read_data(file.path(model,'supp_quota_parameters.csv'))
reference <- read_data(file.path(model,'b12_rate_od_reference_b12_release_uptake_death.csv'))
z <- merge(quota,reference,by=keys)
stopifnot(all(subset(z,parameter=='k_per_od')$estimate>=1e-8*(1-1e-10)),
          all(subset(z,parameter=='r_per_od')$estimate>=exp(-20)*(1-1e-10)),
          all(with(subset(z,parameter=='Qmax_per_od'),estimate*od600_reference>=inputs$quota_blank_sd_pM*(1-1e-10))))
kinetic <- read_data(file.path(model,'supp_quota_kinetic_predictions.csv'))
z <- merge(kinetic,reference,by=keys)
equal(z$adjusted_od600,z$od600_reference)
cat('PASS: full time-course summaries, reused coefficients and mixture predictions, shared quota OD, and retained lower bounds.\n')

# Independently evaluate the saved scalar fits against screened alternatives.
# This catches the previously observed worse local solutions without refitting
# either figure workflow or altering the production coefficients.
data <- read_data('data/processed/b12_release_uptake_death/b12_results_b12_release_uptake_death.csv')
uptake <- subset(read_data(file.path(model,'b12_fits_uptake_b12_release_uptake_death.csv')),strain_id!='blank')
initial <- read_data(file.path(model,'b12_release_initial_b12_estimates_b12_release_uptake_death.csv'))
release <- merge(fits,initial,by=keys)
stopifnot(nrow(uptake)==184L,nrow(release)==184L)
checks <- list()
for (setting in c('uptake','release')) {
  table <- if(setting=='uptake') uptake else release
  for(i in seq_len(nrow(table))) {
    f <- table[i,]
    g <- data
    for(k in keys) g <- g[g[[k]]==f[[k]],]
    g <- subset(g,experiment_type=='release_and_uptake' & is.finite(b12_pm) & b12_pm>0)
    g <- if(setting=='uptake') subset(g,added_b12_pM>0) else subset(g,added_b12_pM==0)
    objective <- function(log_rate) {
      rate <- exp(log_rate)
      pred <- if(setting=='uptake') inputs$blank_offset+g$added_b12_pM*exp(-rate*g$sample_time_minute) else
        inputs$blank_offset+rate/f$k_per_cell+(f$release_initial_b12_pM-rate/f$k_per_cell)*exp(-f$k_per_cell*g$sample_time_minute)
      sum((log(g$b12_pm)-log(pred))^2)
    }
    saved <- objective(log(if(setting=='uptake') f$k_per_cell else f$r_per_cell))
    screened <- screen_scalar_fit(objective,1e-4,if(setting=='uptake') 1e3 else 1e4)$objective
    stopifnot(saved <= screened+1e-6*max(1,screened))
    checks[[length(checks)+1L]] <- data.frame(f[keys],setting,saved_objective=saved,
       screened_objective=screened,objective_gap=saved-screened)
  }
}
write.csv(bind_rows(checks),file.path(model,'canonical_scalar_optimizer_checks.csv'),row.names=FALSE)
cat('PASS: 368 canonical uptake/release fits are no worse than independently screened scalar alternatives.\n')
