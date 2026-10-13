# Run from the repository root after the survey/classification analyses.
source('R/utils/b12_figure_cohort.R')
summary <- read.csv('data/processed/b12_trait_survey/b12_trait_survey_strain_summary.csv')
cohort <- b12_figure_cohort(summary)
tree <- ape::read.tree('data/processed/gtdbtk_tree_analysis/gtdbtk_bac120_study_genomes_bacillati_pseudomonadati_split_rooted_b12_complete_traits.tree')
model <- read.csv('data/processed/b12_trait_classification_analysis/b12_trait_classification_model_input.csv')
stopifnot(nrow(cohort)==202L, !anyDuplicated(cohort$strain_id),
          setequal(cohort$strain_id,tree$tip.label),
          setequal(cohort$strain_id,model$strain_id),
          sum(model$strain_id=='KC611216')==1L,
          !'KC61126' %in% model$strain_id,
          sum(cohort$total_b12_gm>50)==79L,
          any(cohort$n_samples==1L),
          any(!is.finite(cohort$total_b12_gsd)))
# Sample count and unavailable uncertainty do not exclude otherwise valid data.
example <- data.frame(strain_id=c('KC61126','outside','JFK3','sie1401'),
                      total_b12_gm=c(100,100,NA,100),
                      uptake_b12_mean=c(1,1,1,Inf),n_samples=1L,
                      total_b12_gsd=NA_real_)
stopifnot(identical(b12_figure_cohort(example)$strain_id,'KC61126'))
source_dir <- 'data/processed/figure_source_data'
paths <- c('b12_trait_survey/15B_dead_cell_fraction_vs_ec_b12_fraction_trait_group2_source.csv',
           'b12_trait_survey/15C_ec_b12_fraction_over_dead_cell_fraction_vs_uptake_trait_group2_source.csv',
           'b12_trait_survey/12D_trait_group2_boxplot_source.csv')
for(p in paths) {
 data <- read.csv(file.path(source_dir,p))
 stopifnot(nrow(data)==79L,!anyDuplicated(data$strain_id),
           setequal(data$strain_id,cohort$strain_id[cohort$total_b12_gm>50]))
}
folds <- read.csv('data/processed/b12_trait_classification_analysis/random_fold_assignments.csv')
stopifnot(!anyDuplicated(folds[c('trait','strain_id')]))
for(target in unique(folds$trait)) {
 expected <- if(grepl('total_b12_gt50$',target)) cohort$strain_id[cohort$total_b12_gm>50] else cohort$strain_id
 stopifnot(setequal(folds$strain_id[folds$trait==target],expected))
}
cat('Shared cohort, KC identifier, synthesizer sources and fold coverage checks passed.\n')
