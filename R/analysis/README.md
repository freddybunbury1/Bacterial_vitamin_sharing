# Analysis scripts

- `b12_trait_survey_analysis.R`: estimates B12 concentrations, joins survey
  measurements, classifies strain traits, and exports survey figure sources.
- `b12_release_uptake_death_analysis.R`: produces empirical summaries, fits the
  single canonical two-stage kinetic model, and exports empirical/model figure
  sources.
- `growth_dynamics_b12_methionine_analysis.R`: models blank OD, summarizes
  growth curves, and exports the Figure S2 comparison.
- `gtdbtk_tree_analysis.R`: parses GTDB taxonomy and builds the split-rooted
  study-genome and complete-trait trees.
- `b12_trait_classification_analysis.R`: constructs the KO/phenotype input,
  aggregates GWA results, runs single-stage order-block and matched random-fold
  validation, null simulations, and feature-importance summaries.

Classification random folds are constructed once per trait, sampled without
replacement, mutually exclusive, and exhaustive. Their fold count equals the
corresponding number of GTDB order blocks (7 for full-tree traits and 6 for the
total-B12-filtered traits).
