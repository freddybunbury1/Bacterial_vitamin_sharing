# Utilities

`b12_figure_cohort.R` defines the shared QC-passing genome and finite-phenotype
cohort and canonicalizes the legacy KC strain identifier for joins. See
`docs/PIPELINE.md` for the 202-strain cohort and 79-strain synthesizer subset.

Shared pipeline functions belong here. The release/uptake pipeline retains
`release_initial_b12_pM` only for release-plus-uptake fits and predictions;
uptake-only fits start from nominal added B12. The model-structure comparison
utility provides focused first-order and quota fits for Figures S12–S14.
These supplementary fits use explicit OD-specific parameters and include
release in the live–dead predictions. The canonical first-order pipeline
fits culture-level coefficients (despite legacy `_per_cell` column names)
and provides the uptake-only live–dead predictions used in the main figures.
Its OD-reference helper matches initial untreated-source OD to final kinetic
OD (averaged across the 0- and 500-pM wells), then divides each replicate's
coefficients by the mean of initial and final OD. Explicit `k_per_od600` and
`r_per_od600` fields feed the main rate summaries; the legacy fields retain
the culture-level values required for unchanged concentration predictions.

Quota uptake-only trajectories use their analytical solution. Quota release
and live–dead trajectories use deSolve LSODA (relative and absolute tolerances
1e-8); fitted trajectories are checked at 1e-10. Nonnegativity, capacity and
mass conservation are checked on each trajectory. Failed quota L-BFGS-B fits
are retried with bounded nlminb using the same objective and parameter bounds.
The normal analysis exports supplementary fit status, parameters, predictions,
metrics and quota numerical checks under the modeling directory.
