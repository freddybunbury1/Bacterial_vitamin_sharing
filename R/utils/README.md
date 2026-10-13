# Utilities

`b12_figure_cohort.R` defines the shared QC-passing genome and finite-phenotype
cohort and canonicalizes the legacy KC strain identifier for joins. See
`docs/PIPELINE.md` for the 202-strain cohort and 79-strain synthesizer subset.

Shared pipeline functions belong here. The release/uptake pipeline retains
`release_initial_b12_pM` only for release-plus-uptake fits and predictions;
uptake-only fits start from nominal added B12. The model-structure comparison
utility fits the quota comparison and reuses the canonical first-order fits
and predictions for Figures S12–S14; it never refits the first-order comparator.
Both structures use the same mean initial/final OD reference. S13 predicts
uptake only after the initial killed-cell release; S14 retains ongoing release. The canonical first-order pipeline
fits culture-level coefficients (despite legacy `_per_cell` column names)
and provides the uptake-only live–dead predictions used in the main figures.
Its OD-reference helper matches initial untreated-source OD to final kinetic
OD (averaged across the 0- and 500-pM wells), then divides each replicate's
coefficients by the mean of initial and final OD. Explicit `k_per_od600` and
`r_per_od600` fields feed the main rate summaries; the legacy fields retain
the culture-level values used for concentration predictions.
Main 4A, 4C and 5A are strain subsets of the shared first-order summaries in
S12, S13 and S14. Scalar fitting screens candidate minima before refinement,
retaining the original first-order rate floors; accepted solutions must be
no worse than the screened objective within numerical tolerance.

Quota uptake-only trajectories use their analytical solution. Quota release
trajectories use deSolve LSODA (relative and absolute tolerances
1e-8); fitted trajectories are checked at 1e-10. Nonnegativity, capacity and
mass conservation are checked on each trajectory. Failed quota L-BFGS-B fits
are retried with bounded nlminb using the same objective and parameter bounds.
The normal analysis exports supplementary fit status, parameters, predictions,
metrics and quota numerical checks under the modeling directory.
