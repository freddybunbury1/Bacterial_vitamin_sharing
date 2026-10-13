# Pipeline dependency map

`R/compiled_figures/compiled_figures.Rmd` is the public entry point and executes
the following fixed order:

1. `b12_trait_survey_preprocess.R`
2. `b12_release_uptake_death_preprocess.R`
3. `growth_dynamics_b12_methionine_preprocess.R`
4. `b12_trait_survey_analysis.R`
5. `gtdbtk_tree_analysis.R`
6. `growth_dynamics_b12_methionine_analysis.R`
7. `b12_release_uptake_death_analysis.R`
8. `b12_trait_classification_analysis.R`
9. export Tables S1–S3, preprocess figure panels, and assemble Figures 1–5 and S1–S19

## Analysis dependencies

| Analysis | Principal inputs | Principal outputs |
| --- | --- | --- |
| Trait survey | survey intermediates, canonical strain metadata | strain traits, GWA phenotypes, survey figure sources |
| GTDB tree | GTDB summary/tree, strain-to-genome map, trait summary | public taxonomy, study-genome tree, complete-trait tree, tip metadata |
| Growth dynamics | growth intermediate | growth metrics and Figure S3 source |
| Release/uptake/death | release intermediates, trait summary, GTDB taxonomy | empirical tables, uptake/release model fits, compact figure sources |
| Trait classification | KO matrix, traits, GTDB taxonomy/tree, KEGG mapping, raw GWA summaries | matched validation results, nulls, importance, combined GWA table, compact figure sources |

The classification analysis contains the GWA aggregation formerly performed by
`genome_annotations_and_gwa_analysis.R`; the latter is not part of the public
pipeline. GTDB taxonomy and tree processing are both owned by
`gtdbtk_tree_analysis.R`. Empirical and kinetic release analyses are both owned
by `b12_release_uptake_death_analysis.R`.

The canonical first-order uptake model starts at the nominal added-B12
concentration and therefore has no fitted initial-condition offset. The
release-plus-uptake model retains the release-specific
`release_initial_b12_pM`, estimated from the earliest zero-added-B12
observation for each fitting group. Figure 5C–D summaries use only the
zero-added-glucose fits and combine replicate-specific positive parameters by
geometric mean. Figures S12–S14 compare first-order (red) and quota (blue)
model structures. Figure S16 prunes the GTDB tree to the 23 selected
release/uptake strains and displays broad environment, survey total B12,
uptake, and EC-B12 fraction, trait group, and glucose-zero first-order release
and uptake rates.
Each main-figure dynamic rate is normalized per replicate to the mean of
initial and final blank-adjusted OD600, then summarized by geometric mean.

The primary first-order implementation estimates culture-level coefficients,
despite the legacy names `k_per_cell` and `r_per_cell`; its equations can be
written equivalently as K = k x and R = r x for constant culture OD x.
The fitted products are retained for predictions. Figures 4B, 5B–C and S16
display `k_per_od600 = K / x` and `r_per_od600 = R / x`. Initial OD is the
untreated live–dead source (0 min, killed proportion 0); final OD is the
arithmetic mean of the matched 0- and 500-pM kinetic wells at 240 min.
Matching uses strain, experiment and biological replicate, with final OD
also matched by glucose condition. Both time points receive equal weight.
`b12_rate_od_reference_b12_release_uptake_death.csv` records these denominators;
the uptake/release fit exports include the OD-specific coefficients and
conditional SEs (OD measurement uncertainty is not propagated). Missing or
nonpositive non-blank OD stops the calculation; blank controls are not
normalized. The common OD denominator cancels in the Figure 5D k/r ratio.
First-order fits are performed once in the canonical two-stage pipeline;
the supplementary comparison reuses its coefficients and replicate predictions.
Main panels 4A, 4C and 5A select strains from the same summarized first-order
results shown in S12, S13 and S14, respectively. Quota fits also use the mean
initial/final OD reference above, including x times the surviving fraction
for live–dead predictions. S13 omits ongoing release and the release-specific
initial concentration, predicting uptake alone after instantaneous killed-cell
release. S14 retains the independently fitted release coefficient and initial
condition. Quota uptake depends on remaining capacity,
with Q(0) = 0 and the total capacity constrained to at least one standard
deviation of the 500-pM blank-control measurements.

Figures S12 and S13 use the glucose-zero treatment; Figure S14 includes 0 and
10 mM carbon from glucose. Black points and error bars show empirical means
and SDs. Coloured lines join means of separately fitted replicate predictions
at the sampled time points, with alpha 0.8. Facets are ordered blank, Retainers,
Providers, Reclaimers, using the same Table S3 ordinals as Figure S16. All
facets share continuous axes, and Figure S13 has a y-axis upper limit of 1.5.
Canonical scalar first-order fitting screens log-spaced candidate rates and
refines sampled local minima before least-squares refinement. It retains the
existing culture-effective lower bounds of 1e-4 for uptake and release.
Quota release fitting also screens scalar log-parameter minima. The existing
quota uptake bounds and analytical/adaptive solver checks are retained.
Run `Rscript --vanilla tests/test_shared_kinetic_models.R` after rebuilding to
check exact main/supplementary source equality and canonical fit provenance.
Figures S17–S19 contain candidate-gene associations, bioassay comparison, and
survey fluorescence distributions. Manuscript Methods and captions are maintained separately in Overleaf.

All paths are relative to the repository root. The R Markdown entry point finds
that root from its own location, so it can be rendered while the caller's
working directory is elsewhere.

## Optional upstream workflows

Table S1 joins the canonical bacterial strain metadata to GTDB-Tk taxonomy
from phylum through genus,
Table S2 is generated from the canonical media-recipe CSV, and Table S3 joins
the selected release/uptake/death culture wells to GTDB-Tk taxonomy and the
canonical figure-panel strain labels. The workflows
under `upstream/` document how coauthor analyses produced the
fixed genome, GTDB-Tk, and Pyseer inputs under `data/raw/`. They are deliberately
separate from the manuscript R pipeline and are not sourced by
`compiled_figures.Rmd`. Their local `output/` directories are ignored by Git;
replacing a fixed manuscript input is an explicit, reviewed action.

## Shared genome-linked phenotype cohort

The trait panels in Figures 1B--C, 3, S4--S5 and S10 use the same eligibility
as Figure 1A: a genome in the QC-passing study-genome tree and finite total
B12 and uptake estimates. `R/utils/b12_figure_cohort.R` supplies this shared
filter. There is no minimum sample-count or finite-SD requirement for a point;
unavailable uncertainty omits the corresponding error bar. The full cohort
contains 202 unique strains, and the `total_b12_gm > 50` subset contains 79.
Required positive extracellular/death values are additionally checked where
logarithms or ratios are plotted; these checks retain all 79 synthesizers.

The classifier uses the same cohort. The legacy genomic-feature identifier
`KC61126` is matched to `KC611216` before joining, retaining that strain once.
The four tasks in Figure 2 use 202, 202, 79 and 79 strains; the S6 tree panels
use 202 and 79. S7 is recalculated from the resulting random forests.
Supplied pyseer association summaries are not refitted or modified.

## Figure S9 measurements

Panel A reads raw OD600 and raw heat-killed DCS1 fluorescence (`dead_stain`)
from `b12_trait_survey_results.csv`. It retains 1,152 finite, positive pairs,
including 12 blanks, independently of the genome-linked phenotype cohort.
Blanks are drawn last in red; other samples are black. Marginal densities
are fitted in log10 coordinates and normalized separately for the two groups.

Panel B reads the live–dead mixture source table and matches source-culture
OD at time zero and the 0%-killed fluorescence baseline by strain, experiment
and biological replicate. Its x value is killed fraction × blank-adjusted
source OD; its y value is raw mixture DCS1 fluorescence (`live_stain`) minus
the matched 0%-killed value. The pre-filter table retains all differences.
The plotted subset excludes blanks, 0%-killed baselines and nonpositive axis
values, leaving 642 measurements, 92 replicate groups and 23 strains.
Duplicate match keys or missing baselines/source OD stop the calculation.
These transformations and the two-panel assembly are owned by the compiled
notebook; the mixture source table remains the input for baseline subtraction.

## Release checks

The commands under "Check an existing release" in the root README validate
the cohort, OD normalization, immutable inputs, generated inventory and
notebook/source consistency without rerunning the pipeline. Figure layout
settings, including the S1 and S9 panel resizing, are documented in
`R/compiled_figures/README.md`. Historical comparisons and temporary outputs
belong outside the release folder.
