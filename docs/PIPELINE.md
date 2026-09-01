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
9. export Tables S1–S2, preprocess figure panels, and assemble Figures 1–5 and S1–S17

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
