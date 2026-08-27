# Cleaning scripts

The three preprocessors read immutable instrument exports and write only under
`data/intermediate/`:

- `b12_trait_survey_preprocess.R`
- `b12_release_uptake_death_preprocess.R`
- `growth_dynamics_b12_methionine_preprocess.R`

`download_kegg_ko_mapping.R` is the one acquisition step. It retrieves and
validates the stable `data/raw/annotations/kegg_ko_list.tsv` input required by
the classification analysis.
