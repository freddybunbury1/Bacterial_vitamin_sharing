# Cleaning scripts

The three preprocessors read immutable instrument exports and write only under
`data/intermediate/`:

- `b12_trait_survey_preprocess.R`
- `b12_release_uptake_death_preprocess.R`
- `growth_dynamics_b12_methionine_preprocess.R`

OD600 is normalized to a 1 cm path length during preprocessing, before blank
subtraction. Growth dynamics and both initial/final release/uptake readings use
200 uL in 96-well plates (35 mm^2 well area): `od_path_length_cm = 200 / 35 / 10`,
giving a multiplier of 1.75. The survey uses 45 uL in 384-well plates
(11.35 mm^2 well area): `od_path_length_cm = 45 / 11.35 / 10`.
Intermediate OD tables retain the original values as `od600_plate_read` and
record `od_path_length_cm`; existing analysis-facing OD columns contain the
1 cm-normalized values.

`download_kegg_ko_mapping.R` is the one acquisition step. It retrieves and
validates the stable `data/raw/annotations/kegg_ko_list.tsv` input required by
the classification analysis.
