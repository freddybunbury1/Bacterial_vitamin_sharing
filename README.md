# Bacterial vitamin B12 sharing

This repository contains the public raw-data-to-figures pipeline for the
bacterial vitamin B12 sharing manuscript. With one documented KEGG acquisition
step, the repository rebuilds all empirical summaries, kinetic models, genome
classifiers, figure-source tables, and compiled manuscript figures from
`data/raw/`.

## Reproduce all figures

Run these commands from the repository root:

```bash
Rscript --vanilla R/cleaning/download_kegg_ko_mapping.R
Rscript --vanilla -e 'dir.create("results/compiled_figures", recursive = TRUE, showWarnings = FALSE); rmarkdown::render("R/compiled_figures/compiled_figures.Rmd", output_dir = "results/compiled_figures", envir = new.env(parent = globalenv()), clean = TRUE)'
```

The first command retrieves the current KO-to-description table directly from
the KEGG REST API. KEGG asks users to observe its academic-use conditions. The
downloaded files are deliberately ignored by Git. The second command runs the
complete pipeline in dependency order and writes the figures, supplementary
Table S1, Table S2, and Table S3 workbooks, and a self-contained HTML notebook under
`results/compiled_figures/`.

The compiled figure set is Figures 1–5 and S1–S19. Figure S11 compares initial
and four-hour uptake-assay OD600 without added glucose. Figures S12–S14 compare
first-order (red) and quota (blue) kinetic models; Figure S16 presents the
selected-strain phylogeny with first-order release and uptake rates.
Figures 4B, 5B–C and S16 report OD-specific coefficients, normalized per
replicate to the mean of initial and final blank-adjusted OD600 before
geometric averaging. Panels 4A, 4C and 5A reuse the same first-order fits,
predictions and summaries as S12, S13 and S14, respectively. S13 and 4C
include instantaneous killed-cell release followed by uptake, without ongoing
release; S14 and 5A include ongoing release.
Figures S17–S19 contain the candidate-gene associations, bioassay comparison,
and survey fluorescence distributions, respectively. Manuscript Methods and figure captions are maintained separately in Overleaf.

R Markdown requires Pandoc. Running from RStudio uses its bundled Pandoc; a
command-line installation may need Pandoc added to `PATH`.

On macOS, command-line R can use the copy bundled with RStudio by setting:

```bash
export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
```

The path above is for Apple Silicon; on an Intel Mac, replace `aarch64` with
`x86_64`.

A clean end-to-end render took approximately 3–5 minutes on an Apple Silicon
MacBook Pro; runtime will vary with hardware and available cores.

## R dependencies

The pipeline was validated with R 4.5.1. Required CRAN packages are:

```r
install.packages(c(
  "ape", "broom", "cowplot", "deSolve", "dplyr", "drc", "ggExtra", "ggnewscale", "ggplot2",
  "ggrepel", "gtable", "knitr", "lme4", "lubridate", "magick", "maps",
  "minpack.lm", "openxlsx", "patchwork", "phytools", "purrr", "ranger", "rlang",
  "rmarkdown", "scales", "stringr", "tibble", "tidyr"
))
```

## Repository layout

- `data/raw/`: version-controlled instrument exports, curated metadata, and
  primary genome-analysis outputs.
- `assets/figure_source/`: manually authored, version-controlled figure input.
- `R/cleaning/`: the three raw-data preprocessors and the optional KEGG
  acquisition script.
- `R/analysis/`: five canonical analyses that generate processed and
  figure-source tables.
- `R/utils/`: shared assay and kinetic-model functions.
- `R/compiled_figures/compiled_figures.Rmd`: the only production figure
  assembly file and the public end-to-end entry point.
- `upstream/`: optional coauthor workflows documenting how genome annotation,
  GTDB-Tk, and Pyseer outputs under `data/raw/` were generated.
- `data/intermediate/`, `data/processed/`, and `results/`: reproducible,
  Git-ignored outputs, including the journal-ready `TableS1.xlsx`,
  `TableS2.xlsx`, and `TableS3.xlsx` workbooks.

See [docs/PIPELINE.md](docs/PIPELINE.md) for dependencies and
[docs/REPRODUCIBILITY_AUDIT.md](docs/REPRODUCIBILITY_AUDIT.md) for validation.

## Check an existing release

After generating the analysis outputs and notebook, run these checks from the
repository root. They do not refit models or regenerate figures:

```bash
Rscript --vanilla tests/test_figure_cohort.R
Rscript --vanilla tests/test_od_specific_rates.R
Rscript --vanilla tests/test_quota_solver.R
Rscript --vanilla tests/test_shared_kinetic_models.R
python3 tests/check_release_integrity.py
```

The R checks validate cohort membership, OD-specific rates, quota solver
accuracy, conservation and limiting cases, and exact agreement between main
and supplementary kinetic-model results. The shared-model check also exports
scalar optimizer diagnostics without changing the fitted coefficients. The Python 3
check uses only its standard library to verify input-manifest checksums,
the compiled-output inventory, notebook/source consistency, and to flag temporary or exploratory files inside the release.
Current figure layouts and S9 transformations are documented in
[R/compiled_figures/README.md](R/compiled_figures/README.md).

## Licensing

Copyright (c) 2026 Freddy Bunbury and contributors.

Source code is available under the [MIT License](LICENSE). Data,
documentation, and figure assets are available under the
[Creative Commons Attribution 4.0 International License](LICENSE-DATA).
