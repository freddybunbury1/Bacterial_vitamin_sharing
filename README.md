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
complete pipeline in dependency order and writes the figures and a
self-contained HTML notebook under `results/compiled_figures/`.

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
  "ape", "broom", "cowplot", "dplyr", "drc", "ggnewscale", "ggplot2",
  "ggrepel", "knitr", "lme4", "lubridate", "magick", "maps",
  "minpack.lm", "patchwork", "phytools", "purrr", "ranger", "rlang",
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
  Git-ignored outputs.

See [docs/PIPELINE.md](docs/PIPELINE.md) for dependencies and
[docs/REPRODUCIBILITY_AUDIT.md](docs/REPRODUCIBILITY_AUDIT.md) for validation.

## Licensing

Copyright (c) 2026 Freddy Bunbury and contributors.

Source code is available under the [MIT License](LICENSE). Data,
documentation, and figure assets are available under the
[Creative Commons Attribution 4.0 International License](LICENSE-DATA).
