# R pipeline

The public pipeline has one fixed production path:

1. the three experiment-specific scripts in `cleaning/` parse raw instrument
   files;
2. the five scripts in `analysis/` create canonical results and compact
   figure-source tables;
3. `compiled_figures/compiled_figures.Rmd` runs those stages in order and
   assembles every manuscript figure;
4. `utils/` contains shared assay and kinetic-model functions.

There is no production `R/figures/` layer. Figure-source exports are owned by
the analysis that computes them, and final assembly is owned by the R Markdown
notebook.
