# Reproducibility audit

Latest audit date: 2026-10-13. Earlier validation records are retained below.

## Unified-model sister copy, 2026-10-13

This working copy was made in full from `Bacterial_vitamin_sharing_delta_zero`;
all 2,733 copied files were initially byte-identical. The following model
changes are confined to this sister copy; the previous release is the baseline.

- Canonical first-order uptake and release fits are performed once. S12–S14
  reuse their parameters and replicate predictions, and use the same empirical
  and predicted summary values as main panels 4A, 4C and 5A. The shared-model
  test requires exact equality at all 78 main-panel coordinates and checks
  full first-order time-course summaries and replicate prediction provenance.
- Both structures use the same exported mean initial/final OD reference for
  kinetic trajectories and the live fraction in mixture predictions. S13 uses
  instantaneous killed-cell release followed by uptake without ongoing release;
  S14 includes ongoing release with the release-specific initial condition.
- Existing lower bounds are retained. Scalar fitting screens candidate minima
  before local refinement; all 368 non-blank canonical scalar fits are checked
  against independently evaluated screened alternatives. Seven uptake fits and
  two release fits differ by more than 0.1% from the baseline. The existing
  representative-selection rule now selects fb9 as the reclaimer example.
- Figure 4 agreement is CCC 0.775941177848553, MAE 0.0705287794229651 and
  permutation p 1/1001. Figure 5D slope is -0.284587827584429, R-squared
  0.496798574742114 and p 0.000173124884742930. Both models have 184 successful
  non-blank culture/glucose fits; all quota numerical checks pass.
- `tests/test_shared_kinetic_models.R` checks shared outputs, coefficient reuse,
  quota OD references, lower bounds and scalar optimizer objectives. Existing
  cohort, OD-specific-rate and quota-solver regression tests also pass.

### Subsequent housekeeping validation

Removed the unused release-inclusive live–dead prediction function and its
exports, the duplicate live–dead normalization export, four unused plotting
helpers and their styling definitions, and the unused model-ID helper.
Updated README and session wording to describe the shared current workflow;
historical validation entries remain explicitly labelled.

A full notebook rebuild succeeded after cleanup. Exactly three obsolete
processed CSVs were removed; all 1,237 retained raw/intermediate/processed
files are byte-identical to the pre-cleanup baseline. All 26 PNG exports are
pixel-identical, and all eight compiled CSVs are byte-identical. Shared-model,
OD-specific-rate and quota-solver checks pass, including retained bounds and
all 368 scalar optimizer checks. Release integrity passes 1,066 input checks,
the 64-file compiled inventory and notebook consistency; `git diff --check`
passes. No fit, prediction, reported statistic or figure changed in this pass.

### Upload-readiness validation

The upload candidate contains 1,130 tracked/non-ignored source and input files
(10,641,184 bytes; largest file approximately 3.5 MB). Ten tracked files are
modified and `tests/test_shared_kinetic_models.R` is the one new file to include
in the commit. Licenses and the required figure asset are tracked; generated
outputs, KEGG acquisition files and Finder caches are ignored by Git.
Corrected stale 200-strain counts in the current pipeline guide to 202 and
added the shared-model check to the root README.

A source copy containing only the upload-eligible files rebuilt all figures
and tables successfully, starting with no intermediate or processed analysis
outputs. This used the existing installed R packages and the previously
acquired KEGG annotation cache; it was not a new package-installation
environment or a clone of a newly pushed commit. All 106 processed CSVs are
byte-identical to the working copy, and all 26 PNGs are pixel-identical.
All four R validation checks and release integrity pass in that clean source
copy, including 1,066 input checks, the 64-file compiled inventory, exact
main/supplementary agreement and all 368 scalar optimizer checks. No upload
blockers were found. No files were staged, committed or pushed during the audit.

## Previous release validation on 2026-10-13

The complete production notebook was rendered after the quota solver update,
including all eight production scripts and all 42 executable notebook chunks.
Figures 1–5 and S1–S19, supplementary Tables S1–S3, figure-source tables and
the self-contained notebook were regenerated. Figure 5D conditionally displays
`p<0.001` for calculated p-values below that threshold.

- The shared complete-trait cohort contains 202 strains, including 79
  synthesizers. Figure S2 contains the full set of 231 QC genomes.
- Manual DSM24912 H6 and KC611216 H12 exclusions were removed. DSM24912 H6
  fails automatic OD QC; KC611216 pools eight observations across its two
  positions. Figure S1 deduplicates the 277 unique collection isolates before
  sampling-location clustering and taxonomy aggregation.
- Quota uptake, including the uptake-only S13 live–dead predictions, uses an
  analytical solution; S14 release trajectories use adaptive LSODA with relative/absolute tolerances of 1e-8.
  All 184 quota fits converge and all 1,056 fitted-trajectory checks pass
  against tighter integration and physical/conservation constraints.
- The cohort, OD-specific-rate and quota-solver regression tests pass. The
  quota regression test covers 144 parameter combinations plus limiting cases.
- Both reviewed curated-input manifest entries were refreshed. The two curated
  CSVs were normalized to LF line endings with exact preservation of their
  parsed records; their manifest checksums and sizes were refreshed again.
  All 1,066 input-manifest checksum and size checks pass.
- The compiled inventory contains exactly 64 production files. Notebook source
  and cached output sizes agree with the current source.
  Finder files, the default R plot and the obsolete S15 rate-source table were
  removed. Release integrity and `git diff --check` pass.

Manuscript Methods and captions are maintained in Overleaf. The repository
Methods/legend fragments and empty manual-QC flags file were removed, along
with their production dependencies and generated legend copy. Earlier records
below describe historical files and inventories.

This validates the current working tree and its generated outputs, not a newly
cloned committed release. Required new utilities, documentation and tests must
be included in the release commit. Fresh-clone/environment validation follows
publication of that exact commit. The upstream genome-annotation and Pyseer
folders remain provenance templates; the downstream pipeline uses the fixed
reviewed outputs and does not rerun those upstream workflows.

## Full working-tree rebuild (2026-10-01)

The complete public R Markdown entry point was rerun after the folder audit,
using R 4.5.1, rmarkdown 2.30, Pandoc 3.2 and the existing locally acquired
KEGG mapping. The run completed all 42 executable chunks, including all eight
production scripts, in 173 seconds. Logs, temporary rendering files and the
pre-run snapshot were kept outside the release folder.

- Regenerated all 52 PNG/PDF exports (Figures 1–5, S1–S19 and the two additional
  Figure 5 uptake panels), Tables S1–S3, figure-source/statistics tables and the
  self-contained notebook. The compiled inventory remains 61 files.
- All 1,085 raw-directory files are byte-identical to their pre-run versions;
  all 1,066 manifest entries still pass checksum and size verification.
- Intermediate and processed numerical outputs are byte-identical except for
  floating-point rounding in the 12D and 15C survey source tables. Those tables
  retain the same 79 strain IDs and columns, with all numerical differences
  within relative/absolute tolerance `1e-12`.
- The cohort, OD-specific-rate and release-integrity checks all pass. Independent
  S9 baseline and OD matching again confirms 1,152 panel-A measurements and 642
  panel-B measurements across 92 groups and 23 strains. The full cohort remains
  200 strains, including KC611216 once, with 79 synthesizers.
- Visual inspection of regenerated S1 and S9 confirms the requested legend
  clearance and panel resizing. The notebook contains the current source and
  matching output sizes; the exported legends match their source. No unexpected
  or temporary files remain in the release, and `git diff --check` passes.

This was a complete rebuild in the release working tree using its existing
raw inputs, not a separate clean-room copy.

## Release-folder consistency check (2026-10-01)

- Archived two historical phenotype-filter comparison directories outside the
  release (19 files), along with a stray `Rplots.pdf` and 11 Finder metadata
  files. No primary inputs, numerical analysis tables or current figure exports
  were changed by this cleanup. The two additional Figure 5 uptake subpanels
  remain intentional production outputs.
- Verified all 1,066 input-manifest SHA-256 checksums and sizes. The compiled
  directory contains exactly 61 generated files: 25 PNGs, 25 PDFs, six CSVs,
  three XLSX workbooks, one HTML notebook and one LaTeX legend file.
- The notebook's displayed source blocks and downloadable Rmd agree with the
  current production Rmd. Refreshed 12 cached file-size entries left by the
  cohort update. The exported legend file is identical to its source.
- Rechecked the two-panel S9 transformation independently: A contains 1,152
  measurements including 12 blanks; B contains 642 positive pairs across 92
  strain/experiment/biological-replicate groups and 23 strains. Matched
  0%-killed fluorescence subtraction and killed fraction × adjusted initial
  OD agree with the source records; two nonpositive fluorescence increments
  are excluded from B.
- The compiled-figure README records the current bold panel letters and layout
  settings for Figures 3–5, S1 and S9. The pipeline documentation now explains
  S9's two data sources, matching, baseline subtraction and exclusions. The
  direct `gtable` dependency is declared in the package list and preflight.
- Both existing R checks pass: `tests/test_figure_cohort.R` and
  `tests/test_od_specific_rates.R`. All 17 R scripts and 42 Rmd chunks parse.
  `tests/check_release_integrity.py` adds a read-only check of input checksums,
  expected output inventory, temporary/exploratory artifacts, notebook source
  and cached output sizes, and the exported legend copy. The root README
  documents all three commands. `git diff --check` passes.

This was an audit of the existing release and targeted checks, not a new
end-to-end or clean-room render. All 50 PNG/PDF exports, all numerical CSV/XLSX
outputs, and all embedded notebook plot images were preserved during this audit.

## Figure layout updates (2026-10-01)

Panel letters are bold throughout. Figure 3's two data panels are 5% wider;
Figures 4 and 5 transfer 4% of the specified left-panel width allocation to
their right panels. Figure 4B's y-axis title is 4% of the plot height lower.
Figure S1A's size legend is below the map, and S1B's plotting height is reduced
by 30% with its width unchanged. Figure S9B's plotting width is reduced by 8%.
The panel reductions were measured on the rendered plot viewports. Visual
inspection confirmed clear legends and panel labels. Export dimensions and
text sizes were retained; affected PNGs, PDFs and notebook previews were
updated without rerunning unrelated analyses.

## Shared phenotype cohort update (2026-09-30)

- The trait panels now use the QC-passing tree membership and finite total B12
  and uptake estimates used in Figure 1A. All-strain panels contain 200 strains;
  panels restricted to total B12 > 50 pM contain 79. Neither the three-sample
  threshold nor missing uncertainty excludes an otherwise eligible point.
- The classifier normalizes the legacy KC61126 identifier to KC611216 during
  joining. All four classifiers and their null comparisons were rerun. Full-set
  tasks use seven matched folds; synthesizer tasks use six. Figure 2, S6 and S7
  were regenerated from these outputs. Supplied pyseer summaries are unchanged.
- Figures 1B--C, 3, S4--S5 and S10 were regenerated with statistics calculated
  from the updated data. Figure 3A's lower x-axis limit was extended to 0.001
  to include sia0102 in both the points and fitted line. Figure wording and
  annotation formatting are otherwise retained.
- Figure 1A and 1D assets were reused. Only Figures 1, 2, 3, S4, S5, S6, S7
  and S10 PNG/PDF exports were replaced; all other figure exports are unchanged.
  The notebook's affected images, displayed changed code and downloadable Rmd
  were refreshed without rerunning unrelated analyses or figures.
- `Rscript --vanilla tests/test_figure_cohort.R` verifies the 200-strain cohort,
  its 79-strain synthesizer subset, the corrected KC join, retention of valid
  single-observation estimates, source-table membership and fold coverage.

## OD-specific rate update

A complete working-tree render rebuilt the raw-data-to-figures pipeline and
self-contained notebook after normalizing the main rate summaries to OD.

- All 184 non-blank culture/glucose fits (23 strains, four replicates, two
  glucose conditions) have valid matched initial and final OD measurements.
  The initial untreated-source OD is matched by strain, experiment and
  biological replicate. Final OD averages the matched 0- and 500-pM wells
  within glucose condition, and both time points receive equal weight.
- Fits retain their original culture-level coefficients and SEs unchanged
  (numerical tolerance `1e-12`), with explicit OD-specific columns added.
  Normalization precedes geometric averaging in Figures 4B, 5B–C and S16.
  Their common-strain estimates agree, and Figure 5D's k/r ratio is unchanged.
- Figure 5B's exact paired signed-rank p-values were recalculated from the
  normalized strain summaries: Retainers 0.15625, Providers 0.0078125,
  Reclaimers 0.546875. Axis limits and annotations accommodate all points.
- The three updated figures and two dependent additional uptake panels were
  regenerated with OD-specific units. Visual inspection confirmed readable
  labels and no clipping. All 20 other PNG figures are pixel-identical to
  their preceding versions, including Figures S12–S14. Their fitting workflow
  remains unchanged and is distinguished from the main normalization in the
  Methods.
- All processed tables other than the rate exports and normalization metadata
  are byte-identical to their preceding versions. The fitted concentration
  curves, mixture predictions and Figure 5D statistics are unchanged.
- `METHODS_kinetic_models.tex` uses B throughout, includes the blank directly
  in the equations, removes package/function names and historical initial-
  offset commentary, and documents the mean initial/final OD normalization.
  The S16 legend and pipeline READMEs agree. The exported legend file is
  identical to its source.
- `Rscript --vanilla tests/test_od_specific_rates.R` passes independent OD
  matching, coefficient/SE rescaling, missing/duplicate/nonpositive-input
  guards, blank handling and cross-figure agreement checks. R scripts and
  extracted notebook code parse, the Methods fragment parses as LaTeX in
  Pandoc, and `git diff --check` passes.

## Earlier figure and Methods consistency check (before OD normalization)

A complete `rmarkdown::render` rebuilt the raw-data-to-figures pipeline and
the self-contained HTML notebook using Pandoc 3.2 bundled with RStudio.
The current figure set is Figures 1–5 and S1–S19.

- Figures S12–S14 contain exactly two model structures: first-order in red
  (`#D73027`) and quota in blue (`#2166AC`), with line alpha 0.8. Their original
  axis labels are retained, and the continuous Figure S13 y scale ends at 1.5.
- All three comparison figures contain 24 facets, including the blank. Every
  strain's trait group and ordinal agrees with Table S3 and the 23-strain
  phylogeny in Figure S16.
- The retained models' numerical source tables are unchanged from the preceding
  comparison outputs (after restricting those outputs to the two retained
  structures; numerical comparison tolerance `1e-12`).
- The main uptake/release fit exports are numerically unchanged. Pixel
  comparisons confirmed that all five main figures and all four renumbered
  supplementary figures retain their preceding graphical content.
- Figure S16 is the selected-strain rate phylogeny; S17 contains candidate-gene
  associations; S18 compares bioassays; S19 shows survey assay fluorescence.
  The corresponding source-table names and notebook sections agree.
- Source code, documentation, regenerated CSVs, filenames and the rebuilt
  HTML were scanned for obsolete model names, parameter names, and figure
  references. Superseded generated outputs were removed from the repository.
- `METHODS_kinetic_models.tex` removes the shared initial offset from uptake
  and main live–dead predictions while retaining the measured release-only
  initial condition. It expresses rates in OD-specific notation and explicitly
  distinguishes the main fitted culture-level products from the supplementary
  OD-specific fits. It documents the quota equations and their supplementary-
  only use, and the release term in the supplementary live–dead comparison.
- `FIGURE_LEGENDS.tex` supplies the revised S12–S14 and S16 captions. The
  notebook displays these captions and exports an identical copy alongside
  the figures. All modified R code and extracted notebook code parse, the
  Methods fragment parses as LaTeX in Pandoc, and `git diff --check` passes.

At that stage, the recursive inventory under `results/compiled_figures/` contained
61 generated files: 25 PNGs, 25 PDFs, 6 CSVs, 3 XLSX workbooks, 1 HTML notebook,
and 1 LaTeX legend file. The PNG/PDF counts include the two additional Figure 5
subpanels. This was a full working-tree rebuild, not a new clean-room run.

## Inputs and asset

- 1,082 staged files were present under `data/raw/` before the local KEGG KO
  acquisition output was added; no files were zero bytes.
- All 1,066 entries recorded in `docs/manifests/file_manifest.csv` retained
  their recorded SHA-256 checksum after the public retrieval-script record was
  updated to its final stable-filename implementation.
- `assets/figure_source/Figure1D_cartoon.svg` is byte-identical to the approved
  manuscript artwork: SHA-256
  `1633a080a2674ebd5ae76acd23d7f31ef5ded80b0b2766adc3145c906bfec77a`.
- The locally acquired validation copy of `kegg_ko_list.tsv` had SHA-256
  `c471b70cfc78e6b7ed0a086e1a87be883607c6c85334022d5df89b49452be5ed`.
  This file is not distributed and may change when KEGG updates the mapping.

## Clean-room validation

Two independent temporary copies were initially created with zero files under
`data/intermediate/`, `data/processed/`, and `results/`. Each was rendered from
a working directory outside the repository. Both completed the raw-data-to-
figures pipeline and produced identical inventories before supplementary-table
export was added.

A subsequent complete render on 2026-08-27 validated the count-neutral genome
mapping filename and the Table S1–S2 exports. It produced a 59-file result
inventory:

- 25 PNG figures;
- 25 PDF figures;
- 6 CSV source/statistics tables;
- 2 XLSX supplementary tables;
- 1 self-contained HTML notebook.

`TableS1.xlsx` was verified against the canonical metadata as 288 rows by 15
columns and includes Notes documenting its 277 unique isolates. `TableS2.xlsx`
was verified as 210 rows by 2 non-empty recipe columns and includes a provenance
and presentation-processing Notes worksheet.

### Earlier working-tree validation

A complete render on 2026-09-01 validated the addition of GTDB-Tk taxonomy to
Table S1, the new Table S3 export, and promotion of the manuscript Figure 5 to
the sole canonical Figure 5 output. It produced a 58-file result inventory:

- 24 PNG figures;
- 24 PDF figures;
- 6 CSV source/statistics tables;
- 3 XLSX supplementary tables;
- 1 self-contained HTML notebook.

`TableS1.xlsx` was verified as 288 rows by 20 columns. Its five GTDB-Tk
taxonomy columns populate 233 physical-position rows representing 229 unique
isolates; unmatched rows remain blank. `TableS2.xlsx` remained 210 rows by 2
non-empty recipe columns. `TableS3.xlsx` was verified as 24 rows by 10 columns:
23 bacterial strains plus the `blank` control. All non-control rows have
complete GTDB-Tk phylum-to-genus taxonomy, trait group, ordinal, and panel
label. The blank row retains culture well B12 and has no taxonomy, trait,
ordinal, or panel-label values.

This was a full render of the release working tree, not a new independent
clean-room validation.

### Revised kinetic-model validation

The revised release working tree was validated on 2026-09-19 after separating
the uptake and release-plus-uptake initial-condition treatments. The primary
release/uptake/death analysis completed successfully, all modified R files
parsed, and a `knitr::knit` execution completed all 78 notebook chunks. A
self-contained `rmarkdown::render` was not repeated in this command-line
environment because Pandoc was unavailable; the repository README documents
how to use the Pandoc copy bundled with RStudio.

The checks confirmed that uptake-only fit exports contain no release initial
condition, release fits use `release_initial_b12_pM`, Figure 5C and Figure 5D
use only zero-added-glucose parameter estimates, and the focused comparison
contained additional model-structure comparisons. This earlier validation
was superseded by the current first-order/quota comparison described below.

An additional complete notebook execution on 2026-09-21 validated the
branch-length tree rendering, revised GTDB-Tk class palette, and the selected-
strain rate phylogeny (now Figure S16).
The rate phylogeny contains exactly 23 non-blank selected strains, all present in the
GTDB tree, with glucose-zero first-order release and uptake rates summarized by
geometric mean across replicate-specific fits.

At that earlier stage, the compiled-output inventory contained 61 files:

- 26 PNG figures;
- 26 PDF figures;
- 5 CSV source/statistics tables;
- 3 XLSX supplementary tables;
- 1 self-contained HTML notebook.

The HTML notebook in that inventory is the previously rendered artifact; the
78-chunk validation above generated Markdown in a temporary directory and did
not overwrite it.

All intermediate and processed files were byte-identical between the two runs.
All PNG figures were byte-identical except the ImageMagick-composed Figure 1
container; a pixel comparison of the two Figure 1 PNGs reported zero differing
pixels. PDF and notebook containers may embed run-specific metadata.

A third pristine render passed after removal of the unused kinetic benchmark
branch. The subsequent graphics-device cleanup was verified by rerunning the
classifier: its data products remained byte-identical and it left no
`Rplots.pdf` side effect.

The random-fold assignment audit confirmed one row per trait/strain, exhaustive
fold IDs, and no overlapping test assignments. Fold counts were 7 for the two
full-tree classifiers and 6 for the two total-B12-filtered classifiers.

## Comparison with the working manuscript output

Pixel comparisons against the working manuscript output were identical for 18
of 23 compiled PNGs. Expected differences occurred in Figure 2, Figure S7, and
Figure S8 because the public classifier uses the corrected order-matched random
fold design. Figure 3 and Figure S3 had only 274 and 38 changed pixels,
respectively, associated with the curated public metadata/strain identifiers.
All other compiled PNGs had zero pixel differences.
