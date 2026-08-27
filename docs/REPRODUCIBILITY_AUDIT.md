# Reproducibility audit

Audit date: 2026-08-26.

## Inputs and asset

- 1,081 staged files were present under `data/raw/` before the local KEGG KO
  acquisition output was added; no files were zero bytes.
- All 1,065 entries recorded in `docs/manifests/file_manifest.csv` retained
  their recorded SHA-256 checksum after the public retrieval-script record was
  updated to its final stable-filename implementation.
- `assets/figure_source/Figure1D_cartoon.svg` is byte-identical to the approved
  manuscript artwork: SHA-256
  `1633a080a2674ebd5ae76acd23d7f31ef5ded80b0b2766adc3145c906bfec77a`.
- The locally acquired validation copy of `kegg_ko_list.tsv` had SHA-256
  `c471b70cfc78e6b7ed0a086e1a87be883607c6c85334022d5df89b49452be5ed`.
  This file is not distributed and may change when KEGG updates the mapping.

## Clean-room validation

Two independent temporary copies were created with zero files under
`data/intermediate/`, `data/processed/`, and `results/`. Each was rendered from
a working directory outside the repository. Both completed the raw-data-to-
figures pipeline and produced identical 57-file result inventories:

- 25 PNG figures;
- 25 PDF figures;
- 6 CSV source/statistics tables;
- 1 self-contained HTML notebook.

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
fold design. Figure 3 and Figure S2 had only 274 and 38 changed pixels,
respectively, associated with the curated public metadata/strain identifiers.
All other compiled PNGs had zero pixel differences.
