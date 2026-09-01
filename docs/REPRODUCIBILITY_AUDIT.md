# Reproducibility audit

Audit date: 2026-08-27.

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

### Current working-tree validation

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
