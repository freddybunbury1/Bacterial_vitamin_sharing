# Data layout

## Curated metadata

`data/raw/metadata_strains/bacterial_strain_metadata.csv` is the canonical
15-column strain and experiment-layout table. It contains the array placement,
stock growth medium, isolate provenance, and collection metadata used by the
public analyses. GTDB-Tk is the public taxonomy source; legacy internal 16S
`chosen_*` fields are intentionally absent.

Additional mappings and layouts are stored with the analyses that consume them:

- `data/raw/metadata_strains/microtrait_strain_id_to_genome_filename_20260520.csv`
- `data/raw/b12_release_uptake_death/metadata/`

## Instrument and assay inputs

- `data/raw/b12_trait_survey/{bacteria_od,viability,b12_assay}/`
- `data/raw/b12_release_uptake_death/{bacteria_od,viability,b12_assay}/`
- `data/raw/growth_dynamics_b12_methionine/`
- `data/raw/b12_assay_comparison/Ecoli_chlamy_b12_assay_dataset.csv`

Instrument-export filenames and bytes are preserved.

## Genome-derived inputs

- `data/raw/genome_pipeline_nf_core/ko_count_matrix.tsv`
- `data/raw/gwa/`
- `data/raw/gtdbtk/gtdbtk.bac120.summary.tsv`
- `data/raw/gtdbtk/gtdbtk.bac120.study_genomes.tree`
- `data/raw/metadata_strains/microtrait_strain_id_to_genome_filename_20260520.csv`

These are the current study-genome products. The older 193-genome GTDB-Tk
summary and the former `data/pre-raw/` tree are not inputs and are not included.

## KEGG acquisition exception

`R/cleaning/download_kegg_ko_mapping.R` retrieves:

- `data/raw/annotations/kegg_ko_list.tsv`
- `data/raw/annotations/kegg_ko_info.txt`
- `data/raw/annotations/kegg_ko_retrieval_metadata.csv`

These generated acquisition files are ignored by Git. The repository retains
`kegg_ko_info_2026-08-06.txt` as provenance for the KEGG version used during
manuscript development. Academic users should run the acquisition script under
KEGG's stated conditions.

## Generated outputs

`data/intermediate/`, `data/processed/`, and `results/` are generated from the
inputs above and ignored by Git. No production script reads from the legacy
working repository or from `results/`.
