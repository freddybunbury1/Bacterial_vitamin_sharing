# Assembly-to-KO and GTDB-Tk workflow

Contributor: Thomas Janas

## Purpose

Add the exact code and configuration used to process genome assemblies into the
KO count matrix and GTDB-Tk outputs consumed by the manuscript R pipeline. This
workflow is optional for reproducing the figures because the reviewed outputs
are already included under `data/raw/`.

## Files to add

- Pipeline entry point and any helper scripts.
- Configuration or parameter files used for the manuscript run.
- A reproducible software environment: for example `environment.yml`, a
  container definition, or locked workflow versions.
- A completed `assemblies_manifest.tsv` describing every assembly input.
- The exact command used for the manuscript run.

Please record the annotation software versions, GTDB-Tk version, GTDB database
release, KO annotation method/database version, computational parameters, and
any exclusions or manual interventions.

Do not include machine-specific absolute paths. Commands should accept a local
output directory, for example:

```bash
bash run_pipeline.sh --assemblies assemblies_manifest.tsv --outdir output
```

Adapt that example to the actual workflow language rather than adding a wrapper
that does not reflect how the analysis was run.

## Output contract

| Regenerated workflow output | Fixed manuscript input |
| --- | --- |
| `output/ko_count_matrix.tsv` | `data/raw/genome_pipeline_nf_core/ko_count_matrix.tsv` |
| `output/gtdbtk.bac120.summary.tsv` | `data/raw/gtdbtk/gtdbtk.bac120.summary.tsv` |
| `output/gtdbtk.bac120.study_genomes.tree` | `data/raw/gtdbtk/gtdbtk.bac120.study_genomes.tree` |

The workflow documentation should also explain how assembly identifiers are
connected to
`data/raw/metadata_strains/microtrait_strain_id_to_genome_filename_20260520.csv`.

The workflow must write under `output/` by default and must not overwrite the
fixed `data/raw/` files automatically. Differences between regenerated and
fixed outputs should be reviewed before any manuscript input is replaced.

## Manuscript run record

Complete this section when adding the workflow:

- Date run:
- Code version or commit:
- Exact command:
- Assembly source/accessions:
- Annotation pipeline and version:
- GTDB-Tk version:
- GTDB database release:
- KO annotation source/version:
- Parameters differing from software defaults:
- Excluded assemblies or manual interventions:
- Expected numerical or byte-level agreement with the fixed outputs:
