# Pyseer association workflow

Contributor: Jeffrey Zhang

## Purpose

Add the exact code and configuration used to produce the Pyseer association
summaries consumed by `R/analysis/b12_trait_classification_analysis.R`. This
workflow is optional for reproducing the figures because the reviewed summaries
are already included under `data/raw/gwa/`.

## Files to add

- Pipeline entry point and any helper scripts.
- Configuration or parameter files used for the manuscript run.
- A reproducible software environment: for example `environment.yml`, a
  container definition, or pinned Python/package versions.
- A completed `input_manifest.tsv` identifying phenotype, genomic feature,
  population-structure, tree, and sample-mapping inputs.
- The exact commands used for every phenotype/model combination.

Please record the Pyseer version, Python/environment versions, phenotype
construction, sample inclusion/exclusion rules, MAF and correlation settings,
population-structure correction, tree input, and any post-processing used to
create the four summary CSVs.

Do not include machine-specific absolute paths. Commands should accept a local
output directory, for example:

```bash
bash run_pyseer.sh --input-manifest input_manifest.tsv --outdir output
```

Adapt that example to the actual workflow rather than adding a wrapper that
does not reflect how the analysis was run.

## Output contract

Each regenerated file under `output/` corresponds to the identically named
fixed manuscript input under `data/raw/gwa/`:

- `log10_total_b12_gm_MAF_0.01_COR_1_pyseer_mixed_summary.csv`
- `log10_total_b12_gm_MAF_0.01_COR_1_pyseer_mixed_tree_summary.csv`
- `uptake_b12_mean_MAF_0.01_COR_1_pyseer_mixed_summary.csv`
- `uptake_b12_mean_MAF_0.01_COR_1_pyseer_mixed_tree_summary.csv`

The workflow must write under `output/` by default and must not overwrite the
fixed `data/raw/gwa/` files automatically. Differences between regenerated and
fixed outputs should be reviewed before any manuscript input is replaced.

## Manuscript run record

Complete this section when adding the workflow:

- Date run:
- Code version or commit:
- Exact commands:
- Pyseer version:
- Python/environment versions:
- Phenotype inputs and transformations:
- Genomic feature input:
- Population-structure and tree inputs:
- Sample inclusion/exclusion rules:
- MAF, correlation, and other parameters:
- Post-processing steps:
- Expected numerical or byte-level agreement with the fixed outputs:
