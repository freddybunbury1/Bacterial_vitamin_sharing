# External annotation references

Run `Rscript --vanilla R/cleaning/download_kegg_ko_mapping.R` from the
repository root before compiling figures. The script downloads the KO mapping
and release information directly from KEGG, validates KO identifiers, and
records retrieval provenance. Its three outputs are ignored by Git and are not
redistributed in the repository.
