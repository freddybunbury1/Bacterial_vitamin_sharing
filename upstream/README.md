# Optional upstream workflows

This directory documents computational steps performed before the public R
analysis begins:

- `genome_annotation/`: Thomas Janas's assembly-to-KO-matrix and GTDB-Tk
  workflow.
- `pyseer/`: Jeffrey Zhang's workflow for the Pyseer association summaries.

The manuscript R pipeline starts from the reviewed, checksummed outputs already
included under `data/raw/`. It does not automatically run these upstream
workflows. This separation lets readers reproduce the figures without rerunning
computationally intensive genome analyses, while still documenting how those
boundary inputs were created.

Each upstream workflow should:

1. document its input data and exact command;
2. record software, database, and environment versions;
3. avoid absolute paths tied to one computer;
4. accept a configurable output directory;
5. write regenerated files under its local `output/` directory by default; and
6. identify the corresponding fixed files under `data/raw/`.

Local `output/` contents are ignored by Git. A regenerated file should replace a
fixed manuscript input only after its identifiers, schema, values, and
provenance have been reviewed.
