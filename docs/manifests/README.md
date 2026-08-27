# Input transfer manifest

`file_manifest.csv` records the provenance and checksums of the primary inputs
copied while assembling the public repository, plus the KEGG retrieval script.
It is not intended to inventory every README or source file subsequently added
to the repository; Git provides that versioned file inventory.

- `repository_path`: path relative to the staged repository root.
- `working_source_path`: current path relative to the working project root.
- `data_class`: metadata, instrument, computational, external-reference,
  intermediate, processed, result, script, or documentation.
- `origin`: instrument/run/database/person that produced the file.
- `immutable`: `yes` for byte-preserved primary inputs.
- `upload_policy`: git, git-lfs, external, generated, or undecided.
- `size_bytes`: release-file size in bytes.
- `sha256`: checksum of the release file.
- `notes`: curation or provenance details.
