# Strain metadata

Canonical curated strain metadata and strain-to-genome mappings belong here.

`bacterial_strain_metadata.csv` is the canonical 15-column public derivative of
the dated working metadata table. It contains strain identity, physical inventory,
experimental-array layout, and isolation provenance. The following 16S-derived
taxonomy columns were intentionally excluded from the public version:
`species_name`, `species_name_shortened`, `chosen_phylum`, `chosen_class`,
`chosen_order`, `chosen_family`, and `chosen_genus`.

Each row represents one physical strain-bank or experimental-array position,
not necessarily one unique isolate. Repeated `strain_id` values are therefore
permitted.

`stored_adjusted_od600` records the adjusted OD600 associated with the stored
strain entry. The more specific name distinguishes it from OD600 measurements
made during downstream experiments.

`growth_medium` describes the medium used to generate the stored stocks. It is
defined from the stock array: every position in `marine_strainbank_array2` is
assigned `M13`, and every other position is assigned
`TSB_variant_b_0.1x`.

Storage and experiment-layout columns appear first, followed by isolate and
environmental-provenance columns. Rows whose `strain_id` is `blank` or
`unknown` retain their array-derived `growth_medium`, but all isolate and
environmental-provenance fields are empty.
