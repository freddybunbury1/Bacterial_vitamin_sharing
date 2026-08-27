section_name <- "b12_release_uptake_death"

source_strain_layout_path <- "data/raw/b12_release_uptake_death/metadata/strains_selected_b12_release_uptake_20260520.csv"
strain_metadata_path <- "data/raw/metadata_strains/bacterial_strain_metadata.csv"

raw_dir <- file.path("data/raw", section_name)
intermediate_dir <- file.path("data/intermediate", section_name)
helper_path <- file.path("R/utils", "b12_release_uptake_death_helpers.R")

b12_file_pattern <- "FB_Chlorophyll_384_b12_release_uptake_death(1|2)"
od_file_pattern <- "OD600_96well_b12_release_uptake_death(1|2)"
viability_file_pattern <- "nuc_green_dead_b12_release_uptake_death(1|2)"
standard_concs_pm <- c(5000 * 2.5^(0:-10), 0)

file_manifest_path <- file.path(intermediate_dir, "b12_release_uptake_death_input_file_manifest.csv")
strain_layout_path <- file.path(raw_dir, "metadata", basename(source_strain_layout_path))

ru_culture_layout_path <- file.path(intermediate_dir, "release_uptake_culture_layout_b12_release_uptake_death.csv")
ld_culture_layout_path <- file.path(intermediate_dir, "live_dead_mix_culture_layout_b12_release_uptake_death.csv")
live_dead_source_layout_path <- file.path(intermediate_dir, "live_dead_source_culture_layout_b12_release_uptake_death.csv")
culture_samples_path <- file.path(intermediate_dir, "culture_samples_b12_release_uptake_death.csv")
b12_storage_layout_path <- file.path(intermediate_dir, "b12_storage_layout_b12_release_uptake_death.csv")
b12_assay_metadata_path <- file.path(intermediate_dir, "b12_assay_metadata_b12_release_uptake_death.csv")
od_metadata_path <- file.path(intermediate_dir, "od600_metadata_b12_release_uptake_death.csv")
od_raw_path <- file.path(intermediate_dir, "od600_raw_long_b12_release_uptake_death.csv")
od_combined_path <- file.path(intermediate_dir, "od600_combined_b12_release_uptake_death.csv")
viability_metadata_path <- file.path(intermediate_dir, "viability_metadata_b12_release_uptake_death.csv")
viability_raw_path <- file.path(intermediate_dir, "viability_raw_long_b12_release_uptake_death.csv")
viability_combined_path <- file.path(intermediate_dir, "viability_combined_b12_release_uptake_death.csv")
viability_summary_path <- file.path(intermediate_dir, "viability_summary_b12_release_uptake_death.csv")
b12_raw_path <- file.path(intermediate_dir, "b12_assay_raw_long_b12_release_uptake_death.csv")
b12_assay_combined_path <- file.path(intermediate_dir, "b12_assay_combined_b12_release_uptake_death.csv")

source(helper_path)
check_required_packages(c("dplyr", "purrr", "stringr", "tidyr", "tibble"))

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

source_b12_files <- list.files(
  file.path(raw_dir, "b12_assay"),
  pattern = b12_file_pattern,
  full.names = TRUE
)

source_od_files <- list.files(
  file.path(raw_dir, "bacteria_od"),
  pattern = od_file_pattern,
  full.names = TRUE
)

source_viability_files <- list.files(
  file.path(raw_dir, "viability"),
  pattern = viability_file_pattern,
  full.names = TRUE
)

input_file_index <- dplyr::bind_rows(
  tibble::tibble(
    data_type = "strain_layout",
    source_path = source_strain_layout_path,
    local_path = strain_layout_path
  ),
  tibble::tibble(
    data_type = "od600",
    source_path = source_od_files,
    local_path = source_od_files
  ),
  tibble::tibble(
    data_type = "viability_stain",
    source_path = source_viability_files,
    local_path = source_viability_files
  ),
  tibble::tibble(
    data_type = "b12_assay",
    source_path = source_b12_files,
    local_path = source_b12_files
  )
)

required_input_paths <- c(strain_metadata_path, input_file_index$source_path)

if (any(!file.exists(required_input_paths))) {
  stop(
    "Missing source input file(s): ",
    paste(required_input_paths[!file.exists(required_input_paths)], collapse = ", ")
  )
}

raw_file_manifest <- input_file_index |>
  dplyr::mutate(
    size_bytes = file.info(.data$local_path)$size,
    modified_time = as.character(file.info(.data$local_path)$mtime)
  )

strain_metadata_info <- file.info(strain_metadata_path)
raw_file_manifest <- dplyr::bind_rows(
  tibble::tibble(
    data_type = "strain_metadata",
    source_path = strain_metadata_path,
    local_path = strain_metadata_path,
    size_bytes = strain_metadata_info$size,
    modified_time = as.character(strain_metadata_info$mtime)
  ),
  raw_file_manifest
)
utils::write.csv(raw_file_manifest, file_manifest_path, row.names = FALSE, na = "")

strain_metadata <- utils::read.csv(
  strain_metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

strain_layout <- utils::read.csv(
  strain_layout_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

culture_metadata <- build_release_uptake_death_culture_samples(
  strain_layout = strain_layout
)

ru_culture_layout <- culture_metadata$ru_culture_layout
ld_culture_layout <- culture_metadata$ld_culture_layout
live_dead_source_layout <- culture_metadata$live_dead_source_layout
culture_samples <- culture_metadata$culture_samples

b12_storage_layout <- build_b12_storage_layout(culture_samples)
b12_assay_layout <- build_b12_assay_layout(b12_storage_layout)
od_metadata <- build_release_uptake_death_od_metadata(
  ru_culture_layout = ru_culture_layout,
  live_dead_source_layout = live_dead_source_layout
)
viability_metadata <- build_release_uptake_death_viability_metadata(
  ru_culture_layout = ru_culture_layout,
  ld_culture_layout = ld_culture_layout,
  live_dead_source_layout = live_dead_source_layout
)

od_files <- raw_file_manifest |>
  dplyr::filter(data_type == "od600") |>
  dplyr::pull(local_path)

viability_files <- raw_file_manifest |>
  dplyr::filter(data_type == "viability_stain") |>
  dplyr::pull(local_path)

od_raw <- purrr::map_dfr(od_files, process_release_uptake_death_od_file)

od_combined <- od_raw |>
  dplyr::inner_join(
    od_metadata,
    by = c("experiment_name", "assay_plate_number", "assay_well")
  )

viability_raw <- purrr::map_dfr(
  viability_files,
  process_release_uptake_death_viability_file
)

viability_combined <- viability_raw |>
  dplyr::inner_join(
    viability_metadata,
    by = c("experiment_name", "assay_plate_number", "assay_well")
  )

viability_summary <- viability_combined |>
  dplyr::group_by(
    experiment_name,
    assay_plate_number,
    assay_well,
    viability_plate_role,
    source_layout_file,
    source_culture_plate_numbers,
    pre_b12_collapsed,
    technical_replicate,
    viability_status,
    experiment_type,
    culture_plate_number,
    culture_well,
    source_row,
    source_col,
    sample_time_minute,
    sample_time_hour,
    strain_id,
    biological_replicate,
    killed_proportion,
    added_b12_pM,
    added_glucose_mM_C
  ) |>
  dplyr::summarise(
    fluorescence_sd = stats::sd(fluorescence, na.rm = TRUE),
    fluorescence = mean(fluorescence, na.rm = TRUE),
    temporal_replicate_n = dplyr::n(),
    measurement_minutes = paste(sort(unique(measurement_minute)), collapse = ";"),
    raw_files = paste(sort(unique(raw_file)), collapse = ";"),
    .groups = "drop"
  )

b12_files <- raw_file_manifest |>
  dplyr::filter(data_type == "b12_assay") |>
  dplyr::pull(local_path)

b12_data <- purrr::map_dfr(b12_files, process_b12_release_uptake_file)
assay_batches <- sort(unique(b12_data$assay_batch))

b12_standards <- build_b12_release_uptake_standards(
  assay_batches = assay_batches,
  standard_concs_pm = standard_concs_pm
)

if (length(assay_batches) != 1) {
  stop(
    "Expected one B12 assay batch for this section, but found: ",
    paste(assay_batches, collapse = ", "),
    ". Add explicit sample-to-assay-batch mapping before proceeding."
  )
}

b12_full_metadata <- dplyr::bind_rows(
  b12_assay_layout |>
    dplyr::mutate(
      assay_batch = assay_batches[1],
      sample_id = as.integer(sample_id)
    ),
  b12_standards |>
    dplyr::mutate(sample_id = as.integer(sample_id))
)

metadata_duplicates <- b12_full_metadata |>
  dplyr::count(assay_batch, assay_plate_number, assay_well) |>
  dplyr::filter(n > 1)

if (nrow(metadata_duplicates) > 0) {
  print(utils::head(metadata_duplicates, 25))
  stop("Duplicate assay positions in B12 assay metadata.")
}

b12_assay_combined <- b12_data |>
  dplyr::inner_join(
    b12_full_metadata,
    by = c("assay_batch", "assay_plate_number", "assay_well")
  ) |>
  dplyr::mutate(
    assay_dilution = (assay_volume_added + sample_volume_added) /
      sample_volume_added * storage_dilution_factor,
    adjusted_b12_conc_pm = sample_b12_conc_pm / assay_dilution + 0.1,
    log_b12_fm = log(adjusted_b12_conc_pm * 1000),
    log_fluorescence = log(fluorescence)
  )

utils::write.csv(ru_culture_layout, ru_culture_layout_path, row.names = FALSE, na = "")
utils::write.csv(ld_culture_layout, ld_culture_layout_path, row.names = FALSE, na = "")
utils::write.csv(live_dead_source_layout, live_dead_source_layout_path, row.names = FALSE, na = "")
utils::write.csv(culture_samples, culture_samples_path, row.names = FALSE, na = "")
utils::write.csv(b12_storage_layout, b12_storage_layout_path, row.names = FALSE, na = "")
utils::write.csv(b12_full_metadata, b12_assay_metadata_path, row.names = FALSE, na = "")
utils::write.csv(od_metadata, od_metadata_path, row.names = FALSE, na = "")
utils::write.csv(od_raw, od_raw_path, row.names = FALSE, na = "")
utils::write.csv(od_combined, od_combined_path, row.names = FALSE, na = "")
utils::write.csv(viability_metadata, viability_metadata_path, row.names = FALSE, na = "")
utils::write.csv(viability_raw, viability_raw_path, row.names = FALSE, na = "")
utils::write.csv(viability_combined, viability_combined_path, row.names = FALSE, na = "")
utils::write.csv(viability_summary, viability_summary_path, row.names = FALSE, na = "")
utils::write.csv(b12_data, b12_raw_path, row.names = FALSE, na = "")
utils::write.csv(b12_assay_combined, b12_assay_combined_path, row.names = FALSE, na = "")

message("Wrote raw file manifest to: ", file_manifest_path)
message("Read staged raw OD files from: ", file.path(raw_dir, "bacteria_od"))
message("Read staged raw viability files from: ", file.path(raw_dir, "viability"))
message("Read staged raw B12 assay files from: ", file.path(raw_dir, "b12_assay"))
message("Wrote culture and assay metadata into: ", intermediate_dir)
message("Assay batches parsed: ", paste(assay_batches, collapse = ", "))
message("OD files parsed: ", length(od_files), "; raw rows: ", nrow(od_raw))
message("OD combined rows: ", nrow(od_combined))
message("Viability files parsed: ", length(viability_files), "; raw rows: ", nrow(viability_raw))
message("Viability combined rows with current metadata: ", nrow(viability_combined))
message("Viability temporal-summary rows with current metadata: ", nrow(viability_summary))
message("B12 assay files parsed: ", length(b12_files), "; raw rows: ", nrow(b12_data))
message("B12 assay combined rows: ", nrow(b12_assay_combined))
