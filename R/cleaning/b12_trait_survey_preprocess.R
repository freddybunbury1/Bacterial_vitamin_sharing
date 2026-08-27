section_name <- "b12_trait_survey"

strain_metadata_path <- "data/raw/metadata_strains/bacterial_strain_metadata.csv"

raw_dir <- file.path("data/raw", section_name)
intermediate_dir <- file.path("data/intermediate", section_name)
helper_path <- file.path("R/utils", "b12_trait_survey_helpers.R")

experiment_name_list <- c("B12_transfer2", "B12_transfer3")
sample_time_hour_list <- c(20, 44)
standard_concs_pm <- c(5000 * 2.5^(0:-10), 0)
od_path_length_cm <- 45 / 11.35 / 10

od_file_pattern <- "OD600_384well_B12_transfer(2|3)"
viability_file_pattern <- "nuc_green_dead_B12_transfer(2|3)"
b12_file_pattern <- "FB_Chlorophyll_384_B12_transfer(2|3)"

file_manifest_path <- file.path(intermediate_dir, "b12_trait_survey_input_file_manifest.csv")
preprocess_parameters_path <- file.path(intermediate_dir, "preprocess_parameters_b12_trait_survey.csv")
culture_layout_path <- file.path(intermediate_dir, "culture_layout_b12_trait_survey.csv")
culture_samples_path <- file.path(intermediate_dir, "culture_samples_b12_trait_survey.csv")
od_viability_map_path <- file.path(intermediate_dir, "od_viability_assay_map_b12_trait_survey.csv")
b12_storage_map_path <- file.path(intermediate_dir, "b12_storage_map_b12_trait_survey.csv")
b12_assay_metadata_path <- file.path(intermediate_dir, "b12_assay_metadata_b12_trait_survey.csv")
od_raw_path <- file.path(intermediate_dir, "od600_raw_long_b12_trait_survey.csv")
od_summary_path <- file.path(intermediate_dir, "od600_summary_b12_trait_survey.csv")
viability_raw_path <- file.path(intermediate_dir, "viability_raw_long_b12_trait_survey.csv")
viability_summary_path <- file.path(intermediate_dir, "viability_summary_b12_trait_survey.csv")
od_viability_summary_path <- file.path(intermediate_dir, "od600_viability_summary_b12_trait_survey.csv")
b12_raw_path <- file.path(intermediate_dir, "b12_assay_raw_long_b12_trait_survey.csv")
b12_assay_combined_path <- file.path(intermediate_dir, "b12_assay_combined_b12_trait_survey.csv")
pfister_strains_path <- file.path(intermediate_dir, "pfister_lab_strains_used_b12_trait_survey.csv")

source(helper_path)
check_required_packages(c("dplyr", "purrr", "stringr", "tidyr", "tibble"))

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

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

source_b12_files <- list.files(
  file.path(raw_dir, "b12_assay"),
  pattern = b12_file_pattern,
  full.names = TRUE
)

input_file_index <- dplyr::bind_rows(
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

strainbank <- utils::read.csv(
  strain_metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

culture_metadata <- build_culture_metadata(
  strainbank = strainbank,
  experiment_name_list = experiment_name_list,
  sample_time_hour_list = sample_time_hour_list
)

culture_layout <- culture_metadata$culture_layout
culture_samples <- culture_metadata$culture_samples
blank_ids <- culture_samples |>
  dplyr::filter(strain_id == "blank") |>
  dplyr::pull(sample_id)

oddead_map <- build_od_viability_metadata(culture_samples)
b12_storage_map <- build_b12_storage_metadata(culture_samples)
b12_assay_metadata_samples <- build_b12_assay_sample_metadata(b12_storage_map)
b12_standards <- build_b12_standard_metadata(
  experiment_name_list = experiment_name_list,
  standard_concs_pm = standard_concs_pm
)

b12_full_metadata <- dplyr::bind_rows(b12_assay_metadata_samples, b12_standards)

metadata_duplicates <- b12_full_metadata |>
  dplyr::mutate(assay_batch = dplyr::coalesce(assay_batch, experiment_name)) |>
  dplyr::count(assay_batch, assay_plate_number, assay_well) |>
  dplyr::filter(n > 1)

if (nrow(metadata_duplicates) > 0) {
  print(utils::head(metadata_duplicates, 25))
  stop("Duplicate assay positions in B12 assay metadata.")
}

od_files <- raw_file_manifest |>
  dplyr::filter(data_type == "od600") |>
  dplyr::pull(local_path)

viability_files <- raw_file_manifest |>
  dplyr::filter(data_type == "viability_stain") |>
  dplyr::pull(local_path)

b12_files <- raw_file_manifest |>
  dplyr::filter(data_type == "b12_assay") |>
  dplyr::pull(local_path)

od_raw <- purrr::map_dfr(od_files, process_od_grid) |>
  dplyr::mutate(
    od600_plate_read = od600,
    od_path_length_cm = od_path_length_cm,
    od600 = od600_plate_read / od_path_length_cm
  )

od_summary <- od_raw |>
  dplyr::inner_join(
    oddead_map |>
      dplyr::distinct(sample_id, experiment_name, sample_time_hour, assay_plate_number, assay_well),
    by = c("experiment_name", "sample_time_hour", "assay_plate_number", "assay_well")
  ) |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    od600_plate_read_mean = mean(od600_plate_read, na.rm = TRUE),
    od_path_length_cm = dplyr::first(od_path_length_cm),
    raw_od600 = mean(od600, na.rm = TRUE),
    od600_sd = stats::sd(od600, na.rm = TRUE),
    od600_n = sum(!is.na(od600)),
    .groups = "drop"
  )

viability_raw <- purrr::map_dfr(
  viability_files,
  process_vantastar_list,
  experiment_regex = "(?<=dead_).*(?=_plate)",
  experiment_prefix_regex = "^$"
)

viability_summary <- viability_raw |>
  dplyr::inner_join(
    oddead_map |>
      dplyr::distinct(
        sample_id,
        experiment_name,
        sample_time_hour,
        assay_plate_number,
        assay_well,
        viability_status
      ),
    by = c("experiment_name", "sample_time_hour", "assay_plate_number", "assay_well")
  ) |>
  dplyr::group_by(sample_id, viability_status) |>
  dplyr::summarise(
    stain_fluor = mean(fluorescence, na.rm = TRUE),
    stain_sd = stats::sd(fluorescence, na.rm = TRUE),
    stain_n = sum(!is.na(fluorescence)),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from = viability_status,
    values_from = c(stain_fluor, stain_sd, stain_n),
    names_glue = "{viability_status}_{.value}"
  ) |>
  dplyr::rename(
    live_stain = live_stain_fluor,
    dead_stain = dead_stain_fluor,
    live_stain_sd = live_stain_sd,
    dead_stain_sd = dead_stain_sd,
    live_stain_n = live_stain_n,
    dead_stain_n = dead_stain_n
  )

od_viab_summary <- viability_summary |>
  dplyr::full_join(od_summary, by = "sample_id") |>
  dplyr::mutate(
    od600 = raw_od600 - stats::median(raw_od600[sample_id %in% blank_ids], na.rm = TRUE)
  )

b12_data <- purrr::map_dfr(
  b12_files,
  process_vantastar_list,
  experiment_regex = "(?<=FB_Chlorophyll_384_).*(?=_plate)",
  experiment_prefix_regex = "^$"
) |>
  dplyr::select(
    raw_file,
    assay_batch,
    experiment_name,
    assay_plate_number,
    assay_well,
    measure_time_hours,
    fluorescence
  )

b12_assay_combined <- b12_data |>
  dplyr::right_join(
    b12_full_metadata |>
      dplyr::mutate(assay_batch = dplyr::coalesce(assay_batch, experiment_name)),
    by = c("assay_batch", "experiment_name", "assay_plate_number", "assay_well")
  ) |>
  dplyr::mutate(
    assay_dilution = (assay_volume_added + sample_volume_added) /
      sample_volume_added * storage_dilution_factor,
    adjusted_b12_conc_pm = sample_b12_conc_pm / assay_dilution + 0.1,
    log_b12_fm = log(adjusted_b12_conc_pm * 1000),
    log_fluorescence = log(fluorescence)
  )

pfister_strains <- culture_layout |>
  dplyr::filter(isolation_source == "marine") |>
  dplyr::filter(stringr::str_detect(strain_id, "^(TI|TN|DSM)")) |>
  dplyr::distinct(strain_id) |>
  dplyr::select(strain_id)

preprocess_parameters <- tibble::tibble(
  parameter = "od_path_length_cm",
  value = od_path_length_cm,
  definition = "Path length used to convert raw OD600 plate reads to 1 cm path length-adjusted OD600; adjusted OD600 = plate read / od_path_length_cm."
)

utils::write.csv(culture_layout, culture_layout_path, row.names = FALSE, na = "")
utils::write.csv(preprocess_parameters, preprocess_parameters_path, row.names = FALSE, na = "")
utils::write.csv(culture_samples, culture_samples_path, row.names = FALSE, na = "")
utils::write.csv(oddead_map, od_viability_map_path, row.names = FALSE, na = "")
utils::write.csv(b12_storage_map, b12_storage_map_path, row.names = FALSE, na = "")
utils::write.csv(b12_full_metadata, b12_assay_metadata_path, row.names = FALSE, na = "")
utils::write.csv(od_raw, od_raw_path, row.names = FALSE, na = "")
utils::write.csv(od_summary, od_summary_path, row.names = FALSE, na = "")
utils::write.csv(viability_raw, viability_raw_path, row.names = FALSE, na = "")
utils::write.csv(viability_summary, viability_summary_path, row.names = FALSE, na = "")
utils::write.csv(od_viab_summary, od_viability_summary_path, row.names = FALSE, na = "")
utils::write.csv(b12_data, b12_raw_path, row.names = FALSE, na = "")
utils::write.csv(b12_assay_combined, b12_assay_combined_path, row.names = FALSE, na = "")
utils::write.csv(pfister_strains, pfister_strains_path, row.names = FALSE, na = "")

message("Wrote raw file manifest to: ", file_manifest_path)
message("Read staged raw files from: ", raw_dir)
message("Wrote culture and assay metadata into: ", intermediate_dir)
message("OD files parsed: ", length(od_files), "; rows: ", nrow(od_raw))
message("Viability files parsed: ", length(viability_files), "; rows: ", nrow(viability_raw))
message("B12 assay files parsed: ", length(b12_files), "; rows: ", nrow(b12_data))
message("B12 assay combined rows: ", nrow(b12_assay_combined))
