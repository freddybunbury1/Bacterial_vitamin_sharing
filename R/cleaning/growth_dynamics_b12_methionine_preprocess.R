layout_metadata_path <- "data/raw/metadata_strains/bacterial_strain_metadata.csv"
plate_reader_dir <- "data/raw/growth_dynamics_b12_methionine/20251208_bacteria_b12_methionine_2"
output_dir <- "data/intermediate/growth_dynamics_b12_methionine"

metadata_output_path <- file.path(output_dir, "metadata_b12_methionine_2.csv")
od_long_output_path <- file.path(output_dir, "od600_long_b12_methionine_2.csv")
joined_output_path <- file.path(output_dir, "od600_with_metadata_b12_methionine_2.csv")

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  stop("Package 'stringr' is required.")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required.")
}

if (!requireNamespace("lubridate", quietly = TRUE)) {
  stop("Package 'lubridate' is required.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_wells <- 24 * 96

well_96 <- data.frame(
  well_letter = rep(LETTERS[1:8], each = 12),
  well_number = rep(1:12, times = 8),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(culture_well = paste0(well_letter, well_number)) |>
  dplyr::select(culture_well)

strainbank_layout <- utils::read.csv(
  layout_metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

replacement_strains <- data.frame(
  array_plate = rep("sequenced_strainbank_array1", times = 3),
  culture_well = c("H10", "H11", "H12"),
  strain_id = c("Ecoli_wt", "Ecoli_metE", "Ecoli_metE_metH"),
  broad_environment = rep("laboratory control", times = 3),
  stringsAsFactors = FALSE
)

strain_layout <- strainbank_layout |>
  dplyr::filter(
    array_plate %in% c(
      "sequenced_strainbank_array1",
      "sequenced_strainbank_array2",
      "marine_strainbank_array2"
    )
  ) |>
  dplyr::rename(culture_well = array_well) |>
  dplyr::mutate(
    broad_environment = dplyr::coalesce(.data$broad_environment, "Unresolved")
  ) |>
  dplyr::select(
    strain_id,
    array_plate,
    culture_well,
    broad_environment,
    environmental_sample_id,
    sample_location,
    date_collected,
    latitude_collected,
    longitude_collected,
    isolation_description,
    growth_medium
  ) |>
  dplyr::anti_join(replacement_strains, by = c("array_plate", "culture_well")) |>
  dplyr::bind_rows(replacement_strains)

solution_composition <- data.frame(
  solution_number = 1:4,
  divalent_salts_gl = rep(0.1, length.out = 4),
  nacl_gl = rep(1, length.out = 4),
  soytone_gl = rep(0.05, length.out = 4),
  methionine_uM = rep(c(0, 100), each = 2, length.out = 4),
  b12_nM = rep(c(0, 10), length.out = 4),
  stringsAsFactors = FALSE
)

metadata_b12_methionine_2_raw <- data.frame(
  culture_plate = rep(
    c(
      "P0000000015", "P0000000014", "P0000000012", "P0000000011",
      "P0000000010", "P0000000009", "P0000000008", "P0000000007",
      "P0000000006", "P0000000005", "P0000000004", "P0000000003",
      "P0000000028", "P0000000027", "P0000000026", "P0000000025",
      "P0000000024", "P0000000023", "P0000000022", "P0000000021",
      "P0000000020", "P0000000019", "P0000000018", "P0000000017"
    ),
    each = 96,
    length.out = n_wells
  ),
  culture_plate_number = rep(1:24, each = 96, length.out = n_wells),
  culture_well = rep(well_96$culture_well, length.out = n_wells),
  solution_number = rep(1:4, each = 96, length.out = n_wells),
  array_plate = rep(
    c(
      "sequenced_strainbank_array1",
      "sequenced_strainbank_array2",
      "marine_strainbank_array2"
    ),
    each = 96 * 4,
    length.out = n_wells
  ),
  technical_replicate = rep(c(1, 2), each = n_wells / 2),
  stringsAsFactors = FALSE
)

metadata_b12_methionine_2 <- metadata_b12_methionine_2_raw |>
  dplyr::left_join(solution_composition, by = "solution_number") |>
  dplyr::left_join(strain_layout, by = c("array_plate", "culture_well"))

plate_reader_files <- list.files(
  plate_reader_dir,
  pattern = "_FB_OD600scan_96well_P.*\\.csv$",
  full.names = TRUE
)

process_plate_reader_file <- function(file) {
  file_name <- basename(file)

  date_time_str <- stringr::str_extract(file_name, "^[0-9]{8}_[0-9]{6}")
  date_time <- lubridate::ymd_hms(date_time_str)

  culture_plate <- stringr::str_match(
    file_name,
    "_FB_OD600scan_96well_([^_]+)"
  )[, 2]

  df_raw <- utils::read.csv(
    file,
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  names(df_raw) <- paste0("X", seq_along(names(df_raw)))

  df_parsed <- df_raw |>
    dplyr::slice(18:dplyr::n()) |>
    dplyr::filter(!is.na(X1)) |>
    dplyr::filter(stringr::str_detect(X1, ":")) |>
    tidyr::separate(
      X1,
      into = c("culture_well", "OD600nm"),
      sep = ":",
      convert = TRUE
    ) |>
    dplyr::mutate(
      culture_well = stringr::str_replace(
        culture_well,
        "([A-H])0([0-9])",
        "\\1\\2"
      ),
      culture_plate = culture_plate,
      date_time = date_time
    )

  df_parsed
}

od600_long <- dplyr::bind_rows(lapply(plate_reader_files, process_plate_reader_file)) |>
  dplyr::mutate(
    growth_hours = as.numeric(difftime(date_time, min(date_time), units = "hours"))
  )

global_blank_mean_od600nm <- od600_long |>
  dplyr::left_join(
    metadata_b12_methionine_2 |>
      dplyr::select(culture_plate, culture_well, strain_id),
    by = c("culture_plate", "culture_well")
  ) |>
  dplyr::filter(strain_id == "blank") |>
  dplyr::summarise(global_blank_mean_od600nm = mean(OD600nm, na.rm = TRUE)) |>
  dplyr::pull(global_blank_mean_od600nm)

od600_with_metadata <- metadata_b12_methionine_2 |>
  dplyr::left_join(od600_long, by = c("culture_plate", "culture_well")) |>
  dplyr::mutate(
    global_blank_mean_od600nm = global_blank_mean_od600nm,
    adjusted_OD600nm = OD600nm - global_blank_mean_od600nm
  ) |>
  dplyr::group_by(culture_plate) |>
  dplyr::arrange(date_time, .by_group = TRUE) |>
  dplyr::mutate(measurement_number = dplyr::dense_rank(date_time)) |>
  dplyr::ungroup()

utils::write.csv(metadata_b12_methionine_2, metadata_output_path, row.names = FALSE, na = "")
utils::write.csv(od600_long, od_long_output_path, row.names = FALSE, na = "")
utils::write.csv(od600_with_metadata, joined_output_path, row.names = FALSE, na = "")

message("Wrote metadata to: ", metadata_output_path)
message("Wrote parsed OD data to: ", od_long_output_path)
message("Wrote joined analysis input to: ", joined_output_path)
message("Global blank mean OD600nm used for adjustment: ", round(global_blank_mean_od600nm, 6))
message("Rows in joined analysis input: ", nrow(od600_with_metadata))
