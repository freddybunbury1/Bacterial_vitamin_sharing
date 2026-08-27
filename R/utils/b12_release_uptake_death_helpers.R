check_required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Required package(s) are not installed: ", paste(missing, collapse = ", "))
  }
}

fix_well_zeros <- function(well) {
  gsub("(?<=^[A-Z])0+(?=[1-9])", "", well, perl = TRUE)
}

map_well_96_to_384 <- function(well_96, position = "TL") {
  row_let <- stringr::str_extract(well_96, "[A-Z]+")
  col_num <- as.numeric(stringr::str_extract(well_96, "[0-9]+"))
  r96 <- utf8ToInt(row_let) - 64

  r_offset <- ifelse(position %in% c("BL", "BR"), 1, 0)
  c_offset <- ifelse(position %in% c("TR", "BR"), 1, 0)

  r384_num <- (r96 * 2 - 1) + r_offset
  c384_num <- (col_num * 2 - 1) + c_offset

  paste0(intToUtf8(r384_num + 64), c384_num)
}

map_well_96_to_384_v <- Vectorize(map_well_96_to_384)

copy_input_files <- function(file_index, overwrite = TRUE) {
  invisible(lapply(
    unique(dirname(file_index$local_path)),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  same_path <- normalizePath(file_index$source_path, mustWork = FALSE) ==
    normalizePath(file_index$local_path, mustWork = FALSE)
  copied <- same_path
  copied[!same_path] <- file.copy(
    file_index$source_path[!same_path],
    file_index$local_path[!same_path],
    overwrite = overwrite
  )
  if (any(!copied)) {
    stop(
      "Failed to copy input file(s): ",
      paste(file_index$source_path[!copied], collapse = ", ")
    )
  }

  source_info <- file.info(file_index$source_path)
  local_info <- file.info(file_index$local_path)

  file_index |>
    dplyr::mutate(
      source_size_bytes = source_info$size,
      source_modified_time = as.character(source_info$mtime),
      local_size_bytes = local_info$size,
      local_modified_time = as.character(local_info$mtime)
    )
}

process_b12_release_uptake_file <- function(file) {
  fname <- basename(file)
  raw <- utils::read.csv(file, header = FALSE, stringsAsFactors = FALSE)

  start_idx <- which(raw$V1 == "Chromatic: 1")
  if (length(start_idx) == 0) {
    return(NULL)
  }

  raw[(start_idx + 1):nrow(raw), 1, drop = FALSE] |>
    dplyr::filter(V1 != "" & !grepl("Cycle:", V1)) |>
    tidyr::separate(V1, into = c("assay_well", "fluorescence"), sep = ":", convert = TRUE) |>
    dplyr::mutate(
      raw_file = fname,
      assay_batch = stringr::str_extract(fname, "(?<=FB_Chlorophyll_384_).*(?=_plate)"),
      assay_plate_number = as.numeric(stringr::str_extract(fname, "(?<=plate)\\d+")),
      measure_time_hours = as.numeric(stringr::str_extract(fname, "(?<=hour)\\d+")),
      assay_well = fix_well_zeros(assay_well)
    ) |>
    dplyr::select(
      raw_file,
      assay_batch,
      assay_plate_number,
      assay_well,
      measure_time_hours,
      fluorescence
    )
}

process_release_uptake_death_od_file <- function(file) {
  fname <- basename(file)

  utils::read.csv(file, header = FALSE, col.names = 1:12) |>
    dplyr::mutate(well_row = LETTERS[1:8]) |>
    tidyr::pivot_longer(cols = 1:12, names_to = "well_col", values_to = "od600") |>
    dplyr::mutate(
      raw_file = fname,
      experiment_name = stringr::str_extract(fname, "b12_release_uptake_death[0-9]+"),
      assay_plate_number = as.numeric(stringr::str_extract(fname, "(?<=plate)\\d+")),
      measure_time_hours = as.numeric(stringr::str_extract(fname, "(?<=hour)\\d+")),
      measurement_minute = as.numeric(stringr::str_extract(fname, "(?<=minute)\\d+")),
      well_col = as.numeric(stringr::str_extract(well_col, "\\d+")),
      assay_well = paste0(well_row, well_col)
    ) |>
    dplyr::select(
      raw_file,
      experiment_name,
      assay_plate_number,
      assay_well,
      measure_time_hours,
      measurement_minute,
      od600
    )
}

process_release_uptake_death_viability_file <- function(file) {
  fname <- basename(file)
  raw <- utils::read.csv(file, header = FALSE, stringsAsFactors = FALSE)

  start_idx <- which(raw$V1 == "Chromatic: 1")
  if (length(start_idx) == 0) {
    return(NULL)
  }

  raw[(start_idx + 1):nrow(raw), 1, drop = FALSE] |>
    dplyr::filter(V1 != "" & !grepl("Cycle:", V1)) |>
    tidyr::separate(V1, into = c("assay_well", "fluorescence"), sep = ":", convert = TRUE) |>
    dplyr::mutate(
      raw_file = fname,
      experiment_name = stringr::str_extract(fname, "b12_release_uptake_death[0-9]+"),
      assay_plate_number = as.numeric(stringr::str_extract(fname, "(?<=plate)\\d+")),
      measurement_minute = as.numeric(stringr::str_extract(fname, "(?<=minute)\\d+")),
      assay_well = fix_well_zeros(assay_well)
    ) |>
    dplyr::select(
      raw_file,
      experiment_name,
      assay_plate_number,
      assay_well,
      measurement_minute,
      fluorescence
    )
}

build_release_uptake_death_culture_samples <- function(strain_layout) {
  wells_96 <- as.vector(t(outer(LETTERS[1:8], 1:12, paste0)))

  ru_n_wells <- 192
  ru_sample_time_minute_list <- c(0.25, 5, 20, 60, 120, 240)
  ru_experiment_name_list <- c("b12_release_uptake_death1", "b12_release_uptake_death2")

  ld_n_wells <- 384
  ld_sample_time_minute_list <- c(60)
  ld_experiment_name_list <- c("b12_release_uptake_death1", "b12_release_uptake_death2")

  live_dead_source_layout <- data.frame(
    experiment_type = rep("live_dead_source", length(wells_96)),
    culture_plate_number = rep(1L, length(wells_96)),
    added_b12_pM = rep(0, length(wells_96)),
    added_glucose_mM_C = rep(0, length(wells_96)),
    biological_replicate = rep(c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L), each = 12),
    viability_status = rep(c("live", "live", "live", "live", "dead", "dead", "dead", "dead"), each = 12),
    killed_proportion = rep(c(0, 0, 0, 0, 1, 1, 1, 1), each = 12),
    strain_id = rep(
      c(
        strain_layout$strain_id[1:12],
        strain_layout$strain_id[13:24],
        strain_layout$strain_id[1:12],
        strain_layout$strain_id[13:24],
        strain_layout$strain_id[1:12],
        strain_layout$strain_id[13:24],
        strain_layout$strain_id[1:12],
        strain_layout$strain_id[13:24]
      ),
      length.out = length(wells_96)
    ),
    culture_well = wells_96,
    stringsAsFactors = FALSE
  )

  ru_culture_layout <- data.frame(
    experiment_type = rep("release_and_uptake", length.out = ru_n_wells),
    culture_plate_number = rep(1:2, each = 96, length.out = ru_n_wells),
    added_b12_pM = rep(c(0, 500), each = 96, length.out = ru_n_wells),
    added_glucose_mM_C = rep(c(0, 10), each = 12, length.out = ru_n_wells),
    biological_replicate = rep(1:2, each = 48, length.out = ru_n_wells),
    killed_proportion = rep(0, length.out = ru_n_wells),
    strain_id = rep(
      as.vector(apply(matrix(strain_layout$strain_id, nrow = 12), 2, rep, 2)),
      length.out = ru_n_wells
    ),
    culture_well = rep(wells_96, length.out = ru_n_wells),
    stringsAsFactors = FALSE
  )

  ru_culture_samples <- ru_culture_layout |>
    tidyr::expand_grid(
      experiment_name = ru_experiment_name_list,
      sample_time_minute = ru_sample_time_minute_list
    )

  ld_culture_layout <- data.frame(
    experiment_type = rep("live_dead_mix", length.out = ld_n_wells),
    culture_plate_number = rep(1:4, each = 96, length.out = ld_n_wells),
    added_b12_pM = rep(0, length.out = ld_n_wells),
    added_glucose_mM_C = rep(0, length.out = ld_n_wells),
    biological_replicate = rep(1:2, each = 192, length.out = ld_n_wells),
    killed_proportion = rep(
      c(1.0, 0.95, 0.88, 0.79, 0.665, 0.505, 0.285, 0),
      each = 12,
      length.out = ld_n_wells
    ),
    strain_id = rep(
      as.vector(apply(matrix(strain_layout$strain_id, nrow = 12), 2, rep, 8)),
      length.out = ld_n_wells
    ),
    culture_well = rep(wells_96, length.out = ld_n_wells),
    stringsAsFactors = FALSE
  )

  ld_culture_samples <- ld_culture_layout |>
    tidyr::expand_grid(
      experiment_name = ld_experiment_name_list,
      sample_time_minute = ld_sample_time_minute_list
    )

  culture_samples <- dplyr::bind_rows(ru_culture_samples, ld_culture_samples) |>
    dplyr::mutate(sample_id = dplyr::row_number())

  list(
    ru_culture_layout = ru_culture_layout,
    ld_culture_layout = ld_culture_layout,
    live_dead_source_layout = live_dead_source_layout,
    culture_samples = culture_samples
  )
}

build_release_uptake_death_viability_plate1_metadata <- function(live_dead_source_layout) {
  live_dead_source_layout |>
    dplyr::mutate(
      source_culture_plate_numbers = as.character(culture_plate_number),
      pre_b12_collapsed = FALSE
    ) |>
    tidyr::crossing(technical_replicate = c(1L, 2L)) |>
    dplyr::mutate(
      assay_plate_number = 1L,
      viability_plate_role = "live_dead_source_pre_mix",
      source_layout_file = "live_dead_source_culture_layout_b12_release_uptake_death.csv",
      sample_time_minute = 0,
      sample_time_hour = sample_time_minute / 60,
      source_row = stringr::str_extract(culture_well, "^[A-Z]+"),
      source_col = as.numeric(stringr::str_extract(culture_well, "\\d+")),
      source_row_num = match(source_row, LETTERS),
      assay_row_num = (source_row_num * 2L) - 2L + technical_replicate,
      assay_col = (source_col * 2L) - 1L,
      assay_well = paste0(LETTERS[assay_row_num], assay_col)
    ) |>
    tidyr::expand_grid(experiment_name = c("b12_release_uptake_death1", "b12_release_uptake_death2")) |>
    dplyr::select(
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
    )
}

build_release_uptake_death_viability_plate2_metadata <- function(ld_culture_layout) {
  ld_culture_layout |>
    dplyr::mutate(
      assay_plate_number = 2L,
      viability_plate_role = "live_dead_mix_60min",
      source_layout_file = "live_dead_mix_culture_layout_b12_release_uptake_death.csv",
      source_culture_plate_numbers = as.character(culture_plate_number),
      pre_b12_collapsed = FALSE,
      technical_replicate = NA_integer_,
      viability_status = "live",
      source_row = stringr::str_extract(culture_well, "^[A-Z]+"),
      source_col = as.numeric(stringr::str_extract(culture_well, "\\d+")),
      source_row_num = match(source_row, LETTERS),
      assay_row_num = dplyr::case_when(
        culture_plate_number %in% c(1L, 3L) ~ source_row_num,
        culture_plate_number %in% c(2L, 4L) ~ source_row_num + 8L
      ),
      assay_col = dplyr::case_when(
        culture_plate_number %in% c(1L, 2L) ~ (source_col * 2L) - 1L,
        culture_plate_number %in% c(3L, 4L) ~ source_col * 2L
      ),
      assay_well = paste0(LETTERS[assay_row_num], assay_col),
      sample_time_minute = 60,
      sample_time_hour = sample_time_minute / 60
    ) |>
    tidyr::expand_grid(experiment_name = c("b12_release_uptake_death1", "b12_release_uptake_death2")) |>
    dplyr::select(
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
    )
}

build_release_uptake_death_viability_plate3_metadata <- function(ru_culture_layout) {
  ru_culture_layout |>
    dplyr::mutate(
      assay_plate_number = 3L,
      viability_plate_role = "release_uptake_240min",
      source_layout_file = "release_uptake_culture_layout_b12_release_uptake_death.csv",
      source_culture_plate_numbers = as.character(culture_plate_number),
      pre_b12_collapsed = FALSE,
      technical_replicate = NA_integer_,
      source_row = stringr::str_extract(culture_well, "^[A-Z]+"),
      source_col = as.numeric(stringr::str_extract(culture_well, "\\d+")),
      source_row_num = match(source_row, LETTERS),
      assay_row_num = dplyr::case_when(
        culture_plate_number == 1L ~ source_row_num,
        culture_plate_number == 2L ~ source_row_num + 8L
      ),
      sample_time_minute = 240,
      sample_time_hour = sample_time_minute / 60
    ) |>
    tidyr::crossing(viability_status = c("live", "dead")) |>
    dplyr::mutate(
      assay_col = dplyr::case_when(
        viability_status == "live" ~ (source_col * 2L) - 1L,
        viability_status == "dead" ~ source_col * 2L
      ),
      assay_well = paste0(LETTERS[assay_row_num], assay_col)
    ) |>
    tidyr::expand_grid(experiment_name = c("b12_release_uptake_death1", "b12_release_uptake_death2")) |>
    dplyr::select(
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
    )
}

build_release_uptake_death_viability_metadata <- function(
  ru_culture_layout,
  ld_culture_layout,
  live_dead_source_layout
) {
  plate1_metadata <- build_release_uptake_death_viability_plate1_metadata(live_dead_source_layout)
  plate2_metadata <- build_release_uptake_death_viability_plate2_metadata(ld_culture_layout)
  plate3_metadata <- build_release_uptake_death_viability_plate3_metadata(ru_culture_layout)

  dplyr::bind_rows(plate1_metadata, plate2_metadata, plate3_metadata)
}

build_release_uptake_death_od_metadata <- function(ru_culture_layout, live_dead_source_layout) {
  live_dead_source_od_layout <- live_dead_source_layout |>
    dplyr::mutate(
      assay_plate_number = 1L,
      od_plate_role = "live_dead_source_pre_mix",
      source_layout_file = "live_dead_source_culture_layout_b12_release_uptake_death.csv",
      source_culture_plate_numbers = as.character(culture_plate_number),
      pre_b12_collapsed = FALSE,
      sample_time_minute = 0,
      sample_time_hour = sample_time_minute / 60
    )

  release_final_layout <- ru_culture_layout |>
    dplyr::filter(added_b12_pM == 0) |>
    dplyr::mutate(
      assay_plate_number = 2L,
      od_plate_role = "release_final_240min",
      source_layout_file = "release_uptake_culture_layout_b12_release_uptake_death.csv",
      source_culture_plate_numbers = as.character(culture_plate_number),
      pre_b12_collapsed = FALSE,
      viability_status = NA_character_,
      sample_time_minute = 240,
      sample_time_hour = sample_time_minute / 60
    )

  uptake_final_layout <- ru_culture_layout |>
    dplyr::filter(added_b12_pM == 500) |>
    dplyr::mutate(
      assay_plate_number = 3L,
      od_plate_role = "uptake_final_240min",
      source_layout_file = "release_uptake_culture_layout_b12_release_uptake_death.csv",
      source_culture_plate_numbers = as.character(culture_plate_number),
      pre_b12_collapsed = FALSE,
      viability_status = NA_character_,
      sample_time_minute = 240,
      sample_time_hour = sample_time_minute / 60
    )

  dplyr::bind_rows(live_dead_source_od_layout, release_final_layout, uptake_final_layout) |>
    tidyr::expand_grid(
      experiment_name = c("b12_release_uptake_death1", "b12_release_uptake_death2")
    ) |>
    dplyr::mutate(
      expected_measure_time_hours = dplyr::if_else(assay_plate_number == 1L, 20, 20),
      expected_measurement_minute = dplyr::case_when(
        assay_plate_number == 1L ~ 0,
        assay_plate_number %in% c(2L, 3L) ~ 240,
        TRUE ~ NA_real_
      ),
      assay_well = culture_well
    ) |>
    dplyr::select(
      experiment_name,
      assay_plate_number,
      assay_well,
      od_plate_role,
      source_layout_file,
      expected_measure_time_hours,
      expected_measurement_minute,
      source_culture_plate_numbers,
      pre_b12_collapsed,
      viability_status,
      experiment_type,
      culture_plate_number,
      culture_well,
      sample_time_minute,
      sample_time_hour,
      strain_id,
      biological_replicate,
      killed_proportion,
      added_b12_pM,
      added_glucose_mM_C
    )
}

build_b12_storage_layout <- function(culture_samples) {
  culture_samples |>
    dplyr::mutate(
      temp_row = stringr::str_extract(culture_well, "^[A-Z]+"),
      temp_col = as.numeric(stringr::str_extract(culture_well, "\\d+"))
    ) |>
    dplyr::arrange(
      experiment_name,
      dplyr::desc(experiment_type),
      culture_plate_number,
      sample_time_minute,
      temp_row,
      temp_col
    ) |>
    dplyr::mutate(
      plate_96_idx = (dplyr::row_number() - 1) %/% 96,
      storage_plate_number = (plate_96_idx %/% 4) + 1,
      quad_idx = plate_96_idx %% 4,
      storage_quadrant = dplyr::case_when(
        quad_idx == 0 ~ "TL",
        quad_idx == 1 ~ "BL",
        quad_idx == 2 ~ "TR",
        quad_idx == 3 ~ "BR"
      ),
      storage_well = fix_well_zeros(map_well_96_to_384_v(culture_well, storage_quadrant)),
      storage_dilution_factor = dplyr::if_else(experiment_type == "release_and_uptake", 4, 5)
    ) |>
    dplyr::select(-quad_idx, -plate_96_idx, -temp_row, -temp_col)
}

build_b12_assay_layout <- function(b12_storage_layout) {
  assay_quad_map <- list(
    TL = c("TL", "TR"),
    BL = c("BL", "BR"),
    TR = c("TL", "TR"),
    BR = c("BL", "BR")
  )

  b12_storage_layout |>
    dplyr::mutate(
      assay_plate_number = dplyr::case_when(
        storage_quadrant %in% c("TL", "BL") ~ (storage_plate_number * 2) - 1,
        storage_quadrant %in% c("TR", "BR") ~ storage_plate_number * 2
      ),
      assay_quadrant = lapply(storage_quadrant, function(q) assay_quad_map[[q]])
    ) |>
    tidyr::unnest(assay_quadrant) |>
    dplyr::mutate(
      assay_well = fix_well_zeros(map_well_96_to_384_v(culture_well, assay_quadrant)),
      sample_volume_added = dplyr::if_else(assay_quadrant %in% c("TL", "BL"), 10, 30),
      assay_volume_added = 100 - sample_volume_added,
      b12_assay_strain = "metE4",
      assay_type = "sample",
      sample_b12_conc_pm = NA_real_,
      standard_bacteria_strain_id = NA_character_,
      standard_processing = NA_character_,
      sample_time_hour = sample_time_minute / 60,
      sample_type = "filtrate"
    )
}

build_b12_release_uptake_standards <- function(assay_batches, standard_concs_pm) {
  assay_well_list <- expand.grid(
    well_column = 1:24,
    well_row = LETTERS[1:16],
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(assay_well = paste0(well_row, well_column)) |>
    dplyr::arrange(well_row, well_column)

  tidyr::expand_grid(
    assay_batch = assay_batches,
    standard_processing = c("heat_then_dilute", "dilute_then_heat"),
    assay_well = assay_well_list$assay_well
  ) |>
    dplyr::group_by(assay_batch, standard_processing) |>
    dplyr::mutate(
      assay_plate_number = dplyr::if_else(standard_processing == "heat_then_dilute", 17L, 18L),
      assay_type = "standard",
      b12_assay_strain = "metE4",
      storage_dilution_factor = 1,
      sample_type = rep(c("filtrate", "total"), each = 192, length.out = dplyr::n()),
      sample_b12_conc_pm = rep(standard_concs_pm, each = 2, length.out = dplyr::n()),
      sample_volume_added = rep(c(10, 30), length.out = dplyr::n()),
      assay_volume_added = 100 - sample_volume_added,
      standard_bacteria_strain_id = rep(
        c("blank", "HTWP4", "sic1606", "OTU4908"),
        each = 48,
        length.out = dplyr::n()
      ),
      experiment_type = NA_character_,
      isolation_source = NA_character_,
      strain_id = NA_character_,
      sample_id = NA_integer_,
      sample_time_minute = NA_real_,
      sample_time_hour = NA_real_,
      culture_plate_number = NA_integer_,
      culture_well = NA_character_,
      storage_plate_number = NA_integer_,
      storage_well = NA_character_,
      storage_quadrant = NA_character_,
      added_b12_pM = NA_real_,
      added_glucose_mM_C = NA_real_,
      biological_replicate = NA_integer_,
      killed_proportion = NA_real_,
      experiment_name = NA_character_
    ) |>
    dplyr::ungroup()
}
