check_required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Required package(s) are not installed: ", paste(missing, collapse = ", "))
  }
}

extract_numeric <- function(x) {
  as.integer(stringr::str_extract(x, "\\d+"))
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
  copied <- file.copy(file_index$source_path, file_index$local_path, overwrite = overwrite)
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

build_culture_metadata <- function(strainbank, experiment_name_list, sample_time_hour_list) {
  glycerol_array_to_culture <- tibble::tibble(
    array_plate = c(
      "sequenced_strainbank_array1",
      "sequenced_strainbank_array2",
      "marine_strainbank_array2"
    ),
    culture_plate_number = c(1L, 2L, 3L)
  )

  replacement_strains <- tibble::tibble(
    array_plate = rep("sequenced_strainbank_array1", 3),
    culture_plate_number = 1L,
    culture_well = c("H10", "H11", "H12"),
    strain_id = c("Ecoli_wt", "Ecoli_metE", "Ecoli_metE_metH"),
    isolation_source = "other"
  )

  culture_layout <- strainbank |>
    dplyr::filter(
      array_plate %in% c(
        "sequenced_strainbank_array1",
        "sequenced_strainbank_array2",
        "marine_strainbank_array2"
      )
    ) |>
    dplyr::rename(culture_well = array_well) |>
    dplyr::left_join(glycerol_array_to_culture, by = "array_plate") |>
    dplyr::mutate(
      isolation_source = ifelse(stringr::str_detect(array_plate, "marine"), "marine", "soil")
    ) |>
    dplyr::select(
      isolation_source,
      strain_id,
      array_plate,
      culture_plate_number,
      culture_well
    ) |>
    dplyr::anti_join(replacement_strains, by = c("array_plate", "culture_well")) |>
    dplyr::bind_rows(replacement_strains) |>
    dplyr::arrange(array_plate, culture_well) |>
    dplyr::mutate(layout_row_id = dplyr::row_number())

  culture_samples <- culture_layout |>
    tidyr::expand_grid(
      experiment_name = experiment_name_list,
      sample_time_hour = sample_time_hour_list
    ) |>
    dplyr::mutate(sample_id = dplyr::row_number())

  list(
    culture_layout = culture_layout,
    culture_samples = culture_samples,
    replacement_strains = replacement_strains,
    glycerol_array_to_culture = glycerol_array_to_culture
  )
}

build_od_viability_metadata <- function(culture_samples) {
  plate_to_quadrant <- tibble::tibble(
    culture_plate_number = c(1L, 2L, 3L),
    assay_quadrant = c("TL", "BL", "TR")
  )

  culture_samples |>
    dplyr::inner_join(plate_to_quadrant, by = "culture_plate_number") |>
    dplyr::mutate(assay_well = map_well_96_to_384_v(culture_well, assay_quadrant)) |>
    tidyr::crossing(assay_plate_number = c(1L, 2L)) |>
    dplyr::mutate(viability_status = dplyr::if_else(assay_plate_number == 1L, "live", "dead")) |>
    dplyr::select(
      sample_id,
      experiment_name,
      sample_time_hour,
      culture_plate_number,
      culture_well,
      strain_id,
      isolation_source,
      assay_plate_number,
      assay_well,
      viability_status
    )
}

build_b12_storage_metadata <- function(culture_samples) {
  quadrants <- c("TL", "TR", "BL", "BR")

  quadrant_to_sample_type <- tibble::tibble(
    quadrant = quadrants,
    sample_type = c("filtrate", "total", "uptake_1hr", "uptake_0hr")
  )

  culture_and_time_to_storage <- tibble::tibble(
    culture_plate_number = rep(1:3, times = 2),
    sample_time_hour = rep(c(20, 44), each = 3),
    storage_plate_number = 1:6
  )

  culture_samples |>
    tidyr::crossing(quadrant = quadrants) |>
    dplyr::mutate(
      storage_well = map_well_96_to_384_v(culture_well, quadrant),
      storage_dilution_factor = dplyr::if_else(quadrant == "TR", 5, 2)
    ) |>
    dplyr::left_join(quadrant_to_sample_type, by = "quadrant") |>
    dplyr::left_join(
      culture_and_time_to_storage,
      by = c("culture_plate_number", "sample_time_hour")
    ) |>
    dplyr::select(
      sample_id,
      experiment_name,
      sample_time_hour,
      culture_plate_number,
      culture_well,
      strain_id,
      isolation_source,
      quadrant,
      sample_type,
      storage_plate_number,
      storage_well,
      storage_dilution_factor
    )
}

build_b12_assay_sample_metadata <- function(b12_storage_map) {
  assay_quad_map <- list(
    TL = c("TL", "TR"),
    TR = c("BL", "BR"),
    BL = c("TL", "TR"),
    BR = c("BL", "BR")
  )

  b12_storage_map |>
    dplyr::mutate(
      assay_plate_number = dplyr::case_when(
        quadrant %in% c("TL", "TR") ~ (storage_plate_number * 2) - 1,
        quadrant %in% c("BL", "BR") ~ storage_plate_number * 2
      ),
      assay_quadrant = lapply(quadrant, function(q) assay_quad_map[[q]])
    ) |>
    tidyr::unnest(assay_quadrant) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      s_row_num = utf8ToInt(stringr::str_extract(storage_well, "[A-P]")) - 64,
      s_col_num = as.numeric(stringr::str_extract(storage_well, "\\d+")),
      r96 = ceiling(s_row_num / 2),
      c96 = ceiling(s_col_num / 2),
      r_offset = dplyr::if_else(assay_quadrant %in% c("BL", "BR"), 1, 0),
      c_offset = dplyr::if_else(assay_quadrant %in% c("TR", "BR"), 1, 0),
      a_row_num = (r96 * 2 - 1) + r_offset,
      a_col_num = (c96 * 2 - 1) + c_offset,
      assay_well = paste0(intToUtf8(a_row_num + 64), a_col_num)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sample_volume_added = dplyr::case_when(
        assay_quadrant %in% c("TL", "BL") ~ 10,
        assay_quadrant %in% c("TR", "BR") ~ 30
      ),
      assay_volume_added = 100 - sample_volume_added,
      b12_assay_strain = "metE4",
      assay_type = "sample",
      sample_b12_conc_pm = NA_real_,
      standard_bacteria_strain_id = NA_character_
    ) |>
    dplyr::select(
      sample_id,
      experiment_name,
      sample_time_hour,
      sample_type,
      storage_plate_number,
      storage_well,
      storage_dilution_factor,
      assay_plate_number,
      assay_well,
      assay_type,
      sample_volume_added,
      assay_volume_added,
      b12_assay_strain,
      strain_id,
      isolation_source,
      culture_plate_number,
      culture_well,
      sample_b12_conc_pm,
      standard_bacteria_strain_id
    )
}

build_b12_standard_metadata <- function(experiment_name_list, standard_concs_pm) {
  assay_well <- expand.grid(
    well_column = 1:24,
    well_row = LETTERS[1:16],
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(assay_well = paste0(well_row, well_column)) |>
    dplyr::arrange(well_row, well_column)

  tidyr::expand_grid(
    experiment_name = experiment_name_list,
    assay_well = assay_well$assay_well
  ) |>
    dplyr::group_by(experiment_name) |>
    dplyr::mutate(
      assay_batch = experiment_name,
      assay_plate_number = 13L,
      assay_type = "standard",
      b12_assay_strain = "metE4",
      storage_dilution_factor = 1,
      sample_type = rep(c("filtrate", "total"), each = 192, length.out = dplyr::n()),
      sample_b12_conc_pm = rep(standard_concs_pm, each = 2, length.out = dplyr::n()),
      sample_volume_added = rep(c(10L, 30L), length.out = dplyr::n()),
      assay_volume_added = 100L - sample_volume_added,
      standard_bacteria_strain_id = rep(
        c("blank", "HTWP4", "sic1606", "OTU4908"),
        each = 48,
        length.out = dplyr::n()
      ),
      isolation_source = NA_character_,
      strain_id = NA_character_,
      sample_id = NA_integer_,
      sample_time_hour = NA_real_,
      culture_plate_number = NA_integer_,
      culture_well = NA_character_,
      storage_plate_number = NA_integer_,
      storage_well = NA_character_
    ) |>
    dplyr::ungroup()
}

process_od_grid <- function(file) {
  fname <- basename(file)

  utils::read.csv(file, header = FALSE, col.names = 1:24) |>
    dplyr::mutate(well_row = LETTERS[1:16]) |>
    tidyr::pivot_longer(cols = 1:24, names_to = "well_col", values_to = "od600") |>
    dplyr::mutate(
      experiment_name = stringr::str_extract(fname, "(?<=well_).*(?=_plate)"),
      assay_plate_number = as.numeric(stringr::str_extract(fname, "(?<=plate)\\d+")),
      sample_time_hour = as.numeric(stringr::str_extract(fname, "(?<=hour)\\d+")),
      well_col = extract_numeric(well_col),
      assay_well = paste0(well_row, well_col),
      raw_file = fname
    ) |>
    dplyr::select(
      raw_file,
      experiment_name,
      assay_plate_number,
      assay_well,
      sample_time_hour,
      od600
    )
}

process_vantastar_list <- function(file, experiment_regex, experiment_prefix_regex) {
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
      experiment_name = stringr::str_extract(fname, experiment_regex),
      experiment_name = stringr::str_remove(experiment_name, experiment_prefix_regex),
      experiment_name = stringr::str_remove(experiment_name, "_plate$"),
      assay_batch = experiment_name,
      assay_plate_number = as.numeric(stringr::str_extract(fname, "(?<=plate)\\d+")),
      sample_time_hour = as.numeric(stringr::str_extract(fname, "(?<=hour)\\d+")),
      measure_time_hours = sample_time_hour,
      measurement_minute = as.numeric(stringr::str_extract(fname, "(?<=minute)\\d+")),
      assay_well = fix_well_zeros(assay_well),
      raw_file = fname
    ) |>
    dplyr::select(
      raw_file,
      assay_batch,
      experiment_name,
      assay_plate_number,
      assay_well,
      sample_time_hour,
      measure_time_hours,
      measurement_minute,
      fluorescence
    )
}

geo_mean <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) {
    return(NA_real_)
  }
  exp(mean(log(x)))
}

geo_sd <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) {
    return(NA_real_)
  }
  exp(stats::sd(log(x)))
}
