#!/usr/bin/env Rscript

# Single-stage KO-based B12 trait classification. Random validation is matched
# to the GTDB-order-block design: for each trait, one stratified k-fold
# partition is constructed with k equal to the number of order blocks, and each
# mutually exclusive test fold is evaluated exactly once.

ko_matrix_path <- file.path("data/raw/genome_pipeline_nf_core", "ko_count_matrix.tsv")
phenotype_path <- file.path(
  "data/processed/b12_trait_survey", "b12_trait_survey_gwas_phenotypes.csv"
)
trait_summary_path <- file.path(
  "data/processed/b12_trait_survey", "b12_trait_survey_strain_summary.csv"
)
taxonomy_path <- file.path(
  "data/intermediate/gtdbtk", "gtdbtk_bac120_taxonomy_with_strain_id.csv"
)
tree_path <- file.path(
  "data/processed/gtdbtk_tree_analysis",
  "gtdbtk_tree_231_genomes_bacillati_pseudomonadati_split_rooted_b12_complete_traits.tree"
)
kegg_ko_list_path <- file.path("data/raw/annotations", "kegg_ko_list.tsv")
raw_gwa_dir <- file.path("data/raw", "gwa")

processed_dir <- file.path("data/processed", "b12_trait_classification_analysis")
figure_source_dir <- file.path(
  "data/processed/figure_source_data", "b12_trait_classification_analysis"
)
model_input_path <- file.path(processed_dir, "b12_trait_classification_model_input.csv")
group_metadata_path <- file.path(processed_dir, "order_block_group_metadata.csv")
metrics_path <- file.path(processed_dir, "classification_metrics.csv")
predictions_path <- file.path(processed_dir, "classification_predictions.csv")
importance_path <- file.path(processed_dir, "classification_feature_importance.csv")
random_importance_path <- file.path(processed_dir, "random_fold_feature_importance_summary.csv")
null_metrics_path <- file.path(processed_dir, "train_prevalence_null_metrics.csv")
null_summary_path <- file.path(processed_dir, "null_mean_performance.csv")
parameter_path <- file.path(processed_dir, "analysis_parameters.csv")
ko_mapping_path <- file.path(processed_dir, "ko_gene_name_mapping.csv")
gwa_comparison_path <- file.path(
  processed_dir, "gwa_pyseer_mixed_tree_comparison_all_traits.csv"
)

random_seed <- 20260714L
total_b12_threshold <- log10(50)
uptake_b12_threshold <- 250
min_train_feature_prevalence <- 0.02
minimum_order_test_fraction <- 0.05
null_replicates <- 1000L
rf_num_trees <- 1000L
rf_min_node_size <- 1L
rf_sample_fraction <- 0.632
rf_replace <- FALSE
detected_cores <- parallel::detectCores(logical = TRUE)
rf_num_threads <- if (is.na(detected_cores)) 1L else max(1L, min(4L, detected_cores - 1L))
model_engine <- "ranger_random_forest_single_stage"

required_packages <- c("ape", "dplyr", "purrr", "ranger", "stringr", "tibble", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Required packages are missing: ", paste(missing_packages, collapse = ", "))
}

required_files <- c(
  ko_matrix_path, phenotype_path, trait_summary_path, taxonomy_path, tree_path,
  kegg_ko_list_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required classification inputs are missing:\n",
    paste(missing_files, collapse = "\n"),
    if (kegg_ko_list_path %in% missing_files) {
      paste0(
        "\nThe KEGG list is the intentional acquisition exception. Academic users ",
        "must first run: Rscript --vanilla R/cleaning/download_kegg_ko_mapping.R"
      )
    } else ""
  )
}

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}
write_csv <- function(data, path) {
  utils::write.csv(data, path, row.names = FALSE, na = "")
}
normalize_join_key <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace("\\.(fa|fna|fasta)$", "") |>
    stringr::str_replace("_genomic$", "") |>
    stringr::str_replace_all("[^A-Za-z0-9]", "") |>
    stringr::str_to_upper()
}

make_stratified_folds <- function(labels, k) {
  if (k < 2L) stop("At least two folds are required.")
  fold_id <- rep(NA_integer_, length(labels))
  for (level in levels(factor(labels))) {
    index <- sample(which(labels == level))
    fold_id[index] <- rep(seq_len(k), length.out = length(index))
  }
  if (anyNA(fold_id)) stop("Fold assignment left observations unassigned.")
  fold_id
}

safe_balanced_accuracy <- function(truth, prediction) {
  recalls <- vapply(levels(truth), function(level) {
    denominator <- sum(truth == level, na.rm = TRUE)
    if (denominator == 0L) return(NA_real_)
    sum(truth == level & prediction == level, na.rm = TRUE) / denominator
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

simulate_train_prevalence_null <- function(truth, train_y) {
  class_levels <- levels(train_y)
  probabilities <- as.numeric(table(train_y)[class_levels])
  probabilities[is.na(probabilities)] <- 0
  probabilities <- probabilities / sum(probabilities)
  purrr::map_dfr(seq_len(null_replicates), function(replicate_id) {
    prediction <- factor(
      sample(class_levels, length(truth), replace = TRUE, prob = probabilities),
      levels = class_levels
    )
    tibble::tibble(
      null_replicate = replicate_id,
      null_accuracy = mean(prediction == truth),
      null_balanced_accuracy = safe_balanced_accuracy(truth, prediction)
    )
  })
}

make_order_holdout_groups <- function(data, tree_scope) {
  minimum_size <- ceiling(minimum_order_test_fraction * nrow(data))
  order_counts <- data |>
    dplyr::mutate(gtdb_order = dplyr::coalesce(.data$gtdb_order, "Unassigned")) |>
    dplyr::count(.data$gtdb_order, name = "n_genomes")
  large_orders <- order_counts |>
    dplyr::filter(.data$n_genomes >= minimum_size) |>
    dplyr::arrange(dplyr::desc(.data$n_genomes), .data$gtdb_order)
  small_orders <- order_counts |>
    dplyr::filter(.data$n_genomes < minimum_size) |>
    dplyr::arrange(dplyr::desc(.data$n_genomes), .data$gtdb_order)
  if (nrow(large_orders) == 0L) stop("No GTDB order meets the minimum holdout size.")

  membership <- large_orders |>
    dplyr::transmute(
      gtdb_order = .data$gtdb_order,
      split_id = paste0(
        "order_", stringr::str_replace_all(.data$gtdb_order, "[^A-Za-z0-9]+", "_")
      )
    )
  if (nrow(small_orders) > 0L) {
    if (sum(small_orders$n_genomes) >= minimum_size) {
      membership <- dplyr::bind_rows(
        membership,
        dplyr::transmute(small_orders, gtdb_order = .data$gtdb_order,
                         split_id = "pooled_small_orders")
      )
    } else {
      smallest_valid <- large_orders |>
        dplyr::arrange(.data$n_genomes, .data$gtdb_order) |>
        dplyr::slice(1L) |>
        dplyr::pull(.data$gtdb_order)
      target_id <- membership$split_id[membership$gtdb_order == smallest_valid]
      membership <- dplyr::bind_rows(
        membership,
        dplyr::transmute(small_orders, gtdb_order = .data$gtdb_order,
                         split_id = target_id)
      )
    }
  }

  data |>
    dplyr::mutate(gtdb_order = dplyr::coalesce(.data$gtdb_order, "Unassigned")) |>
    dplyr::select("strain_id", "gtdb_order") |>
    dplyr::left_join(membership, by = "gtdb_order") |>
    dplyr::group_by(.data$split_id) |>
    dplyr::summarise(
      tree_scope = tree_scope,
      held_out_orders = paste(sort(unique(.data$gtdb_order)), collapse = ";"),
      n_test = dplyr::n(),
      test_fraction = dplyr::n() / nrow(data),
      test_tip_labels = paste(sort(.data$strain_id), collapse = ";"),
      test_tips = list(sort(.data$strain_id)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$n_test), .data$split_id)
}

fit_classifier <- function(train_data, test_data, target_col, feature_cols) {
  train_y <- factor(train_data[[target_col]], levels = c("below", "above"))
  test_y <- factor(test_data[[target_col]], levels = c("below", "above"))
  prevalence <- colMeans(train_data[, feature_cols, drop = FALSE] > 0, na.rm = TRUE)
  kept <- names(prevalence)[
    prevalence >= min_train_feature_prevalence &
      prevalence <= (1 - min_train_feature_prevalence)
  ]
  if (length(kept) == 0L) kept <- feature_cols
  variable <- vapply(
    train_data[, kept, drop = FALSE],
    function(x) length(unique(x[!is.na(x)])) > 1L,
    logical(1)
  )
  kept <- kept[variable]
  if (length(kept) == 0L) stop("No variable KO features remain after training filtering.")

  mtry <- max(1L, floor(sqrt(length(kept))))
  fit <- ranger::ranger(
    response ~ .,
    data = data.frame(response = train_y, train_data[, kept, drop = FALSE],
                      check.names = FALSE),
    classification = TRUE,
    importance = "permutation",
    num.trees = rf_num_trees,
    mtry = mtry,
    min.node.size = rf_min_node_size,
    sample.fraction = rf_sample_fraction,
    replace = rf_replace,
    num.threads = rf_num_threads,
    seed = random_seed + length(kept) + 1000L
  )
  prediction <- factor(
    stats::predict(
      fit,
      data = data.frame(test_data[, kept, drop = FALSE], check.names = FALSE)
    )$predictions,
    levels = levels(train_y)
  )
  list(
    truth = test_y,
    prediction = prediction,
    n_features = length(kept),
    mtry = mtry,
    importance = tibble::tibble(
      ko = names(fit$variable.importance),
      importance = as.numeric(fit$variable.importance)
    ) |>
      dplyr::filter(is.finite(.data$importance))
  )
}

run_split <- function(data, target_col, split_type, split_id,
                      train_index, test_index, feature_cols) {
  train_data <- data[train_index, , drop = FALSE]
  test_data <- data[test_index, , drop = FALSE]
  train_y <- factor(train_data[[target_col]], levels = c("below", "above"))
  if (dplyr::n_distinct(train_y) < 2L || nrow(test_data) == 0L) return(NULL)
  result <- fit_classifier(train_data, test_data, target_col, feature_cols)
  majority <- names(sort(table(train_y), decreasing = TRUE))[[1]]
  baseline <- factor(rep(majority, nrow(test_data)), levels = levels(train_y))

  list(
    metrics = tibble::tibble(
      trait = target_col, split_type = split_type, split_id = split_id,
      model_engine = model_engine,
      n_train = nrow(train_data), n_test = nrow(test_data),
      n_train_above = sum(train_y == "above"),
      n_test_above = sum(result$truth == "above"),
      n_features = result$n_features, mtry = result$mtry,
      accuracy = mean(result$prediction == result$truth),
      balanced_accuracy = safe_balanced_accuracy(result$truth, result$prediction),
      baseline_accuracy = mean(baseline == result$truth),
      baseline_balanced_accuracy = safe_balanced_accuracy(result$truth, baseline)
    ),
    predictions = tibble::tibble(
      trait = target_col, split_type = split_type, split_id = split_id,
      strain_id = test_data$strain_id,
      truth = as.character(result$truth),
      prediction = as.character(result$prediction),
      correct = result$prediction == result$truth
    ),
    importance = result$importance |>
      dplyr::mutate(
        trait = target_col, split_type = split_type, split_id = split_id,
        model_engine = model_engine
      ) |>
      dplyr::relocate("trait", "split_type", "split_id", "model_engine", "ko"),
    null_metrics = simulate_train_prevalence_null(result$truth, train_y) |>
      dplyr::mutate(
        trait = target_col, split_type = split_type, split_id = split_id,
        n_train = nrow(train_data), n_test = nrow(test_data),
        train_above_probability = mean(train_y == "above"),
        n_test_above = sum(result$truth == "above")
      ) |>
      dplyr::relocate("trait", "split_type", "split_id", "null_replicate")
  )
}

safe_neg_log10 <- function(pvalue) {
  dplyr::case_when(
    is.na(pvalue) | pvalue < 0 ~ NA_real_,
    pvalue == 0 ~ Inf,
    TRUE ~ -log10(pvalue)
  )
}

combine_gwa_pair <- function(mixed_path) {
  source_prefix <- stringr::str_remove(
    basename(mixed_path), "_pyseer_mixed_summary[.]csv$"
  )
  trait <- stringr::str_remove(source_prefix, "_MAF_.*$")
  tree_file <- file.path(
    raw_gwa_dir, paste0(source_prefix, "_pyseer_mixed_tree_summary.csv")
  )
  if (!file.exists(tree_file)) stop("Missing paired tree-aware GWA file: ", tree_file)
  mixed <- read_csv(mixed_path) |>
    dplyr::mutate(variant = as.character(.data$variant), present_without_tree = TRUE)
  mixed_tree <- read_csv(tree_file) |>
    dplyr::mutate(variant = as.character(.data$variant), present_with_tree = TRUE)
  dplyr::full_join(mixed, mixed_tree, by = "variant",
                   suffix = c("_without_tree", "_with_tree")) |>
    dplyr::mutate(
      trait = trait,
      source_prefix = source_prefix,
      source_file_without_tree = basename(mixed_path),
      source_file_with_tree = basename(tree_file),
      present_without_tree = dplyr::coalesce(.data$present_without_tree, FALSE),
      present_with_tree = dplyr::coalesce(.data$present_with_tree, FALSE),
      present_in_both_models = .data$present_without_tree & .data$present_with_tree,
      neg_log10_lrt_pvalue_without_tree = safe_neg_log10(.data$`lrt.pvalue_without_tree`),
      neg_log10_lrt_pvalue_with_tree = safe_neg_log10(.data$`lrt.pvalue_with_tree`)
    ) |>
    dplyr::relocate(
      "trait", "source_prefix", "variant", "present_in_both_models",
      "neg_log10_lrt_pvalue_without_tree", "neg_log10_lrt_pvalue_with_tree"
    ) |>
    dplyr::arrange(.data$trait, .data$variant)
}

# Construct the one canonical KO/phenotype/taxonomy model table.
ko_matrix <- utils::read.delim(
  ko_matrix_path, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
if (!"genome" %in% names(ko_matrix)) stop("KO matrix must contain a genome column.")
ko_features <- names(ko_matrix)[stringr::str_detect(names(ko_matrix), "^K[0-9]{5}$")]
if (length(ko_features) == 0L) stop("KO matrix contains no K##### feature columns.")
ko_matrix <- ko_matrix |>
  dplyr::mutate(
    strain_id = as.character(.data$genome),
    strain_join_key = normalize_join_key(.data$strain_id)
  ) |>
  dplyr::select("strain_id", "strain_join_key", dplyr::all_of(ko_features)) |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(ko_features), ~ as.integer(dplyr::coalesce(as.numeric(.x), 0) > 0))
  ) |>
  dplyr::group_by(.data$strain_join_key) |>
  dplyr::summarise(
    strain_id_ko = dplyr::first(.data$strain_id),
    dplyr::across(dplyr::all_of(ko_features), max),
    .groups = "drop"
  )

phenotypes <- read_csv(phenotype_path) |>
  dplyr::mutate(strain_join_key = normalize_join_key(.data$strain_id))
taxonomy <- read_csv(taxonomy_path) |>
  dplyr::mutate(strain_join_key = normalize_join_key(.data$strain_id)) |>
  dplyr::select("strain_join_key", "gtdb_phylum", "gtdb_class", "gtdb_order") |>
  dplyr::distinct(.data$strain_join_key, .keep_all = TRUE)

model_input <- phenotypes |>
  dplyr::inner_join(taxonomy, by = "strain_join_key") |>
  dplyr::inner_join(ko_matrix, by = "strain_join_key") |>
  dplyr::filter(
    !is.na(.data$gtdb_class),
    is.finite(.data$log10_total_b12_gm),
    is.finite(.data$uptake_b12_mean)
  ) |>
  dplyr::mutate(
    log10_total_b12_gm_high = factor(
      dplyr::if_else(.data$log10_total_b12_gm >= total_b12_threshold,
                     "above", "below"),
      levels = c("below", "above")
    ),
    uptake_b12_mean_high = factor(
      dplyr::if_else(.data$uptake_b12_mean >= uptake_b12_threshold,
                     "above", "below"),
      levels = c("below", "above")
    )
  ) |>
  dplyr::select(
    "strain_id", "gtdb_phylum", "gtdb_class", "gtdb_order",
    "log10_total_b12_gm", "uptake_b12_mean",
    "log10_total_b12_gm_high", "uptake_b12_mean_high",
    dplyr::all_of(ko_features)
  )
if (anyDuplicated(model_input$strain_id)) stop("Model input is not unique by strain_id.")

trait_summary <- read_csv(trait_summary_path) |>
  dplyr::select("strain_id", "total_b12_gm", "filtrate_proportion_gm") |>
  dplyr::distinct(.data$strain_id, .keep_all = TRUE)
tree <- ape::read.tree(tree_path)
model_input <- model_input |>
  dplyr::left_join(trait_summary, by = "strain_id") |>
  dplyr::filter(.data$strain_id %in% tree$tip.label) |>
  dplyr::mutate(
    uptake_b12_mean_high_total_b12_gt50 = factor(
      dplyr::if_else(.data$uptake_b12_mean >= 250, "above", "below"),
      levels = c("below", "above")
    ),
    filtrate_proportion_gm_high_total_b12_gt50 = factor(
      dplyr::if_else(.data$filtrate_proportion_gm >= 0.05, "above", "below"),
      levels = c("below", "above")
    )
  )

full_groups <- make_order_holdout_groups(model_input, "full")
filtered_model_input <- model_input |>
  dplyr::filter(
    is.finite(.data$total_b12_gm), .data$total_b12_gm > 50,
    is.finite(.data$uptake_b12_mean),
    is.finite(.data$filtrate_proportion_gm)
  )
filtered_groups <- make_order_holdout_groups(filtered_model_input, "total_b12_gt50")
all_groups <- dplyr::bind_rows(full_groups, filtered_groups)

task_specs <- list(
  list(target = "log10_total_b12_gm_high", data = model_input, groups = full_groups),
  list(target = "uptake_b12_mean_high", data = model_input, groups = full_groups),
  list(target = "uptake_b12_mean_high_total_b12_gt50",
       data = filtered_model_input, groups = filtered_groups),
  list(target = "filtrate_proportion_gm_high_total_b12_gt50",
       data = filtered_model_input, groups = filtered_groups)
)

set.seed(random_seed)
split_results <- list()
fold_assignment_records <- list()
for (task in task_specs) {
  k <- nrow(task$groups)
  fold_id <- make_stratified_folds(task$data[[task$target]], k)
  # The partition is exhaustive and sampled without replacement: every row has
  # exactly one test-fold ID and each fold is evaluated once.
  if (!identical(sort(unique(fold_id)), seq_len(k))) {
    stop("Random fold partition is not exhaustive for trait ", task$target)
  }
  fold_assignment_records[[task$target]] <- tibble::tibble(
    trait = task$target,
    strain_id = task$data$strain_id,
    random_fold = fold_id,
    n_matched_folds = k
  )
  for (fold in seq_len(k)) {
    split_results[[length(split_results) + 1L]] <- run_split(
      task$data, task$target, "random_order_matched_fold",
      sprintf("random_fold_%02d_of_%02d", fold, k),
      which(fold_id != fold), which(fold_id == fold), ko_features
    )
  }
  for (row in seq_len(nrow(task$groups))) {
    held_out <- task$groups$test_tips[[row]]
    split_results[[length(split_results) + 1L]] <- run_split(
      task$data, task$target, "gtdb_order_block_held_out",
      task$groups$split_id[[row]],
      which(!task$data$strain_id %in% held_out),
      which(task$data$strain_id %in% held_out),
      ko_features
    )
  }
}
split_results <- purrr::compact(split_results)
metrics <- purrr::map_dfr(split_results, "metrics")
predictions <- purrr::map_dfr(split_results, "predictions")
feature_importance <- purrr::map_dfr(split_results, "importance")
null_metrics <- purrr::map_dfr(split_results, "null_metrics")
fold_assignments <- dplyr::bind_rows(fold_assignment_records)

expected_split_counts <- dplyr::bind_rows(lapply(task_specs, function(task) {
  tibble::tibble(trait = task$target, expected = nrow(task$groups))
}))
observed_split_counts <- metrics |>
  dplyr::count(.data$trait, .data$split_type, name = "observed") |>
  dplyr::left_join(expected_split_counts, by = "trait")
if (any(observed_split_counts$observed != observed_split_counts$expected)) {
  stop("Random and order-block split counts are not matched.")
}
if (anyDuplicated(fold_assignments[c("trait", "strain_id")])) {
  stop("A strain was assigned more than once in a trait's random partition.")
}

random_importance <- tidyr::crossing(
  metrics |>
    dplyr::filter(.data$split_type == "random_order_matched_fold") |>
    dplyr::select("trait", "split_type", "split_id") |>
    dplyr::distinct(),
  tibble::tibble(ko = ko_features)
) |>
  dplyr::left_join(
    feature_importance |>
      dplyr::filter(.data$split_type == "random_order_matched_fold") |>
      dplyr::select("trait", "split_type", "split_id", "ko", "importance"),
    by = c("trait", "split_type", "split_id", "ko")
  ) |>
  dplyr::mutate(
    retained_in_split = !is.na(.data$importance),
    importance = dplyr::coalesce(.data$importance, 0)
  ) |>
  dplyr::group_by(.data$trait, .data$ko) |>
  dplyr::summarise(
    mean_permutation_importance = mean(.data$importance),
    sd_permutation_importance = stats::sd(.data$importance),
    n_splits = dplyr::n_distinct(.data$split_id),
    n_splits_retained = sum(.data$retained_in_split),
    .groups = "drop"
  )

null_per_split <- null_metrics |>
  dplyr::group_by(.data$trait, .data$split_type, .data$split_id) |>
  dplyr::summarise(
    null_mean_balanced_accuracy = mean(.data$null_balanced_accuracy),
    null_mean_accuracy = mean(.data$null_accuracy),
    .groups = "drop"
  )
paired_performance <- metrics |>
  dplyr::select(
    "trait", "split_type", "split_id",
    observed_balanced_accuracy = "balanced_accuracy",
    observed_accuracy = "accuracy"
  ) |>
  dplyr::left_join(null_per_split, by = c("trait", "split_type", "split_id")) |>
  dplyr::group_by(.data$trait, .data$split_type) |>
  dplyr::summarise(
    paired_t_test_p_balanced_accuracy = stats::t.test(
      .data$observed_balanced_accuracy, .data$null_mean_balanced_accuracy,
      paired = TRUE, alternative = "greater"
    )$p.value,
    paired_t_test_p_accuracy = stats::t.test(
      .data$observed_accuracy, .data$null_mean_accuracy,
      paired = TRUE, alternative = "greater"
    )$p.value,
    median_observed_accuracy = stats::median(.data$observed_accuracy),
    median_null_accuracy = stats::median(.data$null_mean_accuracy),
    .groups = "drop"
  )
null_summary <- metrics |>
  dplyr::group_by(.data$trait, .data$split_type) |>
  dplyr::summarise(
    n_splits = dplyr::n_distinct(.data$split_id),
    observed_mean_balanced_accuracy = mean(.data$balanced_accuracy),
    observed_mean_accuracy = mean(.data$accuracy),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    null_per_split |>
      dplyr::group_by(.data$trait, .data$split_type) |>
      dplyr::summarise(
        null_mean_balanced_accuracy = mean(.data$null_mean_balanced_accuracy),
        null_mean_accuracy = mean(.data$null_mean_accuracy),
        .groups = "drop"
      ),
    by = c("trait", "split_type")
  ) |>
  dplyr::left_join(paired_performance, by = c("trait", "split_type"))

kegg_ko <- utils::read.delim(
  kegg_ko_list_path, header = FALSE, col.names = c("ko", "definition"),
  quote = "", stringsAsFactors = FALSE
) |>
  dplyr::mutate(
    gene_name = dplyr::if_else(
      grepl(";", .data$definition, fixed = TRUE),
      trimws(sub(";.*$", "", .data$definition)), NA_character_
    ),
    gene_name = gsub(", *", "-", .data$gene_name),
    ko_gene_label = dplyr::if_else(
      is.na(.data$gene_name) | .data$gene_name == "", .data$ko,
      paste0(.data$ko, "_", .data$gene_name)
    )
  )

mixed_files <- list.files(
  raw_gwa_dir, pattern = "_pyseer_mixed_summary[.]csv$", full.names = TRUE
)
if (length(mixed_files) == 0L) stop("No pyseer mixed-model summaries found.")
gwa_comparison <- purrr::map_dfr(sort(mixed_files), combine_gwa_pair)

write_csv(model_input, model_input_path)
write_csv(all_groups |> dplyr::select(-"test_tips"), group_metadata_path)
write_csv(fold_assignments, file.path(processed_dir, "random_fold_assignments.csv"))
write_csv(metrics, metrics_path)
write_csv(predictions, predictions_path)
write_csv(feature_importance, importance_path)
write_csv(random_importance, random_importance_path)
write_csv(null_metrics, null_metrics_path)
write_csv(null_summary, null_summary_path)
write_csv(kegg_ko |> dplyr::select("ko", "gene_name", "ko_gene_label"), ko_mapping_path)
write_csv(gwa_comparison, gwa_comparison_path)
write_csv(
  tibble::tibble(
    parameter = c(
      "random_seed", "validation_design", "n_order_blocks_full",
      "n_order_blocks_total_b12_gt50", "minimum_order_test_fraction",
      "null_replicates", "rf_num_trees", "rf_num_threads",
      "min_train_feature_prevalence", "n_model_strains", "n_ko_features"
    ),
    value = as.character(c(
      random_seed,
      "single-stage; one non-overlapping random partition matched to order-block count",
      nrow(full_groups), nrow(filtered_groups), minimum_order_test_fraction,
      null_replicates, rf_num_trees, rf_num_threads,
      min_train_feature_prevalence, nrow(model_input), length(ko_features)
    ))
  ),
  parameter_path
)

message("Wrote classification model input: ", model_input_path)
message("Matched folds per full-data trait: ", nrow(full_groups))
message("Matched folds per B12-synthesizer trait: ", nrow(filtered_groups))
message("Wrote classification/GWA outputs to: ", processed_dir)

# Compact source-export stage adapted from the approved classification figure analysis.
processed_dir <- file.path("data/processed", "b12_trait_classification_analysis")

figure_source_dir <- file.path(
    "data/processed/figure_source_data", "b12_trait_classification_analysis"
)

metrics_path <- file.path(processed_dir, "classification_metrics.csv")

feature_importance_path <- file.path(processed_dir, "classification_feature_importance.csv")

random_fold_importance_summary_path <- file.path(
    processed_dir, "random_fold_feature_importance_summary.csv"
)

null_metrics_path <- file.path(processed_dir, "train_prevalence_null_metrics.csv")

clade_metadata_path <- file.path(processed_dir, "order_block_group_metadata.csv")

model_input_path <- file.path(processed_dir, "b12_trait_classification_model_input.csv")

ko_gene_name_mapping_path <- file.path(processed_dir, "ko_gene_name_mapping.csv")

gwa_comparison_path <- file.path(
    processed_dir, "gwa_pyseer_mixed_tree_comparison_all_traits.csv"
)

gtdbtk_tree_pruned_path <- file.path("data/processed/gtdbtk_tree_analysis", "gtdbtk_tree_231_genomes_bacillati_pseudomonadati_split_rooted_b12_complete_traits.tree")

gtdbtk_tip_metadata_path <- file.path("data/processed/gtdbtk_tree_analysis", "gtdbtk_tree_231_tip_metadata.csv")

compact_order_accuracy_source_path <- file.path(figure_source_dir, "01F_classification_accuracy_random_order_block_held_out_uncolored_points_source.csv")

all_strains_binary_tree_source_path <- file.path(figure_source_dir, "02C_gtdbtk_tree_binary_traits_tips_source.csv")

filtered_clade_binary_tree_source_path <- file.path(figure_source_dir, "02D_gtdbtk_total_b12_gt50_pruned_tree_binary_traits_tips_source.csv")

top_feature_compact_source_path <- file.path(figure_source_dir, "03B_top_10_mean_permutation_feature_importance_random_order_block_sd_source.csv")

rf_pyseer_combined_source_path <- file.path(figure_source_dir, "04C_rf_feature_importance_vs_pyseer_lrt_pvalue_by_phenotype_source.csv")

read_csv <- function(path) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

write_csv <- function(data, path) {
    utils::write.csv(data, path, row.names = FALSE, na = "")
}

format_accuracy_trait_label <- function(x) {
    dplyr::recode(x, log10_total_b12_gm_high = "total_b12_gm", uptake_b12_mean_high = "uptake_b12_mean", uptake_b12_mean_high_total_b12_gt50 = "uptake_b12_mean | total_b12_gm > 50", 
        filtrate_proportion_gm_high_total_b12_gt50 = "filtrate_proportion_gm | total_b12_gm > 50", .default = x)
}

accuracy_trait_levels <- c("total_b12_gm", "uptake_b12_mean", "uptake_b12_mean | total_b12_gm > 50", "filtrate_proportion_gm | total_b12_gm > 50")

accuracy_trait_labeller <- c(total_b12_gm = "B12 synthesis", uptake_b12_mean = "B12 uptake", `uptake_b12_mean | total_b12_gm > 50` = "B12 uptake\ntotal_b12_gm > 50", 
    `filtrate_proportion_gm | total_b12_gm > 50` = "Extracellular B12 fraction\ntotal_b12_gm > 50")

add_accuracy_facet_labels <- function(data) {
    dplyr::mutate(data, dataset_label = factor(dplyr::if_else(stringr::str_detect(as.character(.data$trait), "total_b12_gt50$"), 
        "B12 synthesizers", "All strains"), levels = c("All strains", "B12 synthesizers")), trait_facet_label = factor(dplyr::case_when(.data$trait == 
        "log10_total_b12_gm_high" ~ "B12 synthesis", stringr::str_detect(as.character(.data$trait), "uptake_b12_mean") ~ 
        "B12 uptake", stringr::str_detect(as.character(.data$trait), "filtrate_proportion_gm") ~ "EC-B12 fraction", TRUE ~ 
        as.character(.data$trait)), levels = c("B12 synthesis", "B12 uptake", "EC-B12 fraction")))
}

format_split_label <- function(x) {
    dplyr::recode(x, random_order_matched_fold = "Matched random folds", gtdb_order_block_held_out = "Order block held out", 
        .default = x)
}

clean_taxonomy <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "Unassigned"
    x
}

metrics <- read_csv(metrics_path)

feature_importance <- read_csv(feature_importance_path)

random_fold_importance_summary <- read_csv(random_fold_importance_summary_path)

null_metrics <- read_csv(null_metrics_path)

clade_metadata <- read_csv(clade_metadata_path)

order_group_tip_map <- purrr::map_dfr(seq_len(nrow(clade_metadata)), function(i) {
    tibble::tibble(tree_scope = clade_metadata$tree_scope[i], split_id = clade_metadata$split_id[i], held_out_group = dplyr::if_else(clade_metadata$split_id[i] == 
        "pooled_small_orders", "Pooled small orders", clade_metadata$held_out_orders[i]), tip_label = strsplit(clade_metadata$test_tip_labels[i], 
        ";", fixed = TRUE)[[1]])
})

model_input <- read_csv(model_input_path)

ko_gene_name_mapping <- read_csv(ko_gene_name_mapping_path)

gwa_comparison <- read_csv(gwa_comparison_path)

gtdbtk_tip_metadata <- read_csv(gtdbtk_tip_metadata_path) |>
    dplyr::select(-dplyr::any_of("filtrate_proportion_gm"))

gtdbtk_tree_pruned <- ape::read.tree(gtdbtk_tree_pruned_path)

gtdbtk_tree_pruned <- ape::drop.tip(gtdbtk_tree_pruned, setdiff(gtdbtk_tree_pruned$tip.label, model_input$strain_id))

n_random_split_repeats <- dplyr::n_distinct(metrics$split_id[metrics$split_type == "random_order_matched_fold"])

clade_comparison_metrics <- dplyr::filter(metrics, .data$split_type %in% c("random_order_matched_fold", "gtdb_order_block_held_out"))

clade_comparison_null_metrics <- dplyr::filter(null_metrics, .data$split_type %in% c("random_order_matched_fold", "gtdb_order_block_held_out"))

null_split_summary <- dplyr::summarise(dplyr::group_by(clade_comparison_null_metrics, .data$trait, .data$split_type, .data$split_id), 
    accuracy = mean(.data$null_accuracy, na.rm = TRUE), balanced_accuracy = mean(.data$null_balanced_accuracy, na.rm = TRUE), 
    .groups = "drop")

clade_split_labels <- dplyr::transmute(dplyr::bind_rows(tidyr::crossing(dplyr::filter(clade_metadata, .data$tree_scope == 
    "full"), trait = c("log10_total_b12_gm_high", "uptake_b12_mean_high")), tidyr::crossing(dplyr::filter(clade_metadata, 
    .data$tree_scope == "total_b12_gt50"), trait = c("uptake_b12_mean_high_total_b12_gt50", "filtrate_proportion_gm_high_total_b12_gt50"))), 
    trait = .data$trait, split_id = .data$split_id, split_point_label = paste0(.data$held_out_orders, " (", .data$n_test, 
        " tips)"))

accuracy_point_data <- dplyr::mutate(dplyr::left_join(dplyr::bind_rows(dplyr::mutate(clade_comparison_metrics, accuracy_trait_label = format_accuracy_trait_label(.data$trait), 
    accuracy_group = format_split_label(.data$split_type), accuracy_model_type = "Classifier"), dplyr::mutate(null_split_summary, 
    accuracy_trait_label = format_accuracy_trait_label(.data$trait), accuracy_group = format_split_label(.data$split_type), 
    accuracy_model_type = "Null")), clade_split_labels, by = c("trait", "split_id")), random_fold_label = stringr::str_replace(.data$split_id, 
    "^repeat_[0-9]+_", ""), split_point_label = dplyr::case_when(.data$split_type == "random_order_matched_fold" ~ .data$split_id,
    TRUE ~ dplyr::coalesce(.data$split_point_label, .data$split_id)), accuracy_group = factor(.data$accuracy_group, 
    levels = c("Matched random folds", "Order block held out")), accuracy_model_type = factor(.data$accuracy_model_type, levels = c("Classifier", 
    "Null")), accuracy_trait_label = factor(.data$accuracy_trait_label, levels = accuracy_trait_levels))

split_point_levels <- c(
    unique(metrics$split_id[metrics$split_type == "random_order_matched_fold"]),
    unique(clade_split_labels$split_point_label)
)

accuracy_point_data <- dplyr::mutate(accuracy_point_data, split_point_label = factor(.data$split_point_label, levels = split_point_levels))

accuracy_point_data <- add_accuracy_facet_labels(accuracy_point_data)

write_csv(accuracy_point_data, compact_order_accuracy_source_path)

get_tree_layout <- function(tree) {
    device_path <- tempfile(fileext = ".pdf")
    grDevices::pdf(device_path)
    on.exit({
        grDevices::dev.off()
        unlink(device_path)
    }, add = TRUE)
    graphics::plot(
        tree, direction = "downwards", show.tip.label = FALSE, plot = FALSE
    )
    get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
}

tree_layout <- get_tree_layout(gtdbtk_tree_pruned)

n_tree_tips <- length(gtdbtk_tree_pruned$tip.label)

clade_tree_tips <- dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::left_join(tibble::tibble(tip_label = gtdbtk_tree_pruned$tip.label, 
    node = seq_len(n_tree_tips), x = tree_layout$xx[seq_len(n_tree_tips)], y = tree_layout$yy[seq_len(n_tree_tips)]), gtdbtk_tip_metadata, 
    by = "tip_label"), dplyr::select(dplyr::filter(order_group_tip_map, .data$tree_scope == "full"), "tip_label", "split_id", 
    "held_out_group"), by = "tip_label"), dplyr::distinct(dplyr::select(model_input, "strain_id", "filtrate_proportion_gm")), 
    by = c(tip_label = "strain_id")), order_taxonomy = clean_taxonomy(.data$gtdb_order), total_b12_gm = as.numeric(.data$total_b12_gm), 
    uptake_b12_mean = as.numeric(.data$uptake_b12_mean), filtrate_proportion_gm = as.numeric(.data$filtrate_proportion_gm), 
    total_b12_binary = factor(dplyr::if_else(.data$total_b12_gm >= 50, "Above", "Below"), levels = c("Below", "Above")), 
    uptake_b12_binary = factor(dplyr::if_else(.data$uptake_b12_mean >= 250, "Above", "Below"), levels = c("Below", "Above")), 
    filtrate_proportion_binary = factor(dplyr::if_else(.data$filtrate_proportion_gm >= 0.05, "Above", "Below"), levels = c("Below", 
        "Above")))

write_csv(clade_tree_tips, all_strains_binary_tree_source_path)

filtered_tip_ids_02b <- dplyr::pull(dplyr::filter(model_input, is.finite(.data$total_b12_gm), .data$total_b12_gm > 50, is.finite(.data$uptake_b12_mean), 
    is.finite(.data$filtrate_proportion_gm)), .data$strain_id)

filtered_tree_02b <- ape::drop.tip(gtdbtk_tree_pruned, setdiff(gtdbtk_tree_pruned$tip.label, filtered_tip_ids_02b))

tree_layout_02b <- get_tree_layout(filtered_tree_02b)

n_tips_02b <- length(filtered_tree_02b$tip.label)

clade_tree_tips_02b <- dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::left_join(tibble::tibble(tip_label = filtered_tree_02b$tip.label, 
    node = seq_len(n_tips_02b), x = tree_layout_02b$xx[seq_len(n_tips_02b)], y = tree_layout_02b$yy[seq_len(n_tips_02b)]), 
    gtdbtk_tip_metadata, by = "tip_label"), dplyr::select(dplyr::filter(order_group_tip_map, .data$tree_scope == "total_b12_gt50"), 
    "tip_label", "split_id", "held_out_group"), by = "tip_label"), dplyr::distinct(dplyr::select(model_input, "strain_id", 
    "filtrate_proportion_gm")), by = c(tip_label = "strain_id")), order_taxonomy = clean_taxonomy(.data$gtdb_order), total_b12_gm = as.numeric(.data$total_b12_gm), 
    uptake_b12_mean = as.numeric(.data$uptake_b12_mean), filtrate_proportion_gm = as.numeric(.data$filtrate_proportion_gm), 
    total_b12_binary = factor(dplyr::if_else(.data$total_b12_gm >= 50, "Above", "Below"), levels = c("Below", "Above")), 
    uptake_b12_binary = factor(dplyr::if_else(.data$uptake_b12_mean >= 250, "Above", "Below"), levels = c("Below", "Above")), 
    filtrate_proportion_binary = factor(dplyr::if_else(.data$filtrate_proportion_gm >= 0.05, "Above", "Below"), levels = c("Below", 
        "Above")))

write_csv(clade_tree_tips_02b, filtered_clade_binary_tree_source_path)

b12_anaerobic_kos <- c("K02302", "K02303", "K13542", "K02304", "K24866", "K02190", "K03795", "K22011", "K03394", "K05934", 
    "K13541", "K21479", "K05936", "K02189", "K02188", "K05895", "K02191", "K03399", "K00595", "K06042", "K02224")

b12_aerobic_kos <- c("K02303", "K13542", "K03394", "K13540", "K02229", "K05934", "K13541", "K05936", "K02228", "K05895", 
    "K00595", "K06042", "K02224", "K02230", "K09882", "K09883")

b12_late_salvage_kos <- c("K00798", "K19221", "K02232", "K02225", "K02227", "K02231", "K00768", "K02226", "K22316", "K02233")

b12_transporter_kos <- c("K16092", "K06858", "K25034", "K06073", "K25027", "K06074", "K25028", "K02471", "K16927", "K01552", 
    "K16785")

b12_pathway_membership_colors <- c(Anaerobic = "#7B3294", Aerobic = "#2F6F73", `Late/salvage` = "#B65F2A", `Anaerobic + Aerobic` = "#4D7FB8", 
    `Anaerobic + Late/salvage` = "#9A5AA6", `Aerobic + Late/salvage` = "#6F9F65", `Anaerobic + Aerobic + Late/salvage` = "#5F5F5F", 
    `Not in KEGG\nB12 module` = "grey40")

b12_transporter_colors <- c(`B12 transporter` = "#2F6F73", `Not B12 transporter` = "grey40")

add_kegg_gene_id <- function(data) {
    dplyr::mutate(dplyr::left_join(data, dplyr::distinct(dplyr::select(ko_gene_name_mapping, "ko", "gene_name")), by = "ko"), 
        gene_id_label = dplyr::coalesce(.data$gene_name, .data$ko))
}

add_b12_pathway_membership <- function(data) {
    dplyr::select(dplyr::mutate(data, in_anaerobic = .data$ko %in% b12_anaerobic_kos, in_aerobic = .data$ko %in% b12_aerobic_kos, 
        in_late_salvage = .data$ko %in% b12_late_salvage_kos, b12_pathway_membership = dplyr::case_when(.data$in_anaerobic & 
            .data$in_aerobic & .data$in_late_salvage ~ "Anaerobic + Aerobic + Late/salvage", .data$in_anaerobic & .data$in_aerobic ~ 
            "Anaerobic + Aerobic", .data$in_anaerobic & .data$in_late_salvage ~ "Anaerobic + Late/salvage", .data$in_aerobic & 
            .data$in_late_salvage ~ "Aerobic + Late/salvage", .data$in_anaerobic ~ "Anaerobic", .data$in_aerobic ~ "Aerobic", 
            .data$in_late_salvage ~ "Late/salvage", TRUE ~ "Not in KEGG\nB12 module"), b12_pathway_membership = factor(.data$b12_pathway_membership, 
            levels = names(b12_pathway_membership_colors))), -"in_anaerobic", -"in_aerobic", -"in_late_salvage")
}

make_rf_pyseer_data <- function(rf_trait, pyseer_trait) {
    add_kegg_gene_id(dplyr::mutate(dplyr::left_join(dplyr::transmute(dplyr::filter(gwa_comparison, .data$trait == pyseer_trait, 
        is.finite(.data$p_adj_BH_without_tree), .data$p_adj_BH_without_tree > 0), ko = .data$variant, pyseer_lrt_pvalue = .data$lrt.pvalue_without_tree, 
        pyseer_p_adj_BH = .data$p_adj_BH_without_tree, pyseer_p_adj_bonferroni = .data$p_adj_bonferroni_without_tree, pyseer_neg_log10_BH_adjusted_lrt_pvalue = -log10(.data$p_adj_BH_without_tree), 
        pyseer_rejected_BH_005 = .data$rejected_BH_005_without_tree, pyseer_rejected_bonferroni_005 = .data$rejected_bonferroni_005_without_tree), 
        dplyr::select(dplyr::filter(random_fold_importance_summary, .data$trait == rf_trait), "ko", "mean_permutation_importance", 
            "sd_permutation_importance", "n_splits", "n_splits_retained"), by = "ko"), mean_permutation_importance = dplyr::coalesce(.data$mean_permutation_importance, 
        0), sd_permutation_importance = dplyr::coalesce(.data$sd_permutation_importance, 0), n_splits = dplyr::coalesce(.data$n_splits, 
        0L), n_splits_retained = dplyr::coalesce(.data$n_splits_retained, 0L)))
}

rf_pyseer_total_b12_data <- add_b12_pathway_membership(make_rf_pyseer_data("log10_total_b12_gm_high", "log10_total_b12_gm"))

rf_pyseer_uptake_b12_data <- dplyr::mutate(make_rf_pyseer_data("uptake_b12_mean_high", "uptake_b12_mean"), b12_transporter = factor(dplyr::if_else(.data$ko %in% 
    b12_transporter_kos, "B12 transporter", "Not B12 transporter"), levels = names(b12_transporter_colors)))

rf_pyseer_combined_data <- dplyr::arrange(dplyr::mutate(dplyr::bind_rows(dplyr::mutate(rf_pyseer_total_b12_data, phenotype = "total_b12_gm", 
    feature_group = dplyr::case_when(as.character(.data$b12_pathway_membership) %in% c("Anaerobic", "Aerobic", "Anaerobic + Aerobic") ~ 
        "Corrin ring", TRUE ~ as.character(.data$b12_pathway_membership))), dplyr::mutate(rf_pyseer_uptake_b12_data, phenotype = "b12_uptake_mean", 
    feature_group = as.character(.data$b12_transporter))), phenotype = factor(.data$phenotype, levels = c("total_b12_gm", 
    "b12_uptake_mean")), feature_group = factor(.data$feature_group, levels = c("Not in KEGG\nB12 module", "Not B12 transporter", 
    "Corrin ring", "Late/salvage", "B12 transporter"))), .data$feature_group)

write_csv(rf_pyseer_combined_data, rf_pyseer_combined_source_path)

importance_split_index <- dplyr::distinct(dplyr::select(metrics, "trait", "split_type", "split_id"))

importance_ko_index <- tibble::tibble(ko = names(model_input)[stringr::str_detect(names(model_input), "^K[0-9]{5}$")])

feature_importance_complete <- dplyr::mutate(dplyr::left_join(tidyr::crossing(importance_split_index, importance_ko_index), 
    dplyr::select(feature_importance, "trait", "split_type", "split_id", "ko", "importance"), by = c("trait", "split_type", 
        "split_id", "ko")), selected_in_split = !is.na(.data$importance), importance = dplyr::coalesce(.data$importance, 
    0))

top_feature_source <- dplyr::mutate(dplyr::arrange(dplyr::mutate(dplyr::left_join(dplyr::ungroup(dplyr::mutate(dplyr::slice_max(dplyr::group_by(dplyr::summarise(dplyr::group_by(feature_importance_complete, 
    .data$trait, .data$split_type, .data$ko), mean_permutation_importance = mean(.data$importance, na.rm = TRUE), sd_permutation_importance = stats::sd(.data$importance, 
    na.rm = TRUE), n_splits = dplyr::n_distinct(.data$split_id), n_selected_splits = sum(.data$selected_in_split), selection_frequency = mean(.data$selected_in_split), 
    .groups = "drop"), .data$trait, .data$split_type), .data$mean_permutation_importance, n = 10, with_ties = FALSE), rank = dplyr::row_number())), 
    dplyr::select(ko_gene_name_mapping, "ko", "ko_gene_label"), by = "ko"), ko_gene_label = dplyr::coalesce(.data$ko_gene_label, 
    .data$ko), trait_label = factor(format_accuracy_trait_label(.data$trait), levels = accuracy_trait_levels), split_label = factor(format_split_label(.data$split_type), 
    levels = c("Matched random folds", "Order block held out")), panel_label = paste(unname(accuracy_trait_labeller[as.character(.data$trait_label)]), 
    as.character(.data$split_label), sep = "\n"), feature_label = paste(.data$trait, .data$split_type, .data$ko_gene_label, 
    sep = "__")), .data$trait_label, .data$split_label, .data$mean_permutation_importance), feature_label = factor(.data$feature_label, 
    levels = unique(.data$feature_label)), panel_label = factor(.data$panel_label, levels = unique(.data$panel_label)))

top_feature_compact_task_levels <- unname(accuracy_trait_labeller[accuracy_trait_levels])

top_feature_compact_split_levels <- c("Matched random folds", "Order block held out")

top_feature_compact_panel_levels <- unlist(lapply(top_feature_compact_split_levels, function(split_label) {
    paste(top_feature_compact_task_levels, split_label, sep = "\n")
}), use.names = FALSE)

top_feature_compact_source <- dplyr::mutate(dplyr::arrange(dplyr::mutate(dplyr::filter(top_feature_source, .data$split_type %in% 
    c("random_order_matched_fold", "gtdb_order_block_held_out")), prediction_task = factor(unname(accuracy_trait_labeller[as.character(.data$trait_label)]), 
    levels = top_feature_compact_task_levels), held_out_test_set = factor(as.character(.data$split_label), levels = top_feature_compact_split_levels), 
    panel_label = factor(paste(.data$prediction_task, .data$held_out_test_set, sep = "\n"), levels = top_feature_compact_panel_levels)), 
    .data$held_out_test_set, .data$prediction_task, .data$mean_permutation_importance), feature_label = factor(paste(.data$trait, 
    .data$split_type, .data$ko_gene_label, sep = "__"), levels = unique(paste(.data$trait, .data$split_type, .data$ko_gene_label, 
    sep = "__"))))

write_csv(top_feature_compact_source, top_feature_compact_source_path)
