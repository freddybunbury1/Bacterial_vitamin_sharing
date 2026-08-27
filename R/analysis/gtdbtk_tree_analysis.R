#!/usr/bin/env Rscript

# Public GTDB-Tk taxonomy and phylogeny analysis. The GTDB-Tk summary is joined
# to strain IDs, the supplied 231-strain tree is rooted at the best-supported
# Bacillati/Pseudomonadati split, and a trait-complete tree is exported for the
# classification and compiled-figure stages.

summary_path <- file.path("data/raw/gtdbtk", "gtdbtk.bac120.summary.tsv")
tree_path <- file.path("data/raw/gtdbtk", "gtdbtk_tree_231_genomes.tree")
genome_map_path <- file.path(
  "data/raw/metadata_strains",
  "microtrait_strain_id_to_genome_filename_20260520.csv"
)
trait_path <- file.path(
  "data/processed/b12_trait_survey",
  "b12_trait_survey_strain_summary.csv"
)

taxonomy_dir <- file.path("data/intermediate", "gtdbtk")
processed_dir <- file.path("data/processed", "gtdbtk_tree_analysis")
taxonomy_path <- file.path(taxonomy_dir, "gtdbtk_bac120_taxonomy_with_strain_id.csv")
split_tree_path <- file.path(
  processed_dir,
  "gtdbtk_tree_231_genomes_bacillati_pseudomonadati_split_rooted.tree"
)
complete_tree_path <- file.path(
  processed_dir,
  "gtdbtk_tree_231_genomes_bacillati_pseudomonadati_split_rooted_b12_complete_traits.tree"
)
tip_metadata_path <- file.path(processed_dir, "gtdbtk_tree_231_tip_metadata.csv")
rooting_summary_path <- file.path(processed_dir, "gtdbtk_tree_231_rooting_summary.csv")

required_packages <- c("ape", "dplyr", "phytools", "stringr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Required packages are missing: ", paste(missing_packages, collapse = ", "))
}

required_files <- c(summary_path, tree_path, genome_map_path, trait_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Required input files are missing:\n", paste(missing_files, collapse = "\n"))
}

dir.create(taxonomy_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  utils::read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A")
  )
}

strip_fasta_extension <- function(x) {
  sub("\\.(fna|fa|fasta)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
}

clean_taxon <- function(x, prefix) {
  value <- stringr::str_remove(x, paste0("^", prefix, "__"))
  value[value == ""] <- NA_character_
  value
}

split_gtdb_taxonomy <- function(classification) {
  ranks <- stringr::str_split_fixed(classification, ";", 7)
  tibble::tibble(
    gtdb_domain = clean_taxon(ranks[, 1], "d"),
    gtdb_phylum = clean_taxon(ranks[, 2], "p"),
    gtdb_class = clean_taxon(ranks[, 3], "c"),
    gtdb_order = clean_taxon(ranks[, 4], "o"),
    gtdb_family = clean_taxon(ranks[, 5], "f"),
    gtdb_genus = clean_taxon(ranks[, 6], "g"),
    gtdb_species = clean_taxon(ranks[, 7], "s")
  )
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0L) NA else x[[1]]
}

get_component_tips <- function(tree, start_node, blocked_node) {
  n_tips <- length(tree$tip.label)
  all_nodes <- seq_len(n_tips + tree$Nnode)
  adjacency <- lapply(all_nodes, function(x) integer())
  for (i in seq_len(nrow(tree$edge))) {
    node_a <- tree$edge[i, 1]
    node_b <- tree$edge[i, 2]
    adjacency[[node_a]] <- c(adjacency[[node_a]], node_b)
    adjacency[[node_b]] <- c(adjacency[[node_b]], node_a)
  }
  seen <- rep(FALSE, length(all_nodes))
  seen[blocked_node] <- TRUE
  queue <- start_node
  while (length(queue) > 0L) {
    node <- queue[[1]]
    queue <- queue[-1]
    if (seen[node]) next
    seen[node] <- TRUE
    queue <- c(queue, adjacency[[node]][!seen[adjacency[[node]]]])
  }
  tree$tip.label[which(seen[seq_len(n_tips)])]
}

score_split_edges <- function(tree, metadata) {
  known <- dplyr::filter(metadata, !is.na(.data$supergroup))
  if (nrow(known) == 0L) stop("No GTDB phyla could be assigned to rooting supergroups.")

  dplyr::bind_rows(lapply(seq_len(nrow(tree$edge)), function(edge_index) {
    node_a <- tree$edge[edge_index, 1]
    node_b <- tree$edge[edge_index, 2]
    side_a_tips <- get_component_tips(tree, node_a, node_b)
    side_a <- known$tip_label %in% side_a_tips
    bacillati_a <- sum(side_a & known$supergroup == "Bacillati")
    pseudomonadati_a <- sum(side_a & known$supergroup == "Pseudomonadati")
    bacillati_b <- sum(!side_a & known$supergroup == "Bacillati")
    pseudomonadati_b <- sum(!side_a & known$supergroup == "Pseudomonadati")
    best_agreement <- max(
      bacillati_a + pseudomonadati_b,
      pseudomonadati_a + bacillati_b
    ) / nrow(known)
    tibble::tibble(
      edge_index = edge_index,
      node_a = node_a,
      node_b = node_b,
      edge_length = tree$edge.length[edge_index],
      split_agreement = best_agreement,
      smaller_known_side = min(sum(side_a), sum(!side_a))
    )
  })) |>
    dplyr::arrange(
      dplyr::desc(.data$split_agreement),
      dplyr::desc(.data$smaller_known_side),
      .data$edge_index
    )
}

gtdb_summary <- utils::read.delim(
  summary_path, check.names = FALSE, stringsAsFactors = FALSE,
  na.strings = c("", "NA", "N/A")
)
genome_map <- read_csv(genome_map_path)
traits <- read_csv(trait_path)

required_summary_columns <- c("user_genome", "classification")
if (!all(required_summary_columns %in% names(gtdb_summary))) {
  stop("GTDB-Tk summary lacks required columns: ", paste(
    setdiff(required_summary_columns, names(gtdb_summary)), collapse = ", "
  ))
}
if (!all(c("strain_id", "fna_filename", "genome_source") %in% names(genome_map))) {
  stop("Genome mapping lacks strain_id, fna_filename, or genome_source.")
}
if (anyDuplicated(genome_map$fna_filename[!is.na(genome_map$fna_filename)])) {
  stop("Genome mapping contains duplicated non-missing fna_filename values.")
}

gtdb_summary <- gtdb_summary |>
  dplyr::bind_cols(split_gtdb_taxonomy(gtdb_summary$classification))
genome_map <- genome_map |>
  dplyr::filter(!is.na(.data$fna_filename), .data$fna_filename != "") |>
  dplyr::mutate(user_genome = strip_fasta_extension(.data$fna_filename))

taxonomy_all <- gtdb_summary |>
  dplyr::left_join(genome_map, by = "user_genome") |>
  dplyr::mutate(
    priority_rerun = stringr::str_detect(
      stringr::str_to_lower(dplyr::coalesce(.data$fna_filename, "")), "rerun"
    ),
    priority_source = dplyr::recode(
      .data$genome_source,
      plasmidsaurus = 1L, zymo = 2L, other = 3L, ncbi = 4L,
      .default = 5L
    )
  )

unmatched_user_genomes <- taxonomy_all$user_genome[is.na(taxonomy_all$strain_id)]
taxonomy <- taxonomy_all |>
  dplyr::filter(!is.na(.data$strain_id)) |>
  dplyr::arrange(
    .data$strain_id,
    dplyr::desc(.data$priority_rerun),
    .data$priority_source,
    .data$fna_filename,
    .data$user_genome
  ) |>
  dplyr::group_by(.data$strain_id) |>
  dplyr::slice(1L) |>
  dplyr::ungroup() |>
  dplyr::select(-"priority_rerun", -"priority_source")

if (anyDuplicated(taxonomy$strain_id)) stop("GTDB taxonomy is not unique by strain_id.")
utils::write.csv(taxonomy, taxonomy_path, row.names = FALSE, na = "")

tree <- ape::read.tree(tree_path)
tree$tip.label[tree$tip.label == "KC61126"] <- "KC611216"
if (anyDuplicated(tree$tip.label)) stop("The 231-genome tree has duplicated tip labels.")

trait_columns <- c(
  "strain_id", "total_b12_gm", "total_b12_per_od_gm", "filtrate_b12_gm",
  "filtrate_proportion_gm", "uptake_b12_mean", "b12_trait_group2", "n_samples"
)
missing_trait_columns <- setdiff(trait_columns, names(traits))
if (length(missing_trait_columns) > 0L) {
  stop("Trait summary lacks required columns: ", paste(missing_trait_columns, collapse = ", "))
}

tip_metadata <- tibble::tibble(tip_label = tree$tip.label) |>
  dplyr::left_join(
    taxonomy |>
      dplyr::select(
        "strain_id", "user_genome", "fna_filename", "genome_source",
        dplyr::starts_with("gtdb_")
      ),
    by = c("tip_label" = "strain_id")
  ) |>
  dplyr::left_join(
    traits |> dplyr::select(dplyr::all_of(trait_columns)),
    by = c("tip_label" = "strain_id")
  ) |>
  dplyr::mutate(
    taxonomy_source = dplyr::if_else(
      !is.na(.data$gtdb_phylum), "GTDB-Tk bac120 summary", NA_character_
    ),
    supergroup = dplyr::case_when(
      .data$gtdb_phylum %in% c("Actinomycetota", "Bacillota") ~ "Bacillati",
      .data$gtdb_phylum %in% c("Bacteroidota", "Pseudomonadota") ~ "Pseudomonadati",
      TRUE ~ NA_character_
    )
  )

edge_scores <- score_split_edges(tree, tip_metadata)
selected_edge <- edge_scores[1, , drop = FALSE]
root_position <- selected_edge$edge_length[[1]] / 2
split_tree <- phytools::reroot(
  tree,
  node.number = selected_edge$node_b[[1]],
  position = root_position
)

complete_tips <- tip_metadata |>
  dplyr::filter(
    is.finite(.data$total_b12_gm),
    is.finite(.data$uptake_b12_mean)
  ) |>
  dplyr::pull(.data$tip_label)
if (length(complete_tips) < 2L) stop("Fewer than two tree tips have complete B12 traits.")
complete_tree <- ape::keep.tip(split_tree, complete_tips)

ape::write.tree(split_tree, file = split_tree_path)
ape::write.tree(complete_tree, file = complete_tree_path)
utils::write.csv(tip_metadata, tip_metadata_path, row.names = FALSE, na = "")
utils::write.csv(
  data.frame(
    rooting_method = "midpoint of edge maximizing Bacillati/Pseudomonadati agreement",
    edge_index = selected_edge$edge_index,
    split_agreement = selected_edge$split_agreement,
    n_tree_tips = length(tree$tip.label),
    n_complete_trait_tips = length(complete_tree$tip.label),
    stringsAsFactors = FALSE
  ),
  rooting_summary_path, row.names = FALSE, na = ""
)

message("Wrote GTDB taxonomy: ", taxonomy_path)
message("Wrote split-rooted tree: ", split_tree_path)
message("Wrote complete-trait tree: ", complete_tree_path)
message("Wrote tree tip metadata: ", tip_metadata_path)
message("Tree tips: ", length(tree$tip.label), "; complete-trait tips: ", length(complete_tree$tip.label))
message("Unmatched GTDB user genomes retained only in the raw summary: ", length(unmatched_user_genomes))
