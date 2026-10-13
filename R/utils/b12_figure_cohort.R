# Shared eligibility for genome-linked trait figures and classifiers.
# The study-genome tree contains the genomes retained after genome QC.
canonical_b12_strain_id <- function(x) {
  x <- as.character(x)
  x[x == "KC61126" & !is.na(x)] <- "KC611216"
  x
}

b12_qc_genome_ids <- function(
    tree_path = "data/raw/gtdbtk/gtdbtk.bac120.study_genomes.tree") {
  unique(canonical_b12_strain_id(ape::read.tree(tree_path)$tip.label))
}

b12_figure_cohort <- function(data, genome_ids = b12_qc_genome_ids()) {
  required <- c("strain_id", "total_b12_gm", "uptake_b12_mean")
  if (!all(required %in% names(data))) stop("Missing phenotype cohort columns.")
  keep <- canonical_b12_strain_id(data$strain_id) %in% genome_ids &
    is.finite(data$total_b12_gm) & is.finite(data$uptake_b12_mean)
  data[keep, , drop = FALSE]
}
