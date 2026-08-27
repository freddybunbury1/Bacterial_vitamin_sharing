#!/usr/bin/env Rscript

# Download the KEGG Orthology label table required for figure annotations.
# KEGG REST access is restricted to academic use by academic users:
# https://www.kegg.jp/kegg/rest/

output_dir <- "data/raw/annotations"
retrieval_date <- format(Sys.Date(), "%Y-%m-%d")
ko_list_path <- file.path(output_dir, "kegg_ko_list.tsv")
ko_info_path <- file.path(output_dir, "kegg_ko_info.txt")
retrieval_metadata_path <- file.path(output_dir, "kegg_ko_retrieval_metadata.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

utils::download.file(
  "https://rest.kegg.jp/list/ko",
  ko_list_path,
  mode = "wb",
  quiet = FALSE
)
utils::download.file(
  "https://rest.kegg.jp/info/ko",
  ko_info_path,
  mode = "wb",
  quiet = FALSE
)

ko_list <- utils::read.delim(
  ko_list_path,
  header = FALSE,
  col.names = c("ko", "definition"),
  quote = "",
  stringsAsFactors = FALSE
)

if (nrow(ko_list) == 0L || anyDuplicated(ko_list$ko) ||
    !all(grepl("^K[0-9]{5}$", ko_list$ko))) {
  stop("Downloaded KEGG KO list failed validation.")
}

utils::write.csv(
  data.frame(
    retrieval_date = retrieval_date,
    list_url = "https://rest.kegg.jp/list/ko",
    info_url = "https://rest.kegg.jp/info/ko",
    academic_use_notice = "https://www.kegg.jp/kegg/rest/",
    stringsAsFactors = FALSE
  ),
  retrieval_metadata_path,
  row.names = FALSE
)

message("Wrote ", nrow(ko_list), " KEGG KO entries: ", ko_list_path)
message("Wrote KEGG KO release metadata: ", ko_info_path)
message("Wrote retrieval provenance: ", retrieval_metadata_path)
