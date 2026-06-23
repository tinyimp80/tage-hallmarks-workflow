#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% args

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "pheatmap", "reticulate"
)

bioc_packages <- c("DESeq2", "GSVA", "limma", "Biobase")

install_if_missing <- function(packages, installer) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(character())
  if (check_only) return(missing)
  installer(missing)
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

missing_cran <- install_if_missing(cran_packages, function(pkgs) {
  install.packages(pkgs, repos = "https://cloud.r-project.org")
})

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  if (check_only) {
    missing_biocmanager <- "BiocManager"
  } else {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
    missing_biocmanager <- character()
  }
} else {
  missing_biocmanager <- character()
}

missing_bioc <- install_if_missing(bioc_packages, function(pkgs) {
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)
})

missing_tage <- if (!requireNamespace("tAge", quietly = TRUE)) "tAge" else character()

all_missing <- unique(c(missing_cran, missing_biocmanager, missing_bioc, missing_tage))

if (length(all_missing) > 0) {
  message("Missing packages: ", paste(all_missing, collapse = ", "))
  if ("tAge" %in% all_missing) {
    message("The tAge R package and trained models are not bundled with this repository.")
    message("Model source: https://zenodo.org/records/18763485")
  }
  quit(status = 1)
}

message("All required R packages are available.")
