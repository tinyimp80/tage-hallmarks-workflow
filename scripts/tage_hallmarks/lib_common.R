options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
})

project_dir <- normalizePath(Sys.getenv("TAGE_PROJECT_DIR", unset = getwd()), mustWork = TRUE)
zenodo_record_url <- "https://zenodo.org/records/18763485"

default_hallmark_gmt <- file.path(
  project_dir,
  "results", "aging_chemical_hallmarks", "hallmarks_aging", "results",
  "open_genes_go_pairwise_termsim", "open_genes_hallmark_GO_input_sets_ensembl.gmt"
)

hallmark_order <- c(
  "Genomic Instability", "Telomere Attrition", "Epigenetic Alterations",
  "Loss of Proteostasis", "Disabled Macroautophagy", "Deregulated Nutrient Sensing",
  "Mitochondrial Dysfunction", "Cellular Senescence", "Stem Cell Exhaustion",
  "Altered Intercellular Communication", "Chronic Inflammation"
)
set_mode_order <- c("aged_up", "aged_down")

usage <- function() {
  cat("
One-line tAge + Open Genes aging Hallmark ssGSEA workflow

Required:
  --metadata FILE
  --sample-id-col COL
  --group-col COL
  --control-group LABEL
  --out-dir DIR
  and exactly one of:
    --counts FILE
    --rsem-dir DIR

Core options:
  --species human|mouse|rat|monkey                 default: human
  --gene-mapping-type Ensembl|Gene.Symbol          default: Ensembl
  --hallmark-gmt FILE                              default: existing Open Genes aged-direction GMT

tAge model options:
  --model-family EN|BR|all                         default: EN
  --target mortality|chronoage|normalizedage|all   default: mortality
  --model-species Multispecies|Mouse|Rodents|auto  default: Multispecies
  --model-tissue Multitissue|Liver|Brain|Kidney|Skeletal_muscle|auto default: Multitissue
  --preprocess scaled_diff|yugene_diff|all         default: scaled_diff

Optional:
  --tissue-col COL
  --tissue-value VALUE
  --low-count-min-total N                          default: 10
  --min-gene-set-size N                            default: 5
  --max-gene-set-size N                            default: 1000
  --config-only                                    validate inputs/model selection only
  --help
")
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected positional argument: ", key)
    key <- sub("^--", "", key)
    if (key %in% c("help", "config-only")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args)) stop("Missing value for --", key)
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

arg_get <- function(opts, key, default = NULL) {
  if (!is.null(opts[[key]])) opts[[key]] else default
}

stop_if_missing <- function(opts, keys) {
  missing <- keys[vapply(keys, function(k) is.null(opts[[k]]) || identical(opts[[k]], ""), logical(1))]
  if (length(missing) > 0) stop("Missing required option(s): --", paste(missing, collapse = ", --"))
}

msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " | ", paste0(..., collapse = ""))
}

safe_name <- function(x) {
  x |>
    str_replace_all("[^A-Za-z0-9._-]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "fg", flag = "#"))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 4, format = "f")))
}

stars_from_p <- function(pvalue) {
  dplyr::case_when(
    is.na(pvalue) ~ "",
    pvalue < 0.001 ~ "***",
    pvalue < 0.01 ~ "**",
    pvalue < 0.05 ~ "*",
    pvalue < 0.1 ~ ".",
    TRUE ~ ""
  )
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sets <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    unique(parts[-c(1, 2)])
  })
  names(sets) <- vapply(strsplit(lines, "\t", fixed = TRUE), `[[`, character(1), 1)
  sets
}

parse_set_metadata <- function(set_name) {
  mode <- dplyr::case_when(
    str_ends(set_name, "__aged_up") ~ "aged_up",
    str_ends(set_name, "__aged_down") ~ "aged_down",
    str_ends(set_name, "__whole") ~ "whole",
    TRUE ~ NA_character_
  )
  hallmark <- set_name |>
    str_remove("__(aged_up|aged_down|whole)$") |>
    str_replace_all("_", " ")
  tibble(input_set_name = set_name, aging_hallmark = hallmark, set_mode = mode)
}

read_counts_csv <- function(path) {
  if (!file.exists(path)) stop("Counts file not found: ", path)
  dt <- fread(path)
  if (ncol(dt) < 2) stop("Counts file must contain gene column and sample columns.")
  genes <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1, with = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (anyNA(mat)) stop("Counts matrix contains NA or non-numeric values.")
  rownames(mat) <- sub("\\.[0-9]+$", "", genes)
  if (any(duplicated(rownames(mat)))) mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  mat
}

read_rsem_counts <- function(rsem_dir, metadata, sample_id_col) {
  if (!dir.exists(rsem_dir)) stop("RSEM directory not found: ", rsem_dir)
  files <- sort(list.files(rsem_dir, pattern = "[.]genes[.]results$", full.names = TRUE))
  if (length(files) == 0) stop("No RSEM *.genes.results files found in ", rsem_dir)
  raw_samples <- sub("[.]genes[.]results$", "", basename(files))
  meta_samples <- as.character(metadata[[sample_id_col]])
  if (anyDuplicated(raw_samples)) stop("Duplicated RSEM sample names after suffix stripping.")
  used <- meta_samples[meta_samples %in% raw_samples]
  if (length(used) == 0) stop("No overlapping RSEM sample IDs and metadata sample IDs.")
  metadata_only <- setdiff(meta_samples, raw_samples)
  if (length(metadata_only) > 0) stop("Metadata samples without RSEM files: ", paste(metadata_only, collapse = ";"))
  selected_files <- files[match(used, raw_samples)]
  read_one <- function(path) {
    x <- fread(path)
    if (!all(c("gene_id", "expected_count") %in% names(x))) stop("RSEM file missing gene_id/expected_count columns: ", path)
    ids <- sub("\\.[0-9]+$", "", as.character(x$gene_id))
    counts <- rowsum(as.numeric(x$expected_count), group = ids, reorder = FALSE)
    setNames(as.numeric(counts[, 1]), rownames(counts))
  }
  vectors <- lapply(selected_files, read_one)
  reference_genes <- names(vectors[[1]])
  if (!all(vapply(vectors, function(x) identical(names(x), reference_genes), logical(1)))) stop("Gene identifiers/order differ across RSEM files.")
  mat <- do.call(cbind, vectors)
  rownames(mat) <- reference_genes
  colnames(mat) <- used
  list(counts = mat, raw_samples = raw_samples, raw_only = setdiff(raw_samples, meta_samples))
}

load_inputs <- function(opts) {
  metadata_path <- arg_get(opts, "metadata")
  sample_id_col <- arg_get(opts, "sample-id-col")
  group_col <- arg_get(opts, "group-col")
  if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path)
  metadata <- fread(metadata_path) |> as_tibble()
  if (!sample_id_col %in% names(metadata)) stop("sample-id-col not found in metadata: ", sample_id_col)
  if (!group_col %in% names(metadata)) stop("group-col not found in metadata: ", group_col)
  metadata[[sample_id_col]] <- as.character(metadata[[sample_id_col]])
  metadata[[group_col]] <- as.character(metadata[[group_col]])
  if (anyDuplicated(metadata[[sample_id_col]])) stop("Metadata contains duplicated sample IDs.")

  tissue_col <- arg_get(opts, "tissue-col")
  tissue_value <- arg_get(opts, "tissue-value")
  if (!is.null(tissue_col) && !tissue_col %in% names(metadata)) stop("tissue-col not found in metadata: ", tissue_col)
  if (!is.null(tissue_value)) {
    if (is.null(tissue_col)) stop("--tissue-value requires --tissue-col")
    metadata <- metadata |> filter(.data[[tissue_col]] == tissue_value)
    if (nrow(metadata) == 0) stop("No metadata rows remain after tissue filter: ", tissue_col, " == ", tissue_value)
  }

  has_counts <- !is.null(opts[["counts"]])
  has_rsem <- !is.null(opts[["rsem-dir"]])
  if (identical(has_counts, has_rsem)) stop("Provide exactly one of --counts or --rsem-dir.")
  if (has_counts) {
    counts <- read_counts_csv(opts[["counts"]])
    raw_only <- character()
  } else {
    r <- read_rsem_counts(opts[["rsem-dir"]], metadata, sample_id_col)
    counts <- r$counts
    raw_only <- r$raw_only
  }
  meta_samples <- metadata[[sample_id_col]]
  if (!setequal(colnames(counts), meta_samples)) {
    stop("Counts/RSEM samples do not match metadata after filtering. Only in counts: ",
         paste(setdiff(colnames(counts), meta_samples), collapse = ";"),
         " Only in metadata: ", paste(setdiff(meta_samples, colnames(counts)), collapse = ";"))
  }
  counts <- counts[, meta_samples, drop = FALSE]
  metadata <- as.data.frame(metadata)
  metadata$sample_id <- metadata[[sample_id_col]]
  metadata$group <- metadata[[group_col]]
  rownames(metadata) <- metadata$sample_id
  list(counts = counts, metadata = metadata, raw_only = raw_only)
}

build_model_registry <- function(models_dir, config_file) {
  cfg <- fread(config_file) |> as_tibble()
  rows <- list()
  ri <- 0L
  for (i in seq_len(nrow(cfg))) {
    for (preprocess in c("scaled_diff", "yugene_diff")) {
      file_col <- if (preprocess == "scaled_diff") "scaled_diff_file" else "yugene_diff_file"
      md5_col <- if (preprocess == "scaled_diff") "scaled_diff_md5" else "yugene_diff_md5"
      fn <- cfg[[file_col]][[i]]
      parts <- strsplit(sub("[.]pkl$", "", fn), "_", fixed = TRUE)[[1]]
      ri <- ri + 1L
      rows[[ri]] <- tibble(
        model_family = cfg$model_family[[i]],
        target = tolower(cfg$target[[i]]),
        output_stem = cfg$output_stem[[i]],
        preprocess = preprocess,
        model_species = parts[[3]],
        model_tissue = paste(parts[4:(length(parts) - 1L)], collapse = "_"),
        file = fn,
        path = file.path(models_dir, fn),
        expected_md5 = cfg[[md5_col]][[i]],
        exists = file.exists(file.path(models_dir, fn))
      )
    }
  }
  bind_rows(rows)
}

expand_choice <- function(x, allowed) {
  if (tolower(x) == "all") allowed else x
}

resolve_auto_model_species <- function(species, registry, requested) {
  if (tolower(requested) != "auto") return(requested)
  preferred <- switch(tolower(species),
    mouse = c("Mouse", "Multispecies"),
    rat = c("Rodents", "Multispecies"),
    monkey = c("Multispecies"),
    human = c("Multispecies"),
    c("Multispecies")
  )
  preferred[preferred %in% registry$model_species][[1]]
}

resolve_auto_model_tissue <- function(tissue_value, registry, requested) {
  if (tolower(requested) != "auto") return(requested)
  if (!is.null(tissue_value)) {
    normalized <- safe_name(tissue_value)
    hit <- registry$model_tissue[tolower(registry$model_tissue) == tolower(normalized)]
    if (length(hit) > 0) return(hit[[1]])
  }
  "Multitissue"
}

select_models <- function(opts, registry) {
  families <- expand_choice(toupper(arg_get(opts, "model-family", "EN")), c("EN", "BR"))
  targets <- expand_choice(tolower(arg_get(opts, "target", "mortality")), c("chronoage", "mortality", "normalizedage"))
  preprocesses <- expand_choice(tolower(arg_get(opts, "preprocess", "scaled_diff")), c("scaled_diff", "yugene_diff"))
  model_species <- resolve_auto_model_species(arg_get(opts, "species", "human"), registry, arg_get(opts, "model-species", "Multispecies"))
  model_tissue <- resolve_auto_model_tissue(arg_get(opts, "tissue-value"), registry, arg_get(opts, "model-tissue", "Multitissue"))
  selected <- registry |>
    filter(.data$model_family %in% families, .data$target %in% targets, .data$preprocess %in% preprocesses,
           .data$model_species == model_species, .data$model_tissue == model_tissue)
  if (nrow(selected) == 0) {
    target_display <- c(chronoage = "Chronoage", mortality = "Mortality", normalizedage = "NormalizedAge")
    expected <- as.vector(outer(families, targets, Vectorize(function(f, t) {
      paste0(f, "_", target_display[[t]], "_", model_species, "_", model_tissue, "_", preprocesses, ".pkl")
    })))
    stop("Requested model combination is not in the local registry. Expected pattern(s): ",
         paste(unique(expected), collapse = "; "),
         ". Models are not auto-downloaded in v1. Source: ", zenodo_record_url)
  }
  missing <- selected |> filter(!.data$exists)
  if (nrow(missing) > 0) stop("Selected model file(s) missing locally: ", paste(missing$file, collapse = "; "))
  observed <- as.character(tools::md5sum(selected$path))
  if (!all(observed == selected$expected_md5)) {
    bad <- selected$file[observed != selected$expected_md5]
    stop("Selected model MD5 mismatch: ", paste(bad, collapse = "; "))
  }
  selected$observed_md5 <- observed
  selected
}

direction_correct <- function(set_mode, delta) {
  dplyr::case_when(set_mode == "aged_up" ~ -delta, set_mode == "aged_down" ~ delta, TRUE ~ NA_real_)
}

make_out_dirs <- function(out_dir) {
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  list(
    base = out_dir,
    prepared = file.path(out_dir, "prepared"),
    qc = file.path(out_dir, "qc"),
    tage = file.path(out_dir, "tage_predictions"),
    hallmark = file.path(out_dir, "hallmark_results"),
    figures = file.path(out_dir, "figures"),
    report = file.path(out_dir, "report"),
    logs = file.path(out_dir, "logs")
  )
}

save_run_config <- function(opts, out_dirs) {
  cfg <- list(opts = opts, out_dirs = out_dirs, project_dir = project_dir, zenodo_record_url = zenodo_record_url)
  saveRDS(cfg, file.path(out_dirs$qc, "run_config.rds"))
  cfg
}

load_run_config <- function(path) {
  readRDS(path)
}

script_path <- function(name) file.path(project_dir, "scripts", "tage_hallmarks", name)
