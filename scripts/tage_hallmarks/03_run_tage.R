#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this_file), "lib_common.R"))
suppressPackageStartupMessages({
  library(tAge)
  library(Biobase)
  library(reticulate)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- sub("^--config=", "", args[startsWith(args, "--config=")][1])
if (is.na(config_file)) stop("Usage: 03_run_tage.R --config=/path/run_config.rds")
cfg <- load_run_config(config_file)
opts <- cfg$opts
out_dirs <- cfg$out_dirs

msg("Step 03 run tAge")
counts <- readRDS(file.path(out_dirs$prepared, "raw_count_matrix.rds"))
metadata <- readRDS(file.path(out_dirs$prepared, "metadata_used.rds"))
selected_models <- readRDS(file.path(out_dirs$qc, "selected_models.rds"))

python_bin <- file.path(Sys.getenv("CONDA_PREFIX"), "bin", "python")
if (!file.exists(python_bin)) stop("Pixi Python not found: ", python_bin)
Sys.setenv(RETICULATE_PYTHON = python_bin)
reticulate::use_python(python_bin, required = TRUE)
reticulate::py_run_string("import joblib, pandas, sklearn")

species <- tolower(arg_get(opts, "species", "human"))
gene_mapping_type <- arg_get(opts, "gene-mapping-type", "Ensembl")
control_group <- arg_get(opts, "control-group")
group_col <- arg_get(opts, "group-col")
sample_id_col <- arg_get(opts, "sample-id-col")

eset <- make_ExpressionSet(counts, metadata, verbose = FALSE)
processed <- tAge_preprocessing(
  eset = eset,
  species = species,
  gene_mapping_type = gene_mapping_type,
  verbose = FALSE,
  control_group_column = group_col,
  control_group_label = control_group
)

prediction_tables <- list()
group_summaries <- list()
test_summaries <- list()
model_keys <- unique(paste(selected_models$model_family, selected_models$target, selected_models$model_species, selected_models$model_tissue, sep = "_"))
for (key in model_keys) {
  sub <- selected_models[paste(selected_models$model_family, selected_models$target, selected_models$model_species, selected_models$model_tissue, sep = "_") == key, ]
  model_family <- sub$model_family[[1]]
  target <- sub$target[[1]]
  paths <- as.list(setNames(sub$path, sub$preprocess))
  pred <- predict_tAge(processed, paths, species = species, mode = model_family)
  pred$model_family <- model_family
  pred$target <- target
  pred$model_species <- sub$model_species[[1]]
  pred$model_tissue <- sub$model_tissue[[1]]
  pred_file <- file.path(out_dirs$tage, paste0("tAge_", model_family, "_", target, "_", sub$model_species[[1]], "_", sub$model_tissue[[1]], "_predictions.csv"))
  fwrite(pred, pred_file)
  prediction_tables[[key]] <- pred
  value_cols <- paste0(sub$preprocess, "_", model_family, "_tAge")
  for (value_col in value_cols) {
    if (!value_col %in% names(pred)) next
    long <- pred |>
      as_tibble() |>
      transmute(
        sample_id = .data[[sample_id_col]],
        group = .data[[group_col]],
        model_family = model_family,
        target = target,
        model_species = sub$model_species[[1]],
        model_tissue = sub$model_tissue[[1]],
        preprocess = sub("_[A-Z]+_tAge$", "", value_col),
        predicted_tAge = .data[[value_col]]
      )
    fwrite(long, file.path(out_dirs$tage, paste0("tAge_", model_family, "_", target, "_", long$preprocess[[1]], "_values.csv")))
    group_summary <- long |>
      group_by(group, model_family, target, model_species, model_tissue, preprocess) |>
      summarise(n = n(), mean_tAge = mean(predicted_tAge, na.rm = TRUE), sd_tAge = sd(predicted_tAge, na.rm = TRUE),
                median_tAge = median(predicted_tAge, na.rm = TRUE), .groups = "drop")
    control_vals <- long$predicted_tAge[long$group == control_group]
    tests <- lapply(setdiff(unique(long$group), control_group), function(g) {
      vals <- long$predicted_tAge[long$group == g]
      tibble(group = g, control = control_group, model_family = model_family, target = target,
             model_species = sub$model_species[[1]], model_tissue = sub$model_tissue[[1]], preprocess = long$preprocess[[1]],
             n_group = length(vals), n_control = length(control_vals),
             mean_group = mean(vals, na.rm = TRUE), mean_control = mean(control_vals, na.rm = TRUE),
             mean_difference = mean(vals, na.rm = TRUE) - mean(control_vals, na.rm = TRUE),
             pvalue = tryCatch(t.test(vals, control_vals)$p.value, error = function(e) NA_real_))
    }) |> bind_rows() |> mutate(padj = p.adjust(pvalue, method = "BH"), significance = stars_from_p(pvalue))
    group_summaries[[paste(key, value_col)]] <- group_summary
    test_summaries[[paste(key, value_col)]] <- tests
  }
}

group_summary_all <- bind_rows(group_summaries)
test_summary_all <- bind_rows(test_summaries)
fwrite(group_summary_all, file.path(out_dirs$tage, "tAge_group_summary.csv"))
fwrite(test_summary_all, file.path(out_dirs$tage, "tAge_vs_control_tests.csv"))
fwrite(tibble(parameter = c("retained_genes_after_tAge_preprocessing"), value = c(nrow(processed$RLE_normalized))),
       file.path(out_dirs$tage, "tAge_preprocessing_qc.tsv"), sep = "\t")
msg("Step 03 completed")
