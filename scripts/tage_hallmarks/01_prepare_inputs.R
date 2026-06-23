#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this_file), "lib_common.R"))

args <- commandArgs(trailingOnly = TRUE)
config_file <- sub("^--config=", "", args[startsWith(args, "--config=")][1])
if (is.na(config_file)) stop("Usage: 01_prepare_inputs.R --config=/path/run_config.rds")
cfg <- load_run_config(config_file)
opts <- cfg$opts
out_dirs <- cfg$out_dirs

msg("Step 01 prepare inputs")
inputs <- load_inputs(opts)
control_group <- arg_get(opts, "control-group")
if (!control_group %in% inputs$metadata$group) stop("control-group not found in metadata group column: ", control_group)

saveRDS(inputs$counts, file.path(out_dirs$prepared, "raw_count_matrix.rds"))
saveRDS(inputs$metadata, file.path(out_dirs$prepared, "metadata_used.rds"))
fwrite(data.frame(gene = rownames(inputs$counts), inputs$counts, check.names = FALSE), file.path(out_dirs$prepared, "raw_count_matrix.csv"))
fwrite(inputs$metadata, file.path(out_dirs$prepared, "metadata_used.csv"))

if (length(inputs$raw_only) > 0) {
  fwrite(tibble(sample_id = inputs$raw_only, reason = "raw_RSEM_without_metadata"), file.path(out_dirs$qc, "excluded_samples.tsv"), sep = "\t")
} else {
  fwrite(tibble(sample_id = character(), reason = character()), file.path(out_dirs$qc, "excluded_samples.tsv"), sep = "\t")
}

input_qc <- tibble(
  input_genes = nrow(inputs$counts),
  input_samples = ncol(inputs$counts),
  metadata_rows = nrow(inputs$metadata),
  control_group = control_group,
  groups = paste(unique(inputs$metadata$group), collapse = ";"),
  min_library_size = min(colSums(inputs$counts, na.rm = TRUE)),
  max_library_size = max(colSums(inputs$counts, na.rm = TRUE)),
  noninteger_values = sum(inputs$counts != round(inputs$counts), na.rm = TRUE),
  negative_values = sum(inputs$counts < 0, na.rm = TRUE)
)
fwrite(input_qc, file.path(out_dirs$qc, "input_qc.tsv"), sep = "\t")
msg("Step 01 completed: ", nrow(inputs$counts), " genes, ", ncol(inputs$counts), " samples")
