#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this_file), "lib_common.R"))

args <- commandArgs(trailingOnly = TRUE)
config_file <- sub("^--config=", "", args[startsWith(args, "--config=")][1])
if (is.na(config_file)) stop("Usage: 02_validate_models.R --config=/path/run_config.rds")
cfg <- load_run_config(config_file)
opts <- cfg$opts
out_dirs <- cfg$out_dirs

msg("Step 02 validate tAge model registry")
registry <- build_model_registry(file.path(project_dir, "models"), file.path(project_dir, "config", "validation_models.csv"))
selected_models <- select_models(opts, registry)
fwrite(registry, file.path(out_dirs$qc, "model_registry.csv"))
fwrite(selected_models, file.path(out_dirs$qc, "selected_model_checksums.csv"))
saveRDS(selected_models, file.path(out_dirs$qc, "selected_models.rds"))
msg("Step 02 completed: ", nrow(selected_models), " model file(s) selected")
