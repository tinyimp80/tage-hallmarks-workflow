#!/usr/bin/env Rscript

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(this_file, mustWork = TRUE))
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
Sys.setenv(TAGE_PROJECT_DIR = project_root)
source(file.path(script_dir, "tage_hallmarks", "lib_common.R"))

run_step <- function(step_script, config_file, log_file) {
  rscript <- file.path(R.home("bin"), "Rscript")
  cmd <- c(step_script, paste0("--config=", config_file))
  msg("Running ", basename(step_script))
  output <- system2(rscript, cmd, stdout = TRUE, stderr = TRUE)
  write(output, file = log_file, append = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    cat(output, sep = "\n")
    stop("Step failed: ", step_script)
  }
  invisible(TRUE)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(opts$help)) {
    usage()
    quit(status = 0)
  }
  stop_if_missing(opts, c("metadata", "sample-id-col", "group-col", "control-group", "out-dir"))
  valid_species <- c("human", "mouse", "rat", "monkey")
  if (!tolower(arg_get(opts, "species", "human")) %in% valid_species) stop("--species must be one of: ", paste(valid_species, collapse = ", "))
  if (!arg_get(opts, "gene-mapping-type", "Ensembl") %in% c("Ensembl", "Gene.Symbol")) stop("--gene-mapping-type must be Ensembl or Gene.Symbol")
  if (identical(!is.null(opts[["counts"]]), !is.null(opts[["rsem-dir"]]))) stop("Provide exactly one of --counts or --rsem-dir.")

  out_dirs <- make_out_dirs(arg_get(opts, "out-dir"))
  invisible(lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  config_file <- file.path(out_dirs$qc, "run_config.rds")
  save_run_config(opts, out_dirs)
  log_file <- file.path(out_dirs$logs, "run_tage_hallmarks_report.log")
  write(paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " | Starting modular tAge + Hallmark workflow"), file = log_file)

  steps <- c(
    script_path("01_prepare_inputs.R"),
    script_path("02_validate_models.R")
  )
  if (!isTRUE(opts[["config-only"]])) {
    steps <- c(steps, script_path("03_run_tage.R"), script_path("04_run_hallmark_ssgsea.R"), script_path("05_make_report.R"))
  }
  for (step in steps) run_step(step, config_file, log_file)
  if (isTRUE(opts[["config-only"]])) {
    msg("Config-only validation completed")
  } else {
    msg("Completed")
    msg("HTML report: ", file.path(out_dirs$report, "tage_hallmarks_report.html"))
    msg("PDF report: ", file.path(out_dirs$report, "tage_hallmarks_report.pdf"))
  }
}

tryCatch(main(), error = function(e) {
  message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " | ERROR: ", conditionMessage(e))
  quit(status = 1)
})
