#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this_file), "lib_common.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(GSVA)
  library(limma)
})

run_ssgsea <- function(expr, gene_sets, min_size, max_size) {
  gene_sets <- lapply(gene_sets, function(x) intersect(unique(x), rownames(expr)))
  gene_sets <- gene_sets[lengths(gene_sets) >= min_size & lengths(gene_sets) <= max_size]
  if (length(gene_sets) == 0) stop("No hallmark gene sets passed size filters after intersecting expression genes.")
  if (exists("ssgseaParam", asNamespace("GSVA"))) {
    param <- GSVA::ssgseaParam(exprData = expr, geneSets = gene_sets, minSize = min_size, maxSize = max_size, normalize = TRUE)
    GSVA::gsva(param, verbose = TRUE)
  } else {
    GSVA::gsva(expr, gene_sets, method = "ssgsea", min.sz = min_size, max.sz = max_size, ssgsea.norm = TRUE, verbose = TRUE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
config_file <- sub("^--config=", "", args[startsWith(args, "--config=")][1])
if (is.na(config_file)) stop("Usage: 04_run_hallmark_ssgsea.R --config=/path/run_config.rds")
cfg <- load_run_config(config_file)
opts <- cfg$opts
out_dirs <- cfg$out_dirs

msg("Step 04 run DESeq2 VST and Open Genes Hallmark ssGSEA")
counts <- readRDS(file.path(out_dirs$prepared, "raw_count_matrix.rds"))
metadata <- readRDS(file.path(out_dirs$prepared, "metadata_used.rds"))
control_group <- arg_get(opts, "control-group")
deg_results <- list()

count_round <- round(counts)
storage.mode(count_round) <- "integer"
keep <- rowSums(count_round) >= as.integer(arg_get(opts, "low-count-min-total", 10))
if (sum(keep) == 0) stop("No genes passed low-count filter.")
metadata$group <- factor(metadata$group)
if (!control_group %in% levels(metadata$group)) stop("Control group absent from metadata group column: ", control_group)
metadata$group <- relevel(metadata$group, ref = control_group)
dds <- DESeqDataSetFromMatrix(countData = count_round[keep, , drop = FALSE], colData = metadata, design = ~ group)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)
saveRDS(expr, file.path(out_dirs$prepared, "vst_normalized_expression.rds"))
fwrite(as.data.frame(expr) |> rownames_to_column("gene"), file.path(out_dirs$prepared, "vst_normalized_expression.csv"))
for (g in setdiff(levels(metadata$group), control_group)) {
  res <- results(dds, contrast = c("group", g, control_group)) |> as.data.frame() |> rownames_to_column("gene")
  deg_results[[g]] <- as_tibble(res)
  fwrite(res, file.path(out_dirs$hallmark, paste0("DESeq2_", safe_name(g), "_vs_", safe_name(control_group), ".csv")))
}

gmt_file <- arg_get(opts, "hallmark-gmt", default_hallmark_gmt)
if (!file.exists(gmt_file)) stop("Hallmark GMT not found: ", gmt_file)
all_sets <- read_gmt(gmt_file)
set_meta <- bind_rows(lapply(names(all_sets), parse_set_metadata)) |>
  filter(set_mode %in% set_mode_order, aging_hallmark %in% hallmark_order)
gene_sets <- all_sets[set_meta$input_set_name]
names(gene_sets) <- set_meta$input_set_name

scores <- run_ssgsea(expr, gene_sets, as.integer(arg_get(opts, "min-gene-set-size", 5)), as.integer(arg_get(opts, "max-gene-set-size", 1000)))
set_meta <- set_meta |> filter(input_set_name %in% rownames(scores))
directed_gene_membership <- bind_rows(lapply(seq_len(nrow(set_meta)), function(i) {
  tibble(
    input_set_name = set_meta$input_set_name[[i]],
    aging_hallmark = set_meta$aging_hallmark[[i]],
    set_mode = set_meta$set_mode[[i]],
    gene = intersect(unique(all_sets[[set_meta$input_set_name[[i]]]]), rownames(expr))
  )
})) |>
  distinct(input_set_name, aging_hallmark, set_mode, gene)
fwrite(directed_gene_membership, file.path(out_dirs$hallmark, "open_genes_aged_direction_gene_membership.csv"))
fwrite(set_meta, file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_gene_set_metadata.csv"))
fwrite(as.data.frame(scores) |> rownames_to_column("input_set_name"), file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_score_matrix.csv"))

long <- as.data.frame(scores) |>
  rownames_to_column("input_set_name") |>
  as_tibble() |>
  pivot_longer(-input_set_name, names_to = "sample_id", values_to = "ssgsea_score") |>
  left_join(set_meta, by = "input_set_name") |>
  left_join(as_tibble(metadata) |> select(sample_id, group), by = "sample_id")
fwrite(long, file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_scores_long.csv"))

metadata$group <- factor(metadata$group)
design <- model.matrix(~ 0 + group, data = metadata)
original_levels <- sub("^group", "", colnames(design))
safe_levels <- make.names(original_levels)
colnames(design) <- safe_levels
group_key <- setNames(safe_levels, original_levels)
contrast_groups <- setdiff(original_levels, control_group)
contrast_matrix <- makeContrasts(contrasts = paste0(group_key[contrast_groups], "-", group_key[[control_group]]), levels = design)
fit <- eBayes(contrasts.fit(lmFit(scores, design), contrast_matrix))
limma_tbl <- bind_rows(lapply(seq_along(contrast_groups), function(i) {
  g <- contrast_groups[[i]]
  topTable(fit, coef = i, number = Inf, sort.by = "none") |>
    rownames_to_column("input_set_name") |>
    as_tibble() |>
    transmute(input_set_name, intervention = g, control = control_group,
              ssgsea_delta_vs_control = logFC, average_ssgsea_score = AveExpr,
              t_statistic = t, pvalue = P.Value, padj = adj.P.Val, B = B)
})) |>
  left_join(set_meta, by = "input_set_name") |>
  mutate(direction_corrected_ssgsea_reversal = direction_correct(set_mode, ssgsea_delta_vs_control),
         reversal_call = case_when(
           pvalue < 0.05 & direction_corrected_ssgsea_reversal > 0 ~ "pvalue_significant_reversal",
           pvalue < 0.05 & direction_corrected_ssgsea_reversal < 0 ~ "pvalue_significant_aging_direction",
           direction_corrected_ssgsea_reversal > 0 ~ "trend_reversal",
           direction_corrected_ssgsea_reversal < 0 ~ "trend_aging_direction",
           TRUE ~ "neutral"
         ),
         pvalue_significance = stars_from_p(pvalue)) |>
  arrange(factor(aging_hallmark, levels = hallmark_order), factor(set_mode, levels = set_mode_order), intervention)

summary <- limma_tbl |>
  group_by(intervention) |>
  summarise(n_direction_sets_tested = n(),
            n_reversal_trend = sum(direction_corrected_ssgsea_reversal > 0, na.rm = TRUE),
            n_aging_direction_trend = sum(direction_corrected_ssgsea_reversal < 0, na.rm = TRUE),
            fraction_reversal_trend = n_reversal_trend / n_direction_sets_tested,
            n_pvalue_significant_reversal = sum(reversal_call == "pvalue_significant_reversal", na.rm = TRUE),
            n_pvalue_significant_aging_direction = sum(reversal_call == "pvalue_significant_aging_direction", na.rm = TRUE),
            mean_direction_corrected_ssgsea_reversal = mean(direction_corrected_ssgsea_reversal, na.rm = TRUE),
            median_direction_corrected_ssgsea_reversal = median(direction_corrected_ssgsea_reversal, na.rm = TRUE), .groups = "drop") |>
  mutate(net_pvalue_significant_reversal_sets = n_pvalue_significant_reversal - n_pvalue_significant_aging_direction)
hallmark_summary <- limma_tbl |>
  group_by(aging_hallmark, intervention) |>
  summarise(n_direction_sets_tested = n(), n_reversal_trend = sum(direction_corrected_ssgsea_reversal > 0, na.rm = TRUE),
            n_pvalue_significant_reversal = sum(reversal_call == "pvalue_significant_reversal", na.rm = TRUE),
            mean_direction_corrected_ssgsea_reversal = mean(direction_corrected_ssgsea_reversal, na.rm = TRUE), .groups = "drop")

deg_reversal_long <- bind_rows(lapply(names(deg_results), function(g) {
  deg_results[[g]] |>
    select(gene, log2FoldChange, pvalue, padj) |>
    mutate(intervention = g, control = control_group)
})) |>
  inner_join(directed_gene_membership, by = "gene") |>
  filter(!is.na(log2FoldChange)) |>
  mutate(
    direction_corrected_log2FC = direction_correct(set_mode, log2FoldChange),
    improved = direction_corrected_log2FC > 0,
    aging_direction = direction_corrected_log2FC < 0
  )

directed_background <- deg_reversal_long |>
  group_by(intervention, gene, set_mode) |>
  summarise(
    direction_corrected_log2FC = mean(direction_corrected_log2FC, na.rm = TRUE),
    improved = direction_corrected_log2FC > 0,
    .groups = "drop"
  )

fisher_rows <- lapply(split(deg_reversal_long, list(deg_reversal_long$intervention, deg_reversal_long$aging_hallmark), drop = TRUE), function(x) {
  intervention <- unique(x$intervention)[[1]]
  hallmark <- unique(x$aging_hallmark)[[1]]
  bg <- directed_background[directed_background$intervention == intervention, , drop = FALSE]
  hallmark_keys <- x |> distinct(gene, set_mode) |> mutate(key = paste(gene, set_mode, sep = "\r")) |> pull(key)
  bg <- bg |> mutate(key = paste(gene, set_mode, sep = "\r"), in_hallmark = key %in% hallmark_keys)
  a <- sum(bg$in_hallmark & bg$improved, na.rm = TRUE)
  b <- sum(bg$in_hallmark & !bg$improved, na.rm = TRUE)
  c <- sum(!bg$in_hallmark & bg$improved, na.rm = TRUE)
  d <- sum(!bg$in_hallmark & !bg$improved, na.rm = TRUE)
  ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater"), error = function(e) NULL)
  tibble(
    intervention = intervention,
    aging_hallmark = hallmark,
    fisher_improved_in_hallmark = a,
    fisher_not_improved_in_hallmark = b,
    fisher_improved_background = c,
    fisher_not_improved_background = d,
    fisher_odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
    fisher_pvalue = if (is.null(ft)) NA_real_ else ft$p.value
  )
}) |>
  bind_rows() |>
  group_by(intervention) |>
  mutate(fisher_padj = p.adjust(fisher_pvalue, method = "BH")) |>
  ungroup()

gene_reversal_summary <- deg_reversal_long |>
  group_by(intervention, aging_hallmark) |>
  summarise(
    n_detected_directional_gene_entries = n(),
    n_unique_detected_genes = n_distinct(gene),
    n_improved_gene_entries = sum(improved, na.rm = TRUE),
    n_aging_direction_gene_entries = sum(aging_direction, na.rm = TRUE),
    fraction_improved = n_improved_gene_entries / n_detected_directional_gene_entries,
    percent_improved = 100 * fraction_improved,
    reversal_balance = 2 * fraction_improved - 1,
    mean_direction_corrected_log2FC = mean(direction_corrected_log2FC, na.rm = TRUE),
    median_direction_corrected_log2FC = median(direction_corrected_log2FC, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(fisher_rows, by = c("intervention", "aging_hallmark")) |>
  mutate(
    # Bounded gene-level score:
    #   magnitude_component = direction-corrected DEG strength, compressed to [-1, 1]
    #   percentage_component = percent improved centered at 50%, also in [-1, 1]
    # Positive combined_rejuvenation_score means the Hallmark moves opposite to aging.
    magnitude_component = tanh(mean_direction_corrected_log2FC),
    percentage_component = reversal_balance,
    combined_rejuvenation_score = 0.5 * magnitude_component + 0.5 * percentage_component,
    fisher_signed_log10_pvalue = if_else(
      is.na(fisher_pvalue) | is.na(fisher_odds_ratio),
      NA_real_,
      -log10(pmax(fisher_pvalue, .Machine$double.xmin)) * sign(log2(fisher_odds_ratio))
    ),
    score_interpretation = case_when(
      combined_rejuvenation_score > 0 ~ "rejuvenation_direction",
      combined_rejuvenation_score < 0 ~ "aging_direction",
      TRUE ~ "neutral"
    )
  ) |>
  arrange(factor(aging_hallmark, levels = hallmark_order), intervention)

integrated_hallmark_summary <- hallmark_summary |>
  left_join(gene_reversal_summary, by = c("aging_hallmark", "intervention")) |>
  mutate(
    # Final Hallmark score combines sample-level ssGSEA reversal and gene-level
    # reversal percentage/strength. Positive values indicate rejuvenation direction.
    integrated_rejuvenation_score = 0.5 * tanh(mean_direction_corrected_ssgsea_reversal) + 0.5 * combined_rejuvenation_score,
    integrated_score_interpretation = case_when(
      integrated_rejuvenation_score > 0 ~ "rejuvenation_direction",
      integrated_rejuvenation_score < 0 ~ "aging_direction",
      TRUE ~ "neutral"
    )
  )

integrated_intervention_summary <- integrated_hallmark_summary |>
  group_by(intervention) |>
  summarise(
    n_hallmarks_tested = n(),
    n_rejuvenation_hallmarks = sum(integrated_rejuvenation_score > 0, na.rm = TRUE),
    fraction_rejuvenation_hallmarks = n_rejuvenation_hallmarks / n_hallmarks_tested,
    mean_integrated_rejuvenation_score = mean(integrated_rejuvenation_score, na.rm = TRUE),
    median_integrated_rejuvenation_score = median(integrated_rejuvenation_score, na.rm = TRUE),
    mean_gene_combined_rejuvenation_score = mean(combined_rejuvenation_score, na.rm = TRUE),
    mean_ssgsea_rejuvenation_score = mean(mean_direction_corrected_ssgsea_reversal, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(limma_tbl, file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_limma_vs_control.csv"))
fwrite(summary, file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_intervention_rejuvenation_summary.csv"))
fwrite(hallmark_summary, file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_hallmark_rejuvenation_summary.csv"))
fwrite(deg_reversal_long, file.path(out_dirs$hallmark, "open_genes_aged_direction_gene_level_reversal.csv"))
fwrite(gene_reversal_summary, file.path(out_dirs$hallmark, "open_genes_hallmark_gene_reversal_percentage_fisher_summary.csv"))
fwrite(integrated_hallmark_summary, file.path(out_dirs$hallmark, "open_genes_hallmark_integrated_rejuvenation_summary.csv"))
fwrite(integrated_intervention_summary, file.path(out_dirs$hallmark, "open_genes_intervention_integrated_rejuvenation_summary.csv"))
msg("Step 04 completed")
