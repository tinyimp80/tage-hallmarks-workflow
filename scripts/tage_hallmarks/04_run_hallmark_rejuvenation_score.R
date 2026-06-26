#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this_file), "lib_common.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
config_file <- sub("^--config=", "", args[startsWith(args, "--config=")][1])
if (is.na(config_file)) stop("Usage: 04_run_hallmark_rejuvenation_score.R --config=/path/run_config.rds")
cfg <- load_run_config(config_file)
opts <- cfg$opts
out_dirs <- cfg$out_dirs

msg("Step 04 run DESeq2 VST and Open Genes directional Hallmark scoring")
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

directed_gene_membership <- bind_rows(lapply(seq_len(nrow(set_meta)), function(i) {
  tibble(
    input_set_name = set_meta$input_set_name[[i]],
    aging_hallmark = set_meta$aging_hallmark[[i]],
    set_mode = set_meta$set_mode[[i]],
    gene = unique(all_sets[[set_meta$input_set_name[[i]]]])
  )
})) |>
  distinct(input_set_name, aging_hallmark, set_mode, gene)
fwrite(directed_gene_membership, file.path(out_dirs$hallmark, "open_genes_aged_direction_gene_membership.csv"))

directional_gene_membership <- directed_gene_membership |>
  filter(gene %in% rownames(expr)) |>
  mutate(direction_weight = case_when(
    set_mode == "aged_up" ~ 1,
    set_mode == "aged_down" ~ -1,
    TRUE ~ NA_real_
  )) |>
  filter(!is.na(direction_weight))

min_size <- as.integer(arg_get(opts, "min-gene-set-size", 5))
max_size <- as.integer(arg_get(opts, "max-gene-set-size", 1000))
hallmark_gene_summary <- directional_gene_membership |>
  group_by(aging_hallmark) |>
  summarise(
    n_detected_directional_gene_entries = n(),
    n_detected_aged_up_entries = sum(set_mode == "aged_up"),
    n_detected_aged_down_entries = sum(set_mode == "aged_down"),
    n_unique_detected_genes = n_distinct(gene),
    .groups = "drop"
  ) |>
  filter(n_detected_directional_gene_entries >= min_size, n_detected_directional_gene_entries <= max_size)
if (nrow(hallmark_gene_summary) == 0) stop("No Hallmark passed directional gene count filters.")
directional_gene_membership <- directional_gene_membership |>
  filter(aging_hallmark %in% hallmark_gene_summary$aging_hallmark)
fwrite(directional_gene_membership, file.path(out_dirs$hallmark, "open_genes_directional_hallmark_gene_membership_detected.csv"))
fwrite(hallmark_gene_summary, file.path(out_dirs$hallmark, "open_genes_directional_hallmark_gene_set_summary.csv"))

expr_z <- t(scale(t(expr)))
expr_z[!is.finite(expr_z)] <- 0
score_rows <- lapply(hallmark_order[hallmark_order %in% hallmark_gene_summary$aging_hallmark], function(hallmark) {
  m <- directional_gene_membership |> filter(aging_hallmark == hallmark)
  weighted_expr <- expr_z[m$gene, , drop = FALSE] * m$direction_weight
  colMeans(weighted_expr, na.rm = TRUE)
})
scores <- do.call(rbind, score_rows)
rownames(scores) <- hallmark_order[hallmark_order %in% hallmark_gene_summary$aging_hallmark]
fwrite(as.data.frame(scores) |> rownames_to_column("aging_hallmark"), file.path(out_dirs$hallmark, "open_genes_directional_hallmark_activity_score_matrix.csv"))

long <- as.data.frame(scores) |>
  rownames_to_column("aging_hallmark") |>
  as_tibble() |>
  pivot_longer(-aging_hallmark, names_to = "sample_id", values_to = "directional_hallmark_aging_activity_score") |>
  left_join(as_tibble(metadata) |> select(sample_id, group), by = "sample_id")
fwrite(long, file.path(out_dirs$hallmark, "open_genes_directional_hallmark_activity_scores_long.csv"))

metadata$group <- factor(metadata$group)
design <- model.matrix(~ 0 + group, data = metadata)
original_levels <- sub("^group", "", colnames(design))
safe_levels <- make.names(original_levels)
colnames(design) <- safe_levels
group_key <- setNames(safe_levels, original_levels)
contrast_groups <- setdiff(original_levels, control_group)
contrast_matrix <- makeContrasts(contrasts = paste0(group_key[contrast_groups], "-", group_key[[control_group]]), levels = design)
fit <- eBayes(contrasts.fit(lmFit(scores, design), contrast_matrix))
hallmark_activity_limma <- bind_rows(lapply(seq_along(contrast_groups), function(i) {
  g <- contrast_groups[[i]]
  topTable(fit, coef = i, number = Inf, sort.by = "none") |>
    rownames_to_column("aging_hallmark") |>
    as_tibble() |>
    transmute(aging_hallmark, intervention = g, control = control_group,
              aging_activity_delta_vs_control = logFC,
              directional_hallmark_rejuvenation_score = -logFC,
              average_directional_hallmark_aging_activity_score = AveExpr,
              t_statistic = t, pvalue = P.Value, padj = adj.P.Val, B = B)
})) |>
  left_join(hallmark_gene_summary, by = "aging_hallmark") |>
  mutate(
    reversal_call = case_when(
      pvalue < 0.05 & directional_hallmark_rejuvenation_score > 0 ~ "pvalue_significant_reversal",
      pvalue < 0.05 & directional_hallmark_rejuvenation_score < 0 ~ "pvalue_significant_aging_direction",
      directional_hallmark_rejuvenation_score > 0 ~ "trend_reversal",
      directional_hallmark_rejuvenation_score < 0 ~ "trend_aging_direction",
      TRUE ~ "neutral"
    ),
    pvalue_significance = stars_from_p(pvalue)
  ) |>
  arrange(factor(aging_hallmark, levels = hallmark_order), intervention)

hallmark_activity_intervention_summary <- hallmark_activity_limma |>
  group_by(intervention) |>
  summarise(n_hallmarks_tested = n(),
            n_reversal_trend = sum(directional_hallmark_rejuvenation_score > 0, na.rm = TRUE),
            n_aging_direction_trend = sum(directional_hallmark_rejuvenation_score < 0, na.rm = TRUE),
            fraction_reversal_trend = n_reversal_trend / n_hallmarks_tested,
            n_pvalue_significant_reversal = sum(reversal_call == "pvalue_significant_reversal", na.rm = TRUE),
            n_pvalue_significant_aging_direction = sum(reversal_call == "pvalue_significant_aging_direction", na.rm = TRUE),
            mean_directional_hallmark_rejuvenation_score = mean(directional_hallmark_rejuvenation_score, na.rm = TRUE),
            median_directional_hallmark_rejuvenation_score = median(directional_hallmark_rejuvenation_score, na.rm = TRUE),
            .groups = "drop") |>
  mutate(net_pvalue_significant_reversal_hallmarks = n_pvalue_significant_reversal - n_pvalue_significant_aging_direction)

deg_reversal_long <- bind_rows(lapply(names(deg_results), function(g) {
  deg_results[[g]] |>
    select(gene, log2FoldChange, pvalue, padj) |>
    mutate(intervention = g, control = control_group)
})) |>
  inner_join(directional_gene_membership, by = "gene") |>
  filter(!is.na(log2FoldChange)) |>
  mutate(
    direction_corrected_log2FC = direction_correct(set_mode, log2FoldChange),
    improved = direction_corrected_log2FC > 0,
    aging_direction = direction_corrected_log2FC < 0
  )

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
  mutate(
    fraction_reversed_genes = fraction_improved,
    percent_reversed_genes = percent_improved,
    mean_rejuvenation_log2FC = mean_direction_corrected_log2FC,
    median_rejuvenation_log2FC = median_direction_corrected_log2FC,
    # Bounded gene-level score:
    #   magnitude_component = direction-corrected DEG strength, compressed to [-1, 1]
    #   percentage_component = percent improved centered at 50%, also in [-1, 1]
    # Positive combined_rejuvenation_score means the Hallmark moves opposite to aging.
    magnitude_component = tanh(mean_direction_corrected_log2FC),
    percentage_component = reversal_balance,
    combined_rejuvenation_score = 0.5 * magnitude_component + 0.5 * percentage_component,
    score_interpretation = case_when(
      combined_rejuvenation_score > 0 ~ "rejuvenation_direction",
      combined_rejuvenation_score < 0 ~ "aging_direction",
      TRUE ~ "neutral"
    )
  ) |>
  arrange(factor(aging_hallmark, levels = hallmark_order), intervention)

gene_reversal_intervention_summary <- gene_reversal_summary |>
  group_by(intervention) |>
  summarise(
    mean_fraction_reversed_genes = mean(fraction_reversed_genes, na.rm = TRUE),
    median_fraction_reversed_genes = median(fraction_reversed_genes, na.rm = TRUE),
    mean_rejuvenation_log2FC = mean(mean_rejuvenation_log2FC, na.rm = TRUE),
    median_rejuvenation_log2FC = median(median_rejuvenation_log2FC, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(hallmark_activity_limma, file.path(out_dirs$hallmark, "open_genes_directional_hallmark_activity_limma_vs_control.csv"))
fwrite(hallmark_activity_intervention_summary, file.path(out_dirs$hallmark, "open_genes_directional_hallmark_activity_intervention_rejuvenation_summary.csv"))
fwrite(deg_reversal_long, file.path(out_dirs$hallmark, "open_genes_aged_direction_gene_level_reversal.csv"))
fwrite(gene_reversal_summary, file.path(out_dirs$hallmark, "open_genes_hallmark_gene_reversal_percentage_summary.csv"))
fwrite(gene_reversal_intervention_summary, file.path(out_dirs$hallmark, "open_genes_hallmark_gene_reversal_intervention_summary.csv"))
msg("Step 04 completed")
