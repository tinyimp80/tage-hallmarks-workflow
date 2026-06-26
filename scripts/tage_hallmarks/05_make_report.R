#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(this_file), "lib_common.R"))

args <- commandArgs(trailingOnly = TRUE)
config_file <- sub("^--config=", "", args[startsWith(args, "--config=")][1])
if (is.na(config_file)) stop("Usage: 05_make_report.R --config=/path/run_config.rds")
cfg <- load_run_config(config_file)
opts <- cfg$opts
out_dirs <- cfg$out_dirs

msg("Step 05 make figures and report")
rejuvenation_green <- "#2A9D8F"
aging_red <- "#E76F51"
neutral_white <- "white"
metadata <- readRDS(file.path(out_dirs$prepared, "metadata_used.rds"))
control_group <- arg_get(opts, "control-group")
group_levels <- c(control_group, setdiff(unique(metadata$group), control_group))
candidate_intervention_groups <- setdiff(group_levels, control_group)
auto_young_reference_groups <- candidate_intervention_groups[
  str_detect(str_to_lower(candidate_intervention_groups), "^(young|young[_ .-]*con|young[_ .-]*control|quiescence|quiescent)$")
]
explicit_young_reference_groups <- arg_get(opts, "young-control-group")
if (is.null(explicit_young_reference_groups)) explicit_young_reference_groups <- arg_get(opts, "young-reference-group")
explicit_young_reference_groups <- if (is.null(explicit_young_reference_groups)) character(0) else str_split(explicit_young_reference_groups, ",", simplify = FALSE)[[1]] |> str_trim()
explicit_exclude_groups <- arg_get(opts, "hallmark-plot-exclude-groups")
explicit_exclude_groups <- if (is.null(explicit_exclude_groups)) character(0) else str_split(explicit_exclude_groups, ",", simplify = FALSE)[[1]] |> str_trim()
hallmark_plot_exclude_groups <- unique(c(auto_young_reference_groups, explicit_young_reference_groups, explicit_exclude_groups))
hallmark_display_groups <- setdiff(candidate_intervention_groups, hallmark_plot_exclude_groups)
if (length(hallmark_display_groups) == 0) {
  warning("No Hallmark display groups remain after excluding young/reference groups; using all non-control groups.")
  hallmark_display_groups <- candidate_intervention_groups
}
msg("Hallmark plots exclude reference groups: ", ifelse(length(hallmark_plot_exclude_groups) == 0, "none", paste(hallmark_plot_exclude_groups, collapse = ", ")))
tests <- fread(file.path(out_dirs$tage, "tAge_vs_control_tests.csv")) |> as_tibble()
activity_tbl <- fread(file.path(out_dirs$hallmark, "open_genes_directional_hallmark_activity_limma_vs_control.csv")) |> as_tibble()
activity_summary <- fread(file.path(out_dirs$hallmark, "open_genes_directional_hallmark_activity_intervention_rejuvenation_summary.csv")) |> as_tibble()
gene_rev <- fread(file.path(out_dirs$hallmark, "open_genes_hallmark_gene_reversal_percentage_summary.csv")) |> as_tibble()
gene_support_intervention <- fread(file.path(out_dirs$hallmark, "open_genes_hallmark_gene_reversal_intervention_summary.csv")) |> as_tibble()

primary <- tests |> filter(model_family == toupper(arg_get(opts, "model-family", "EN")) | tolower(arg_get(opts, "model-family", "EN")) == "all")
primary <- primary |> filter(target == tolower(arg_get(opts, "target", "mortality")) | tolower(arg_get(opts, "target", "mortality")) == "all")
if (nrow(primary) == 0) primary <- tests
primary_preprocess <- if ("scaled_diff" %in% primary$preprocess) "scaled_diff" else primary$preprocess[[1]]
primary <- primary |> filter(preprocess == primary_preprocess)
values_file <- list.files(out_dirs$tage, pattern = paste0("_", primary_preprocess, "_values[.]csv$"), full.names = TRUE)[[1]]
values <- fread(values_file) |> as_tibble()
values$group <- factor(values$group, levels = group_levels[group_levels %in% values$group])

p_box <- ggplot(values, aes(x = group, y = predicted_tAge, fill = group)) +
  geom_hline(yintercept = mean(values$predicted_tAge[values$group == control_group], na.rm = TRUE), linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.shape = NA, alpha = 0.78) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.85) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Predicted tAge", title = "tAge clock", subtitle = "Primary selected model/preprocessing")
ggsave(file.path(out_dirs$figures, "tAge_boxplot.png"), p_box, width = 10.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dirs$figures, "tAge_boxplot.pdf"), p_box, width = 10.5, height = 6.2, bg = "white")

primary$group <- factor(primary$group, levels = setdiff(group_levels, control_group))
p_delta <- ggplot(primary, aes(x = group, y = mean_difference, fill = mean_difference < 0)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(width = 0.72) +
  geom_text(aes(label = significance), vjust = ifelse(primary$mean_difference >= 0, -0.25, 1.15), size = 4) +
  scale_fill_manual(values = c(`TRUE` = rejuvenation_green, `FALSE` = aging_red), guide = "none") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = "Mean tAge difference vs control", title = "tAge delta vs control", subtitle = "Negative values indicate lower predicted age than control")
ggsave(file.path(out_dirs$figures, "tAge_delta_vs_control.png"), p_delta, width = 10.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dirs$figures, "tAge_delta_vs_control.pdf"), p_delta, width = 10.5, height = 6.2, bg = "white")

activity_summary_plot <- activity_summary |>
  filter(intervention %in% hallmark_display_groups) |>
  mutate(intervention = factor(intervention, levels = hallmark_display_groups))
p_activity <- ggplot(activity_summary_plot, aes(x = intervention, y = mean_directional_hallmark_rejuvenation_score)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(aes(fill = mean_directional_hallmark_rejuvenation_score > 0), width = 0.72) +
  geom_text(aes(label = paste0("p<0.05 rev: ", n_pvalue_significant_reversal, "\ntrend: ", n_reversal_trend, "/", n_hallmarks_tested)),
            vjust = ifelse(activity_summary_plot$mean_directional_hallmark_rejuvenation_score >= 0, -0.15, 1.1), size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = c(`TRUE` = rejuvenation_green, `FALSE` = aging_red), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.28, 0.28))) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.margin = margin(12, 24, 28, 18)) +
  labs(x = NULL, y = "Mean directional Hallmark rejuvenation score", title = "Directional Hallmark activity summary", subtitle = "One signed score per Hallmark; positive values indicate movement opposite to aging direction")
ggsave(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_summary.png"), p_activity, width = 10.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_summary.pdf"), p_activity, width = 10.5, height = 6.2, bg = "white")

dot_heatmap_df <- activity_tbl |>
  select(aging_hallmark, intervention, directional_hallmark_rejuvenation_score, pvalue, padj) |>
  left_join(
    gene_rev |>
      select(aging_hallmark, intervention, percent_reversed_genes, fraction_reversed_genes),
    by = c("aging_hallmark", "intervention")
  ) |>
  filter(intervention %in% hallmark_display_groups) |>
  mutate(
    intervention = factor(intervention, levels = hallmark_display_groups),
    aging_hallmark = factor(aging_hallmark, levels = rev(hallmark_order)),
    activity_pvalue_label = stars_from_p(pvalue),
    pvalue_display = case_when(
      is.na(pvalue) ~ "NA",
      pvalue < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", pvalue)
    )
  )
fwrite(dot_heatmap_df, file.path(out_dirs$hallmark, "open_genes_directional_hallmark_rejuvenation_dotplot_values.csv"))
fwrite(dot_heatmap_df, file.path(out_dirs$hallmark, "open_genes_directional_hallmark_rejuvenation_dot_heatmap_values.csv"))

score_lim <- max(abs(dot_heatmap_df$directional_hallmark_rejuvenation_score), na.rm = TRUE)
if (!is.finite(score_lim) || score_lim == 0) score_lim <- 1
size_range <- range(dot_heatmap_df$percent_reversed_genes, na.rm = TRUE)
if (!all(is.finite(size_range))) size_range <- c(0, 100)
size_limits <- c(max(0, floor(size_range[[1]] / 10) * 10), min(100, ceiling(size_range[[2]] / 10) * 10))
if (diff(size_limits) < 20) {
  mid <- mean(size_limits)
  size_limits <- c(max(0, mid - 10), min(100, mid + 10))
}
dot_heatmap_df <- dot_heatmap_df |>
  mutate(percent_reversed_genes_for_size = pmin(pmax(percent_reversed_genes, size_limits[[1]]), size_limits[[2]]))
p_dot_heatmap <- ggplot(dot_heatmap_df, aes(x = intervention, y = aging_hallmark)) +
  geom_point(aes(size = percent_reversed_genes_for_size, fill = directional_hallmark_rejuvenation_score),
             shape = 21, color = "#1f2933", stroke = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%", percent_reversed_genes)), color = "#1f2933", size = 2.35, vjust = 0.95) +
  geom_text(aes(label = activity_pvalue_label), color = "#1f2933", size = 2.9, vjust = -0.45) +
  scale_fill_gradient2(
    low = aging_red,
    mid = neutral_white,
    high = rejuvenation_green,
    midpoint = 0,
    limits = c(-score_lim, score_lim),
    name = "Rejuvenation\nscore"
  ) +
  scale_size_continuous(
    range = c(1.8, 24),
    limits = size_limits,
    breaks = pretty(size_limits, n = 4),
    name = "% reversed\ngenes"
  ) +
  coord_fixed(ratio = 0.72) +
  theme_bw(base_size = 10.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.35),
    legend.position = "right",
    plot.title.position = "plot"
  ) +
  labs(
    title = "Directional Hallmark rejuvenation dot plot",
    subtitle = "Color = VST-based rejuvenation score. Size and label = % genes reversed. Stars = activity p-value. Young/reference controls omitted."
  )
ggsave(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_dotplot.png"), p_dot_heatmap, width = 12.8, height = 7.8, dpi = 240, bg = "white")
ggsave(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_dotplot.pdf"), p_dot_heatmap, width = 12.8, height = 7.8, bg = "white")
# Backward-compatible filenames for older report links.
ggsave(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_dot_heatmap.png"), p_dot_heatmap, width = 12.8, height = 7.8, dpi = 240, bg = "white")
ggsave(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_dot_heatmap.pdf"), p_dot_heatmap, width = 12.8, height = 7.8, bg = "white")

mat <- activity_tbl |>
  filter(intervention %in% hallmark_display_groups) |>
  select(aging_hallmark, intervention, directional_hallmark_rejuvenation_score) |>
  pivot_wider(names_from = intervention, values_from = directional_hallmark_rejuvenation_score) |>
  column_to_rownames("aging_hallmark") |>
  as.matrix()
mat <- mat[hallmark_order[hallmark_order %in% rownames(mat)], hallmark_display_groups[hallmark_display_groups %in% colnames(mat)], drop = FALSE]
lim <- max(abs(mat), na.rm = TRUE)
if (!is.finite(lim) || lim == 0) lim <- 1
png(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_heatmap.png"), width = 2200, height = 1550, res = 220)
pheatmap(mat, color = colorRampPalette(c(aging_red, neutral_white, rejuvenation_green))(100), breaks = seq(-lim, lim, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE, main = "Directional Hallmark activity rejuvenation score", fontsize_row = 8.5, fontsize_col = 10, border_color = NA)
dev.off()
pdf(file.path(out_dirs$figures, "hallmark_directional_activity_rejuvenation_heatmap.pdf"), width = 9.5, height = 7)
pheatmap(mat, color = colorRampPalette(c(aging_red, neutral_white, rejuvenation_green))(100), breaks = seq(-lim, lim, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE, main = "Directional Hallmark activity rejuvenation score", fontsize_row = 8.5, fontsize_col = 10, border_color = NA)
dev.off()

merged <- tests |> group_by(group) |> slice(1) |> ungroup() |>
  filter(group %in% hallmark_display_groups) |>
  left_join(
    activity_summary |> select(intervention, mean_directional_hallmark_rejuvenation_score, n_reversal_trend, n_hallmarks_tested),
    by = c("group" = "intervention")
  ) |>
  left_join(gene_support_intervention, by = c("group" = "intervention"))
table_rows <- paste(apply(merged, 1, function(r) {
  paste0("<tr><td><strong>", html_escape(r[["group"]]), "</strong></td><td class='num'>", fmt_num(as.numeric(r[["mean_difference"]]), 4),
         "</td><td class='num'>", fmt_p(as.numeric(r[["pvalue"]])),
         "</td><td class='num'>", fmt_num(as.numeric(r[["mean_directional_hallmark_rejuvenation_score"]]), 4),
         "</td><td class='num'>", fmt_num(as.numeric(r[["mean_fraction_reversed_genes"]]) * 100, 3),
         "</td><td class='num'>", fmt_num(as.numeric(r[["mean_rejuvenation_log2FC"]]), 4),
         "</td><td class='num'>", r[["n_reversal_trend"]], "/", r[["n_hallmarks_tested"]], "</td></tr>")
}), collapse = "\n")
html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>tAge + Aging Hallmark Report</title>",
  "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;color:#1d1d1f}section,header{padding:56px 7vw}.hero{background:#f5f5f7}.wrap{max-width:1180px;margin:auto}h1{font-size:48px;line-height:1.08;margin:0 0 16px}h2{font-size:32px;margin:0 0 12px}.lead{font-size:21px;color:#333;max-width:900px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:22px}.fig img{width:100%;border:1px solid #e5e5ea}.caption{font-size:13px;color:#6e6e73}table{border-collapse:collapse;width:100%;margin-top:24px}th,td{padding:10px 12px;border-bottom:1px solid #e5e5ea;text-align:left;font-size:14px}.num{text-align:right;font-variant-numeric:tabular-nums}.pill{display:inline-block;border:1px solid #d2d2d7;border-radius:999px;padding:7px 12px;margin:4px;background:white;font-size:13px}@media(max-width:900px){.grid{grid-template-columns:1fr}h1{font-size:36px}}</style>",
  "</head><body><header class='hero'><div class='wrap'><h1>tAge + Open Genes Hallmark Report</h1>",
  "<div class='lead'>One-line RNA-seq report from raw counts/RSEM expected counts. Lower tAge and positive directional Hallmark activity reversal scores indicate rejuvenation-like movement versus control.</div>",
  "<p><span class='pill'>Species: ", html_escape(arg_get(opts, "species", "human")), "</span><span class='pill'>Control: ", html_escape(control_group), "</span><span class='pill'>Model: ", html_escape(arg_get(opts, "model-family", "EN")), " ", html_escape(arg_get(opts, "target", "mortality")), "</span><span class='pill'>Tissue model: ", html_escape(arg_get(opts, "model-tissue", "Multitissue")), "</span><span class='pill'>Preprocess: ", html_escape(arg_get(opts, "preprocess", "scaled_diff")), "</span></p>",
  "</div></header>",
  "<section><div class='wrap'><h2>tAge clock</h2><div class='grid'><div class='fig'><img src='../figures/tAge_boxplot.png'><div class='caption'>Predicted tAge by group.</div></div><div class='fig'><img src='../figures/tAge_delta_vs_control.png'><div class='caption'>Mean tAge difference versus control.</div></div></div></div></section>",
  "<section><div class='wrap'><h2>Directional Hallmark activity</h2><div class='grid'><div class='fig'><img src='../figures/hallmark_directional_activity_rejuvenation_summary.png'><div class='caption'>Mean Activity_Rejuvenation_Score across intervention groups only. Green means rejuvenation direction; red means aging-like direction.</div></div><div class='fig'><img src='../figures/hallmark_directional_activity_rejuvenation_dotplot.png'><div class='caption'>Primary Hallmark summary. Dot color shows VST-based rejuvenation score; dot size and inner percentage show percent of Hallmark genes reversed toward young. Text stars show activity-score p-values. Young/reference controls are omitted.</div></div></div></div></section>",
  "<section><div class='wrap'><h2>Summary</h2><table><thead><tr><th>Group</th><th class='num'>tAge delta</th><th class='num'>tAge p</th><th class='num'>Activity score</th><th class='num'>Mean reversed genes %</th><th class='num'>Mean rejuv log2FC</th><th class='num'>Rev Hallmarks</th></tr></thead><tbody>",
  table_rows,
  "</tbody></table><p class='caption'>Activity score is Activity_Rejuvenation_Score = -1 * (mean activity_intervention - mean activity_aged_control). Positive values indicate rejuvenation-like reversal. Gene-level metrics are supporting evidence from DESeq2 intervention versus aged control.</p></div></section>",
  "<section><div class='wrap'><h2>Methods</h2><p>tAge models are selected from local model files and validated by MD5 checksums. Model inventory follows Zenodo record ",
  html_escape(zenodo_record_url), ". Hallmark activity scores use DESeq2 VST expression and Open Genes aged_up/down genes in one signed score per Hallmark: aged_up genes have positive weights and aged_down genes have negative weights after gene-wise z-scoring. Limma tests directional Hallmark activity versus control. Gene-level reversal percentages are reported only as supporting count-based evidence.</p></div></section>",
  "</body></html>"
)
html_file <- file.path(out_dirs$report, "tage_hallmarks_report.html")
writeLines(html, html_file)
chrome <- Sys.which("google-chrome")
pdf_file <- file.path(out_dirs$report, "tage_hallmarks_report.pdf")
if (nzchar(chrome)) {
  invisible(system2(chrome, c("--headless", "--no-sandbox", "--disable-gpu", paste0("--print-to-pdf=", pdf_file), paste0("file://", normalizePath(html_file, mustWork = TRUE))), stdout = TRUE, stderr = TRUE))
}
manifest <- tibble(parameter = c("status", "project_dir", "out_dir", "html_report", "pdf_report", "species", "gene_mapping_type", "control_group", "zenodo_record"),
                   value = c("success", project_dir, out_dirs$base, html_file, pdf_file, arg_get(opts, "species", "human"), arg_get(opts, "gene-mapping-type", "Ensembl"), control_group, zenodo_record_url))
fwrite(manifest, file.path(out_dirs$qc, "run_manifest.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out_dirs$logs, "session_info.txt"))
msg("Step 05 completed: ", html_file)
