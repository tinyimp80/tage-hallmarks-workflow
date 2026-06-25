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
rejuvenation_green <- "#1B9E77"
aging_red <- "#D95F02"
neutral_white <- "white"
metadata <- readRDS(file.path(out_dirs$prepared, "metadata_used.rds"))
control_group <- arg_get(opts, "control-group")
group_levels <- c(control_group, setdiff(unique(metadata$group), control_group))
tests <- fread(file.path(out_dirs$tage, "tAge_vs_control_tests.csv")) |> as_tibble()
limma_tbl <- fread(file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_limma_vs_control.csv")) |> as_tibble()
hs <- fread(file.path(out_dirs$hallmark, "open_genes_aged_direction_ssgsea_intervention_rejuvenation_summary.csv")) |> as_tibble()
gene_rev <- fread(file.path(out_dirs$hallmark, "open_genes_hallmark_gene_reversal_percentage_fisher_summary.csv")) |> as_tibble()
integrated_hallmark <- fread(file.path(out_dirs$hallmark, "open_genes_hallmark_integrated_rejuvenation_summary.csv")) |> as_tibble()
integrated_intervention <- fread(file.path(out_dirs$hallmark, "open_genes_intervention_integrated_rejuvenation_summary.csv")) |> as_tibble()

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

hs$intervention <- factor(hs$intervention, levels = setdiff(group_levels, control_group))
p_hs <- ggplot(hs, aes(x = intervention, y = mean_direction_corrected_ssgsea_reversal)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(aes(fill = mean_direction_corrected_ssgsea_reversal > 0), width = 0.72) +
  geom_text(aes(label = paste0("p<0.05 rev: ", n_pvalue_significant_reversal, "\ntrend: ", n_reversal_trend, "/", n_direction_sets_tested)),
            vjust = ifelse(hs$mean_direction_corrected_ssgsea_reversal >= 0, -0.15, 1.1), size = 3.0, lineheight = 0.9) +
  scale_fill_manual(values = c(`TRUE` = rejuvenation_green, `FALSE` = aging_red), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.28, 0.28))) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.margin = margin(12, 24, 28, 18)) +
  labs(x = NULL, y = "Mean direction-corrected ssGSEA delta", title = "Hallmark ssGSEA rejuvenation summary", subtitle = "Positive values indicate movement opposite to aging direction")
ggsave(file.path(out_dirs$figures, "hallmark_ssgsea_rejuvenation_summary.png"), p_hs, width = 10.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dirs$figures, "hallmark_ssgsea_rejuvenation_summary.pdf"), p_hs, width = 10.5, height = 6.2, bg = "white")

integrated_intervention$intervention <- factor(integrated_intervention$intervention, levels = setdiff(group_levels, control_group))
p_integrated <- ggplot(integrated_intervention, aes(x = intervention, y = mean_integrated_rejuvenation_score)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(aes(fill = mean_integrated_rejuvenation_score > 0), width = 0.72) +
  geom_text(
    aes(label = paste0("rev hallmarks: ", n_rejuvenation_hallmarks, "/", n_hallmarks_tested)),
    vjust = ifelse(integrated_intervention$mean_integrated_rejuvenation_score >= 0, -0.2, 1.15),
    size = 3.2
  ) +
  scale_fill_manual(values = c(`TRUE` = rejuvenation_green, `FALSE` = aging_red), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.22, 0.22))) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.margin = margin(12, 24, 28, 18)) +
  labs(
    x = NULL,
    y = "Mean integrated rejuvenation score",
    title = "Integrated Hallmark rejuvenation score",
    subtitle = "Combines ssGSEA reversal magnitude with gene-level percent improved"
  )
ggsave(file.path(out_dirs$figures, "hallmark_integrated_rejuvenation_summary.png"), p_integrated, width = 10.5, height = 6.2, dpi = 220, bg = "white")
ggsave(file.path(out_dirs$figures, "hallmark_integrated_rejuvenation_summary.pdf"), p_integrated, width = 10.5, height = 6.2, bg = "white")

p_percent <- gene_rev |>
  mutate(
    intervention = factor(intervention, levels = setdiff(group_levels, control_group)),
    aging_hallmark = factor(aging_hallmark, levels = rev(hallmark_order))
  ) |>
  ggplot(aes(x = intervention, y = aging_hallmark)) +
  geom_point(aes(size = percent_improved, color = fisher_pvalue), alpha = 0.9) +
  scale_color_viridis_c(option = "magma", trans = "reverse", na.value = "grey75") +
  scale_size_continuous(range = c(1.5, 8), limits = c(0, 100)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.minor = element_blank()) +
  labs(
    x = NULL,
    y = NULL,
    size = "% improved",
    color = "Fisher p",
    title = "Gene-level Hallmark reversal percentage",
    subtitle = "Improved means aged_up genes decreased and aged_down genes increased versus control"
  )
ggsave(file.path(out_dirs$figures, "hallmark_gene_percent_improved_dotplot.png"), p_percent, width = 11, height = 7.5, dpi = 220, bg = "white")
ggsave(file.path(out_dirs$figures, "hallmark_gene_percent_improved_dotplot.pdf"), p_percent, width = 11, height = 7.5, bg = "white")

mat <- limma_tbl |>
  mutate(row_id = paste(aging_hallmark, set_mode, sep = " | ")) |>
  select(row_id, intervention, direction_corrected_ssgsea_reversal) |>
  pivot_wider(names_from = intervention, values_from = direction_corrected_ssgsea_reversal) |>
  column_to_rownames("row_id") |>
  as.matrix()
lim <- max(abs(mat), na.rm = TRUE)
if (!is.finite(lim) || lim == 0) lim <- 1
png(file.path(out_dirs$figures, "hallmark_ssgsea_rejuvenation_heatmap.png"), width = 2200, height = 1850, res = 220)
pheatmap(mat, color = colorRampPalette(c(aging_red, neutral_white, rejuvenation_green))(100), breaks = seq(-lim, lim, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE, main = "Direction-corrected ssGSEA reversal", fontsize_row = 8, fontsize_col = 10, border_color = NA)
dev.off()
pdf(file.path(out_dirs$figures, "hallmark_ssgsea_rejuvenation_heatmap.pdf"), width = 9.5, height = 8)
pheatmap(mat, color = colorRampPalette(c(aging_red, neutral_white, rejuvenation_green))(100), breaks = seq(-lim, lim, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE, main = "Direction-corrected ssGSEA reversal", fontsize_row = 8, fontsize_col = 10, border_color = NA)
dev.off()

integrated_mat <- integrated_hallmark |>
  select(aging_hallmark, intervention, integrated_rejuvenation_score) |>
  pivot_wider(names_from = intervention, values_from = integrated_rejuvenation_score) |>
  column_to_rownames("aging_hallmark") |>
  as.matrix()
integrated_mat <- integrated_mat[hallmark_order[hallmark_order %in% rownames(integrated_mat)], setdiff(group_levels, control_group)[setdiff(group_levels, control_group) %in% colnames(integrated_mat)], drop = FALSE]
ilim <- max(abs(integrated_mat), na.rm = TRUE)
if (!is.finite(ilim) || ilim == 0) ilim <- 1
png(file.path(out_dirs$figures, "hallmark_integrated_rejuvenation_heatmap.png"), width = 2200, height = 1550, res = 220)
pheatmap(integrated_mat, color = colorRampPalette(c(aging_red, neutral_white, rejuvenation_green))(100), breaks = seq(-ilim, ilim, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE, main = "Integrated Hallmark rejuvenation score", fontsize_row = 8.5, fontsize_col = 10, border_color = NA)
dev.off()
pdf(file.path(out_dirs$figures, "hallmark_integrated_rejuvenation_heatmap.pdf"), width = 9.5, height = 7)
pheatmap(integrated_mat, color = colorRampPalette(c(aging_red, neutral_white, rejuvenation_green))(100), breaks = seq(-ilim, ilim, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE, main = "Integrated Hallmark rejuvenation score", fontsize_row = 8.5, fontsize_col = 10, border_color = NA)
dev.off()

merged <- tests |> group_by(group) |> slice(1) |> ungroup() |> left_join(hs, by = c("group" = "intervention"))
merged <- merged |> left_join(integrated_intervention, by = c("group" = "intervention"))
table_rows <- paste(apply(merged, 1, function(r) {
  paste0("<tr><td><strong>", html_escape(r[["group"]]), "</strong></td><td class='num'>", fmt_num(as.numeric(r[["mean_difference"]]), 4),
         "</td><td class='num'>", fmt_p(as.numeric(r[["pvalue"]])),
         "</td><td class='num'>", fmt_num(as.numeric(r[["mean_integrated_rejuvenation_score"]]), 4),
         "</td><td class='num'>", fmt_num(as.numeric(r[["mean_direction_corrected_ssgsea_reversal"]]), 4),
         "</td><td class='num'>", fmt_num(as.numeric(r[["mean_gene_combined_rejuvenation_score"]]), 4),
         "</td><td class='num'>", r[["n_rejuvenation_hallmarks"]], "/", r[["n_hallmarks_tested"]], "</td></tr>")
}), collapse = "\n")
html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>tAge + Aging Hallmark Report</title>",
  "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;color:#1d1d1f}section,header{padding:56px 7vw}.hero{background:#f5f5f7}.wrap{max-width:1180px;margin:auto}h1{font-size:48px;line-height:1.08;margin:0 0 16px}h2{font-size:32px;margin:0 0 12px}.lead{font-size:21px;color:#333;max-width:900px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:22px}.fig img{width:100%;border:1px solid #e5e5ea}.caption{font-size:13px;color:#6e6e73}table{border-collapse:collapse;width:100%;margin-top:24px}th,td{padding:10px 12px;border-bottom:1px solid #e5e5ea;text-align:left;font-size:14px}.num{text-align:right;font-variant-numeric:tabular-nums}.pill{display:inline-block;border:1px solid #d2d2d7;border-radius:999px;padding:7px 12px;margin:4px;background:white;font-size:13px}@media(max-width:900px){.grid{grid-template-columns:1fr}h1{font-size:36px}}</style>",
  "</head><body><header class='hero'><div class='wrap'><h1>tAge + Open Genes Hallmark ssGSEA Report</h1>",
  "<div class='lead'>One-line RNA-seq report from raw counts/RSEM expected counts. Lower tAge and positive integrated Hallmark scores indicate rejuvenation-like movement versus control.</div>",
  "<p><span class='pill'>Species: ", html_escape(arg_get(opts, "species", "human")), "</span><span class='pill'>Control: ", html_escape(control_group), "</span><span class='pill'>Model: ", html_escape(arg_get(opts, "model-family", "EN")), " ", html_escape(arg_get(opts, "target", "mortality")), "</span><span class='pill'>Tissue model: ", html_escape(arg_get(opts, "model-tissue", "Multitissue")), "</span><span class='pill'>Preprocess: ", html_escape(arg_get(opts, "preprocess", "scaled_diff")), "</span></p>",
  "</div></header>",
  "<section><div class='wrap'><h2>tAge clock</h2><div class='grid'><div class='fig'><img src='../figures/tAge_boxplot.png'><div class='caption'>Predicted tAge by group.</div></div><div class='fig'><img src='../figures/tAge_delta_vs_control.png'><div class='caption'>Mean tAge difference versus control.</div></div></div></div></section>",
  "<section><div class='wrap'><h2>Integrated Hallmark score</h2><div class='grid'><div class='fig'><img src='../figures/hallmark_integrated_rejuvenation_summary.png'><div class='caption'>Mean integrated Hallmark score. Green means rejuvenation direction; red means aging-like direction.</div></div><div class='fig'><img src='../figures/hallmark_integrated_rejuvenation_heatmap.png'><div class='caption'>Integrated score by Hallmark and intervention. Green means rejuvenation direction; red means aging-like direction.</div></div></div></div></section>",
  "<section><div class='wrap'><h2>Supporting Hallmark evidence</h2><div class='grid'><div class='fig'><img src='../figures/hallmark_ssgsea_rejuvenation_summary.png'><div class='caption'>ssGSEA reversal component.</div></div><div class='fig'><img src='../figures/hallmark_gene_percent_improved_dotplot.png'><div class='caption'>Gene-level percentage component and Fisher enrichment p-value.</div></div></div></div></section>",
  "<section><div class='wrap'><h2>Integrated summary</h2><table><thead><tr><th>Group</th><th class='num'>tAge delta</th><th class='num'>tAge p</th><th class='num'>Integrated score</th><th class='num'>ssGSEA component</th><th class='num'>Gene component</th><th class='num'>Rev Hallmarks</th></tr></thead><tbody>",
  table_rows,
  "</tbody></table><p class='caption'>Integrated score = 0.5*tanh(ssGSEA reversal) + 0.5*gene-level combined score. Gene-level combined score = 0.5*tanh(mean direction-corrected log2FC) + 0.5*(2*fraction improved - 1). Positive values indicate rejuvenation-like reversal.</p></div></section>",
  "<section><div class='wrap'><h2>Methods</h2><p>tAge models are selected from local model files and validated by MD5 checksums. Model inventory follows Zenodo record ",
  html_escape(zenodo_record_url), ". Hallmark scores use DESeq2 VST expression, Open Genes aged_up/down gene sets, GSVA ssGSEA, limma contrasts versus control, and DESeq2 gene-level reversal percentages. Fisher exact tests use detected Open Genes aged-direction entries as the directional background.</p></div></section>",
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
