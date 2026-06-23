# tAge Hallmarks RNA-seq Workflow

Reusable one-line R workflow for transcriptomic aging analysis from raw RNA-seq counts or RSEM `*.genes.results` files.

The workflow runs:

1. Gladyshev `tAge` prediction.
2. DESeq2 VST normalization from raw/RSEM expected counts.
3. Open Genes aged-up/down Hallmarks of Aging ssGSEA.
4. limma intervention-vs-control tests on ssGSEA scores.
5. Static HTML and PDF reports.

## Requirements

This repository contains workflow code only. It does not include raw data, trained clock models, or generated results.

Install the R/Python runtime with pixi:

```bash
pixi install
pixi run Rscript scripts/install_requirements.R --check-only
```

Or create the same runtime with mamba or conda:

```bash
mamba env create -f environment.yml
mamba activate tage-hallmarks-workflow
Rscript scripts/install_requirements.R --check-only
```

The equivalent conda commands are:

```bash
conda env create -f environment.yml
conda activate tage-hallmarks-workflow
Rscript scripts/install_requirements.R --check-only
```

The workflow expects:

- the `tAge` R package to be installed in the runtime environment
- local tAge model files under `models/`
- a model registry at `config/validation_models.csv`

Missing tAge models are not downloaded automatically in v1.

## Sources

- Gladyshev tAge code: [Gladyshev-Lab/tAge](https://github.com/Gladyshev-Lab/tAge)
- Gladyshev tAge publication: [https://doi.org/10.1038/s41586-026-10542-3](https://doi.org/10.1038/s41586-026-10542-3)
- Open Genes genes: [Open Genes genes](https://open-genes.com/genes)
- Open Genes API documentation: [Open Genes API docs](https://open-genes.com/api/docs)
- Open Genes publication: [https://doi.org/10.1093/nar/gkad712](https://doi.org/10.1093/nar/gkad712)
- Hallmarks of Aging framework: [https://doi.org/10.1016/j.cell.2022.11.001](https://doi.org/10.1016/j.cell.2022.11.001)
- Reference data used for aged-up/down direction assignment: [https://doi.org/10.18632/aging.204896](https://doi.org/10.18632/aging.204896)

## Included Gene Sets

This repository includes an Ensembl GMT file:

```text
gene_sets/open_genes_hallmark_aged_up_down_ensembl.gmt
```

The gene sets were built from Open Genes genes assigned to Hallmarks of Aging categories following the Hallmarks framework described in [https://doi.org/10.1016/j.cell.2022.11.001](https://doi.org/10.1016/j.cell.2022.11.001). Genes were then split into direction-specific `aged_up` and `aged_down` sets using a reference differential-expression contrast from the dataset associated with [https://doi.org/10.18632/aging.204896](https://doi.org/10.18632/aging.204896):

- reference contrast: `Quiescence` versus `Senescence`
- `aged_log2FC = -1 * log2FoldChange`
- `aged_up`: Open Genes Hallmark genes with `aged_log2FC > 0`
- `aged_down`: Open Genes Hallmark genes with `aged_log2FC < 0`
- no adjusted-p-value filter was applied for the direction split; fold-change sign alone defines direction
- genes without a mapped Ensembl ID or without a nonzero reference fold-change direction were not included in direction-specific sets

The resulting GMT contains 22 gene sets: 11 Open Genes Hallmark categories, each represented as `aged_up` and `aged_down`. `Dysbiosis` is not included because it does not have a direct Open Genes Hallmark gene mapping in this workflow.

The companion summary table is:

```text
gene_sets/open_genes_hallmark_aged_up_down_gene_set_summary.csv
```

## Quick Start

Run from a project directory containing this repository's `scripts/`, local tAge `models/`, and `config/validation_models.csv`:

```bash
pixi run Rscript scripts/run_tage_hallmarks_report.R \
  --rsem-dir /path/to/rsem_gene_results \
  --metadata /path/to/metadata.csv \
  --sample-id-col sample_id \
  --group-col condition \
  --control-group Control \
  --species human \
  --gene-mapping-type Ensembl \
  --model-family EN \
  --target mortality \
  --model-species Multispecies \
  --model-tissue Multitissue \
  --preprocess scaled_diff \
  --out-dir results/example_tage_hallmarks
```

If using conda or mamba, activate the environment first and call `Rscript` directly:

```bash
mamba activate tage-hallmarks-workflow

Rscript scripts/run_tage_hallmarks_report.R \
  --rsem-dir /path/to/rsem_gene_results \
  --metadata /path/to/metadata.csv \
  --sample-id-col sample_id \
  --group-col condition \
  --control-group Control \
  --species human \
  --gene-mapping-type Ensembl \
  --model-family EN \
  --target mortality \
  --model-species Multispecies \
  --model-tissue Multitissue \
  --preprocess scaled_diff \
  --out-dir results/example_tage_hallmarks
```

For a raw count matrix instead of RSEM files:

```bash
pixi run Rscript scripts/run_tage_hallmarks_report.R \
  --counts /path/to/raw_counts.csv \
  --metadata /path/to/metadata.csv \
  --sample-id-col sample_id \
  --group-col condition \
  --control-group Control \
  --species human \
  --gene-mapping-type Ensembl \
  --model-family EN \
  --target mortality \
  --model-species Multispecies \
  --model-tissue Multitissue \
  --preprocess scaled_diff \
  --out-dir results/example_tage_hallmarks
```

## Input Formats

### RSEM Directory

Use `--rsem-dir` for a directory containing files named:

```text
<sample_id>.genes.results
```

Each file must contain:

- `gene_id`
- `expected_count`

### Count Matrix

Use `--counts` for a CSV matrix:

- first column: gene ID or gene symbol
- remaining columns: sample IDs
- values: raw counts or RSEM expected counts

### Metadata

Required columns:

- sample ID column passed by `--sample-id-col`
- group/condition column passed by `--group-col`

Optional:

- tissue column passed by `--tissue-col`
- subset value passed by `--tissue-value`

## Species and Model Options

Accepted analysis species:

```text
human, mouse, rat, monkey
```

tAge model options:

```text
--model-family EN|BR|all
--target mortality|chronoage|normalizedage|all
--model-species Multispecies|Mouse|Rodents|auto
--model-tissue Multitissue|Liver|Brain|Kidney|Skeletal_muscle|auto
--preprocess scaled_diff|yugene_diff|all
```

Defaults:

```text
EN, mortality, Multispecies, Multitissue, scaled_diff
```

If `--model-tissue auto` is used, the workflow attempts to match `--tissue-value` to a locally available model tissue and otherwise falls back to `Multitissue`.

## Output Layout

The workflow writes under `--out-dir`:

```text
prepared/
qc/
tage_predictions/
hallmark_results/
figures/
report/
logs/
```

Main outputs:

- `report/tage_hallmarks_report.html`
- `report/tage_hallmarks_report.pdf`
- `tage_predictions/tAge_vs_control_tests.csv`
- `hallmark_results/open_genes_aged_direction_ssgsea_limma_vs_control.csv`
- `hallmark_results/open_genes_aged_direction_ssgsea_intervention_rejuvenation_summary.csv`

## Interpretation

For Hallmark ssGSEA:

- `aged_up`: rejuvenation direction is lower score in intervention versus control.
- `aged_down`: rejuvenation direction is higher score in intervention versus control.
- `direction_corrected_ssgsea_reversal > 0` means the intervention moved opposite to the aging direction.

Raw counts are not used directly for pathway scoring. The workflow uses DESeq2 VST-normalized expression for GSVA/ssGSEA.
