# tAge Hallmarks RNA-seq Workflow

Reusable R workflow for transcriptomic aging analysis from raw RNA-seq counts or RSEM `*.genes.results` files.

It runs:

- Gladyshev [`tAge`](https://github.com/Gladyshev-Lab/tAge) prediction, following the Gladyshev tAge publication ([Nature, 2026](https://doi.org/10.1038/s41586-026-10542-3))
- DESeq2 VST normalization
- [Open Genes](https://open-genes.com/genes) Hallmarks of Aging aged-up/down ssGSEA, using gene metadata available through the [Open Genes API](https://open-genes.com/api/docs) and Open Genes publication ([NAR, 2024](https://doi.org/10.1093/nar/gkad712))
- limma intervention-vs-control tests
- HTML and PDF report generation

## Requirements

This repository contains workflow code and Open Genes Hallmark gene sets. It requires local Gladyshev tAge model files under `models/` and a model registry at `config/validation_models.csv`.

PDF report rendering uses headless Chrome with `google-chrome --headless --print-to-pdf`.

Prepare local tAge model files:

```bash
mkdir -p models config
cp config/validation_models.example.csv config/validation_models.csv
md5sum models/*.pkl
```

The example registry lists the expected Gladyshev tAge model filenames and checksums for the default multispecies multitissue models. If you use a different model set, edit `config/validation_models.csv` so each `*_file` entry matches a file in `models/` and each `*_md5` entry matches the corresponding checksum. The workflow uses this registry to select models and validate file integrity before running tAge.

Install with pixi:

```bash
pixi install
pixi run install-tage
pixi run check-requirements
```

Or with mamba:

```bash
mamba env create -f environment.yml
mamba activate tage-hallmarks-workflow
Rscript -e 'remotes::install_github("Gladyshev-Lab/tAge", upgrade = "never", dependencies = c("Depends", "Imports", "LinkingTo"), build_vignettes = FALSE)'
Rscript -e 'pkgs <- c("data.table", "dplyr", "tidyr", "tibble", "stringr", "ggplot2", "pheatmap", "reticulate", "remotes", "DESeq2", "GSVA", "limma", "Biobase", "edgeR", "tAge"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", ")); message("All required R packages are available.")'
```

For conda, replace `mamba` with `conda` in the commands above.

## Quick Start

RSEM input:

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

Raw count matrix input:

```bash
pixi run Rscript scripts/run_tage_hallmarks_report.R \
  --counts /path/to/raw_counts.csv \
  --metadata /path/to/metadata.csv \
  --sample-id-col sample_id \
  --group-col condition \
  --control-group Control \
  --species human \
  --gene-mapping-type Ensembl \
  --out-dir results/example_tage_hallmarks
```

## Inputs

- RSEM directory: files named `<sample_id>.genes.results` containing `gene_id` and `expected_count`
- Count matrix: CSV with genes in the first column and samples in the remaining columns
- Metadata: CSV containing the sample ID column and group/condition column

Accepted analysis species: `human`, `mouse`, `rat`, `monkey`.

Main model options:

```text
--model-family EN|BR|all
--target mortality|chronoage|normalizedage|all
--model-species Multispecies|Mouse|Rodents|auto
--model-tissue Multitissue|Liver|Brain|Kidney|Skeletal_muscle|auto
--preprocess scaled_diff|yugene_diff|all
```

Defaults: `EN`, `mortality`, `Multispecies`, `Multitissue`, `scaled_diff`.

## Gene Sets

Included files:

```text
gene_sets/open_genes_hallmark_aged_up_down_ensembl.gmt
gene_sets/open_genes_hallmark_aged_up_down_gene_set_summary.csv
```

The GMT was built from [Open Genes](https://open-genes.com/genes) genes assigned to Hallmarks of Aging categories following the [Hallmarks framework](https://doi.org/10.1016/j.cell.2022.11.001). Genes were split into `aged_up` and `aged_down` using a reference `Quiescence` versus `Senescence` differential-expression contrast from the dataset associated with [Aging, 2022](https://doi.org/10.18632/aging.204896):

- `aged_log2FC = -1 * log2FoldChange`
- `aged_up`: `aged_log2FC > 0`
- `aged_down`: `aged_log2FC < 0`
- no adjusted-p-value filter was applied for this direction split

`Dysbiosis` is excluded because host bulk RNA-seq alone cannot directly measure dysbiosis; microbiome or related microbial-composition data would be required.

## Outputs

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

Key files:

- `report/tage_hallmarks_report.html`
- `report/tage_hallmarks_report.pdf`
- `tage_predictions/tAge_vs_control_tests.csv`
- `hallmark_results/open_genes_aged_direction_ssgsea_limma_vs_control.csv`
- `hallmark_results/open_genes_aged_direction_ssgsea_intervention_rejuvenation_summary.csv`

## Interpretation

For Hallmark ssGSEA, `direction_corrected_ssgsea_reversal > 0` means the intervention moved opposite to the aging direction:

- `aged_up`: lower intervention score than control is interpreted as rejuvenation direction
- `aged_down`: higher intervention score than control is interpreted as rejuvenation direction

Raw counts are not scored directly. The workflow uses DESeq2 VST-normalized expression for GSVA/ssGSEA.
