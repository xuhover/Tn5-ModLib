# Tn5-ModLib PBMC single-cell RNA-seq analysis

This directory contains a reproducible R workflow for the eight PBMC libraries used for the downstream single-cell analysis in the Tn5-ModLib manuscript. The workflow starts from the processed Cell Ranger filtered count matrices and cell-level annotation tables deposited in GigaDB.

The script deliberately separates two concepts that were mixed in the original working notebook:

1. `celltype` contains the marker-based manual annotations assigned after each library was processed and clustered independently. The reviewed labels are imported from the deposited annotation tables by barcode so that the manuscript table and composition figure can be reproduced exactly even if a future Seurat version changes cluster numbering.
2. `integrated_cluster` is generated after all eight libraries are merged and integrated with Harmony. It is used for the joint UMAP and technical assessment of integration, but it does not overwrite the independently assigned `celltype`.

Only `celltype` is exported as the biological cell-type annotation. The original notebook's `celltype_fine` field is not used.

## Input files

Place these eight GigaDB files in one directory without renaming them:

- `PBMC_3prime-1_10x.tar.gz`
- `PBMC_3prime-1_Tn5.tar.gz`
- `PBMC_3prime-2_10x.tar.gz`
- `PBMC_3prime-2_Tn5.tar.gz`
- `PBMC_5prime-1_10x.tar.gz`
- `PBMC_5prime-1_Tn5.tar.gz`
- `PBMC_5prime-2_10x.tar.gz`
- `PBMC_5prime-2_Tn5.tar.gz`

Place the eight corresponding annotation files in the same directory:

- `PBMC_3prime-1_10x_cell_annotations.tsv.gz`
- `PBMC_3prime-1_Tn5_cell_annotations.tsv.gz`
- `PBMC_3prime-2_10x_cell_annotations.tsv.gz`
- `PBMC_3prime-2_Tn5_cell_annotations.tsv.gz`
- `PBMC_5prime-1_10x_cell_annotations.tsv.gz`
- `PBMC_5prime-1_Tn5_cell_annotations.tsv.gz`
- `PBMC_5prime-2_10x_cell_annotations.tsv.gz`
- `PBMC_5prime-2_Tn5_cell_annotations.tsv.gz`

A validated local copy of these eight reviewed tables is included in `input_cell_annotations/`. For a local run, copy or symlink them into the same input directory as the eight matrix archives.

Each archive must contain `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`. The script extracts each archive into a temporary directory and does not alter the deposited files.

## Software environment

The cleaned workflow was prepared against the following environment recorded on the original analysis server:

| Software | Version |
|---|---:|
| R | 4.5.2 |
| Seurat | 5.5.0 |
| SeuratObject | 5.4.0 |
| sctransform | 0.4.3 |
| glmGamPoi | 1.22.0 |
| harmony | 2.0.2 |
| dplyr | 1.1.4 |
| ggplot2 | 4.0.3 |
| patchwork | 1.3.2 |

Exact cluster numbering can change between package versions and even when random-number state differs. Therefore, reviewed cell types are joined from the deposited annotation files rather than being reassigned from newly generated cluster numbers. Marker plots should still be inspected when reproducing the analysis.

## Run the complete workflow

```bash
Rscript tn5_pbmc_analysis.R \
  --input-dir /path/to/gigadb_matrices \
  --output-dir /path/to/tn5_pbmc_results \
  --save-objects false
```

Set `--save-objects true` to save the eight independently analyzed Seurat objects and the Harmony-integrated object as RDS files. These objects are large and are not required to generate the tabular annotations and figures.

The script uses sequential execution for reproducibility when the `future` package is available. A machine with sufficient RAM for approximately 70,000 cells is required; the exact requirement depends on the sparse matrices and installed Seurat version.

## Analysis parameters

For each library independently, the script performs:

- Seurat object creation from the Cell Ranger filtered matrix;
- mitochondrial fraction calculation using genes beginning with `MT-`;
- recalculation of the original filters (`nFeature_RNA > 200`, `nFeature_RNA < 5000`, and `percent.mt < 10`) and verification of every reviewed QC-passed cell against those thresholds;
- selection of the audited `passed_seurat_qc = TRUE` barcode set from the deposited annotation table, so that the exact manuscript cell set is retained even when a matrix-format or software-version boundary changes by one cell;
- `SCTransform` with `method = "glmGamPoi"`, `vst.flavor = "v2"`, and mitochondrial fraction regressed out;
- PCA with 50 components;
- neighbor graph construction using PCs 1–30;
- graph clustering at resolution 0.5;
- UMAP using PCs 1–10;
- cluster-to-cell-type mapping based on the marker-gene inspection recorded in the original analysis notebook.

For the joint analysis, QC-passed cells are merged and processed with RNA normalization, 2,000 variable features, scaling with mitochondrial fraction regression, 50-component PCA, and Harmony integration. PCs 1–24 of the Harmony embedding are used for the neighbor graph, clustering at resolution 0.2, and the integrated UMAP.

## Output structure

```text
tn5_pbmc_results/
├── annotations/
│   └── eight *_cell_annotations.tsv.gz files
├── figures/
│   ├── per-library cluster, cell-type, and marker plots
│   ├── integrated UMAP plots
│   ├── integrated marker dot plot
│   └── cell-type composition plot
├── tables/
│   ├── input_manifest.tsv
│   ├── reported_cluster_to_celltype_mapping.tsv
│   └── celltype_counts_and_percentages.tsv
├── logs/
│   ├── run_summary.txt
│   └── software_versions.tsv
└── objects/
    └── optional RDS files when --save-objects true
```

Each annotation file contains every barcode from its corresponding submitted matrix. Cells retained after Seurat QC have `passed_seurat_qc = TRUE` and receive `celltype`, `seurat_cluster`, and `integrated_cluster` values. Cells not retained for downstream analysis remain in the table with `passed_seurat_qc = FALSE` and `NA` analysis fields.

## Cell-type definitions

The exported `celltype` values are:

- B
- CD14+ Mono
- CD8+ T
- DC
- FCGR3A+ Mono
- Memory CD4+ T
- NK
- Naive CD4+ T
- Platelet

The original cluster-to-cell-type mappings are written to `tables/reported_cluster_to_celltype_mapping.tsv` for auditability. Because cluster numbering is version-sensitive, the script does not automatically apply these historical mappings to newly calculated clusters; it imports the reviewed per-cell labels from the deposited annotation files instead.

## QC boundary audit

Validation identified five cells whose QC status changes when the metrics are recalculated from the submitted matrix archives rather than the original Cell Ranger HDF5 files. All are exact boundary cases:

| Sample | Barcode | Original HDF5 | Submitted archive | Reviewed action |
|---|---|---:|---:|---|
| `3prime_1_10x` | `GCAACCGCATGCACTA-1` | 5,002 features | 4,999 features | Exclude |
| `3prime_2_10x` | `CATCCGTCAGATTCGT-1` | 5,000 features | 4,999 features | Exclude |
| `3prime_2_Tn5` | `GGTGTCGCACGCGCTA-1` | 5,001 features | 4,999 features | Exclude |
| `5prime_1_10x` | `GCTCCTACATTCACTT-1` | 5,001 features | 4,999 features | Exclude |
| `5prime_2_10x` | `CTAGTGACAGTTAACC-1` | 9.991% mitochondrial reads | 10.000% mitochondrial reads | Include |

The reviewed action preserves the exact cells analyzed for the manuscript under the stated strict filters (`nFeature_RNA < 5000` and `percent.mt < 10`). The script reports these barcodes as warnings so that this decision is visible rather than silently changing the analysis population.

## Reproducibility notes

- Random seeds are explicit and recorded in `logs/run_summary.txt`. To preserve the original settings recorded in the Seurat command log, SCTransform uses seed 1448145, PCA and UMAP use seed 42, and graph clustering uses seed 0. No additional global seed is set because doing so changed the original SCTransform/glmGamPoi clustering result.
- The R and package versions are written to `logs/software_versions.tsv`. Hostnames, usernames, platform paths, environment variables, and locale settings are not recorded.
- Cell-type percentages are calculated directly from the analyzed objects; no percentages are manually entered.
- Marker-based manual cell types are joined by sample and barcode from the deposited annotation tables, while clustering and integration are recalculated from the matrices.
- In the original exploratory notebook, a 16-cluster field was created from a PCA neighbor graph on the merged object, while Harmony was used for the joint UMAP. The cleaned workflow keeps these concepts distinct: `integrated_cluster` is recalculated from the Harmony embedding (14 clusters in the tested environment), and the manuscript cell-type composition is reproduced from the reviewed per-sample `celltype` labels.
- The script stops if any reviewed QC-passed barcode is absent from its matrix. It prints explicit barcode-level warnings for cells that change sides of a QC threshold when recalculated from the submitted matrix: reviewed manuscript cells remain included, and additional threshold-passing cells not present in the reviewed manuscript set remain excluded.
- Server-specific absolute paths, interactive object reassignment, and unrelated exploratory code from the original notebook have been removed.
- The script does not record input/output absolute paths, login details, hostnames, or user identifiers in its output logs.
- The script never writes to the input directory.

Before public release, the generated marker plots, cell counts, and Table 3 percentages should be compared with the final accepted manuscript. If a cluster changes because a different package version was used, its marker expression must be reviewed rather than silently reusing an old cluster number.
