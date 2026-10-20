#!/usr/bin/env Rscript

# Reproducible PBMC analysis for the Tn5-ModLib manuscript.
#
# This script starts from the eight submitted Cell Ranger filtered-matrix
# archives and their deposited cell-annotation tables. It performs independent
# per-library analysis followed by Harmony integration for joint clustering and
# visualization. The manually reviewed cell types are joined by barcode from
# the annotation tables; integrated clusters are stored in a separate field.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript tn5_pbmc_analysis.R --input-dir DIR --output-dir DIR ",
      "[--save-objects true|false]\n\n",
      "Required:\n",
      "  --input-dir      Directory containing the eight submitted .tar.gz matrices\n",
      "  --output-dir     Directory in which results will be written\n\n",
      "Optional:\n",
      "  --save-objects   Save per-sample and integrated Seurat RDS files (default: false)\n"
    )
  )
}

read_option <- function(name, default = NULL) {
  exact <- match(name, args)
  if (!is.na(exact)) {
    if (exact == length(args)) stop("Missing value after ", name)
    return(args[[exact + 1L]])
  }
  prefix <- paste0(name, "=")
  inline <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(inline)) return(sub(paste0("^", prefix), "", inline[[1L]]))
  default
}

if (any(args %in% c("-h", "--help"))) {
  usage()
  quit(save = "no", status = 0)
}

input_dir_arg <- read_option("--input-dir")
output_dir_arg <- read_option("--output-dir")
save_objects_arg <- tolower(read_option("--save-objects", "false"))

if (is.null(input_dir_arg) || is.null(output_dir_arg)) {
  usage()
  stop("Both --input-dir and --output-dir are required.")
}
if (!save_objects_arg %in% c("true", "false")) {
  stop("--save-objects must be true or false.")
}
save_objects <- identical(save_objects_arg, "true")

# Explicitly preserve the Seurat defaults used by the original notebook.
# These are kept separate because SCTransform, FindClusters, RunPCA, and
# RunUMAP use different defaults.
sct_seed <- 1448145L
pca_seed <- 42L
cluster_seed <- 0L
umap_seed <- 42L

required_packages <- c(
  "Seurat", "SeuratObject", "sctransform", "glmGamPoi", "harmony",
  "dplyr", "ggplot2", "patchwork"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". See README.md for the tested versions."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}
options(future.globals.maxSize = 8 * 1024^3)

input_dir <- normalizePath(input_dir_arg, mustWork = TRUE)
dir.create(output_dir_arg, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir_arg, mustWork = TRUE)

directories <- c(
  annotations = file.path(output_dir, "annotations"),
  figures = file.path(output_dir, "figures"),
  tables = file.path(output_dir, "tables"),
  logs = file.path(output_dir, "logs"),
  objects = file.path(output_dir, "objects")
)
invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))

manifest <- data.frame(
  sample = c(
    "3prime_1_10x", "3prime_1_Tn5", "3prime_2_10x", "3prime_2_Tn5",
    "5prime_1_10x", "5prime_1_Tn5", "5prime_2_10x", "5prime_2_Tn5"
  ),
  file = c(
    "PBMC_3prime-1_10x.tar.gz", "PBMC_3prime-1_Tn5.tar.gz",
    "PBMC_3prime-2_10x.tar.gz", "PBMC_3prime-2_Tn5.tar.gz",
    "PBMC_5prime-1_10x.tar.gz", "PBMC_5prime-1_Tn5.tar.gz",
    "PBMC_5prime-2_10x.tar.gz", "PBMC_5prime-2_Tn5.tar.gz"
  ),
  annotation_file = c(
    "PBMC_3prime-1_10x_cell_annotations.tsv.gz",
    "PBMC_3prime-1_Tn5_cell_annotations.tsv.gz",
    "PBMC_3prime-2_10x_cell_annotations.tsv.gz",
    "PBMC_3prime-2_Tn5_cell_annotations.tsv.gz",
    "PBMC_5prime-1_10x_cell_annotations.tsv.gz",
    "PBMC_5prime-1_Tn5_cell_annotations.tsv.gz",
    "PBMC_5prime-2_10x_cell_annotations.tsv.gz",
    "PBMC_5prime-2_Tn5_cell_annotations.tsv.gz"
  ),
  gigadb_sample_id = c(
    "PBMC3prime1", "PBMC3prime1", "PBMC3prime2", "PBMC3prime2",
    "PBMC5prime1", "PBMC5prime1", "PBMC5prime2", "PBMC5prime2"
  ),
  assay = c("3prime", "3prime", "3prime", "3prime", "5prime", "5prime", "5prime", "5prime"),
  library_method = c("10x", "Tn5-ModLib", "10x", "Tn5-ModLib", "10x", "Tn5-ModLib", "10x", "Tn5-ModLib"),
  stringsAsFactors = FALSE
)

cluster_to_celltype <- list(
  "3prime_1_10x" = c(
    "0" = "Naive CD4+ T", "1" = "Memory CD4+ T", "2" = "CD14+ Mono",
    "3" = "CD8+ T", "4" = "B", "5" = "Naive CD4+ T", "6" = "NK",
    "7" = "CD8+ T", "8" = "Naive CD4+ T", "9" = "Platelet",
    "10" = "FCGR3A+ Mono", "11" = "DC", "12" = "DC"
  ),
  "3prime_1_Tn5" = c(
    "0" = "Memory CD4+ T", "1" = "Naive CD4+ T", "2" = "CD14+ Mono",
    "3" = "CD8+ T", "4" = "B", "5" = "Naive CD4+ T", "6" = "NK",
    "7" = "CD8+ T", "8" = "Naive CD4+ T", "9" = "Naive CD4+ T",
    "10" = "Memory CD4+ T", "11" = "Platelet", "12" = "FCGR3A+ Mono",
    "13" = "DC", "14" = "DC"
  ),
  "3prime_2_10x" = c(
    "0" = "Naive CD4+ T", "1" = "Memory CD4+ T", "2" = "CD14+ Mono",
    "3" = "CD8+ T", "4" = "Naive CD4+ T", "5" = "CD8+ T", "6" = "NK",
    "7" = "Naive CD4+ T", "8" = "Memory CD4+ T", "9" = "B", "10" = "B",
    "11" = "Memory CD4+ T", "12" = "Platelet", "13" = "FCGR3A+ Mono",
    "14" = "DC", "15" = "DC"
  ),
  "3prime_2_Tn5" = c(
    "0" = "Naive CD4+ T", "1" = "Memory CD4+ T", "2" = "CD14+ Mono",
    "3" = "CD8+ T", "4" = "Naive CD4+ T", "5" = "CD8+ T", "6" = "NK",
    "7" = "Naive CD4+ T", "8" = "B", "9" = "Memory CD4+ T", "10" = "B",
    "11" = "Memory CD4+ T", "12" = "Platelet", "13" = "FCGR3A+ Mono",
    "14" = "NK", "15" = "DC"
  ),
  "5prime_1_10x" = c(
    "0" = "CD14+ Mono", "1" = "Naive CD4+ T", "2" = "Memory CD4+ T",
    "3" = "CD8+ T", "4" = "NK", "5" = "Naive CD4+ T", "6" = "Naive CD4+ T",
    "7" = "B", "8" = "B", "9" = "Naive CD4+ T", "10" = "Memory CD4+ T",
    "11" = "Naive CD4+ T", "12" = "CD8+ T", "13" = "Naive CD4+ T",
    "14" = "Platelet", "15" = "Memory CD4+ T", "16" = "FCGR3A+ Mono",
    "17" = "DC"
  ),
  "5prime_1_Tn5" = c(
    "0" = "CD14+ Mono", "1" = "Naive CD4+ T", "2" = "Memory CD4+ T",
    "3" = "CD8+ T", "4" = "NK", "5" = "Naive CD4+ T", "6" = "Naive CD4+ T",
    "7" = "B", "8" = "B", "9" = "Memory CD4+ T", "10" = "Naive CD4+ T",
    "11" = "Naive CD4+ T", "12" = "CD8+ T", "13" = "Naive CD4+ T",
    "14" = "Platelet", "15" = "Naive CD4+ T", "16" = "FCGR3A+ Mono",
    "17" = "DC"
  ),
  "5prime_2_10x" = c(
    "0" = "CD14+ Mono", "1" = "Memory CD4+ T", "2" = "Naive CD4+ T",
    "3" = "CD8+ T", "4" = "Naive CD4+ T", "5" = "NK", "6" = "Naive CD4+ T",
    "7" = "B", "8" = "B", "9" = "Naive CD4+ T", "10" = "Naive CD4+ T",
    "11" = "Naive CD4+ T", "12" = "Platelet", "13" = "CD8+ T",
    "14" = "Memory CD4+ T", "15" = "Naive CD4+ T", "16" = "FCGR3A+ Mono",
    "17" = "Naive CD4+ T", "18" = "DC"
  ),
  "5prime_2_Tn5" = c(
    "0" = "CD14+ Mono", "1" = "Memory CD4+ T", "2" = "Naive CD4+ T",
    "3" = "CD8+ T", "4" = "NK", "5" = "Naive CD4+ T", "6" = "Naive CD4+ T",
    "7" = "B", "8" = "Naive CD4+ T", "9" = "B", "10" = "Naive CD4+ T",
    "11" = "Naive CD4+ T", "12" = "CD8+ T", "13" = "Naive CD4+ T",
    "14" = "Memory CD4+ T", "15" = "Platelet", "16" = "Naive CD4+ T",
    "17" = "FCGR3A+ Mono", "18" = "DC", "19" = "DC"
  )
)

marker_genes <- c(
  "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ",
  "PPBP", "CD8A", "IL7R", "S100A4", "CCR7", "MS4A7", "NKG7", "CST3"
)

write_tsv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt", compression = 9)
  on.exit(close(connection), add = TRUE)
  write.table(
    data, connection, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA"
  )
}

read_submitted_matrix <- function(archive, sample_name) {
  extract_dir <- tempfile(pattern = paste0("tn5_", sample_name, "_"))
  dir.create(extract_dir, recursive = TRUE)
  on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(archive, exdir = extract_dir)

  matrix_files <- list.files(
    extract_dir, pattern = "^matrix[.]mtx([.]gz)?$",
    recursive = TRUE, full.names = TRUE
  )
  if (length(matrix_files) != 1L) {
    stop(
      "Expected one matrix.mtx(.gz) in ", basename(archive),
      "; found ", length(matrix_files)
    )
  }
  matrix_dir <- dirname(matrix_files[[1L]])
  counts <- Read10X(data.dir = matrix_dir, unique.features = TRUE)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      warning("Multiple modalities found for ", sample_name, "; using the first.")
      counts <- counts[[1L]]
    }
  }
  counts
}

read_reviewed_annotations <- function(annotation_path) {
  annotation <- read.delim(
    gzfile(annotation_path),
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  celltype_column <- if ("celltype" %in% names(annotation)) {
    "celltype"
  } else if ("cell_type" %in% names(annotation)) {
    "cell_type"
  } else {
    stop("No celltype or cell_type column in ", basename(annotation_path))
  }
  required <- c("barcode", "passed_seurat_qc", celltype_column)
  missing <- setdiff(required, names(annotation))
  if (length(missing)) {
    stop(
      "Missing annotation columns in ", basename(annotation_path), ": ",
      paste(missing, collapse = ", ")
    )
  }
  if (anyDuplicated(annotation$barcode)) {
    stop("Duplicate barcodes in ", basename(annotation_path))
  }
  retained <- annotation$passed_seurat_qc %in% c(TRUE, "TRUE", "True", "true", 1, "1")
  annotation[retained, , drop = FALSE]
}

attach_celltypes <- function(object, annotation_path, sample_name) {
  annotation <- read_reviewed_annotations(annotation_path)
  celltype_column <- if ("celltype" %in% names(annotation)) "celltype" else "cell_type"
  matched <- match(object$barcode, annotation$barcode)
  if (anyNA(matched)) {
    missing_barcodes <- object$barcode[is.na(matched)]
    stop(
      "Annotation table does not cover all QC-passed cells for ", sample_name,
      ": ", length(missing_barcodes), " missing barcodes: ",
      paste(utils::head(missing_barcodes, 10L), collapse = ", ")
    )
  }
  labels <- annotation[[celltype_column]][matched]
  if (anyNA(labels) || any(labels == "")) {
    stop(
      "Missing celltype values for QC-passed cells in ",
      basename(annotation_path)
    )
  }
  object$celltype <- labels
  object
}

archive_paths <- file.path(input_dir, manifest$file)
annotation_paths <- file.path(input_dir, manifest$annotation_file)
missing_archives <- manifest$file[!file.exists(archive_paths)]
if (length(missing_archives)) {
  stop("Missing input archives: ", paste(missing_archives, collapse = ", "))
}
missing_annotations <- manifest$annotation_file[!file.exists(annotation_paths)]
if (length(missing_annotations)) {
  stop("Missing input annotation files: ", paste(missing_annotations, collapse = ", "))
}

write.table(
  manifest, file.path(directories[["tables"]], "input_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

mapping_table <- do.call(
  rbind,
  lapply(names(cluster_to_celltype), function(sample_name) {
    mapping <- cluster_to_celltype[[sample_name]]
    data.frame(
      sample = sample_name,
      seurat_cluster = names(mapping),
      celltype = unname(mapping),
      stringsAsFactors = FALSE
    )
  })
)
write.table(
  mapping_table, file.path(directories[["tables"]], "reported_cluster_to_celltype_mapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

analyzed_objects <- list()
raw_barcodes <- list()

for (i in seq_len(nrow(manifest))) {
  sample_name <- manifest$sample[[i]]
  archive <- archive_paths[[i]]
  message("[", i, "/", nrow(manifest), "] Reading and analyzing ", sample_name)

  counts <- read_submitted_matrix(archive, sample_name)
  raw_barcodes[[sample_name]] <- colnames(counts)

  object <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  object$barcode <- colnames(object)
  object$sample <- sample_name
  object$gigadb_sample_id <- manifest$gigadb_sample_id[[i]]
  object$assay_type <- manifest$assay[[i]]
  object$library_method <- manifest$library_method[[i]]
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

  threshold_barcodes <- colnames(object)[
    object$nFeature_RNA > 200 & object$nFeature_RNA < 5000 & object$percent.mt < 10
  ]
  reviewed_annotations <- read_reviewed_annotations(annotation_paths[[i]])
  reviewed_barcodes <- reviewed_annotations$barcode
  missing_from_matrix <- setdiff(reviewed_barcodes, colnames(object))
  if (length(missing_from_matrix)) {
    stop(
      "Reviewed QC barcodes are absent from the submitted matrix for ", sample_name,
      ": ", paste(utils::head(missing_from_matrix, 10L), collapse = ", ")
    )
  }
  reviewed_outside_thresholds <- setdiff(reviewed_barcodes, threshold_barcodes)
  if (length(reviewed_outside_thresholds)) {
    warning(
      sample_name, ": ", length(reviewed_outside_thresholds),
      " reviewed manuscript cell(s) fall outside the thresholds when recalculated from the submitted matrix but were retained: ",
      paste(utils::head(reviewed_outside_thresholds, 10L), collapse = ", ")
    )
  }
  threshold_only <- setdiff(threshold_barcodes, reviewed_barcodes)
  if (length(threshold_only)) {
    warning(
      sample_name, ": ", length(threshold_only),
      " threshold-passing cell(s) were not present in the reviewed manuscript cell set and were excluded: ",
      paste(utils::head(threshold_only, 10L), collapse = ", ")
    )
  }
  object <- subset(object, cells = reviewed_barcodes)

  object <- SCTransform(
    object,
    method = "glmGamPoi",
    vars.to.regress = "percent.mt",
    vst.flavor = "v2",
    seed.use = sct_seed,
    verbose = FALSE
  )
  object <- RunPCA(object, npcs = 50, verbose = FALSE, seed.use = pca_seed)
  object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
  object <- FindClusters(
    object, resolution = 0.5, random.seed = cluster_seed, verbose = FALSE
  )
  object <- RunUMAP(
    object, dims = 1:10, seed.use = umap_seed,
    umap.method = "uwot", metric = "cosine", verbose = FALSE
  )
  object <- attach_celltypes(object, annotation_paths[[i]], sample_name)
  Idents(object) <- "celltype"
  object <- RenameCells(object, add.cell.id = sample_name)
  analyzed_objects[[sample_name]] <- object

  cluster_plot <- DimPlot(
    object, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE, raster = FALSE
  ) + ggtitle(paste(sample_name, "clusters"))
  celltype_plot <- DimPlot(
    object, reduction = "umap", group.by = "celltype", label = TRUE,
    repel = TRUE, raster = FALSE
  ) + ggtitle(paste(sample_name, "cell types"))
  ggsave(
    file.path(directories[["figures"]], paste0(sample_name, "_umap.pdf")),
    cluster_plot + celltype_plot, width = 13, height = 6.5, units = "in"
  )

  present_markers <- intersect(marker_genes, rownames(object))
  marker_plot <- DotPlot(
    object, features = present_markers, group.by = "seurat_clusters"
  ) + RotatedAxis() + ggtitle(paste(sample_name, "marker expression"))
  ggsave(
    file.path(directories[["figures"]], paste0(sample_name, "_marker_dotplot.pdf")),
    marker_plot, width = 12, height = 7, units = "in"
  )

  if (save_objects) {
    saveRDS(
      object,
      file.path(directories[["objects"]], paste0(sample_name, "_seurat.rds")),
      compress = "gzip"
    )
  }

  rm(counts, object, cluster_plot, celltype_plot, marker_plot)
  gc()
}

message("Merging QC-passed cells and performing Harmony integration")
integration_inputs <- lapply(analyzed_objects, function(object) {
  DefaultAssay(object) <- "RNA"
  object[["SCT"]] <- NULL
  object@reductions <- list()
  object@graphs <- list()
  object
})

integrated <- integration_inputs[[1L]]
if (length(integration_inputs) > 1L) {
  integrated <- merge(
    x = integrated,
    y = integration_inputs[-1L],
    merge.data = FALSE,
    merge.dr = FALSE
  )
}
DefaultAssay(integrated) <- "RNA"
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

integrated <- NormalizeData(integrated, verbose = FALSE)
integrated <- FindVariableFeatures(
  integrated, selection.method = "vst", nfeatures = 2000, verbose = FALSE
)
integrated <- ScaleData(integrated, vars.to.regress = "percent.mt", verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE, seed.use = pca_seed)
integrated[["RNA"]] <- split(integrated[["RNA"]], f = integrated$sample)
integrated <- IntegrateLayers(
  object = integrated,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  verbose = FALSE
)
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

integrated <- FindNeighbors(
  integrated,
  reduction = "integrated.harmony",
  dims = 1:24,
  graph.name = c("harmony_nn", "harmony_snn"),
  verbose = FALSE
)
integrated <- FindClusters(
  integrated,
  graph.name = "harmony_snn",
  resolution = 0.2,
  cluster.name = "integrated_cluster",
  random.seed = cluster_seed,
  verbose = FALSE
)
integrated <- RunUMAP(
  integrated,
  reduction = "integrated.harmony",
  dims = 1:24,
  reduction.name = "umap.harmony",
  seed.use = umap_seed,
  umap.method = "uwot",
  metric = "cosine",
  verbose = FALSE
)

n_integrated_clusters <- length(unique(as.character(integrated$integrated_cluster)))

integrated_cluster_plot <- DimPlot(
  integrated,
  reduction = "umap.harmony",
  group.by = "integrated_cluster",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + ggtitle("Harmony-integrated clusters")
ggsave(
  file.path(directories[["figures"]], "integrated_umap_clusters.pdf"),
  integrated_cluster_plot, width = 8, height = 7, units = "in"
)

integrated_sample_plot <- DimPlot(
  integrated,
  reduction = "umap.harmony",
  group.by = "sample",
  raster = FALSE
) + ggtitle("Harmony-integrated cells by library")
ggsave(
  file.path(directories[["figures"]], "integrated_umap_samples.pdf"),
  integrated_sample_plot, width = 9, height = 7, units = "in"
)

integrated_split_plot <- DimPlot(
  integrated,
  reduction = "umap.harmony",
  group.by = "integrated_cluster",
  split.by = "sample",
  ncol = 4,
  raster = FALSE
)
ggsave(
  file.path(directories[["figures"]], "integrated_umap_split_by_sample.pdf"),
  integrated_split_plot, width = 18, height = 11, units = "in"
)

present_markers <- intersect(marker_genes, rownames(integrated))
integrated_marker_plot <- DotPlot(
  integrated,
  features = present_markers,
  group.by = "integrated_cluster",
  assay = "RNA"
) + RotatedAxis() + ggtitle("Marker expression by integrated cluster")
ggsave(
  file.path(directories[["figures"]], "integrated_marker_dotplot.pdf"),
  integrated_marker_plot, width = 13, height = 8, units = "in"
)

metadata <- integrated[[]]
metadata$cell_id <- rownames(metadata)

celltype_counts <- as.data.frame(
  table(sample = metadata$sample, celltype = metadata$celltype),
  stringsAsFactors = FALSE
)
names(celltype_counts)[[3L]] <- "cell_count"
celltype_counts <- celltype_counts |>
  group_by(sample) |>
  mutate(percentage = 100 * cell_count / sum(cell_count)) |>
  ungroup()
write.table(
  celltype_counts,
  file.path(directories[["tables"]], "celltype_counts_and_percentages.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

composition_plot <- ggplot(
  celltype_counts,
  aes(x = sample, y = percentage, fill = celltype)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  labs(x = NULL, y = "Cells (%)", fill = "Cell type") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  file.path(directories[["figures"]], "celltype_composition.pdf"),
  composition_plot, width = 11, height = 6.5, units = "in"
)

for (i in seq_len(nrow(manifest))) {
  sample_name <- manifest$sample[[i]]
  sample_metadata <- metadata[metadata$sample == sample_name, , drop = FALSE]
  analyzed_index <- match(raw_barcodes[[sample_name]], sample_metadata$barcode)

  annotation <- data.frame(
    barcode = raw_barcodes[[sample_name]],
    sample = sample_name,
    gigadb_sample_id = manifest$gigadb_sample_id[[i]],
    assay = manifest$assay[[i]],
    library_method = manifest$library_method[[i]],
    passed_seurat_qc = !is.na(analyzed_index),
    cell_id = sample_metadata$cell_id[analyzed_index],
    seurat_cluster = as.character(sample_metadata$seurat_clusters[analyzed_index]),
    integrated_cluster = as.character(sample_metadata$integrated_cluster[analyzed_index]),
    celltype = as.character(sample_metadata$celltype[analyzed_index]),
    stringsAsFactors = FALSE
  )

  write_tsv_gz(
    annotation,
    file.path(directories[["annotations"]], manifest$annotation_file[[i]])
  )
}

if (save_objects) {
  saveRDS(
    integrated,
    file.path(directories[["objects"]], "PBMC_8libraries_harmony_integrated_seurat.rds"),
    compress = "gzip"
  )
}

software_versions <- data.frame(
  software = c("R", required_packages),
  version = c(
    as.character(getRversion()),
    vapply(
      required_packages,
      function(package) as.character(utils::packageVersion(package)),
      character(1)
    )
  ),
  stringsAsFactors = FALSE
)
write.table(
  software_versions,
  file.path(directories[["logs"]], "software_versions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  c(
    paste0("Completed: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("SCTransform seed: ", sct_seed),
    paste0("PCA seed: ", pca_seed),
    paste0("Clustering seed: ", cluster_seed),
    paste0("UMAP seed: ", umap_seed),
    paste0("QC-passed cells: ", ncol(integrated)),
    paste0("Integrated clusters: ", n_integrated_clusters),
    paste0("Objects saved: ", save_objects)
  ),
  file.path(directories[["logs"]], "run_summary.txt")
)

message("Analysis completed successfully.")
