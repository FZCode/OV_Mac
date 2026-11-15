# =============================================================================
# Pre-processing helpers
# =============================================================================

#' Standard Seurat workflow with Harmony integration
#'
#' Normalization -> highly variable genes -> scaling -> PCA -> Harmony ->
#' neighbour graph -> clustering at several resolutions -> UMAP.
ProcessSingleCellData <- function(scobj,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000,
                                  selection.method = "vst",
                                  nfeatures = 2000,
                                  dims.use = 1:20,
                                  resolution.seq = seq(0.1, 0.6, 0.2),
                                  group.by.vars = "orig.ident") {
  scobj <- NormalizeData(scobj, normalization.method = normalization.method,
                         scale.factor = scale.factor)
  scobj <- FindVariableFeatures(scobj, selection.method = selection.method, nfeatures = nfeatures)
  scobj <- ScaleData(scobj)
  scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj))
  scobj <- harmony::RunHarmony(scobj, reduction = "pca",
                               group.by.vars = group.by.vars,
                               reduction.save = "harmony",
                               dims.use = dims.use)
  scobj <- FindNeighbors(scobj, reduction = "harmony", dims = dims.use)
  scobj <- FindClusters(scobj, resolution = resolution.seq)
  scobj <- RunUMAP(scobj, reduction = "harmony", dims = dims.use, reduction.name = "umap")
  scobj
}

#' Find positive markers for several clustering resolutions in parallel
#'
#' Large clusters are down-sampled to `rare_threshold` cells before
#' FindAllMarkers. One CSV is written per resolution:
#' <output_dir>/<output_prefix><resolution column>.csv
FindClusterMarkersParallel <- function(scobj,
                                       output_prefix = "marker_",
                                       output_dir = "Findmarker",
                                       resolution_columns = paste0("RNA_snn_res.0.", 1:8),
                                       rare_threshold = 5000,
                                       logfc_threshold = 0.25,
                                       min_pct = 0.2,
                                       workers = 8,
                                       seed = 2025) {
  suppressPackageStartupMessages({
    library(foreach)
    library(doParallel)
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  valid_columns <- intersect(resolution_columns, colnames(scobj@meta.data))
  if (length(valid_columns) == 0) stop("None of the resolution columns exist in meta.data")
  if (length(valid_columns) < length(resolution_columns)) {
    warning("Missing resolution columns: ",
            paste(setdiff(resolution_columns, valid_columns), collapse = ", "))
  }
  message("Finding markers for ", length(valid_columns), " resolution(s)...")

  cl <- makeCluster(workers)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  results <- foreach(i = seq_along(valid_columns), .combine = c, .packages = "Seurat",
                     .export = c("rare_threshold", "logfc_threshold", "min_pct",
                                 "output_dir", "output_prefix", "seed")) %dopar% {
    column <- valid_columns[i]
    tryCatch({
      set.seed(seed + i)
      Idents(scobj) <- scobj@meta.data[[column]]
      # Down-sample large clusters to speed up FindAllMarkers
      keep_cells <- unlist(lapply(unique(Idents(scobj)), function(cluster) {
        cells <- WhichCells(scobj, idents = cluster)
        if (length(cells) <= rare_threshold) cells else sample(cells, rare_threshold)
      }))
      sub_obj <- subset(scobj, cells = keep_cells)
      Idents(sub_obj) <- column
      markers <- FindAllMarkers(sub_obj, logfc.threshold = logfc_threshold,
                                min.pct = min_pct, only.pos = TRUE, verbose = FALSE)
      write.csv(markers, file.path(output_dir, paste0(output_prefix, column, ".csv")),
                row.names = FALSE)
      paste("Done:", column)
    }, error = function(e) {
      paste("Failed:", column, "-", e$message)
    })
  }

  cat(results, sep = "\n")
  message("Marker tables saved to: ", output_dir)
  invisible(results)
}
