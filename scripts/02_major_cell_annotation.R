# =============================================================================
# 02  Major cell type annotation
#
# Harmony integration, clustering, marker detection and manual annotation of
# the major cell types.
#
# Input : results/objects/scobj_merged.rds
# Output: results/objects/scobj_annotated.rds, UMAP figures, marker tables
# =============================================================================

source("R/config.R")
source("R/preprocessing.R")
source("R/plotting.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

scobj <- readRDS(file.path(object_dir, "scobj_merged.rds"))

# -----------------------------------------------------------------------------
# 1. Normalisation, Harmony integration, clustering and UMAP
# -----------------------------------------------------------------------------
scobj <- ProcessSingleCellData(scobj, group.by.vars = "orig.ident",
                               resolution.seq = seq(0.1, 0.6, 0.1))

# -----------------------------------------------------------------------------
# 2. Marker genes for every resolution (used for manual annotation)
# -----------------------------------------------------------------------------
FindClusterMarkersParallel(scobj, output_prefix = "all_", output_dir = marker_dir,
                           workers = n_workers)

# Optional: cluster stability across resolutions
# clustree::clustree(scobj, prefix = "RNA_snn_res.")

# -----------------------------------------------------------------------------
# 3. Manual annotation at resolution 0.2
# NOTE: cluster IDs refer to the original run. Re-check them against
# results/markers/all_RNA_snn_res.0.2.csv whenever the clustering is re-run.
# -----------------------------------------------------------------------------
Idents(scobj) <- "RNA_snn_res.0.2"
scobj <- RenameIdents(scobj,
                      "0"  = "T",
                      "1"  = "CAF",
                      "2"  = "Epi",
                      "3"  = "Mac",
                      "4"  = "End",
                      "5"  = "CAF",
                      "6"  = "Proliferating_cluster",
                      "7"  = "Epi",
                      "8"  = "CAF",
                      "9"  = "B",
                      "10" = "mCAF",
                      "11" = "Removed",
                      "12" = "Mac",
                      "13" = "B",
                      "14" = "Removed")

# Clusters 11 and 14 could not be assigned and are removed
scobj <- subset(scobj, idents = "Removed", invert = TRUE)
scobj$main_cluster <- droplevels(Idents(scobj))
print(table(scobj$main_cluster, scobj$group))

saveRDS(scobj, file.path(object_dir, "scobj_annotated.rds"))

# -----------------------------------------------------------------------------
# 4. UMAP figures
# -----------------------------------------------------------------------------
p <- plot_umap_labeled(scobj, "main_cluster", pal_major, split_by = "group",
                       title = "UMAP plot by Group", legend_text_size = 16)
save_pdf(p, file.path(figure_dir, "umap_allcluster_split_group.pdf"), width = 15, height = 8)

p <- plot_umap_labeled(scobj, "main_cluster", pal_major,
                       title = "UMAP plot", legend_text_size = 18)
save_pdf(p, file.path(figure_dir, "umap_allcluster.pdf"), width = 10, height = 6)
