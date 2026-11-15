# =============================================================================
# 03  Epithelial (tumour) cell analysis
#
# Sub-clustering of epithelial cells, cholesterol (steroid) biosynthesis
# activity with scMetabolism, and PCSK9 expression.
#
# Input : results/objects/scobj_annotated.rds
# Output: results/objects/Epi_OV2.rds, figures, marker tables
# =============================================================================

source("R/config.R")
source("R/preprocessing.R")
source("R/plotting.R")
source("R/scmetabolism_patch.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(scMetabolism)
})

scobj <- readRDS(file.path(object_dir, "scobj_annotated.rds"))

# -----------------------------------------------------------------------------
# 1. Sub-cluster epithelial cells
# -----------------------------------------------------------------------------
epi <- subset(scobj, idents = "Epi")
rm(scobj); gc()
epi <- ProcessSingleCellData(epi, group.by.vars = "orig.ident", resolution.seq = seq(0.1, 0.6, 0.1))
FindClusterMarkersParallel(epi, output_prefix = "Epi_", output_dir = marker_dir, workers = n_workers,
                           resolution_columns = paste0("RNA_snn_res.0.", 1:4))

# Remove contaminating clusters (resolution 0.1: 5, 6, 7; check marker tables after re-running)
epi <- epi[, !epi$RNA_snn_res.0.1 %in% c("5", "6", "7")]

# Re-cluster the remaining tumour cells and name the subclusters
epi <- ProcessSingleCellData(epi, group.by.vars = "orig.ident", resolution.seq = seq(0.1, 0.6, 0.1))
Idents(epi) <- "RNA_snn_res.0.2"
epi <- RenameIdents(epi, "0" = "Epi1", "1" = "Epi2", "2" = "Epi3", "3" = "Epi4", "4" = "Epi5")
epi$celltype <- Idents(epi)
FindClusterMarkersParallel(epi, output_prefix = "Epi_sub_", output_dir = marker_dir, workers = n_workers,
                           resolution_columns = paste0("RNA_snn_res.0.", 1:4))
saveRDS(epi, file.path(object_dir, "Epi_OV2.rds"))

# -----------------------------------------------------------------------------
# 2. UMAP figures
# -----------------------------------------------------------------------------
p <- plot_umap_labeled(epi, "celltype", pal_epi, split_by = "group",
                       title = "UMAP plot by Group", legend_text_size = 16)
save_pdf(p, file.path(figure_dir, "Epi_umap_split_group.pdf"), width = 10, height = 6)

p <- plot_umap_labeled(epi, "celltype", pal_epi, title = "UMAP plot", legend_text_size = 18)
save_pdf(p, file.path(figure_dir, "Epi_umap.pdf"), width = 10, height = 6)

# -----------------------------------------------------------------------------
# 3. Steroid (cholesterol) biosynthesis activity
# -----------------------------------------------------------------------------
epi_metab <- sc.metabolism.Seurat(obj = epi, method = "VISION", imputation = FALSE,
                                  ncores = n_cores, metabolism.type = "KEGG")
input_pathway <- "Steroid biosynthesis"

for (pheno in c("celltype", "group")) {
  p <- DotPlot.metabolism(obj = epi_metab, pathway = input_pathway, phenotype = pheno, norm = "y",
                          axis.text.x.size = 16, axis.text.y.size = 16,
                          legend.text.size = 14, legend.title.size = 14)
  ggsave(file.path(figure_dir, paste0("Epi_metabolism_dotplot_by_", pheno, ".pdf")),
         p, width = 10, height = 6, units = "in")
}

# -----------------------------------------------------------------------------
# 4. PCSK9 expression
# -----------------------------------------------------------------------------
p <- plot_gene_violin(epi, "PCSK9", "celltype", colors = pal_violin_epi, x_label = "Cluster",
                      title_size = 27, axis_title_size = 26, axis_text_size = 24)
save_pdf(p, file.path(figure_dir, "Epi_PCSK9_by_celltype.pdf"), width = 10, height = 8)

p <- plot_gene_violin(epi, "PCSK9", "group", colors = pal_violin, add_stats = TRUE,
                      stat_method = "t.test", x_label = "Group",
                      title_size = 27, axis_title_size = 26, axis_text_size = 24)
save_pdf(p, file.path(figure_dir, "Epi_PCSK9_by_group.pdf"), width = 10, height = 8)
