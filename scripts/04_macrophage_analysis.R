# =============================================================================
# 04  Macrophage analysis
#
# Sub-clustering of macrophages, cholesterol metabolism, M1/M2/MDSC signature
# scores, polarization gene heatmap and per-gene feature plots.
#
# Input : results/objects/scobj_annotated.rds
#         data/immune_signature.csv (optional, for signature scores)
# Output: results/objects/Mac_OV2.rds, figures, marker / DEG tables
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
  library(ComplexHeatmap)
  library(circlize)
})

scobj <- readRDS(file.path(object_dir, "scobj_annotated.rds"))

# -----------------------------------------------------------------------------
# 1. Sub-cluster macrophages
# -----------------------------------------------------------------------------
mac <- subset(scobj, idents = "Mac")
rm(scobj); gc()
mac <- ProcessSingleCellData(mac, group.by.vars = "orig.ident", resolution.seq = seq(0.1, 0.6, 0.1))
FindClusterMarkersParallel(mac, output_prefix = "Mac_", output_dir = marker_dir, workers = n_workers,
                           resolution_columns = "RNA_snn_res.0.1")

# Remove cluster 4 (resolution 0.1; check marker tables after re-running) and re-cluster
mac <- mac[, mac$RNA_snn_res.0.1 != "4"]
mac <- ProcessSingleCellData(mac, group.by.vars = "orig.ident", resolution.seq = 0.1)

# Name the subclusters. NOTE: verify this mapping with the marker table and the
# TAM marker dot plot below whenever the clustering is re-run.
Idents(mac) <- "RNA_snn_res.0.1"
mac <- RenameIdents(mac, "0" = "APOE_Mac", "1" = "MT_Mac", "2" = "Inflam_TAM", "3" = "KI67_TAM")
mac$celltype <- Idents(mac)
FindClusterMarkersParallel(mac, output_prefix = "Mac_sub_", output_dir = marker_dir, workers = n_workers,
                           resolution_columns = "RNA_snn_res.0.1")

# TAM subset markers per subcluster
p <- DotPlot(mac, features = intersect(tam_markers, rownames(mac)), group.by = "celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figure_dir, "Mac_TAM_markers_dotplot.pdf"), p, width = 12, height = 4)

# -----------------------------------------------------------------------------
# 2. UMAP figures
# -----------------------------------------------------------------------------
p <- plot_umap_labeled(mac, "celltype", pal_mac, split_by = "group",
                       title = "UMAP plot by Group", legend_text_size = 16)
save_pdf(p, file.path(figure_dir, "Mac_umap_split_group.pdf"), width = 10, height = 6)

p <- plot_umap_labeled(mac, "celltype", pal_mac, title = "UMAP plot", legend_text_size = 18)
save_pdf(p, file.path(figure_dir, "Mac_umap.pdf"), width = 10, height = 6)

# -----------------------------------------------------------------------------
# 3. Steroid (cholesterol) biosynthesis activity
# -----------------------------------------------------------------------------
mac_metab <- sc.metabolism.Seurat(obj = mac, method = "VISION", imputation = FALSE,
                                  ncores = n_cores, metabolism.type = "KEGG")
input_pathway <- "Steroid biosynthesis"
saveRDS(mac_metab, file.path(object_dir, "Mac_metabolism.rds"))

for (pheno in c("celltype", "group")) {
  p <- DotPlot.metabolism(obj = mac_metab, pathway = input_pathway, phenotype = pheno, norm = "y",
                          axis.text.x.size = 16, axis.text.y.size = 16,
                          legend.text.size = 14, legend.title.size = 14)
  ggsave(file.path(figure_dir, paste0("Mac_metabolism_dotplot_by_", pheno, ".pdf")),
         p, width = 10, height = 6, units = "in")
}
p <- BoxPlot.metabolism(obj = mac_metab, pathway = input_pathway, phenotype = "group", ncol = 1)
ggsave(file.path(figure_dir, "Mac_metabolism_boxplot_by_group.pdf"), p, width = 6, height = 6)

# -----------------------------------------------------------------------------
# 4. Cholesterol metabolism genes
# -----------------------------------------------------------------------------
p <- plot_gene_violin(mac, "PCSK9", "group", colors = pal_violin, add_stats = TRUE, x_label = "Group",
                      title_size = 27, axis_title_size = 26, axis_text_size = 24)
save_pdf(p, file.path(figure_dir, "Mac_PCSK9_by_group.pdf"), width = 10, height = 8)

lipid_summary <- plot_genes_by_groups(
  mac, genes = cholesterol_genes, grouping_vars = c("group", "celltype"),
  output_dir = file.path(figure_dir, "Mac_cholesterol_genes"), output_prefix = "lipid_genes",
  colors = pal_mac, title_size = 26, axis_title_size = 25, axis_text_size = 24
)
print(lipid_summary)

# -----------------------------------------------------------------------------
# 5. Immune signature scores (M1 / M2 macrophage, MDSC)
# -----------------------------------------------------------------------------
if (file.exists(signature_file)) {
  sig <- read.csv(signature_file, stringsAsFactors = FALSE)
  names(sig)[1:3] <- c("gene", "group", "species")
  sig$group[sig$group == "M-MDSC (top100)"] <- "M-MDSC"
  sig$group[sig$group == "G-MDSC (top100)"] <- "G-MDSC"
  sig$gene <- toupper(sig$gene)
  sig <- sig[sig$gene %in% rownames(mac), ]
  signature_list <- split(sig$gene, sig$group)

  if (all(c("M1 Macrophage", "M2 Macrophage") %in% names(signature_list))) {
    sig_m1m2 <- setNames(signature_list[c("M1 Macrophage", "M2 Macrophage")], c("M1", "M2"))
    saveRDS(sig_m1m2, file.path(object_dir, "Mac_M1M2_signatures.rds"))
  }

  # AddModuleScore names the columns SigScore1, SigScore2, ...; rename them to the signature names
  mac <- AddModuleScore(mac, features = signature_list, name = "SigScore")
  score_cols <- paste0("SigScore", seq_along(signature_list))
  colnames(mac@meta.data)[match(score_cols, colnames(mac@meta.data))] <- names(signature_list)

  # All signatures by subcluster (one page per signature)
  pdf(file.path(figure_dir, "Mac_signature_scores_by_celltype.pdf"), width = 8, height = 6)
  for (s in names(signature_list)) {
    print(VlnPlot(mac, features = s, group.by = "celltype", cols = pal_mac, pt.size = 0) +
            geom_hline(yintercept = mean(mac@meta.data[[s]]), lty = 2))
  }
  dev.off()

  # M2 score by group
  if ("M2 Macrophage" %in% names(signature_list)) {
    p <- VlnPlot(mac, features = "M2 Macrophage", group.by = "group", cols = pal_mac, pt.size = 0) +
      geom_hline(yintercept = mean(mac$`M2 Macrophage`), lty = 2) +
      theme(plot.title   = element_text(size = 20, face = "bold"),
            axis.title   = element_text(size = 20),
            axis.text.x  = element_text(size = 18, angle = 45, hjust = 1),
            axis.text.y  = element_text(size = 18),
            legend.title = element_text(size = 18),
            legend.text  = element_text(size = 18))
    ggsave(file.path(figure_dir, "Mac_M2_group.pdf"), p, width = 10, height = 6, units = "in")
  }
} else {
  message("Signature file not found (", signature_file, "); skipping signature scores")
}

saveRDS(mac, file.path(object_dir, "Mac_OV2.rds"))

# -----------------------------------------------------------------------------
# 6. Polarization gene heatmap (mean expression per subcluster, row Z-score)
# -----------------------------------------------------------------------------
heat_genes <- intersect(polarization_genes, rownames(mac))
message(length(heat_genes), " of ", length(polarization_genes), " polarization genes found; missing: ",
        paste(setdiff(polarization_genes, heat_genes), collapse = ", "))

avg_expr <- group_average(GetAssayData(mac, assay = "RNA", layer = "data")[heat_genes, , drop = FALSE],
                          mac$celltype)
gene_var <- apply(avg_expr, 1, var, na.rm = TRUE)
avg_expr <- avg_expr[!is.na(gene_var) & gene_var > 0, , drop = FALSE]
scaled_expr <- t(scale(t(avg_expr)))

ht <- Heatmap(scaled_expr,
              name = "Z-score",
              col = colorRamp2(c(-2, 0, 2), c("#007ABA", "white", "#ac6894")),
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D2",
              width = unit(4, "cm"),
              height = unit(15, "cm"),
              rect_gp = gpar(col = "black", lwd = 0.5),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 10))
pdf(file.path(figure_dir, "Mac_polarization_heatmap.pdf"), width = 4, height = 7)
draw(ht)
dev.off()

# -----------------------------------------------------------------------------
# 7. One density plot and one FeaturePlot per polarization gene
# -----------------------------------------------------------------------------
plot_genes_per_file(mac, heat_genes, file.path(figure_dir, "Mac_gene_density"), type = "density")
plot_genes_per_file(mac, heat_genes, file.path(figure_dir, "Mac_gene_featureplot"), type = "feature")

# -----------------------------------------------------------------------------
# 8. Subcluster DEGs
# -----------------------------------------------------------------------------
degs <- FindAllMarkers(mac)
write.csv(degs, file.path(results_dir, "Mac_celltype_DEGs.csv"))
