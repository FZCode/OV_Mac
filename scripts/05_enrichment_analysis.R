# =============================================================================
# 05  Enrichment analysis of epithelial cells (SCP)
#
# GSEA between epithelial subclusters and over-representation analysis of
# Cancer vs Normal DEGs.
#
# Input : results/objects/Epi_OV2.rds
# Output: results/objects/Epi_OV2_enrichment.rds, figures
# =============================================================================

source("R/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SCP)
  library(BiocParallel)
  library(ggplot2)
})

epi <- readRDS(file.path(object_dir, "Epi_OV2.rds"))

# Parallel back-end (single SOCK worker on Windows to avoid serialization issues)
bpp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 1, type = "SOCK", stop.on.error = FALSE)
} else {
  MulticoreParam(workers = min(8, parallel::detectCores()))
}

# -----------------------------------------------------------------------------
# 1. GSEA between epithelial subclusters
# -----------------------------------------------------------------------------
epi <- RunDEtest(srt = epi, group_by = "celltype", fc.threshold = 1, only.pos = FALSE)
epi <- RunGSEA(srt = epi, group_by = "celltype", db = c("KEGG", "GO_BP"),
               species = "Homo_sapiens", DE_threshold = "p_val_adj < 0.1", BPPARAM = bpp)

# GO:0006695 = cholesterol biosynthetic process
pdf(file.path(figure_dir, "Epi3_GSEA_cholesterol_biosynthesis.pdf"), width = 8, height = 6)
print(GSEAPlot(srt = epi, group_by = "celltype", db = "GO_BP",
               group_use = "Epi3", id_use = "GO:0006695"))
dev.off()

# -----------------------------------------------------------------------------
# 2. Cancer vs Normal: DEGs and KEGG over-representation
# -----------------------------------------------------------------------------
epi <- RunDEtest(srt = epi, group_by = "group", fc.threshold = 1, only.pos = FALSE)

pdf(file.path(figure_dir, "Epi_group_volcano.pdf"), width = 8, height = 6)
print(VolcanoPlot(srt = epi, group_by = "group"))
dev.off()

epi <- RunEnrichment(srt = epi, group_by = "group", db = "KEGG", species = "Homo_sapiens",
                     DE_threshold = "avg_log2FC > log2(0.25) & p_val_adj < 0.05")

pdf(file.path(figure_dir, "Epi_Cancer_KEGG_enrichment.pdf"), width = 10, height = 8)
print(EnrichmentPlot(srt = epi, group_by = "group", group_use = "Cancer",
                     topTerm = 20, plot_type = "bar"))
dev.off()

saveRDS(epi, file.path(object_dir, "Epi_OV2_enrichment.rds"))
