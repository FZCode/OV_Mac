# =============================================================================
# 01  Data preparation
#
# Extract the GSE184880 raw files, organise them into one 10x folder per
# sample, build one Seurat object per sample and merge them.
#
# Input : data/GSE184880_RAW.tar
# Output: results/objects/scobj_merged.rds
# =============================================================================

source("R/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
})

# -----------------------------------------------------------------------------
# 1. Extract the GEO archive
# -----------------------------------------------------------------------------
if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE)
  untar(raw_tar, exdir = raw_dir)
}

# -----------------------------------------------------------------------------
# 2. One folder per sample with standard 10x file names, e.g.
#    GSM5599220_Norm1.barcodes.tsv.gz -> GSM5599220_Norm1/barcodes.tsv.gz
#    genes.tsv.gz is renamed to features.tsv.gz
# -----------------------------------------------------------------------------
raw_files  <- list.files(raw_dir, pattern = "\\.(tsv|mtx)\\.gz$")
sample_ids <- sub("\\..*$", "", raw_files)

for (sid in unique(sample_ids)) {
  sample_dir <- file.path(raw_dir, sid)
  dir.create(sample_dir, showWarnings = FALSE)
  for (f in raw_files[sample_ids == sid]) {
    new_name <- sub("^[^.]*\\.", "", f)
    new_name <- sub("^genes\\.tsv\\.gz$", "features.tsv.gz", new_name)
    file.rename(file.path(raw_dir, f), file.path(sample_dir, new_name))
  }
}

# -----------------------------------------------------------------------------
# 3. Create and merge Seurat objects
# -----------------------------------------------------------------------------
sample_dirs <- list.dirs(raw_dir, full.names = TRUE, recursive = FALSE)
scobj_list <- lapply(sample_dirs, function(d) {
  CreateSeuratObject(counts = Read10X(data.dir = d), project = basename(d),
                     min.cells = 3, min.features = 300)
})
names(scobj_list) <- basename(sample_dirs)

scobj_merged <- merge(x = scobj_list[[1]], y = scobj_list[-1], add.cell.ids = names(scobj_list))
scobj_merged[["RNA"]] <- JoinLayers(scobj_merged[["RNA"]])

# Sample name (e.g. Norm1 / Cancer7) and group (Normal / Cancer)
scobj_merged$sample <- str_extract(scobj_merged$orig.ident, "(?<=_).*")
scobj_merged$group  <- ifelse(grepl("Norm", scobj_merged$orig.ident), "Normal", "Cancer")
print(table(scobj_merged$sample, scobj_merged$group))

saveRDS(scobj_merged, file.path(object_dir, "scobj_merged.rds"))
