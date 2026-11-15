# Cholesterol metabolism of tumour cells and macrophages in ovarian cancer (scRNA-seq)

Single-cell RNA-seq analysis of normal ovary and ovarian cancer samples from
GEO [GSE184880](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184880).
The workflow annotates the major cell types, sub-clusters epithelial (tumour)
cells and macrophages, and compares cholesterol metabolism (KEGG *Steroid
biosynthesis*, cholesterol synthesis / uptake / efflux genes, PCSK9) between
subclusters and between Normal and Cancer samples.

## Repository structure

```
.
├── R/
│   ├── config.R                 # paths, palettes, gene sets
│   ├── preprocessing.R          # Seurat + Harmony workflow, parallel marker detection
│   ├── plotting.R               # UMAP, violin, heatmap and feature-plot helpers
│   └── scmetabolism_patch.R     # Seurat v5 compatible scMetabolism functions
├── scripts/
│   ├── 01_data_preparation.R    # extract GEO files, build and merge Seurat objects
│   ├── 02_major_cell_annotation.R
│   ├── 03_epithelial_analysis.R
│   ├── 04_macrophage_analysis.R
│   └── 05_enrichment_analysis.R # SCP GSEA / enrichment of epithelial cells
├── data/                        # raw data (not tracked)
└── results/                     # objects, marker tables, figures (not tracked)
```

## Data

1. Download `GSE184880_RAW.tar` from GEO and place it in `data/`.
2. (Optional) Put an immune signature table in `data/immune_signature.csv`
   with three columns: `gene`, `group`, `species`. Groups named
   `M1 Macrophage`, `M2 Macrophage`, `M-MDSC (top100)` and `G-MDSC (top100)`
   are used in script 04. Without this file the signature-score step is skipped.

## Requirements

R >= 4.3 and the following packages:

- CRAN: `Seurat` (>= 5.0), `harmony`, `ggplot2`, `ggh4x`, `ggpubr`, `dplyr`,
  `stringr`, `foreach`, `doParallel`, `wesanderson`, `circlize`, `viridis`,
  `clustree` (optional)
- Bioconductor: `ComplexHeatmap`, `BiocParallel`, `AUCell`, `GSEABase`, `GSVA`
  (only for the non-default scMetabolism methods), `Nebulosa`
- GitHub: [`scMetabolism`](https://github.com/wu-yc/scMetabolism),
  [`VISION`](https://github.com/YosefLab/VISION),
  [`SCP`](https://github.com/zhanghao-njmu/SCP)

## Usage

Run the scripts in order from the repository root:

```r
source("scripts/01_data_preparation.R")
source("scripts/02_major_cell_annotation.R")
source("scripts/03_epithelial_analysis.R")
source("scripts/04_macrophage_analysis.R")
source("scripts/05_enrichment_analysis.R")
```

Each script reads the object saved by the previous step from
`results/objects/`, so steps can be re-run independently. Paths and the number
of parallel workers can be changed in `R/config.R`.

## Notes

- Cluster annotations (`RenameIdents`) refer to the cluster numbers of the
  original run. Leiden/Louvain cluster numbers can change between package
  versions, so check the marker tables in `results/markers/` before reusing
  them.
- No mitochondrial / doublet filtering is applied beyond
  `min.cells = 3` and `min.features = 300`.

## Citation

If you use the metabolism scoring, please cite scMetabolism:
Wu Y, et al. *Cancer Discovery* 2021. https://pubmed.ncbi.nlm.nih.gov/34417225/
