# =============================================================================
# Shared configuration: paths, palettes and gene sets.
# All scripts assume that the working directory is the repository root.
# =============================================================================

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
data_dir    <- "data"
raw_tar     <- file.path(data_dir, "GSE184880_RAW.tar")   # downloaded from GEO
raw_dir     <- file.path(data_dir, "GSE184880")           # extracted 10x files
results_dir <- "results"
object_dir  <- file.path(results_dir, "objects")
marker_dir  <- file.path(results_dir, "markers")
figure_dir  <- file.path(results_dir, "figures")
for (d in c(object_dir, marker_dir, figure_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Optional immune signature table used in script 04
# (three columns: gene, group, species; e.g. "M1 Macrophage", "M2 Macrophage", "M-MDSC (top100)")
signature_file <- file.path(data_dir, "immune_signature.csv")

# -----------------------------------------------------------------------------
# Computing resources
# -----------------------------------------------------------------------------
n_workers <- 20   # parallel workers for marker detection
n_cores   <- 15   # cores for scMetabolism (VISION)

# -----------------------------------------------------------------------------
# Colour palettes
# -----------------------------------------------------------------------------
pal_major      <- c("#d4de9c", "#94c58f", "#86c7b4", "#9cd2ed", "#bb82b1", "#ea9994", "#f2c396", "#bf812d")
pal_epi        <- c("#4b6aa8", "#3ca0cf", "#83ab8e", "#e6b884", "#cc7f73")
pal_mac        <- c("#E7CB94", "#8d689d", "#537eb7", "#d8a0c0")
pal_violin     <- c("#d4de9c", "#94c58f", "#86c7b4", "#9cd2ed", "#a992c0", "#ea9994", "#f2c396", "#bb82b1")
pal_violin_epi <- c("#d4de9c", "#86c7b4", "#9cd2ed", "#a992c0", "#ea9994", "#f2c396", "#bb82b1")

# -----------------------------------------------------------------------------
# Gene sets
# -----------------------------------------------------------------------------
# Cholesterol metabolism
cholesterol_genes <- c(
  "HMGCR", "SQLE",                      # synthesis
  "LDLR", "SCARB1",                     # uptake
  "ABCA1", "ABCG1", "ABCG5", "ABCG8",   # efflux
  "ACAT1", "ACAT2",                     # esterification
  "SREBF2", "NRF1",                     # transcriptional regulation
  "PCSK9"                               # LDLR degradation
)

# Tumour-associated macrophage (TAM) subset markers
tam_markers <- c(
  "CXCL10", "CD274", "ISG15", "IDO1",                              # IFN-TAMs
  "ARG1", "MRC1", "CX3CR1",                                        # Reg-TAMs
  "IL1B", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CCL3", "CCL3L1",    # Inflam-TAMs
  "VEGFA", "SPP1",                                                 # Angio-TAMs
  "APOC1", "APOE", "ACP5", "FABP5",                                # LA-TAMs
  "LYVE1", "HES1", "FOLR2",                                        # RTM-TAMs
  "MKI67", "TOP2A", "PCNA"                                         # Prolif-TAMs
)

# Macrophage polarization / effector genes (heatmap and per-gene feature plots).
# Entries that are not official gene symbols (e.g. MHCII, CD206) are dropped automatically.
polarization_genes <- unique(c(
  "IL12", "IL23", "TNF", "IL6", "CD86", "MHCII", "IL1B", "MARCO", "NOS2", "CD64", "CD80",
  "CXCR10", "CXCL9", "CXCL10", "CXCL11", "IL1A", "TNFa", "CCL5", "IRF5", "IRF1", "CD40",
  "IDO1", "KYNU", "CCR7", "ARG1", "ARG2", "IL10", "CD32", "CD163", "CD23", "CD200R1",
  "PDCD1LG2", "CD274", "CSF1R", "CD206", "IL1RN", "IL1R2", "IL4R", "CCL4", "CCL13",
  "CCL20", "CCL17", "CCL18", "CCL22", "CCL24", "LYVE1", "VEGFA", "VEGFB", "VEGFC",
  "VEGFD", "EGF", "CTSA", "CTSB", "CSTC", "CTSD", "TGFB1", "TGFB2", "TGFB3", "MMP14",
  "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASL", "TNFSF12", "TNFSF8", "CD276", "VTCN1",
  "MSR1", "FN1", "IRF4"
))

set.seed(2025)
