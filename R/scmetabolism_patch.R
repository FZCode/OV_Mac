# =============================================================================
# scMetabolism helpers
#
# Patched versions of two functions from the scMetabolism package
# (Wu et al., Cancer Discovery 2021; https://pubmed.ncbi.nlm.nih.gov/34417225/):
#   - sc.metabolism.Seurat: reads counts with `layer =` (Seurat v5 compatible)
#   - DotPlot.metabolism  : adds control over font and dot sizes
# Defined in the global environment, so they take precedence over the package
# versions even after library(scMetabolism).
# =============================================================================

#' Score KEGG / REACTOME metabolic pathways for every cell
#'
#' Scores are stored in obj@assays$METABOLISM$score (pathways x cells).
sc.metabolism.Seurat <- function(obj, method = "VISION", imputation = FALSE, ncores = 2,
                                 metabolism.type = "KEGG") {
  countexp <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "counts"))

  gmt_file <- switch(
    metabolism.type,
    KEGG     = system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism"),
    REACTOME = system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolism"),
    stop("metabolism.type must be 'KEGG' or 'REACTOME'")
  )
  message("Pathway database: ", metabolism.type)

  if (imputation) {
    # ALRA imputation (Linderman et al., bioRxiv 397588); alra() comes from the ALRA scripts
    message("Running ALRA imputation...")
    imputed <- alra(countexp)[[3]]
    rownames(imputed) <- rownames(countexp)
    countexp <- imputed
  }

  message("Quantifying metabolic activity with ", method, "...")
  if (method == "VISION") {
    n_umi <- colSums(countexp)
    scaled_counts <- t(t(countexp) / n_umi) * median(n_umi)
    vis <- VISION::Vision(scaled_counts, signatures = gmt_file)
    options(mc.cores = ncores)
    vis <- VISION::analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))
  } else if (method == "AUCell") {
    rankings <- AUCell::AUCell_buildRankings(countexp, nCores = ncores, plotStats = FALSE)
    cells_auc <- AUCell::AUCell_calcAUC(GSEABase::getGmt(gmt_file), rankings)
    signature_exp <- data.frame(AUCell::getAUC(cells_auc))
  } else if (method %in% c("ssGSEA", "GSVA")) {
    # Uses the classic gsva() interface (GSVA < 1.50); newer versions require gsvaParam()/ssgseaParam()
    es <- GSVA::gsva(countexp, GSEABase::getGmt(gmt_file), method = tolower(method),
                     kcdf = "Poisson", parallel.sz = ncores)
    signature_exp <- data.frame(es)
  } else {
    stop("method must be one of: VISION, AUCell, ssGSEA, GSVA")
  }

  obj@assays$METABOLISM$score <- signature_exp
  obj
}

#' Dot plot of the median pathway score per group
#'
#' @param norm "y" = scale each pathway to [0, 1] across groups,
#'             "x" = scale each group to [0, 1] across pathways, "na" = raw medians
DotPlot.metabolism <- function(obj, pathway, phenotype, norm = "y",
                               axis.text.x.size = 12, axis.text.y.size = 12,
                               legend.text.size = 10, legend.title.size = 12,
                               dot.min.size = 1, dot.max.size = 8) {
  meta  <- obj@meta.data
  score <- obj@assays$METABOLISM$score

  df <- data.frame(
    group   = rep(as.character(meta[[phenotype]]), times = length(pathway)),
    pathway = rep(pathway, each = nrow(meta)),
    value   = unlist(lapply(pathway, function(pw) as.numeric(score[pw, ])))
  )
  med <- aggregate(value ~ group + pathway, data = df, FUN = median)

  range01 <- function(x) (x - min(x)) / (max(x) - min(x))
  if (norm == "y") med$value <- ave(med$value, med$pathway, FUN = range01)
  if (norm == "x") med$value <- ave(med$value, med$group, FUN = range01)

  pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
  ggplot(med, aes(x = group, y = pathway, color = value, size = value)) +
    geom_point() +
    xlab(phenotype) + ylab("Metabolic Pathway") +
    scale_size(range = c(dot.min.size, dot.max.size)) +
    scale_color_gradientn(colours = pal) +
    labs(color = "Value", size = "Value") +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = axis.text.x.size),
      axis.text.y  = element_text(size = axis.text.y.size),
      legend.text  = element_text(size = legend.text.size),
      legend.title = element_text(size = legend.title.size),
      panel.grid   = element_blank()
    )
}
