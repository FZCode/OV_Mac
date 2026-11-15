# =============================================================================
# Plotting helpers
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Save a ggplot as PDF (cairo_pdf when available, otherwise the default device)
save_pdf <- function(plot, file, width, height) {
  tryCatch(
    ggsave(file, plot = plot, width = width, height = height, units = "in", device = cairo_pdf),
    error = function(e) {
      warning("cairo_pdf not available, falling back to the default PDF device")
      ggsave(file, plot = plot, width = width, height = height, units = "in", device = "pdf")
    }
  )
  invisible(file)
}

#' Minimal UMAP theme (arrowed, truncated axes; no ticks)
umap_theme <- function(legend_text_size = 16) {
  theme_classic() +
    theme(
      legend.key.size  = unit(1, "cm"),
      legend.text      = element_text(size = legend_text_size),
      legend.title     = element_blank(),
      axis.line        = element_line(arrow = arrow(type = "closed")),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title.x     = element_text(hjust = 0, vjust = -0.5, size = 14),
      axis.title.y     = element_text(hjust = 0, vjust = 1.5, size = 14),
      plot.title       = element_text(size = 16, face = "plain", hjust = 0.5),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text       = element_text(face = "plain", size = 12),
      panel.spacing    = unit(0.5, "cm")
    )
}

#' UMAP coloured by a metadata column, with labels at the cluster medians
#'
#' @param label_col metadata column used for colours and labels
#' @param split_by  optional metadata column for facetting (e.g. "group")
plot_umap_labeled <- function(obj, label_col, colors, split_by = NULL,
                              title = "UMAP plot", legend_text_size = 16,
                              point_size = 0.4, label_size = 6) {
  emb <- Embeddings(obj, reduction = "umap")
  df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                   label = obj@meta.data[[label_col]])
  group_vars <- "label"
  if (!is.null(split_by)) {
    df$split <- obj@meta.data[[split_by]]
    group_vars <- c("split", "label")
  }
  centers <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")

  axis_trunc <- ggh4x::guide_axis_truncated(trunc_lower = unit(0, "npc"),
                                            trunc_upper = unit(3, "cm"))
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = factor(label))) +
    geom_point(size = point_size, alpha = 0.8) +
    geom_text(data = centers, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              color = "black", size = label_size, fontface = "plain") +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(title = NULL, override.aes = list(shape = 15, size = 3)),
           x = axis_trunc, y = axis_trunc) +
    umap_theme(legend_text_size) +
    labs(title = title, x = "UMAP_1", y = "UMAP_2")
  if (!is.null(split_by)) p <- p + facet_wrap(~ split, ncol = 2)
  p
}

#' Violin plot of one gene across the levels of a metadata column
#'
#' When `add_stats = TRUE`, all pairwise comparisons with >= 3 cells per group
#' are tested with ggpubr::stat_compare_means.
plot_gene_violin <- function(obj, gene, group_col, colors = pal_violin,
                             add_stats = FALSE, stat_method = "wilcox.test",
                             p_adjust_method = "bonferroni", hide_ns = TRUE,
                             p_label_format = "p.format", p_label_size = NULL,
                             x_label = group_col, x_angle = 0,
                             title_size = 26, axis_title_size = 25, axis_text_size = 23,
                             font_family = "sans") {
  if (!gene %in% rownames(obj)) stop("Gene not found: ", gene)
  if (!group_col %in% colnames(obj@meta.data)) stop("Column not found in meta.data: ", group_col)

  df <- data.frame(
    Expression = GetAssayData(obj, assay = "RNA", layer = "data")[gene, ],
    Group      = obj@meta.data[[group_col]]
  )
  df <- df[!is.na(df$Expression) & !is.na(df$Group), ]
  df$Group <- factor(df$Group)
  n_groups <- nlevels(df$Group)

  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
    geom_violin(trim = TRUE, alpha = 0.8) +
    geom_jitter(size = 0.4, alpha = 0.5, width = 0.2, color = "black") +
    labs(title = gene, x = x_label, y = "Expression Level") +
    scale_fill_manual(values = rep(colors, length.out = n_groups)) +
    coord_cartesian(ylim = c(0, NA)) +
    theme_classic(base_family = font_family) +
    theme(
      legend.position = "none",
      plot.title   = element_text(hjust = 0.5, size = title_size, face = "bold"),
      axis.title   = element_text(size = axis_title_size),
      axis.text.x  = element_text(size = axis_text_size, angle = x_angle,
                                  hjust = if (x_angle == 0) 0.5 else 1),
      axis.text.y  = element_text(size = axis_text_size)
    )

  if (add_stats && n_groups >= 2) {
    n_per_group <- table(df$Group)
    comparisons <- combn(levels(df$Group), 2, simplify = FALSE)
    comparisons <- Filter(function(cc) all(n_per_group[cc] >= 3), comparisons)
    if (length(comparisons) > 0) {
      p <- p + ggpubr::stat_compare_means(
        comparisons     = comparisons,
        method          = stat_method,
        p.adjust.method = p_adjust_method,
        hide.ns         = hide_ns,
        step.increase   = 0.1,
        size            = if (is.null(p_label_size)) max(2.5, axis_text_size * 0.25) else p_label_size,
        label           = p_label_format,
        vjust           = -0.2
      )
    }
  }
  p
}

#' Violin plots for several genes x several grouping variables
#'
#' Statistical tests are added for every grouping except cell-type groupings.
#' Writes one PDF per gene/grouping and a summary CSV with the smallest raw
#' pairwise p-value per plot. Extra arguments are passed to plot_gene_violin().
plot_genes_by_groups <- function(obj, genes, grouping_vars = c("group", "celltype"),
                                 output_dir, output_prefix = "gene_expression",
                                 colors = pal_violin, stat_method = "wilcox.test",
                                 width = 10, height = 8, ...) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  genes_ok <- intersect(genes, rownames(obj))
  if (length(genes_ok) == 0) stop("None of the genes are present in the object")
  if (length(genes_ok) < length(genes)) {
    warning("Genes not found: ", paste(setdiff(genes, genes_ok), collapse = ", "))
  }
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")

  rows <- list()
  for (gene in genes_ok) {
    for (var in grouping_vars) {
      is_celltype <- grepl("celltype", var, ignore.case = TRUE)
      p <- plot_gene_violin(obj, gene, var, colors = colors, add_stats = !is_celltype,
                            stat_method = stat_method, x_angle = 45, ...)
      g <- factor(obj@meta.data[[var]])
      w <- if (is_celltype && nlevels(g) > 8) width * 1.5 else width
      save_pdf(p, file.path(output_dir, paste0(output_prefix, "_", gene, "_by_", var, ".pdf")), w, height)

      # Smallest raw pairwise p-value (for the summary table only)
      min_p <- NA_real_; best <- ""
      if (!is_celltype && nlevels(g) >= 2) {
        x <- expr[gene, ]
        for (cc in combn(levels(g), 2, simplify = FALSE)) {
          a <- x[which(g == cc[1])]; b <- x[which(g == cc[2])]
          if (length(a) < 3 || length(b) < 3) next
          pv <- tryCatch(
            if (stat_method == "t.test") t.test(a, b)$p.value else wilcox.test(a, b)$p.value,
            error = function(e) NA_real_)
          if (!is.na(pv) && (is.na(min_p) || pv < min_p)) {
            min_p <- pv; best <- paste(cc[1], "vs", cc[2])
          }
        }
      }
      rows[[length(rows) + 1]] <- data.frame(
        Gene = gene, Grouping = var, NumGroups = nlevels(g),
        HasSignificantDiff = !is.na(min_p) && min_p < 0.05,
        MostSignificantComparison = best, MostSignificantPvalue = min_p)
    }
  }
  summary_tab <- do.call(rbind, rows)
  write.csv(summary_tab, file.path(output_dir, paste0(output_prefix, "_summary.csv")), row.names = FALSE)
  summary_tab
}

#' Mean expression of every feature (rows) within each group of cells (columns)
group_average <- function(matrix_data, group_vector) {
  matrix_data <- as.matrix(matrix_data)
  if (length(group_vector) != ncol(matrix_data)) {
    stop("Length of group_vector must equal the number of columns of matrix_data")
  }
  idx <- split(seq_along(group_vector), group_vector)
  idx <- idx[lengths(idx) > 0]
  avg <- sapply(idx, function(i) rowMeans(matrix_data[, i, drop = FALSE], na.rm = TRUE))
  rownames(avg) <- rownames(matrix_data)
  avg
}

#' One PDF per gene: Nebulosa density plot (falls back to FeaturePlot) or FeaturePlot
plot_genes_per_file <- function(obj, genes, output_dir, type = c("density", "feature"),
                                min_cells = 3, width = 6, height = 4) {
  type <- match.arg(type)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  genes <- intersect(genes, rownames(obj))
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")

  for (gene in genes) {
    n_pos <- sum(expr[gene, ] > 0)
    if (n_pos < min_cells) {
      message("Skip ", gene, ": expressed in ", n_pos, " cell(s)")
      next
    }
    p <- NULL
    if (type == "density") {
      p <- tryCatch(Nebulosa::plot_density(obj, features = gene, pal = "cividis"),
                    error = function(e) { message(gene, ": plot_density failed - ", e$message); NULL })
    }
    if (is.null(p)) {
      p <- tryCatch(
        if (type == "density") {
          FeaturePlot(obj, features = gene, cols = viridis::magma(100), pt.size = 0.5, order = TRUE)
        } else {
          FeaturePlot(obj, features = gene, order = TRUE)
        },
        error = function(e) { message(gene, ": FeaturePlot failed - ", e$message); NULL })
    }
    if (is.null(p)) next
    ggsave(file.path(output_dir, paste0(gene, "_", type, ".pdf")), p, width = width, height = height)
  }
  message("Plots saved to: ", output_dir)
}
