#' Bar Plots for Gene Expression
#'
#' Create bar and box plots for gene expression data from BulkRNAseq objects.
#'
#' @param object A BulkRNAseq object
#' @param data Character, data source for expression data or name of gsva result
#' @param genes Character vector, gene names to plot
#' @param group Character, column name in metadata for grouping
#' @param palette Character vector, colors for groups
#' @param method Character, statistical test method for comparison
#' @importFrom methods setGeneric setMethod
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot geom_jitter scale_color_manual
#' @importFrom ggplot2 scale_fill_manual labs theme_bw facet_wrap element_text
#' @importFrom ggplot2 geom_hline geom_vline theme element_blank element_line
#' @importFrom stats dist hclust model.matrix
#' @importFrom utils data
#' @import dplyr tidyr
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' data(counts)
#' data(metadata)
#' bulk_obj <- BulkRNAseq(counts, metadata)
#' bar(bulk_obj, genes = "CD36")
#' bar(bulk_obj, genes = c("CD36", "GAPDH"))
#' }
setGeneric("bar", function(object,
                           data = "data",
                           genes = "CD36",
                           group = "group",
                           palette = c("#00AFBB", "#E7B800"),
                           method = "t.test") {
  standardGeneric("bar")
})

#' @rdname bar
setMethod("bar", "BulkRNAseq", function(object,
                                        data = "data",
                                        genes = "CD36",
                                        group = "group",
                                        palette = c("#00AFBB", "#E7B800"),
                                        method = "t.test") {

  # Check required packages
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required for bar plots. Please install it.", call. = FALSE)
  }

  # Validate inputs
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }

  if (!group %in% colnames(object@metadata)) {
    stop(paste0("Group variable '", group, "' not found in metadata"), call. = FALSE)
  }

  # Prepare expression data
  if (data == "data") {
    if (!all(genes %in% rownames(object@data))) {
      missing_genes <- genes[!genes %in% rownames(object@data)]
      stop(paste0("Genes not found in expression data: ",
                 paste(missing_genes, collapse = ", ")), call. = FALSE)
    }
    expr_data <- data.frame(
      group = object@metadata[[group]],
      t(object@data[genes, , drop = FALSE])
    )
  } else if (data %in% names(object@gsva)) {
    if (!all(genes %in% rownames(object@gsva[[data]]))) {
      missing_genes <- genes[!genes %in% rownames(object@gsva[[data]])]
      stop(paste0("Genes not found in GSVA data '", data, "': ",
                 paste(missing_genes, collapse = ", ")), call. = FALSE)
    }
    expr_data <- data.frame(
      group = object@metadata[[group]],
      t(object@gsva[[data]][genes, , drop = FALSE])
    )
  } else {
    stop("data must be 'data' or a name from gsva results", call. = FALSE)
  }

  # Ensure group is factor
  if (!is.factor(expr_data$group)) {
    expr_data$group <- as.factor(expr_data$group)
  }

  # Plot based on number of genes
  if (length(genes) == 1) {
    .plot_single_gene(expr_data, genes, palette, method)
  } else {
    .plot_multiple_genes(expr_data, genes, palette, method)
  }
})

# Internal function for single gene plotting
.plot_single_gene <- function(expr_data, gene_name, palette, method) {
  comparisons <- list(rev(levels(expr_data$group)))

  gggpubr::ggboxplot(
    expr_data,
    x = "group",
    y = gene_name,
    color = "group",
    palette = palette,
    add = "jitter"
  ) +
    gggpubr::stat_compare_means(
      comparisons = comparisons,
      method = method
    ) +
    ggplot2::ggtitle(paste("Expression of", gene_name))
}

# Internal function for multiple genes plotting
.plot_multiple_genes <- function(expr_data, gene_list, palette, method) {
  plot_data <- expr_data %>%
    dplyr::select(group, dplyr::all_of(gene_list)) %>%
    tidyr::pivot_longer(
      cols = -group,
      names_to = "gene",
      values_to = "expression"
    )

  my_comparisons <- list(rev(levels(expr_data$group)))

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = group, y = expression, fill = group)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~ gene, scales = "free_y") +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons,
      label = "p.signif",
      method = method
    ) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(
      title = "Expression of multiple genes",
      x = "Group",
      y = "Expression"
    )
}

#' Heatmap Visualization
#'
#' Create heatmaps for gene expression data from BulkRNAseq objects.
#'
#' @param object A BulkRNAseq object
#' @param data Character, data source: "data" for expression data or name of gsva result
#' @param genes Character vector of gene names or numeric for top differential genes
#' @param group Character, column name in metadata for annotation
#' @param show_rownames Logical, whether to show gene names
#' @param scale Character, scaling method ("row", "column", "none")
#' @param color_palette Color palette for heatmap
#'
#' @return A pheatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' hmap(bulk_obj, genes = c("CD36", "GAPDH"))
#' hmap(bulk_obj, genes = 10)  # Top 10 differential genes
#' }
setGeneric("hmap", function(object,
                           data = "data",
                           genes = "CD36",
                           group = "group",
                           show_rownames = TRUE,
                           scale = "row",
                           color_palette = NULL) {
  standardGeneric("hmap")
})

#' @rdname hmap
setMethod("hmap", "BulkRNAseq", function(object,
                                        data = "data",
                                        genes = "CD36",
                                        group = "group",
                                        show_rownames = TRUE,
                                        scale = "row",
                                        color_palette = NULL) {

  # Check required packages
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required for heatmaps. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for color palettes. Please install it.", call. = FALSE)
  }

  # Set default color palette
  if (is.null(color_palette)) {
    color_palette <- viridis::viridis(10, alpha = 1, begin = 0.5, end = 1, direction = 1)
  }

  # Validate group variable
  if (!group %in% colnames(object@metadata)) {
    stop(paste0("Group variable '", group, "' not found in metadata"), call. = FALSE)
  }

  # Prepare expression data based on genes input
  if (is.numeric(genes)) {
    expr_data <- .get_top_diff_genes(object, data, genes)
  } else if (is.character(genes)) {
    expr_data <- .get_specific_genes(object, data, genes)
  } else {
    stop("genes must be character vector or numeric", call. = FALSE)
  }

  # Create annotation
  annotation_col <- data.frame(group = object@metadata[[group]]) %>%
    `rownames<-`(colnames(expr_data))

  # Create heatmap
  pheatmap::pheatmap(
    expr_data,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_legend = TRUE,
    show_rownames = show_rownames,
    scale = scale,
    color = color_palette,
    cellwidth = 30,
    cellheight = 12,
    fontsize = 10
  )
})

# Internal function to get top differential genes
.get_top_diff_genes <- function(object, data, n_genes) {
  if (data == "data") {
    if (length(object@allDiff) == 0 || is.null(object@allDiff$allDiff)) {
      stop("No differential analysis results found. Run differential analysis first.", call. = FALSE)
    }
    diff_data <- object@allDiff$allDiff
    expr_matrix <- object@data
  } else if (data %in% names(object@gsva)) {
    if (!data %in% names(object@allDiff)) {
      stop(paste0("No differential analysis results found for '", data, "'"), call. = FALSE)
    }
    diff_data <- object@allDiff[[data]]
    expr_matrix <- object@gsva[[data]]
  } else {
    stop("data must be 'data' or a name from gsva results", call. = FALSE)
  }

  # Get top differential genes
  top_genes <- diff_data %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::mutate(up_down = ifelse(logFC > 0, "up", "down")) %>%
    dplyr::group_by(up_down) %>%
    dplyr::arrange(desc(abs(logFC))) %>%
    dplyr::slice_head(n = n_genes) %>%
    dplyr::pull(gene_id)

  if (length(top_genes) == 0) {
    stop("No significant differential genes found", call. = FALSE)
  }

  return(expr_matrix[top_genes, , drop = FALSE])
}

# Internal function to get specific genes
.get_specific_genes <- function(object, data, genes) {
  if (data == "data") {
    expr_matrix <- object@data
  } else if (data %in% names(object@gsva)) {
    expr_matrix <- object@gsva[[data]]
  } else {
    stop("data must be 'data' or a name from gsva results", call. = FALSE)
  }

  # Check if genes exist
  missing_genes <- genes[!genes %in% rownames(expr_matrix)]
  if (length(missing_genes) > 0) {
    stop(paste0("Genes not found: ", paste(missing_genes, collapse = ", ")), call. = FALSE)
  }

  return(expr_matrix[genes, , drop = FALSE])
}

#' Volcano Plot
#'
#' Create volcano plots for differential expression results.
#'
#' @param object A BulkRNAseq object
#' @param data Character, name of differential analysis result
#' @param fc_threshold Numeric, fold change threshold
#' @param highlight_fc_threshold Numeric, fold change threshold for highlighting genes
#' @param adj.P_value_threshold Numeric, adjusted p-value threshold
#' @param point_size Numeric, size of points
#' @param highlight_size Numeric, size of highlighted points
#' @param colors Character vector, colors for up, non-significant, down, and highlighted genes
#' @param symmetry Logical, whether to make x-axis symmetric
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' volcano(bulk_obj, data = "allDiff")
#' }
setGeneric("volcano", function(object,
                              data = "allDiff",
                              fc_threshold = 1,
                              highlight_fc_threshold = 3,
                              adj.P_value_threshold = 0.05,
                              point_size = 1.2,
                              highlight_size = 3,
                              colors = c("darkred", "black", "royalblue4", "green2"),
                              symmetry = FALSE) {
  standardGeneric("volcano")
})

#' @rdname volcano
setMethod("volcano", "BulkRNAseq", function(object,
                                           data = "allDiff",
                                           fc_threshold = 1,
                                           highlight_fc_threshold = 3,
                                           adj.P_value_threshold = 0.05,
                                           point_size = 1.2,
                                           highlight_size = 3,
                                           colors = c("darkred", "black", "royalblue4", "green2"),
                                           symmetry = FALSE) {

  # Check required packages
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for volcano plots. Please install it.", call. = FALSE)
  }

  # Get differential analysis data
  if (data == "allDiff" && "allDiff" %in% names(object@allDiff)) {
    diff_data <- object@allDiff$allDiff
  } else if (data %in% names(object@allDiff)) {
    diff_data <- object@allDiff[[data]]
  } else {
    stop(paste0("Differential analysis result '", data, "' not found"), call. = FALSE)
  }

  .plot_volcano(diff_data, fc_threshold, highlight_fc_threshold,
                adj.P_value_threshold, point_size, highlight_size,
                colors, symmetry)
})

# Internal volcano plot function
.plot_volcano <- function(diff_data, fc_threshold, highlight_fc_threshold,
                         adj.P_value_threshold, point_size, highlight_size,
                         colors, symmetry) {

  # Ensure gene_id column exists
  if (!"gene_id" %in% colnames(diff_data)) {
    diff_data$gene_id <- rownames(diff_data)
  }

  # Classify genes
  diff_data <- diff_data %>%
    dplyr::mutate(
      group = dplyr::case_when(
        adj.P.Val < adj.P_value_threshold & logFC > fc_threshold ~ "Up",
        adj.P.Val < adj.P_value_threshold & logFC < -fc_threshold ~ "Down",
        TRUE ~ "NS"
      ) %>% factor(levels = c("Up", "NS", "Down"))
    )

  # Get genes to highlight
  highlight_genes <- diff_data %>%
    dplyr::filter(abs(logFC) >= highlight_fc_threshold & adj.P.Val < adj.P_value_threshold)

  # Count differential genes
  gene_counts <- diff_data %>%
    dplyr::count(group) %>%
    dplyr::filter(group != "NS")

  plot_title <- sprintf("%d down, %d up",
                       sum(gene_counts$n[gene_counts$group == "Down"], na.rm = TRUE),
                       sum(gene_counts$n[gene_counts$group == "Up"], na.rm = TRUE))

  # Create plot
  p <- ggplot2::ggplot(diff_data, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
    ggplot2::geom_point(alpha = 0.8, size = point_size) +
    ggplot2::scale_color_manual(values = colors[1:3]) +
    ggplot2::labs(
      x = "log2 (fold change)",
      y = "-log10 (adj.P.Val)",
      title = plot_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "top"
    ) +
    ggplot2::geom_hline(yintercept = -log10(adj.P_value_threshold),
                       lty = 4, lwd = 0.6, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                       lty = 4, lwd = 0.6, alpha = 0.8)

  # Add highlighted genes
  if (nrow(highlight_genes) > 0) {
    p <- p +
      ggplot2::geom_point(data = highlight_genes,
                         alpha = 0.8, size = highlight_size,
                         color = colors[4]) +
      ggrepel::geom_text_repel(data = highlight_genes,
                              ggplot2::aes(label = gene_id),
                              color = "black", alpha = 0.8)
  }

  # Symmetric x-axis
  if (symmetry) {
    ce <- ceiling(max(abs(diff_data$logFC), na.rm = TRUE))
    p <- p + ggplot2::scale_x_continuous(limits = c(-ce, ce), expand = c(0, 0))
  }

  return(p)
}

#' MA Plot
#'
#' Create MA plots for differential expression results.
#'
#' @param object A BulkRNAseq object
#' @param data Character, name of differential analysis result
#' @param fc_threshold Numeric, fold change threshold
#' @param adj.P_value_threshold Numeric, adjusted p-value threshold
#' @param highlight_fc_threshold Numeric, fold change threshold for highlighting
#' @param point_size Numeric, size of points
#' @param highlight_size Numeric, size of highlighted points
#' @param colors Character vector, colors for different groups
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' MA(bulk_obj, data = "allDiff")
#' }
setGeneric("MA", function(object,
                          data = "allDiff",
                          fc_threshold = 1,
                          adj.P_value_threshold = 0.05,
                          highlight_fc_threshold = 4,
                          point_size = 1.2,
                          highlight_size = 3,
                          colors = c("darkred", "grey50", "royalblue4", "green4")) {
  standardGeneric("MA")
})

#' @rdname MA
setMethod("MA", "BulkRNAseq", function(object,
                                       data = "allDiff",
                                       fc_threshold = 1,
                                       adj.P_value_threshold = 0.05,
                                       highlight_fc_threshold = 4,
                                       point_size = 1.2,
                                       highlight_size = 3,
                                       colors = c("darkred", "grey50", "royalblue4", "green4")) {

  # Check required packages
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for MA plots. Please install it.", call. = FALSE)
  }

  # Get differential analysis data
  if (data == "allDiff" && "allDiff" %in% names(object@allDiff)) {
    diff_data <- object@allDiff$allDiff
  } else if (data %in% names(object@allDiff)) {
    diff_data <- object@allDiff[[data]]
  } else {
    stop(paste0("Differential analysis result '", data, "' not found"), call. = FALSE)
  }

  .plot_MA(diff_data, fc_threshold, adj.P_value_threshold,
           highlight_fc_threshold, point_size, highlight_size, colors)
})

# Internal MA plot function
.plot_MA <- function(diff_data, fc_threshold, adj.P_value_threshold,
                    highlight_fc_threshold, point_size, highlight_size, colors) {

  # Ensure required columns exist
  if (!"gene_id" %in% colnames(diff_data)) {
    diff_data$gene_id <- rownames(diff_data)
  }

  # Create AveExpr column if not present
  if (!"AveExpr" %in% colnames(diff_data)) {
    if ("baseMean" %in% colnames(diff_data)) {
      diff_data$AveExpr <- log2(diff_data$baseMean + 1)
    } else if ("logCPM" %in% colnames(diff_data)) {
      diff_data$AveExpr <- diff_data$logCPM
    } else {
      stop("Cannot find average expression column (baseMean, logCPM, or AveExpr)", call. = FALSE)
    }
  }

  # Classify genes
  diff_data <- diff_data %>%
    dplyr::mutate(
      group = dplyr::case_when(
        adj.P.Val < adj.P_value_threshold & logFC > fc_threshold ~ "Up",
        adj.P.Val < adj.P_value_threshold & logFC < -fc_threshold ~ "Down",
        TRUE ~ "NS"
      ) %>% factor(levels = c("Up", "NS", "Down"))
    )

  # Get genes to highlight
  highlight_genes <- diff_data %>%
    dplyr::filter(abs(logFC) >= highlight_fc_threshold & adj.P.Val < adj.P_value_threshold)

  # Count genes
  up_count <- sum(diff_data$group == "Up", na.rm = TRUE)
  down_count <- sum(diff_data$group == "Down", na.rm = TRUE)
  plot_title <- paste0(down_count, " down, ", up_count, " up")

  # Create plot
  p <- ggplot2::ggplot(diff_data, ggplot2::aes(x = AveExpr, y = logFC, color = group)) +
    ggplot2::geom_point(alpha = 0.8, size = point_size) +
    ggplot2::scale_color_manual(values = colors[1:3]) +
    ggplot2::labs(
      y = "log2 (Fold Change)",
      x = "log2 (Average Expression)",
      title = plot_title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.position = "top"
    ) +
    ggplot2::geom_hline(yintercept = c(fc_threshold, -fc_threshold),
                       lty = 2, lwd = 1) +
    ggplot2::geom_hline(yintercept = 0, lwd = 1.2)

  # Add highlighted genes
  if (nrow(highlight_genes) > 0) {
    p <- p +
      ggplot2::geom_point(data = highlight_genes,
                         alpha = 0.8, size = highlight_size,
                         color = colors[4]) +
      ggrepel::geom_text_repel(data = highlight_genes,
                              ggplot2::aes(label = gene_id),
                              color = "black", alpha = 0.8)
  }

  return(p)
}

#' PCA and Hierarchical Clustering Plot
#'
#' Create combined PCA and hierarchical clustering plots.
#'
#' @param object A BulkRNAseq object
#' @param dist_method Distance method for clustering
#' @param hclust_method Hierarchical clustering method
#' @param k Number of clusters
#' @param horiz Logical, horizontal dendrogram
#' @param k_colors Color palette for clusters
#' @param rect Logical, draw rectangles around clusters
#' @param rect_fill Logical, fill rectangles
#' @param pca_style PCA plot style
#' @param add_ellipses Logical, add confidence ellipses
#' @param repel Logical, use repel for labels
#'
#' @return A combined plot object
#' @export
#'
#' @examples
#' \dontrun{
#' pca_hc(bulk_obj)
#' }
setGeneric("pca_hc", function(object,
                              dist_method = c("euclidean", "manhattan"),
                              hclust_method = c("ward.D2", "complete", "average"),
                              k = 2,
                              horiz = FALSE,
                              k_colors = "lancet",
                              rect = TRUE,
                              rect_fill = TRUE,
                              pca_style = "default",
                              add_ellipses = TRUE,
                              repel = TRUE) {
  standardGeneric("pca_hc")
})

#' @rdname pca_hc
setMethod("pca_hc", "BulkRNAseq", function(object,
                                           dist_method = c("euclidean", "manhattan"),
                                           hclust_method = c("ward.D2", "complete", "average"),
                                           k = 2,
                                           horiz = FALSE,
                                           k_colors = "lancet",
                                           rect = TRUE,
                                           rect_fill = TRUE,
                                           pca_style = "default",
                                           add_ellipses = TRUE,
                                           repel = TRUE) {

  # Check required packages
  required_packages <- c("factoextra", "tinyarray", "ggsci", "patchwork")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required. Please install it."), call. = FALSE)
    }
  }

  # Check if group column exists
  if (!"group" %in% colnames(object@metadata)) {
    stop("'group' column not found in metadata", call. = FALSE)
  }

  # Match arguments
  dist_method <- match.arg(dist_method)
  hclust_method <- match.arg(hclust_method)

  # Get expression data and group information
  expr_data <- object@data
  group <- object@metadata$group

  # Create hierarchical clustering plot
  p1 <- .plot_hclust_dendrogram(
    expr_data = expr_data,
    dist_method = dist_method,
    hclust_method = hclust_method,
    k = k,
    horiz = horiz,
    k_colors = k_colors,
    rect = rect,
    rect_fill = rect_fill
  )

  # Create PCA plot
  p2 <- tinyarray::draw_pca(
    exp = expr_data,
    group_list = group,
    color = ggsci::pal_lancet(palette = "lanonc")(9),
    style = pca_style,
    addEllipses = add_ellipses,
    repel = repel
  )

  # Combine plots
  p1 | p2
})

# Internal function for hierarchical clustering dendrogram
.plot_hclust_dendrogram <- function(expr_data, dist_method, hclust_method,
                                   k, horiz, k_colors, rect, rect_fill) {

  # Calculate distance matrix
  dist_mat <- dist(t(expr_data), method = dist_method)

  # Hierarchical clustering
  hc <- hclust(dist_mat, method = hclust_method)

  # Create dendrogram
  factoextra::fviz_dend(
    hc,
    k = k,
    cex = 1,
    horiz = horiz,
    type = "rectangle",
    k_colors = k_colors,
    color_labels_by_k = TRUE,
    rect = rect,
    rect_fill = rect_fill,
    rect_border = k_colors
  )
}

