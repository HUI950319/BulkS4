#' @importFrom stats setNames
#' @importFrom utils head tail

#' @title Run Differential Expression Analysis on BulkRNAseq Object
#' @description Perform differential expression analysis on BulkRNAseq object using DESeq2, edgeR, or limma-voom
#'
#' @param object A BulkRNAseq object
#' @param method Character, analysis method: "DESeq2", "edgeR", or "limma_voom"
#' @param group_col Character, column name in metadata for grouping, default "group"
#' @param treat_level Character, treatment level name, default "treat"
#' @param control_level Character, control level name, default "con"
#' @param fc_cutoff Numeric, log2 fold change threshold, default 1
#' @param padj_cutoff Numeric, adjusted p-value threshold, default 0.05
#' @param min_counts Numeric, minimum count threshold (DESeq2 only), default 1
#' @param normalize_method Character, normalization method (limma-voom only), default "quantile"
#' @param result_name Character, name for storing results in object, default NULL (auto-generated)
#'
#' @return Updated BulkRNAseq object with analysis results added
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example object
#' counts_mat <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' rownames(counts_mat) <- paste0("Gene", 1:100)
#' colnames(counts_mat) <- paste0("Sample", 1:10)
#'
#' metadata <- data.frame(
#'   group = rep(c("con", "treat"), each = 5),
#'   row.names = colnames(counts_mat)
#' )
#'
#' bulk_obj <- BulkRNAseq(counts_mat, metadata)
#'
#' # Run different analyses
#' bulk_obj <- runDiffAnalysis(bulk_obj, method = "DESeq2")
#' bulk_obj <- runDiffAnalysis(bulk_obj, method = "edgeR")
#' bulk_obj <- runDiffAnalysis(bulk_obj, method = "limma_voom")
#' }
runDiffAnalysis <- function(object,
                           method = c("DESeq2", "edgeR", "limma_voom"),
                           group_col = "group",
                           treat_level = "treat",
                           control_level = "con",
                           fc_cutoff = 1,
                           padj_cutoff = 0.05,
                           min_counts = 1,
                           normalize_method = "quantile",
                           result_name = NULL) {

  # Input validation
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }

  method <- match.arg(method)

  # Check if group column exists
  if (!group_col %in% colnames(object@metadata)) {
    stop("Column '", group_col, "' not found in metadata", call. = FALSE)
  }

  # Get counts and metadata
  counts <- getCounts(object)
  metadata <- getMetadata(object)

  # Prepare metadata for analysis
  if (group_col != "group") {
    metadata$group <- metadata[[group_col]]
  }

  # Check group levels
  unique_groups <- unique(metadata$group)
  if (!all(c(treat_level, control_level) %in% unique_groups)) {
    stop("Treatment level '", treat_level, "' and control level '", control_level,
         "' must be present in ", group_col, " column", call. = FALSE)
  }

  # Recode group levels to standard format
  metadata$group <- ifelse(metadata$group == treat_level, "treat", "con")

  # Generate result name if not provided
  if (is.null(result_name)) {
    result_name <- paste0(method, "_", treat_level, "_vs_", control_level)
  }

  # Run analysis based on method
  if (method == "DESeq2") {
    result <- run_DESeq2(counts, metadata, fc_cutoff, padj_cutoff, min_counts)
  } else if (method == "edgeR") {
    result <- run_edgeR(counts, metadata, fc_cutoff, padj_cutoff)
  } else if (method == "limma_voom") {
    result <- run_limma_voom(counts, metadata, normalize_method, fc_cutoff, padj_cutoff)
  }

  # Update object with results
  object@allDiff[[result_name]] <- result$allDiff

  # Update normalized data if this is the first analysis or if requested
  if (length(object@allDiff) == 1) {
    object@data <- as.matrix(result$exprSet)
  }

  # Print summary
  n_total <- nrow(result$allDiff)
  n_sig <- nrow(result$fltDiff)
  n_up <- sum(result$fltDiff$logFC > 0, na.rm = TRUE)
  n_down <- sum(result$fltDiff$logFC < 0, na.rm = TRUE)

  message("\n=====================================")
  message(" Differential Expression Analysis")
  message("=====================================")
  message(" Method: ", method)
  message(" Comparison: ", treat_level, " vs ", control_level)
  message("-------------------------------------")
  message(" Total genes tested: ", n_total)
  message(" Significant genes: ", n_sig)
  message(" Up-regulated: ", n_up)
  message(" Down-regulated: ", n_down)
  message(" FC cutoff: ", fc_cutoff)
  message(" Adj. p-value cutoff: ", padj_cutoff)
  message("=====================================\n")

  return(object)
}

#' @title Compare Multiple Differential Expression Methods
#' @description Run multiple differential expression methods and compare results
#'
#' @param object A BulkRNAseq object
#' @param methods Character vector, methods to compare: "DESeq2", "edgeR", "limma_voom"
#' @param group_col Character, column name in metadata for grouping
#' @param treat_level Character, treatment level name
#' @param control_level Character, control level name
#' @param fc_cutoff Numeric, log2 fold change threshold
#' @param padj_cutoff Numeric, adjusted p-value threshold
#' @param plot_comparison Logical, whether to create comparison plots
#'
#' @return Updated BulkRNAseq object with all analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' bulk_obj <- compareDiffMethods(bulk_obj,
#'                               methods = c("DESeq2", "edgeR", "limma_voom"))
#' }
compareDiffMethods <- function(object,
                              methods = c("DESeq2", "edgeR", "limma_voom"),
                              group_col = "group",
                              treat_level = "treat",
                              control_level = "con",
                              fc_cutoff = 1,
                              padj_cutoff = 0.05,
                              plot_comparison = TRUE) {

  # Run each method
  for (method in methods) {
    message("Running ", method, " analysis...")
    object <- runDiffAnalysis(
      object = object,
      method = method,
      group_col = group_col,
      treat_level = treat_level,
      control_level = control_level,
      fc_cutoff = fc_cutoff,
      padj_cutoff = padj_cutoff
    )
  }

  # Compare results
  if (length(methods) > 1) {
    message("\n=====================================")
    message(" Method Comparison Summary")
    message("=====================================")

    result_names <- paste0(methods, "_", treat_level, "_vs_", control_level)

    for (i in seq_along(methods)) {
      result <- object@allDiff[[result_names[i]]]
      n_sig <- sum(result$adj.P.Val < padj_cutoff & abs(result$logFC) > fc_cutoff, na.rm = TRUE)
      message(" ", methods[i], ": ", n_sig, " significant genes")
    }

    # Find common significant genes
    if (length(methods) >= 2) {
      sig_genes_list <- list()
      for (i in seq_along(methods)) {
        result <- object@allDiff[[result_names[i]]]
        sig_genes <- result$gene_id[result$adj.P.Val < padj_cutoff &
                                   abs(result$logFC) > fc_cutoff]
        sig_genes_list[[methods[i]]] <- sig_genes[!is.na(sig_genes)]
      }

      common_genes <- Reduce(intersect, sig_genes_list)
      message("-------------------------------------")
      message(" Common significant genes: ", length(common_genes))
      message("=====================================\n")
    }
  }

  return(object)
}

#' @title Extract Differential Expression Results
#' @description Extract differential expression results from BulkRNAseq object
#'
#' @param object A BulkRNAseq object
#' @param analysis_name Character, name of the analysis result to extract
#' @param significant_only Logical, whether to return only significant genes
#' @param fc_cutoff Numeric, fold change cutoff for significance
#' @param padj_cutoff Numeric, adjusted p-value cutoff for significance
#'
#' @return Data frame with differential expression results
#' @export
#'
#' @examples
#' \dontrun{
#' # Get all results
#' all_results <- extractDiffResults(bulk_obj, "DESeq2_treat_vs_con")
#'
#' # Get only significant results
#' sig_results <- extractDiffResults(bulk_obj, "DESeq2_treat_vs_con",
#'                                  significant_only = TRUE)
#' }
extractDiffResults <- function(object,
                              analysis_name,
                              significant_only = FALSE,
                              fc_cutoff = 1,
                              padj_cutoff = 0.05) {

  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }

  if (!analysis_name %in% names(object@allDiff)) {
    stop("Analysis '", analysis_name, "' not found. Available analyses: ",
         paste(names(object@allDiff), collapse = ", "), call. = FALSE)
  }

  result <- object@allDiff[[analysis_name]]

  if (significant_only) {
    result <- result[!is.na(result$adj.P.Val) &
                    result$adj.P.Val < padj_cutoff &
                    abs(result$logFC) > fc_cutoff, ]
  }

  return(result)
}

#' @title Get Analysis Summary
#' @description Get summary of all differential expression analyses in BulkRNAseq object
#'
#' @param object A BulkRNAseq object
#' @param fc_cutoff Numeric, fold change cutoff for counting significant genes
#' @param padj_cutoff Numeric, adjusted p-value cutoff for counting significant genes
#'
#' @return Data frame with analysis summary
#' @export
#'
#' @examples
#' \dontrun{
#' summary_df <- getAnalysisSummary(bulk_obj)
#' print(summary_df)
#' }
getAnalysisSummary <- function(object, fc_cutoff = 1, padj_cutoff = 0.05) {

  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }

  if (length(object@allDiff) == 0) {
    message("No differential expression analyses found in object")
    return(data.frame())
  }

  summary_list <- list()

  for (analysis_name in names(object@allDiff)) {
    result <- object@allDiff[[analysis_name]]

    n_total <- nrow(result)
    n_sig <- sum(!is.na(result$adj.P.Val) &
                result$adj.P.Val < padj_cutoff &
                abs(result$logFC) > fc_cutoff)
    n_up <- sum(!is.na(result$adj.P.Val) &
               result$adj.P.Val < padj_cutoff &
               result$logFC > fc_cutoff)
    n_down <- sum(!is.na(result$adj.P.Val) &
                 result$adj.P.Val < padj_cutoff &
                 result$logFC < -fc_cutoff)

    summary_list[[analysis_name]] <- data.frame(
      Analysis = analysis_name,
      Total_Genes = n_total,
      Significant_Genes = n_sig,
      Up_Regulated = n_up,
      Down_Regulated = n_down,
      FC_Cutoff = fc_cutoff,
      Padj_Cutoff = padj_cutoff,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  return(summary_df)
}
 