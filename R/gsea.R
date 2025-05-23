

#' @title Gene Set Enrichment Analysis (GSEA) for BulkRNAseq Objects
#' @description Perform Gene Set Enrichment Analysis using clusterProfiler
#' 
#' @param object A BulkRNAseq object containing differential expression results
#' @param geneSet_name Character vector, names of gene sets to analyze. Default "H" (Hallmark)
#' @param analysis_name Character, name of differential analysis results to use. 
#'   If NULL, uses the first available analysis
#' @param logfc_col Character, column name for log fold change values. Default "logFC"
#' @param geneid_col Character, column name for gene identifiers. Default "gene_id"
#' @param pvalueCutoff Numeric, p-value cutoff for significance. Default 0.05
#' @param pAdjustMethod Character, p-value adjustment method. Default "BH"
#' @param minGSSize Integer, minimum gene set size. Default 15
#' @param maxGSSize Integer, maximum gene set size. Default 500
#' @param verbose Logical, whether to print progress messages. Default TRUE
#' 
#' @return Updated BulkRNAseq object with GSEA results stored in gsea slot
#' 
#' @examples
#' \dontrun{
#' # Basic GSEA analysis
#' bulk_obj <- gsea(bulk_obj, geneSet_name = "H")
#' 
#' # Multiple gene sets
#' bulk_obj <- gsea(bulk_obj, geneSet_name = c("H", "C2_CP_KEGG"))
#' 
#' # Custom parameters
#' bulk_obj <- gsea(bulk_obj, 
#'                  geneSet_name = "H",
#'                  pvalueCutoff = 0.01,
#'                  minGSSize = 10)
#' }
#' 
#' @export
setGeneric("gsea", function(object, 
                           geneSet_name = "H", 
                           analysis_name = NULL,
                           logfc_col = "logFC", 
                           geneid_col = "gene_id",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           minGSSize = 15,
                           maxGSSize = 500,
                           verbose = TRUE)
  standardGeneric("gsea"))

#' @rdname gsea
#' @importFrom methods setGeneric setMethod
#' @importFrom stats setNames complete.cases
#' @import clusterProfiler  
#' @export
setMethod("gsea", "BulkRNAseq", function(object, 
                                         geneSet_name = "H",
                                         analysis_name = NULL,
                                         logfc_col = "logFC", 
                                         geneid_col = "gene_id",
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH",
                                         minGSSize = 15,
                                         maxGSSize = 500,
                                         verbose = TRUE) {
  
  # Input validation
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  
  # Check required packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required for GSEA analysis. ",
         "Install with: BiocManager::install('clusterProfiler')", call. = FALSE)
  }
  
  # Validate parameters
  if (!is.numeric(pvalueCutoff) || pvalueCutoff <= 0 || pvalueCutoff > 1) {
    stop("pvalueCutoff must be a number between 0 and 1", call. = FALSE)
  }
  
  if (!is.numeric(minGSSize) || minGSSize < 1) {
    stop("minGSSize must be a positive integer", call. = FALSE)
  }
  
  if (!is.numeric(maxGSSize) || maxGSSize < minGSSize) {
    stop("maxGSSize must be greater than or equal to minGSSize", call. = FALSE)
  }
  
  # Handle multiple gene sets
  if (length(geneSet_name) > 1) {
    if (verbose) {
      message("Processing ", length(geneSet_name), " gene set collections...")
    }
    
    for (gs_name in geneSet_name) {
      if (verbose) {
        message("Processing gene set collection: ", gs_name)
      }
      object <- gsea(object, gs_name, analysis_name, logfc_col, geneid_col,
                    pvalueCutoff, pAdjustMethod, minGSSize, maxGSSize, verbose)
    }
    return(object)
  }
  
  # Get differential expression results
  diff_results <- get_diff_results(object, analysis_name, verbose)
  
  # Validate required columns
  if (!logfc_col %in% colnames(diff_results)) {
    stop("Column '", logfc_col, "' not found in differential expression results", call. = FALSE)
  }
  
  if (!geneid_col %in% colnames(diff_results)) {
    stop("Column '", geneid_col, "' not found in differential expression results", call. = FALSE)
  }
  
  # Prepare gene ranking list
  gene_order <- prepare_gene_ranking(diff_results, logfc_col, geneid_col, verbose)
  
  # Get gene sets
  geneSet <- get_gene_sets(object, geneSet_name, verbose)
  
  # Perform GSEA analysis
  if (verbose) {
    message("Running GSEA analysis for gene set collection: ", geneSet_name)
  }
  
  gsea_result <- tryCatch({
    clusterProfiler::GSEA(
      geneList = gene_order,
      TERM2GENE = geneSet,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      verbose = verbose
    )
  }, error = function(e) {
    stop("GSEA analysis failed: ", e$message, call. = FALSE)
  })
  
  # Store results
  object@gsea[[geneSet_name]] <- gsea_result
  
  if (verbose) {
    n_significant <- sum(gsea_result@result$p.adjust < pvalueCutoff, na.rm = TRUE)
    message("GSEA completed. Found ", n_significant, " significantly enriched gene sets")
  }
  
  return(object)
})

# Helper function to get differential expression results
get_diff_results <- function(object, analysis_name, verbose) {
  if (length(object@allDiff) == 0) {
    stop("No differential expression results found. Please run differential analysis first.", call. = FALSE)
  }
  
  if (is.null(analysis_name)) {
    analysis_name <- names(object@allDiff)[1]
    if (verbose) {
      message("Using differential analysis results: ", analysis_name)
    }
  } else {
    if (!analysis_name %in% names(object@allDiff)) {
      stop("Analysis '", analysis_name, "' not found. Available analyses: ",
           paste(names(object@allDiff), collapse = ", "), call. = FALSE)
    }
  }
  
  return(object@allDiff[[analysis_name]])
}

# Helper function to prepare gene ranking
prepare_gene_ranking <- function(diff_results, logfc_col, geneid_col, verbose) {
  # Remove rows with missing values
  complete_rows <- complete.cases(diff_results[, c(logfc_col, geneid_col)])
  if (sum(!complete_rows) > 0) {
    if (verbose) {
      message("Removing ", sum(!complete_rows), " genes with missing logFC or gene ID values")
    }
    diff_results <- diff_results[complete_rows, ]
  }
  
  # Create named vector and sort by logFC
  gene_order <- stats::setNames(diff_results[[logfc_col]], diff_results[[geneid_col]])
  gene_order <- sort(gene_order, decreasing = TRUE)
  
  # Remove duplicated gene names (keep the one with highest logFC)
  gene_order <- gene_order[!duplicated(names(gene_order))]
  
  if (verbose) {
    message("Prepared gene ranking list with ", length(gene_order), " genes")
  }
  
  return(gene_order)
}

# Helper function to get gene sets
get_gene_sets <- function(object, geneSet_name, verbose) {
  if (!geneSet_name %in% names(object@geneSet)) {
    available_sets <- names(object@geneSet)
    stop("Gene set collection '", geneSet_name, "' not found. ",
         "Available collections: ", paste(available_sets, collapse = ", "), call. = FALSE)
  }
  
  geneSet <- object@geneSet[[geneSet_name]]
  
  # Validate gene set structure
  required_cols <- c("gs_name", "gene_symbol")
  if (!all(required_cols %in% colnames(geneSet))) {
    stop("Gene set must contain columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }
  
  # Remove rows with missing values
  complete_rows <- complete.cases(geneSet[, required_cols])
  if (sum(!complete_rows) > 0) {
    if (verbose) {
      message("Removing ", sum(!complete_rows), " gene set entries with missing values")
    }
    geneSet <- geneSet[complete_rows, ]
  }
  
  if (verbose) {
    n_sets <- length(unique(geneSet$gs_name))
    n_genes <- length(unique(geneSet$gene_symbol))
    message("Using gene set collection '", geneSet_name, "' with ", n_sets, 
            " gene sets and ", n_genes, " unique genes")
  }
  
  return(geneSet[, required_cols])
}