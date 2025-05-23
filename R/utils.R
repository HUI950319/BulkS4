#' @importFrom methods is

#' @title Get Counts Matrix
#' @description Extract the raw counts matrix from BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' 
#' @return Matrix of raw counts
#' @export
getCounts <- function(object) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  return(object@counts)
}

#' @title Get Normalized Data Matrix
#' @description Extract the normalized data matrix from BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' 
#' @return Matrix of normalized data
#' @export
getData <- function(object) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  return(object@data)
}

#' @title Get Metadata
#' @description Extract the metadata from BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' 
#' @return Data.frame of sample metadata
#' @export
getMetadata <- function(object) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  return(object@metadata)
}

#' @title Get Differential Analysis Results
#' @description Extract differential analysis results from BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' @param name Character. Name of the differential analysis result to extract.
#'   If NULL, returns all results.
#' 
#' @return List or data.frame of differential analysis results
#' @export
getDiffResults <- function(object, name = NULL) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  
  if (is.null(name)) {
    return(object@allDiff)
  } else {
    if (!name %in% names(object@allDiff)) {
      stop("Differential analysis result '", name, "' not found", call. = FALSE)
    }
    return(object@allDiff[[name]])
  }
}

#' @title Get Gene Sets
#' @description Extract gene sets from BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' @param name Character. Name of the gene set to extract.
#'   If NULL, returns all gene sets.
#' 
#' @return List or data.frame of gene sets
#' @export
getGeneSets <- function(object, name = NULL) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  
  if (is.null(name)) {
    return(object@geneSet)
  } else {
    if (!name %in% names(object@geneSet)) {
      stop("Gene set '", name, "' not found", call. = FALSE)
    }
    return(object@geneSet[[name]])
  }
}

#' @title Set Normalized Data
#' @description Set the normalized data matrix in BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' @param data Matrix of normalized data
#' 
#' @return Updated BulkRNAseq object
#' @export
setData <- function(object, data) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  # Validate dimension consistency
  if (!identical(dim(data), dim(object@counts))) {
    stop("New data matrix dimensions must match original counts matrix", call. = FALSE)
  }
  
  # Validate row and column names consistency
  if (!identical(rownames(data), rownames(object@counts)) ||
      !identical(colnames(data), colnames(object@counts))) {
    stop("New data matrix row and column names must match original counts matrix", call. = FALSE)
  }
  
  object@data <- data
  return(object)
}

#' @title Add Differential Analysis Results
#' @description Add differential analysis results to BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' @param results Data.frame or list of differential analysis results
#' @param name Character. Name for the differential analysis results
#' 
#' @return Updated BulkRNAseq object
#' @export
addDiffResults <- function(object, results, name) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  
  if (missing(name) || is.null(name) || name == "") {
    stop("Result name must be provided", call. = FALSE)
  }
  
  object@allDiff[[name]] <- results
  return(object)
}

#' @title Summary of BulkRNAseq Object
#' @description Provide a detailed summary of BulkRNAseq object contents
#' 
#' @param object A BulkRNAseq object
#' @param ... Additional arguments (for S3 method compatibility)
#' 
#' @return Invisibly returns a list with summary information
#' @export
summary.BulkRNAseq <- function(object, ...) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  
  # Collect summary information
  summary_info <- list(
    n_genes = nrow(object@counts),
    n_samples = ncol(object@counts),
    metadata_vars = colnames(object@metadata),
    diff_analyses = names(object@allDiff),
    gene_sets = names(object@geneSet),
    gsea_results = names(object@gsea),
    gsva_results = names(object@gsva)
  )
  
  # Print summary
  cat("BulkRNAseq Object Summary:\n")
  cat("========================\n")
  cat("Genes:", summary_info$n_genes, "\n")
  cat("Samples:", summary_info$n_samples, "\n")
  cat("Metadata variables:", length(summary_info$metadata_vars), "\n")
  cat("Differential analyses:", length(summary_info$diff_analyses), "\n")
  cat("Gene sets:", length(summary_info$gene_sets), "\n")
  cat("GSEA results:", length(summary_info$gsea_results), "\n")
  cat("GSVA results:", length(summary_info$gsva_results), "\n")
  
  invisible(summary_info)
} 