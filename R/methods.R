#' @importFrom methods setMethod new
#' @importFrom cli col_blue

#' @title Extract Metadata Columns
#' @description Extract columns from metadata using $ operator
#' 
#' @param x A BulkRNAseq object
#' @param name Character. Column name in metadata
#' 
#' @return Vector of values from the specified metadata column
#' @export
#' 
#' @examples
#' # Assuming bulk_obj is a BulkRNAseq object with metadata
#' # bulk_obj$condition  # Extract condition column
setMethod("$", "BulkRNAseq", function(x, name) {
  if (!name %in% colnames(x@metadata)) {
    stop("Column '", name, "' not found in metadata", call. = FALSE)
  }
  return(x@metadata[[name]])
})

#' @title Show BulkRNAseq Object
#' @description Display method for BulkRNAseq objects
#' 
#' @param object A BulkRNAseq object
#' 
#' @return Invisibly returns the object (for method chaining)
#' @export
setMethod("show", "BulkRNAseq", function(object) {
  # Create separator line
  sep_line <- paste(rep("-", 50), collapse = "")
  
  # Print title and basic information
  cat("\n", cli::col_blue("====== BulkRNAseq Object ======"), "\n\n")
  
  # Print basic statistics
  cat(sprintf("%-15s %d\n", "Genes:", nrow(object@counts)))
  cat(sprintf("%-15s %d\n", "Samples:", ncol(object@counts)))
  
  # Print metadata information
  meta_vars <- if (ncol(object@metadata) > 0) {
    paste(colnames(object@metadata), collapse = ", ")
  } else {
    "None"
  }
  cat(sprintf("%-15s %s\n", "Metadata vars:", meta_vars))
  
  # Print analysis results information
  cat("\n", cli::col_blue("------ Analysis Results ------"), "\n\n")
  
  # Differential analysis results
  diff_count <- length(object@allDiff)
  diff_names <- if (diff_count > 0) {
    paste(names(object@allDiff), collapse = ", ")
  } else {
    "None"
  }
  cat(sprintf("%-15s %d: %s\n", "Diff analyses:", diff_count, diff_names))
  
  # Gene set information
  geneset_count <- length(object@geneSet)
  geneset_names <- if (geneset_count > 0) {
    paste(names(object@geneSet), collapse = ", ")
  } else {
    "None"
  }
  cat(sprintf("%-15s %d: %s\n", "Gene sets:", geneset_count, geneset_names))
  
  # Enrichment analysis results
  gsea_count <- length(object@gsea)
  gsea_names <- if (gsea_count > 0) {
    paste(names(object@gsea), collapse = ", ")
  } else {
    "None"
  }
  cat(sprintf("%-15s %d: %s\n", "GSEA results:", gsea_count, gsea_names))
  
  # GSVA results
  gsva_count <- length(object@gsva)
  gsva_names <- if (gsva_count > 0) {
    paste(names(object@gsva), collapse = ", ")
  } else {
    "None"
  }
  cat(sprintf("%-15s %d: %s\n", "GSVA results:", gsva_count, gsva_names))
  
  # End separator line
  cat("\n", sep_line, "\n")
  
  invisible(object)
})

#' @title Get Dimensions
#' @description Get dimensions of the count matrix
#' 
#' @param x A BulkRNAseq object
#' 
#' @return Integer vector of length 2 (nrow, ncol)
#' @export
setMethod("dim", "BulkRNAseq", function(x) {
  dim(x@counts)
})

#' @title Get Row Names
#' @description Get gene names (row names of count matrix)
#' 
#' @param x A BulkRNAseq object
#' 
#' @return Character vector of gene names
#' @export
setMethod("rownames", "BulkRNAseq", function(x) {
  rownames(x@counts)
})

#' @title Get Column Names  
#' @description Get sample names (column names of count matrix)
#' 
#' @param x A BulkRNAseq object
#' 
#' @return Character vector of sample names
#' @export
setMethod("colnames", "BulkRNAseq", function(x) {
  colnames(x@counts)
})

#' @title Subset BulkRNAseq Object
#' @description Subset BulkRNAseq object by genes and/or samples
#' 
#' @param x A BulkRNAseq object
#' @param i Row indices (genes) to subset
#' @param j Column indices (samples) to subset
#' @param drop Logical. Whether to drop dimensions (not used, for compatibility)
#' 
#' @return A subsetted BulkRNAseq object
#' @export
setMethod("[", "BulkRNAseq", function(x, i, j, drop = FALSE) {
  # Handle missing indices
  if (missing(i)) i <- seq_len(nrow(x@counts))
  if (missing(j)) j <- seq_len(ncol(x@counts))
  
  # Subset data
  new_counts <- x@counts[i, j, drop = FALSE]
  new_data <- x@data[i, j, drop = FALSE]
  new_metadata <- x@metadata[j, , drop = FALSE]
  
  # Create new object
  new("BulkRNAseq",
      counts = new_counts,
      data = new_data,
      metadata = new_metadata,
      allDiff = x@allDiff,  # Keep original analysis results
      geneSet = x@geneSet,
      gsea = x@gsea,
      gsva = x@gsva)
}) 