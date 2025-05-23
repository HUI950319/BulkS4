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
    stop("\u5217\u540d '", name, "' \u5728metadata\u4e2d\u4e0d\u5b58\u5728", call. = FALSE)
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
  # \u521b\u5efa\u5206\u9694\u7ebf
  sep_line <- paste(rep("-", 50), collapse = "")
  
  # \u6253\u5370\u6807\u9898\u548c\u57fa\u672c\u4fe1\u606f
  cat("\n", cli::col_blue("====== BulkRNAseq\u5bf9\u8c61 ======"), "\n\n")
  
  # \u6253\u5370\u57fa\u672c\u7edc\u8ba1\u4fe1\u606f
  cat(sprintf("%-15s %d\n", "\u57fa\u56e0\u6570\u91cf:", nrow(object@counts)))
  cat(sprintf("%-15s %d\n", "\u6837\u672c\u6570\u91cf:", ncol(object@counts)))
  
  # \u6253\u5370\u5143\u6570\u636e\u4fe1\u606f
  meta_vars <- if (ncol(object@metadata) > 0) {
    paste(colnames(object@metadata), collapse = ", ")
  } else {
    "\u65e0"
  }
  cat(sprintf("%-15s %s\n", "\u5143\u6570\u636e\u53d8\u91cf:", meta_vars))
  
  # \u6253\u5370\u5206\u6790\u7ed3\u679c\u4fe1\u606f
  cat("\n", cli::col_blue("------ \u5206\u6790\u7ed3\u679c ------"), "\n\n")
  
  # \u5dee\u5f02\u5206\u6790\u7ed3\u679c
  diff_count <- length(object@allDiff)
  diff_names <- if (diff_count > 0) {
    paste(names(object@allDiff), collapse = ", ")
  } else {
    "\u65e0"
  }
  cat(sprintf("%-15s %d\u4e2a: %s\n", "\u5dee\u5f02\u5206\u6790\u7ed3\u679c:", diff_count, diff_names))
  
  # \u57fa\u56e0\u96c6\u4fe1\u606f
  geneset_count <- length(object@geneSet)
  geneset_names <- if (geneset_count > 0) {
    paste(names(object@geneSet), collapse = ", ")
  } else {
    "\u65e0"
  }
  cat(sprintf("%-15s %d\u4e2a: %s\n", "\u57fa\u56e0\u96c6:", geneset_count, geneset_names))
  
  # \u5bcc\u96c6\u5206\u6790\u7ed3\u679c
  gsea_count <- length(object@gsea)
  gsea_names <- if (gsea_count > 0) {
    paste(names(object@gsea), collapse = ", ")
  } else {
    "\u65e0"
  }
  cat(sprintf("%-15s %d\u4e2a: %s\n", "GSEA\u7ed3\u679c:", gsea_count, gsea_names))
  
  # GSVA\u7ed3\u679c
  gsva_count <- length(object@gsva)
  gsva_names <- if (gsva_count > 0) {
    paste(names(object@gsva), collapse = ", ")
  } else {
    "\u65e0"
  }
  cat(sprintf("%-15s %d\u4e2a: %s\n", "GSVA\u7ed3\u679c:", gsva_count, gsva_names))
  
  # \u7ed3\u675f\u5206\u9694\u7ebf
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
  # \u5904\u7406\u7f3a\u5931\u7684\u7d22\u5f15
  if (missing(i)) i <- seq_len(nrow(x@counts))
  if (missing(j)) j <- seq_len(ncol(x@counts))
  
  # \u5b50\u96c6\u5316\u6570\u636e
  new_counts <- x@counts[i, j, drop = FALSE]
  new_data <- x@data[i, j, drop = FALSE]
  new_metadata <- x@metadata[j, , drop = FALSE]
  
  # \u521b\u5efa\u65b0\u5bf9\u8c61
  new("BulkRNAseq",
      counts = new_counts,
      data = new_data,
      metadata = new_metadata,
      allDiff = x@allDiff,  # \u4fdd\u7559\u539f\u6709\u5206\u6790\u7ed3\u679c
      geneSet = x@geneSet,
      gsea = x@gsea,
      gsva = x@gsva)
}) 