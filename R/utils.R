#' @title Get Counts Matrix
#' @description Extract the raw counts matrix from BulkRNAseq object
#' 
#' @param object A BulkRNAseq object
#' 
#' @return Matrix of raw counts
#' @export
getCounts <- function(object) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object必须是BulkRNAseq类型", call. = FALSE)
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
    stop("object必须是BulkRNAseq类型", call. = FALSE)
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
    stop("object必须是BulkRNAseq类型", call. = FALSE)
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
    stop("object必须是BulkRNAseq类型", call. = FALSE)
  }
  
  if (is.null(name)) {
    return(object@allDiff)
  } else {
    if (!name %in% names(object@allDiff)) {
      stop("差异分析结果中不存在名为 '", name, "' 的结果", call. = FALSE)
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
    stop("object必须是BulkRNAseq类型", call. = FALSE)
  }
  
  if (is.null(name)) {
    return(object@geneSet)
  } else {
    if (!name %in% names(object@geneSet)) {
      stop("基因集中不存在名为 '", name, "' 的基因集", call. = FALSE)
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
    stop("object必须是BulkRNAseq类型", call. = FALSE)
  }
  
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  # 验证维度一致性
  if (!identical(dim(data), dim(object@counts))) {
    stop("新数据矩阵的维度必须与原始counts矩阵一致", call. = FALSE)
  }
  
  # 验证行名和列名一致性
  if (!identical(rownames(data), rownames(object@counts)) ||
      !identical(colnames(data), colnames(object@counts))) {
    stop("新数据矩阵的行名和列名必须与原始counts矩阵一致", call. = FALSE)
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
    stop("object必须是BulkRNAseq类型", call. = FALSE)
  }
  
  if (missing(name) || is.null(name) || name == "") {
    stop("必须提供结果名称", call. = FALSE)
  }
  
  object@allDiff[[name]] <- results
  return(object)
}

#' @title Summary of BulkRNAseq Object
#' @description Provide a detailed summary of BulkRNAseq object contents
#' 
#' @param object A BulkRNAseq object
#' 
#' @return Invisibly returns a list with summary information
#' @export
summary.BulkRNAseq <- function(object) {
  if (!inherits(object, "BulkRNAseq")) {
    stop("object必须是BulkRNAseq类型", call. = FALSE)
  }
  
  # 收集摘要信息
  summary_info <- list(
    n_genes = nrow(object@counts),
    n_samples = ncol(object@counts),
    metadata_vars = colnames(object@metadata),
    diff_analyses = names(object@allDiff),
    gene_sets = names(object@geneSet),
    gsea_results = names(object@gsea),
    gsva_results = names(object@gsva)
  )
  
  # 打印摘要
  cat("BulkRNAseq对象摘要:\n")
  cat("==================\n")
  cat("基因数量:", summary_info$n_genes, "\n")
  cat("样本数量:", summary_info$n_samples, "\n")
  cat("元数据变量:", length(summary_info$metadata_vars), "个\n")
  cat("差异分析结果:", length(summary_info$diff_analyses), "个\n")
  cat("基因集:", length(summary_info$gene_sets), "个\n")
  cat("GSEA结果:", length(summary_info$gsea_results), "个\n")
  cat("GSVA结果:", length(summary_info$gsva_results), "个\n")
  
  invisible(summary_info)
} 