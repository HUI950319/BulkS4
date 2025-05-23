#' @title BulkRNAseq S4 Class
#' @description An S4 class for storing and analyzing bulk RNA-seq data
#' 
#' @slot counts matrix. Raw count matrix (genes x samples)
#' @slot data matrix. Normalized data matrix (genes x samples)  
#' @slot metadata data.frame. Sample metadata (rownames are sample names)
#' @slot allDiff list. Differential analysis results (list of data.frames)
#' @slot geneSet list. Gene sets for GSEA/GSVA analysis (multiple data.frames)
#' @slot gsea list. GSEA/enrichment analysis results (data.frames)
#' @slot gsva list. GSVA/ssGSEA analysis results (pathways x samples data.frames)
#' 
#' @exportClass BulkRNAseq
setClass(
  "BulkRNAseq",
  slots = list(
    counts = "matrix",           # 原始计数矩阵（基因 x 样本）
    data = "matrix",             # 标准化后的数据矩阵（基因 x 样本）
    metadata = "data.frame",     # 样本元数据（行名为样本名）
    allDiff = "list",            # 差异分析的基因集（list of data.frame）
    geneSet = "list",            # GSEA/GSVA的基因集（多个data.frame）
    gsea = "list",               # enricher/GSEA的结果（data.frame）
    gsva = "list"                # GSVA/ssGSEA的结果（通路 x 样本 data.frame）
  ),
  validity = function(object) {
    errors <- character()
    
    # 检查counts和data矩阵维度一致性
    if (!identical(dim(object@counts), dim(object@data))) {
      errors <- c(errors, "counts和data矩阵维度必须一致")
    }
    
    # 检查metadata行数与矩阵列数一致性
    if (nrow(object@metadata) != ncol(object@counts)) {
      errors <- c(errors, "metadata行数必须与counts列数匹配")
    }
    
    # 检查样本名一致性
    if (!identical(rownames(object@metadata), colnames(object@counts))) {
      errors <- c(errors, "metadata行名必须与counts列名一致")
    }
    
    if (length(errors) == 0) TRUE else errors
  }
)

#' @title Create BulkRNAseq Object
#' @description Constructor function for BulkRNAseq S4 objects
#' 
#' @param counts A matrix or data.frame of raw counts, or a list containing 
#'   required fields (counts, metadata, exprSet, allDiff)
#' @param metadata A data.frame with sample metadata. If NULL, an empty 
#'   data.frame will be created
#' @param add_msig.geneSet Logical. Whether to add MSigDB gene sets to the object
#' 
#' @return A BulkRNAseq S4 object
#' @export
#' 
#' @examples
#' # Create from matrix
#' counts_mat <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' rownames(counts_mat) <- paste0("Gene", 1:100)
#' colnames(counts_mat) <- paste0("Sample", 1:10)
#' 
#' metadata <- data.frame(
#'   condition = rep(c("Control", "Treatment"), each = 5),
#'   row.names = colnames(counts_mat)
#' )
#' 
#' bulk_obj <- BulkRNAseq(counts_mat, metadata)
BulkRNAseq <- function(counts, metadata = NULL, add_msig.geneSet = TRUE) {
  
  # 输入验证辅助函数
  .validate_matrix <- function(x, name) {
    if (!is.matrix(x) && !is.data.frame(x)) {
      stop(sprintf("%s必须是矩阵或数据框", name), call. = FALSE)
    }
    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    return(x)
  }
  
  .validate_metadata <- function(meta, n_samples, sample_names) {
    if (is.null(meta)) {
      meta <- data.frame(row.names = sample_names)
    } else {
      if (!is.data.frame(meta)) {
        stop("metadata必须为数据框", call. = FALSE)
      }
      if (nrow(meta) != n_samples) {
        stop(sprintf("metadata行数(%d)必须与样本数(%d)匹配", 
                    nrow(meta), n_samples), call. = FALSE)
      }
    }
    return(meta)
  }
  
  # 处理列表输入
  if (inherits(counts, "list")) {
    required_fields <- c("counts", "metadata", "exprSet", "allDiff")
    missing_fields <- setdiff(required_fields, names(counts))
    
    if (length(missing_fields) > 0) {
      stop("输入列表缺少必要字段: ", paste(missing_fields, collapse = ", "), 
           call. = FALSE)
    }
    
    # 提取并验证字段
    counts_mat <- .validate_matrix(counts$counts, "counts")
    expr_set <- .validate_matrix(counts$exprSet, "exprSet")
    meta_data <- .validate_metadata(counts$metadata, ncol(counts_mat), 
                                   colnames(counts_mat))
    all_diff <- counts$allDiff
    
    # 验证exprSet与counts的一致性
    if (!identical(dim(expr_set), dim(counts_mat))) {
      stop("exprSet必须与counts有相同的维度", call. = FALSE)
    }
    
    if (!identical(colnames(expr_set), colnames(counts_mat))) {
      stop("exprSet必须与counts有相同的列名", call. = FALSE)
    }
    
    # 确保行名一致
    rownames(meta_data) <- colnames(counts_mat)
    
    # 创建对象
    object <- new("BulkRNAseq",
                  counts = counts_mat,
                  data = expr_set,
                  metadata = meta_data,
                  allDiff = list(allDiff = all_diff),
                  geneSet = list(),
                  gsea = list(),
                  gsva = list())
    
  } else {
    # 处理矩阵/数据框输入
    counts_mat <- .validate_matrix(counts, "counts")
    meta_data <- .validate_metadata(metadata, ncol(counts_mat), 
                                   colnames(counts_mat))
    
    # 确保行名一致
    rownames(meta_data) <- colnames(counts_mat)
    
    # 创建对象
    object <- new("BulkRNAseq",
                  counts = counts_mat,
                  data = counts_mat,  # 初始时data与counts相同
                  metadata = meta_data,
                  allDiff = list(),
                  geneSet = list(),
                  gsea = list(),
                  gsva = list())
  }
  
  # 添加MSigDB基因集
  if (add_msig.geneSet) {
    object <- .add_msigdb_genesets(object)
  }
  
  return(object)
}

#' @title Add MSigDB Gene Sets
#' @description Internal function to add MSigDB gene sets to BulkRNAseq object
#' @param object A BulkRNAseq object
#' @return BulkRNAseq object with gene sets added
#' @keywords internal
.add_msigdb_genesets <- function(object) {
  # 根据操作系统确定路径
  geneSet_path <- if (Sys.info()["sysname"] == "Linux") {
    "/mnt/e/0.scRNA/guozi/0.resource/geneSet_msig.RData"
  } else {
    "E:/0.scRNA/guozi/0.resource/geneSet_msig.RData"
  }
  
  # 检查文件是否存在并加载
  if (file.exists(geneSet_path)) {
    tryCatch({
      load(geneSet_path, envir = environment())
      if (exists("geneSet_msig", envir = environment())) {
        object@geneSet <- get("geneSet_msig", envir = environment())
        message("成功加载MSigDB基因集")
      } else {
        warning("基因集文件中未找到geneSet_msig对象", call. = FALSE)
      }
    }, error = function(e) {
      warning("加载基因集文件时出错: ", e$message, call. = FALSE)
    })
  } else {
    warning("无法找到基因集文件: ", geneSet_path, call. = FALSE)
  }
  
  return(object)
} 