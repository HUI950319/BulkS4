#' @importFrom methods new setClass
#' @importFrom cli col_blue

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
    counts = "matrix",           # Raw count matrix (genes x samples)
    data = "matrix",             # Normalized data matrix (genes x samples)
    metadata = "data.frame",     # Sample metadata (rownames are sample names)
    allDiff = "list",            # Differential analysis gene sets (list of data.frame)
    geneSet = "list",            # Gene sets for GSEA/GSVA (multiple data.frame)
    gsea = "list",               # enricher/GSEA results (data.frame)
    gsva = "list"                # GSVA/ssGSEA results (pathways x samples data.frame)
  ),
  validity = function(object) {
    errors <- character()

    # # Check counts and data matrix dimension consistency
    # if (!identical(dim(object@counts), dim(object@data))) {
    #   errors <- c(errors, "counts and data matrix dimensions must be consistent")
    # }

    # Check metadata rows match matrix columns
    if (nrow(object@metadata) != ncol(object@counts)) {
      errors <- c(errors, "metadata rows must match counts columns")
    }

    # Check sample name consistency
    if (!identical(rownames(object@metadata), colnames(object@counts))) {
      errors <- c(errors, "metadata rownames must match counts column names")
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

  # Input validation helper functions
  .validate_matrix <- function(x, name) {
    if (!is.matrix(x) && !is.data.frame(x)) {
      stop(sprintf("%s must be a matrix or data frame", name), call. = FALSE)
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
        stop("metadata must be a data frame", call. = FALSE)
      }
      if (nrow(meta) != n_samples) {
        stop(sprintf("metadata rows (%d) must match sample count (%d)",
                    nrow(meta), n_samples), call. = FALSE)
      }
    }
    return(meta)
  }

  # Handle list input
  if (inherits(counts, "list")) {
    required_fields <- c("counts", "metadata", "exprSet", "allDiff")
    missing_fields <- setdiff(required_fields, names(counts))

    if (length(missing_fields) > 0) {
      stop("Input list missing required fields: ", paste(missing_fields, collapse = ", "),
           call. = FALSE)
    }

    # Extract and validate fields
    counts_mat <- .validate_matrix(counts$counts, "counts")
    expr_set <- .validate_matrix(counts$exprSet, "exprSet")
    meta_data <- .validate_metadata(counts$metadata, ncol(counts_mat),
                                   colnames(counts_mat))
    all_diff <- counts$allDiff

    # # Validate exprSet and counts consistency---allow inconsistency, gene deduplication
    # if (!identical(dim(expr_set), dim(counts_mat))) {
    #   stop("exprSet must have same dimensions as counts", call. = FALSE)
    # }

    if (!identical(colnames(expr_set), colnames(counts_mat))) {
      stop("exprSet must have same column names as counts", call. = FALSE)
    }

    # Ensure consistent rownames
    rownames(meta_data) <- colnames(counts_mat)

    # Create object
    object <- new("BulkRNAseq",
                  counts = counts_mat,
                  data = expr_set,
                  metadata = meta_data,
                  allDiff = list(allDiff = all_diff),
                  geneSet = list(),
                  gsea = list(),
                  gsva = list())

  } else {
    # Handle matrix/data.frame input
    counts_mat <- .validate_matrix(counts, "counts")
    meta_data <- .validate_metadata(metadata, ncol(counts_mat),
                                   colnames(counts_mat))

    # Ensure consistent rownames
    rownames(meta_data) <- colnames(counts_mat)

    # Create object
    object <- new("BulkRNAseq",
                  counts = counts_mat,
                  data = counts_mat,  # Initially data same as counts
                  metadata = meta_data,
                  allDiff = list(),
                  geneSet = list(),
                  gsea = list(),
                  gsva = list())
  }

  # Add MSigDB gene sets
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
  tryCatch({
    # Load the geneSet_msig dataset from the package data
    # This will load the geneSet_msig object into the current environment
    data("geneSet_msig", envir = environment())

    # Check if the object was loaded successfully
    if (exists("geneSet_msig", envir = environment())) {
      object@geneSet <- get("geneSet_msig", envir = environment())
      message("Successfully loaded MSigDB gene sets from package data")
    } else {
      warning("geneSet_msig object not found in package data", call. = FALSE)
    }
  }, error = function(e) {
    warning("Error loading geneSet_msig from package data: ", e$message, call. = FALSE)
})

  return(object)
}
