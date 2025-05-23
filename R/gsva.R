#' @title Gene Set Variation Analysis (GSVA) for BulkRNAseq Objects
#' @description Perform Gene Set Variation Analysis using GSVA package
#' 
#' @param object A BulkRNAseq object containing expression data
#' @param geneSet_name Character vector, names of gene sets to analyze. Default "H" (Hallmark)
#' @param gsva_method Character, GSVA method to use. Options: "gsva", "ssgsea". Default "gsva"
#' @param kcdf Character, kernel to use for density estimation. 
#'   Options: "Gaussian", "Poisson", "none". Default "Gaussian"
#' @param min_sz Integer, minimum gene set size. Default 1
#' @param max_sz Integer, maximum gene set size. Default Inf
#' @param mx_diff Logical, whether to use max difference. Default TRUE
#' @param tau Numeric, tau parameter for ssGSEA. Default 1
#' @param ssgsea_norm Logical, whether to normalize ssGSEA scores. Default TRUE
#' @param run_limma Logical, whether to run limma differential analysis. Default TRUE
#' @param verbose Logical, whether to print progress messages. Default TRUE
#' 
#' @return Updated BulkRNAseq object with GSVA results stored in gsva slot
#' 
#' @examples
#' \dontrun{
#' # Basic GSVA analysis
#' bulk_obj <- gsva(bulk_obj, geneSet_name = "H")
#' 
#' # ssGSEA analysis
#' bulk_obj <- gsva(bulk_obj, geneSet_name = "H", gsva_method = "ssgsea")
#' 
#' # Multiple gene sets
#' bulk_obj <- gsva(bulk_obj, geneSet_name = c("H", "C2_CP_KEGG"))
#' 
#' # Custom parameters
#' bulk_obj <- gsva(bulk_obj, 
#'                  geneSet_name = "H",
#'                  kcdf = "Poisson",
#'                  min_sz = 10,
#'                  max_sz = 200)
#' }
#' 
#' @export
setGeneric("gsva", function(object, 
                           geneSet_name = "H", 
                           gsva_method = c("gsva", "ssgsea"),
                           kcdf = "Gaussian",
                           min_sz = 1,
                           max_sz = Inf,
                           mx_diff = TRUE,
                           tau = 1,
                           ssgsea_norm = TRUE,
                           run_limma = TRUE,
                           verbose = TRUE)
  standardGeneric("gsva"))

#' @rdname gsva
#' @import GSVA
#' @import limma
#' @importFrom stats complete.cases
#' @export
setMethod("gsva", "BulkRNAseq", function(object, 
                                         geneSet_name = "H", 
                                         gsva_method = c("gsva", "ssgsea"),
                                         kcdf = "Gaussian",
                                         min_sz = 1,
                                         max_sz = Inf,
                                         mx_diff = TRUE,
                                         tau = 1,
                                         ssgsea_norm = TRUE,
                                         run_limma = TRUE,
                                         verbose = TRUE) {
  
  # Input validation
  if (!inherits(object, "BulkRNAseq")) {
    stop("object must be a BulkRNAseq object", call. = FALSE)
  }
  
  # Check required packages
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package 'GSVA' is required for GSVA analysis. ",
         "Install with: BiocManager::install('GSVA')", call. = FALSE)
  }
  
  if (run_limma && !requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for differential analysis. ",
         "Install with: BiocManager::install('limma')", call. = FALSE)
  }
  
  # Validate method
  gsva_method <- match.arg(gsva_method)
  
  # Validate kcdf
  valid_kcdf <- c("Gaussian", "Poisson", "none")
  if (!kcdf %in% valid_kcdf) {
    stop("kcdf must be one of: ", paste(valid_kcdf, collapse = ", "), call. = FALSE)
  }
  
  # Validate parameters
  if (!is.numeric(min_sz) || min_sz < 1) {
    stop("min_sz must be a positive integer", call. = FALSE)
  }
  
  if (!is.numeric(max_sz) || max_sz < min_sz) {
    stop("max_sz must be greater than or equal to min_sz", call. = FALSE)
  }
  
  if (!is.numeric(tau) || tau <= 0) {
    stop("tau must be a positive number", call. = FALSE)
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
      object <- gsva(object, gs_name, gsva_method, kcdf, min_sz, max_sz, 
                    mx_diff, tau, ssgsea_norm, run_limma, verbose)
    }
    return(object)
  }
  
  # Get expression data
  expr_data <- get_expression_data(object, verbose)
  
  # Get and process gene sets
  geneSet_list <- get_and_process_gene_sets(object, geneSet_name, verbose)
  
  # Perform GSVA analysis
  if (verbose) {
    message("Running ", toupper(gsva_method), " analysis for gene set collection: ", geneSet_name)
  }
  
  gsva_scores <- perform_gsva_analysis(expr_data, geneSet_list, gsva_method, 
                                      kcdf, min_sz, max_sz, mx_diff, tau, 
                                      ssgsea_norm, verbose)
  
  # Store GSVA results
  object@gsva[[geneSet_name]] <- as.data.frame(gsva_scores)
  
  # Perform differential analysis if requested
  if (run_limma) {
    diff_results <- perform_limma_analysis(object, gsva_scores, geneSet_name, verbose)
    if (!is.null(diff_results)) {
      object@allDiff[[geneSet_name]] <- diff_results
    }
  }
  
  if (verbose) {
    message("GSVA analysis completed for gene set collection: ", geneSet_name)
  }
  
  return(object)
})

# Helper function to get expression data
get_expression_data <- function(object, verbose) {
  if (is.null(object@data) || ncol(object@data) == 0) {
    stop("No expression data found in object. Please ensure data slot is populated.", call. = FALSE)
  }
  
  expr_data <- object@data
  
  # Check for missing values
  if (any(is.na(expr_data))) {
    n_missing <- sum(is.na(expr_data))
    if (verbose) {
      message("Warning: Found ", n_missing, " missing values in expression data")
    }
  }
  
  if (verbose) {
    message("Using expression data with ", nrow(expr_data), " genes and ", 
            ncol(expr_data), " samples")
  }
  
  return(expr_data)
}

# Helper function to get and process gene sets
get_and_process_gene_sets <- function(object, geneSet_name, verbose) {
  # Check if gene set exists
  if (!geneSet_name %in% names(object@geneSet)) {
    available_sets <- names(object@geneSet)
    stop("Gene set collection '", geneSet_name, "' not found. ",
         "Available collections: ", paste(available_sets, collapse = ", "), call. = FALSE)
  }
  
  geneSet_df <- object@geneSet[[geneSet_name]]
  
  # Validate gene set structure
  required_cols <- c("gs_name", "gene_symbol")
  if (!all(required_cols %in% colnames(geneSet_df))) {
    stop("Gene set must contain columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }
  
  # Remove rows with missing values
  complete_rows <- complete.cases(geneSet_df[, required_cols])
  if (sum(!complete_rows) > 0) {
    if (verbose) {
      message("Removing ", sum(!complete_rows), " gene set entries with missing values")
    }
    geneSet_df <- geneSet_df[complete_rows, ]
  }
  
  # Convert to list format using base R
  geneSet_list <- split(geneSet_df$gene_symbol, geneSet_df$gs_name)
  
  # Remove empty gene sets
  geneSet_list <- geneSet_list[lengths(geneSet_list) > 0]
  
  if (verbose) {
    n_sets <- length(geneSet_list)
    n_genes <- length(unique(unlist(geneSet_list)))
    message("Processed gene set collection '", geneSet_name, "' with ", n_sets, 
            " gene sets and ", n_genes, " unique genes")
  }
  
  return(geneSet_list)
}

# Helper function to perform GSVA analysis
perform_gsva_analysis <- function(expr_data, geneSet_list, gsva_method, kcdf, 
                                 min_sz, max_sz, mx_diff, tau, ssgsea_norm, verbose) {
  
  # Create GSVA parameter object
  if (gsva_method == "gsva") {
    param <- GSVA::gsvaParam(
      exprData = expr_data,
      geneSets = geneSet_list,
      kcdf = kcdf,
      minSize = min_sz,
      maxSize = max_sz,
      maxDiff = mx_diff
    )
  } else if (gsva_method == "ssgsea") {
    param <- GSVA::ssgseaParam(
      exprData = expr_data,
      geneSets = geneSet_list,
      minSize = min_sz,
      maxSize = max_sz,
      tau = tau,
      normalize = ssgsea_norm
    )
  }
  
  # Run GSVA analysis
  gsva_result <- tryCatch({
    GSVA::gsva(param = param, verbose = verbose)
  }, error = function(e) {
    stop("GSVA analysis failed: ", e$message, call. = FALSE)
  })
  
  if (verbose) {
    message("GSVA analysis completed. Generated scores for ", nrow(gsva_result), " gene sets")
  }
  
  return(gsva_result)
}

# Helper function to perform limma differential analysis
perform_limma_analysis <- function(object, gsva_scores, geneSet_name, verbose) {
  # Check if metadata has group column
  if (!"group" %in% colnames(object@metadata)) {
    if (verbose) {
      message("No 'group' column found in metadata. Skipping differential analysis.")
    }
    return(NULL)
  }
  
  # Check sample consistency
  if (ncol(gsva_scores) != nrow(object@metadata)) {
    warning("Sample number mismatch between GSVA scores and metadata. Skipping differential analysis.")
    return(NULL)
  }
  
  if (!all(colnames(gsva_scores) == rownames(object@metadata))) {
    warning("Sample names don't match between GSVA scores and metadata. Skipping differential analysis.")
    return(NULL)
  }
  
  if (verbose) {
    message("Running limma differential analysis for GSVA scores...")
  }
  
  # Create design matrix
  group <- factor(object@metadata$group)
  if (length(levels(group)) < 2) {
    if (verbose) {
      message("Only one group found in metadata. Skipping differential analysis.")
    }
    return(NULL)
  }
  
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", paste0("group", levels(group)[2]))
  
  # Fit linear model and perform differential analysis
  tryCatch({
    fit <- limma::lmFit(gsva_scores, design)
    fit <- limma::eBayes(fit)
    
    # Extract results
    results <- limma::topTable(fit, adjust = 'fdr', coef = 2, number = Inf)
    
    # Add gene set names as a column
    results$gene_id <- rownames(results)
    
    if (verbose) {
      n_significant <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
      message("Limma analysis completed. Found ", n_significant, 
              " significantly differential gene sets (FDR < 0.05)")
    }
    
    return(results)
    
  }, error = function(e) {
    warning("Limma analysis failed: ", e$message)
    return(NULL)
  })
}
