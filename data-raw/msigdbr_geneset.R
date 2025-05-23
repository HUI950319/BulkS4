# Generate MSigDB gene sets for BulkS4 package
# This script downloads and processes MSigDB gene sets using msigdbr package
# The resulting data will be saved as internal package data

# Load required packages
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("Package 'msigdbr' is required. Install with: install.packages('msigdbr')")
}

if (!requireNamespace("pbapply", quietly = TRUE)) {
  message("Package 'pbapply' not found. Using base apply functions.")
  use_pbapply <- FALSE
} else {
  use_pbapply <- TRUE
}

library(msigdbr)

# Function to safely get MSigDB collections
get_msigdb_collections <- function(species = "Homo sapiens") {
  tryCatch({
    collections <- msigdbr::msigdbr_collections()
    message("Successfully retrieved ", nrow(collections), " MSigDB collections")
    return(collections)
  }, error = function(e) {
    stop("Failed to retrieve MSigDB collections: ", e$message, call. = FALSE)
  })
}

# Function to create parameter list for msigdbr
create_msigdbr_args <- function(collections, species = "Homo sapiens") {
  args_list <- vector("list", nrow(collections))
  
  for (i in seq_len(nrow(collections))) {
    collection <- collections$gs_collection[i]
    subcollection <- collections$gs_subcollection[i]
    
    # Handle empty subcollections
    if (is.na(subcollection) || subcollection == "") {
      subcollection <- NULL
    }
    
    args_list[[i]] <- list(
      species = species,
      category = collection,
      subcategory = subcollection
    )
  }
  
  # Create meaningful names for the list
  names(args_list) <- paste0(
    collections$gs_collection,
    ifelse(is.na(collections$gs_subcollection) | collections$gs_subcollection == "", 
           "", paste0("_", collections$gs_subcollection))
  )
  
  return(args_list)
}

# Function to download and process gene sets
download_gene_sets <- function(args_list, use_parallel = TRUE, n_cores = 2) {
  message("Downloading gene sets from MSigDB...")
  message("This may take several minutes...")
  
  # Choose apply function based on availability and preference
  apply_fun <- if (use_pbapply && use_parallel) {
    function(x, fun) pbapply::pblapply(x, fun, cl = n_cores)
  } else if (use_pbapply) {
    pbapply::pblapply
  } else {
    lapply
  }
  
  # Download and process gene sets
  gene_sets <- tryCatch({
    apply_fun(args_list, function(args) {
      # Download gene set
      gs_data <- do.call(msigdbr::msigdbr, args)
      
      # Process and return only necessary columns
      if (nrow(gs_data) > 0) {
        unique_sets <- unique(gs_data[, c("gs_name", "gene_symbol")])
        return(unique_sets)
      } else {
        return(data.frame(gs_name = character(0), gene_symbol = character(0)))
      }
    })
  }, error = function(e) {
    stop("Failed to download gene sets: ", e$message, call. = FALSE)
  })
  
  # Remove empty gene sets
  gene_sets <- gene_sets[sapply(gene_sets, nrow) > 0]
  
  message("Successfully downloaded ", length(gene_sets), " gene set collections")
  return(gene_sets)
}

# Function to validate gene sets
validate_gene_sets <- function(gene_sets) {
  message("Validating gene sets...")
  
  # Check structure
  for (i in seq_along(gene_sets)) {
    gs <- gene_sets[[i]]
    if (!all(c("gs_name", "gene_symbol") %in% colnames(gs))) {
      stop("Invalid gene set structure in collection: ", names(gene_sets)[i], call. = FALSE)
    }
  }
  
  # Calculate statistics
  total_sets <- sum(sapply(gene_sets, function(x) length(unique(x$gs_name))))
  total_genes <- length(unique(unlist(lapply(gene_sets, function(x) x$gene_symbol))))
  object_size_mb <- round(as.numeric(object.size(gene_sets)) / (1024^2), 2)
  
  message("Validation complete:")
  message("  - Total gene sets: ", total_sets)
  message("  - Unique genes: ", total_genes)
  message("  - Object size: ", object_size_mb, " MB")
  
  return(gene_sets)
}

# Main execution
main <- function() {
  message("Starting MSigDB gene set generation for BulkS4 package")
  message(paste(rep("=", 60), collapse = ""))
  
  # Get MSigDB collections
  collections <- get_msigdb_collections()
  
  # Create parameter list
  args_list <- create_msigdbr_args(collections)
  message("Created parameter list for ", length(args_list), " collections")
  
  # Download gene sets
  # Note: Reduce cores if you have memory limitations
  n_cores <- min(4, parallel::detectCores() - 1)
  geneSet_msig <- download_gene_sets(args_list, use_parallel = TRUE, n_cores = n_cores)
  
  # Validate results
  geneSet_msig <- validate_gene_sets(geneSet_msig)
  
  # Save to package data
  message("Saving gene sets to package data...")
  usethis::use_data(geneSet_msig, overwrite = TRUE, internal = FALSE)
  
  message("MSigDB gene set generation completed successfully!")
  message("Data saved as 'geneSet_msig' in package data")
  
  return(invisible(geneSet_msig))
}

# Execute main function
if (!interactive()) {
  main()
} else {
  message("Run main() to execute the gene set generation process")
}

