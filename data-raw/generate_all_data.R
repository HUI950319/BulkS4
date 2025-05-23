# Master script to generate all package data for BulkS4
# This script coordinates the generation of all package datasets

# Set up environment
message("BulkS4 Package Data Generation")
message(paste(rep("=", 50), collapse = ""))
message("This script will generate all package datasets including:")
message("  1. MSigDB gene sets")
message("  2. Example GSE150392 dataset")
message("")

# Check if we're in the correct directory
if (!file.exists("DESCRIPTION") || !file.exists("data-raw")) {
  stop("Please run this script from the package root directory", call. = FALSE)
}

# Function to check and install required packages
check_and_install_packages <- function() {
  message("Checking required packages...")
  
  # Required packages for data generation
  required_cran <- c("usethis", "devtools")
  required_bioc <- c("GEOquery", "msigdbr")
  optional_packages <- c("pbapply", "parallel")
  
  # Check CRAN packages
  missing_cran <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
  if (length(missing_cran) > 0) {
    message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran)
  }
  
  # Check Bioconductor packages
  missing_bioc <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc)
  }
  
  # Check optional packages
  missing_optional <- optional_packages[!sapply(optional_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_optional) > 0) {
    message("Optional packages not found (will use alternatives): ", paste(missing_optional, collapse = ", "))
  }
  
  message("Package check completed")
}

# Function to generate MSigDB gene sets
generate_msigdb_data <- function(skip_if_exists = TRUE) {
  message("\n--- Generating MSigDB Gene Sets ---")
  
  # Check if data already exists
  if (skip_if_exists && file.exists("data/geneSet_msig.rda")) {
    message("MSigDB data already exists, skipping generation")
    message("Set skip_if_exists = FALSE to regenerate")
    return(TRUE)
  }
  
  # Source and run the MSigDB script
  if (file.exists("data-raw/msigdbr_geneset.R")) {
    message("Running MSigDB gene set generation...")
    tryCatch({
      source("data-raw/msigdbr_geneset.R", local = TRUE)
      message("MSigDB gene sets generated successfully")
      return(TRUE)
    }, error = function(e) {
      warning("Failed to generate MSigDB data: ", e$message)
      return(FALSE)
    })
  } else {
    warning("MSigDB generation script not found: data-raw/msigdbr_geneset.R")
    return(FALSE)
  }
}

# Function to generate example data
generate_example_data <- function(skip_if_exists = TRUE, subset_size = 1000) {
  message("\n--- Generating Example GSE150392 Data ---")
  
  # Check if data already exists
  if (skip_if_exists && file.exists("data/example_counts.rda")) {
    message("Example data already exists, skipping generation")
    message("Set skip_if_exists = FALSE to regenerate")
    return(TRUE)
  }
  
  # Source and run the example data script
  if (file.exists("data-raw/example_GSE150392.R")) {
    message("Running example data generation...")
    tryCatch({
      source("data-raw/example_GSE150392.R", local = TRUE)
      message("Example data generated successfully")
      return(TRUE)
    }, error = function(e) {
      warning("Failed to generate example data: ", e$message)
      return(FALSE)
    })
  } else {
    warning("Example data generation script not found: data-raw/example_GSE150392.R")
    return(FALSE)
  }
}

# Function to validate generated data
validate_generated_data <- function() {
  message("\n--- Validating Generated Data ---")
  
  data_files <- c(
    "data/geneSet_msig.rda",
    "data/example_counts.rda",
    "data/example_metadata.rda"
  )
  
  validation_results <- list()
  
  for (file in data_files) {
    if (file.exists(file)) {
      file_size <- file.size(file)
      validation_results[[basename(file)]] <- list(
        exists = TRUE,
        size_mb = round(file_size / (1024^2), 2)
      )
      message("✓ ", basename(file), " (", round(file_size / (1024^2), 2), " MB)")
    } else {
      validation_results[[basename(file)]] <- list(exists = FALSE, size_mb = 0)
      message("✗ ", basename(file), " (missing)")
    }
  }
  
  return(validation_results)
}

# Function to create data documentation
create_data_documentation <- function() {
  message("\n--- Creating Data Documentation ---")
  
  # Create R documentation file for datasets
  doc_file <- "R/data.R"
  
  doc_content <- '#\' MSigDB Gene Sets
#\' 
#\' A comprehensive collection of gene sets from the Molecular Signatures Database (MSigDB).
#\' This dataset contains gene sets organized by collection and subcollection.
#\' 
#\' @format A named list where each element corresponds to an MSigDB collection:
#\' \\describe{
#\'   \\item{collection_name}{Data frame with columns gs_name (gene set name) and gene_symbol}
#\' }
#\' 
#\' @source \\url{https://www.gsea-msigdb.org/gsea/msigdb/}
#\' @references Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015). 
#\'   The molecular signatures database hallmark gene set collection. Cell systems, 1(6), 417-425.
"geneSet_msig"

#\' Example RNA-seq Count Data from GSE150392
#\' 
#\' A subset of RNA-seq count data from GSE150392 study comparing COVID-19 patients with controls.
#\' This dataset contains 1000 genes selected for having reasonable expression levels and variability.
#\' 
#\' @format A matrix with genes as rows and samples as columns:
#\' \\describe{
#\'   \\item{rows}{Gene symbols}
#\'   \\item{columns}{Sample identifiers}
#\'   \\item{values}{Raw count values}
#\' }
#\' 
#\' @source GEO accession GSE150392
#\' @references Overmyer, K. A., et al. (2021). Large-scale multi-omic analysis of COVID-19 severity. 
#\'   Cell systems, 12(1), 23-40.
"example_counts"

#\' Example Metadata for GSE150392 Dataset
#\' 
#\' Sample metadata corresponding to the example_counts dataset.
#\' 
#\' @format A data frame with samples as rows and metadata as columns:
#\' \\describe{
#\'   \\item{sample_id}{Sample identifier}
#\'   \\item{group}{Treatment group (con = control, treat = COVID-19)}
#\'   \\item{condition}{Descriptive condition name}
#\'   \\item{batch}{Batch information}
#\' }
#\' 
#\' @source GEO accession GSE150392
"example_metadata"
'
  
  # Write documentation
  writeLines(doc_content, doc_file)
  message("Data documentation created: ", doc_file)
}

# Main execution function
main <- function(skip_existing = TRUE, 
                 generate_msigdb = TRUE, 
                 generate_examples = TRUE,
                 subset_size = 1000,
                 create_docs = TRUE) {
  
  start_time <- Sys.time()
  
  message("Starting data generation process...")
  message("Options:")
  message("  - Skip existing data: ", skip_existing)
  message("  - Generate MSigDB: ", generate_msigdb)
  message("  - Generate examples: ", generate_examples)
  message("  - Example subset size: ", subset_size)
  message("  - Create documentation: ", create_docs)
  message("")
  
  # Check and install packages
  check_and_install_packages()
  
  # Generate data
  results <- list()
  
  if (generate_msigdb) {
    results$msigdb <- generate_msigdb_data(skip_if_exists = skip_existing)
  }
  
  if (generate_examples) {
    results$examples <- generate_example_data(skip_if_exists = skip_existing, subset_size = subset_size)
  }
  
  # Validate results
  validation <- validate_generated_data()
  
  # Create documentation
  if (create_docs) {
    create_data_documentation()
  }
  
  # Summary
  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  
  message("\n", paste(rep("=", 50), collapse = ""))
  message("Data Generation Summary")
  message(paste(rep("=", 50), collapse = ""))
  message("Duration: ", duration, " minutes")
  message("Results:")
  
  if (generate_msigdb) {
    message("  MSigDB: ", ifelse(results$msigdb, "SUCCESS", "FAILED"))
  }
  if (generate_examples) {
    message("  Examples: ", ifelse(results$examples, "SUCCESS", "FAILED"))
  }
  
  message("\nGenerated files:")
  for (name in names(validation)) {
    status <- ifelse(validation[[name]]$exists, "✓", "✗")
    size <- ifelse(validation[[name]]$exists, paste0(" (", validation[[name]]$size_mb, " MB)"), "")
    message("  ", status, " ", name, size)
  }
  
  message("\nData generation completed!")
  
  return(list(results = results, validation = validation, duration = duration))
}

# Execute if run as script
if (!interactive()) {
  result <- main()
} else {
  message("Run main() to start the data generation process")
  message("\nOptions:")
  message("  main()                                    # Generate all data")
  message("  main(skip_existing = FALSE)              # Regenerate existing data")
  message("  main(generate_msigdb = FALSE)            # Skip MSigDB generation")
  message("  main(generate_examples = FALSE)          # Skip example data")
  message("  main(subset_size = 500)                  # Smaller example dataset")
} 