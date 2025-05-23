# Generate example dataset GSE150392 for BulkS4 package
# This script downloads and processes GSE150392 data from GEO
# The data will be used as example data for package demonstrations

# Load required packages
required_packages <- c("GEOquery", "devtools")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Required packages missing: ", paste(missing_packages, collapse = ", "), 
       "\nInstall with: BiocManager::install(c('", paste(missing_packages, collapse = "', '"), "'))",
       call. = FALSE)
}

# Load the BulkS4 package functions (assuming we're in development)
if (file.exists("R/differential_analysis.R")) {
  source("R/differential_analysis.R")
} else {
  # Try to load from installed package
  if (!requireNamespace("BulkS4", quietly = TRUE)) {
    stop("BulkS4 package functions not available. Make sure you're in the package directory or have BulkS4 installed.",
         call. = FALSE)
  }
}

library(GEOquery)

# Function to download GEO supplementary files
download_geo_data <- function(geo_accession = "GSE150392", destdir = tempdir()) {
  message("Downloading GEO supplementary files for ", geo_accession, "...")
  
  tryCatch({
    gset_sup <- GEOquery::getGEOSuppFiles(geo_accession, baseDir = destdir)
    message("Successfully downloaded ", nrow(gset_sup), " supplementary files")
    return(gset_sup)
  }, error = function(e) {
    stop("Failed to download GEO data: ", e$message, 
         "\nPlease check your internet connection and GEO accession number.", call. = FALSE)
  })
}

# Function to read and process count data
process_count_data <- function(file_path, gene_col = "gene_id") {
  message("Reading count data from: ", basename(file_path))
  
  # Read the data
  tryCatch({
    counts_raw <- read.delim(file_path, header = TRUE, sep = ",", row.names = 1, 
                            stringsAsFactors = FALSE, check.names = FALSE)
    message("Read ", nrow(counts_raw), " genes and ", ncol(counts_raw), " samples")
  }, error = function(e) {
    stop("Failed to read count data: ", e$message, call. = FALSE)
  })
  
  # Convert to data frame and add gene_id column
  counts_df <- data.frame(
    gene_id = rownames(counts_raw),
    counts_raw,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Clean gene IDs (remove version numbers if present)
  message("Cleaning gene IDs...")
  counts_df$gene_id <- ifelse(
    grepl("_", counts_df$gene_id), 
    sub(".*_", "", counts_df$gene_id), 
    counts_df$gene_id
  )
  
  # Remove duplicates using the dedup_genes function
  message("Removing duplicate genes...")
  if (exists("dedup_genes")) {
    counts_clean <- dedup_genes(counts_df, gene_col = gene_col, method = "first", gene_to_rownames = TRUE)
  } else {
    # Fallback method if dedup_genes is not available
    warning("dedup_genes function not available, using simple deduplication")
    counts_clean <- counts_df[!duplicated(counts_df[[gene_col]]), ]
    rownames(counts_clean) <- counts_clean[[gene_col]]
    counts_clean[[gene_col]] <- NULL
  }
  
  message("Final dataset: ", nrow(counts_clean), " genes and ", ncol(counts_clean), " samples")
  return(counts_clean)
}

# Function to create metadata
create_metadata <- function(sample_names, treatment_pattern = "Cov") {
  message("Creating metadata for ", length(sample_names), " samples...")
  
  # Determine treatment groups based on sample names
  group <- ifelse(grepl(treatment_pattern, sample_names, ignore.case = TRUE), "treat", "con")
  group <- factor(group, levels = c("con", "treat"))
  
  # Create additional metadata
  metadata <- data.frame(
    sample_id = sample_names,
    group = group,
    condition = ifelse(group == "treat", "COVID-19", "Control"),
    batch = paste0("Batch", rep(1:ceiling(length(sample_names)/6), each = 6)[1:length(sample_names)]),
    row.names = sample_names,
    stringsAsFactors = FALSE
  )
  
  # Print summary
  message("Metadata summary:")
  message("  - Control samples: ", sum(group == "con"))
  message("  - Treatment samples: ", sum(group == "treat"))
  message("  - Batches: ", length(unique(metadata$batch)))
  
  return(metadata)
}

# Function to validate the dataset
validate_dataset <- function(counts, metadata) {
  message("Validating dataset...")
  
  # Check dimensions
  if (ncol(counts) != nrow(metadata)) {
    stop("Dimension mismatch: counts has ", ncol(counts), " samples but metadata has ", nrow(metadata), " rows")
  }
  
  # Check sample names
  if (!all(colnames(counts) == rownames(metadata))) {
    stop("Sample names don't match between counts and metadata")
  }
  
  # Check for missing values
  na_counts <- sum(is.na(counts))
  if (na_counts > 0) {
    warning("Found ", na_counts, " missing values in count data")
  }
  
  # Check count data type
  if (!all(counts >= 0, na.rm = TRUE)) {
    warning("Found negative values in count data")
  }
  
  # Calculate basic statistics
  total_counts <- sum(counts, na.rm = TRUE)
  median_lib_size <- median(colSums(counts, na.rm = TRUE))
  
  message("Validation complete:")
  message("  - Genes: ", nrow(counts))
  message("  - Samples: ", ncol(counts))
  message("  - Total counts: ", format(total_counts, big.mark = ","))
  message("  - Median library size: ", format(median_lib_size, big.mark = ","))
  
  return(TRUE)
}

# Function to create a smaller example dataset
create_example_subset <- function(counts, metadata, n_genes = 1000, seed = 123) {
  message("Creating example subset with ", n_genes, " genes...")
  
  set.seed(seed)
  
  # Calculate gene statistics for selection
  gene_means <- rowMeans(counts, na.rm = TRUE)
  gene_vars <- apply(counts, 1, var, na.rm = TRUE)
  
  # Select genes with reasonable expression levels
  expressed_genes <- gene_means > quantile(gene_means, 0.1, na.rm = TRUE)
  variable_genes <- gene_vars > quantile(gene_vars, 0.2, na.rm = TRUE)
  
  # Combine criteria and sample
  good_genes <- which(expressed_genes & variable_genes)
  if (length(good_genes) > n_genes) {
    selected_genes <- sample(good_genes, n_genes)
  } else {
    selected_genes <- good_genes
    message("Only ", length(good_genes), " genes meet criteria, using all of them")
  }
  
  # Create subset
  counts_subset <- counts[selected_genes, ]
  
  message("Example subset created with ", nrow(counts_subset), " genes")
  return(list(counts = counts_subset, metadata = metadata))
}

# Main execution function
main <- function(create_subset = TRUE, subset_size = 1000) {
  message("Starting GSE150392 example data generation for BulkS4 package")
  message(paste(rep("=", 60), collapse = ""))
  
  # Download data
  gset_sup <- download_geo_data("GSE150392")
  
  # Process count data (use the first supplementary file)
  count_file <- rownames(gset_sup)[1]
  counts <- process_count_data(count_file)
  
  # Create metadata
  metadata <- create_metadata(colnames(counts), treatment_pattern = "Cov")
  
  # Validate dataset
  validate_dataset(counts, metadata)
  
  # Create example subset if requested
  if (create_subset) {
    example_data <- create_example_subset(counts, metadata, n_genes = subset_size)
    example_counts <- example_data$counts
    example_metadata <- example_data$metadata
    
    # Save example data
    message("Saving example data to package...")
    usethis::use_data(example_counts, example_metadata, overwrite = TRUE)
    
    message("Example data saved as 'example_counts' and 'example_metadata'")
  }
  
  # Save full dataset (optional, commented out due to size)
  # usethis::use_data(counts, metadata, overwrite = TRUE, internal = TRUE)
  
  message("GSE150392 data processing completed successfully!")
  
  # Return data for inspection
  if (create_subset) {
    return(list(
      full_counts = counts,
      full_metadata = metadata,
      example_counts = example_counts,
      example_metadata = example_metadata
    ))
  } else {
    return(list(counts = counts, metadata = metadata))
  }
}

# Execute main function
if (!interactive()) {
  result <- main(create_subset = TRUE, subset_size = 1000)
} else {
  message("Run main() to execute the data generation process")
  message("Options:")
  message("  main(create_subset = TRUE, subset_size = 1000)  # Create 1000-gene example")
  message("  main(create_subset = FALSE)                     # Process full dataset only")
}




