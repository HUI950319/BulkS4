# Test script for data generation functions
# This script tests the data generation capabilities of BulkS4 package

# Set up environment
cat("Testing BulkS4 Data Generation Functions\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Check if we're in the correct directory
if (!file.exists("DESCRIPTION") || !file.exists("data-raw")) {
  stop("Please run this script from the package root directory")
}

# Test helper functions
test_passed <- function(test_name) {
  cat("✓", test_name, "PASSED\n")
}

test_failed <- function(test_name, error) {
  cat("✗", test_name, "FAILED:", error, "\n")
}

test_skipped <- function(test_name, reason) {
  cat("⚠", test_name, "SKIPPED:", reason, "\n")
}

# Test 1: Check data-raw scripts exist
cat("\n=== Test 1: Script Availability ===\n")
required_scripts <- c(
  "data-raw/generate_all_data.R",
  "data-raw/msigdbr_geneset.R", 
  "data-raw/example_GSE150392.R"
)

for (script in required_scripts) {
  if (file.exists(script)) {
    test_passed(paste("Script exists:", basename(script)))
  } else {
    test_failed(paste("Script missing:", basename(script)), "File not found")
  }
}

# Test 2: Check package dependencies
cat("\n=== Test 2: Package Dependencies ===\n")
required_packages <- c("usethis", "devtools")
optional_packages <- c("msigdbr", "GEOquery", "pbapply", "parallel")

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    test_passed(paste("Required package:", pkg))
  } else {
    test_failed(paste("Required package:", pkg), "Not installed")
  }
}

for (pkg in optional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    test_passed(paste("Optional package:", pkg))
  } else {
    test_skipped(paste("Optional package:", pkg), "Not installed")
  }
}

# Test 3: Test individual script functions
cat("\n=== Test 3: Script Function Testing ===\n")

# Test msigdbr script functions
if (file.exists("data-raw/msigdbr_geneset.R")) {
  tryCatch({
    source("data-raw/msigdbr_geneset.R", local = TRUE)
    
    # Test get_msigdb_collections function
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      collections <- get_msigdb_collections()
      if (is.data.frame(collections) && nrow(collections) > 0) {
        test_passed("get_msigdb_collections() function")
      } else {
        test_failed("get_msigdb_collections() function", "Invalid return value")
      }
      
      # Test create_msigdbr_args function
      args_list <- create_msigdbr_args(collections[1:3, ])  # Test with first 3 collections
      if (is.list(args_list) && length(args_list) == 3) {
        test_passed("create_msigdbr_args() function")
      } else {
        test_failed("create_msigdbr_args() function", "Invalid return value")
      }
    } else {
      test_skipped("MSigDB functions", "msigdbr package not available")
    }
  }, error = function(e) {
    test_failed("MSigDB script loading", e$message)
  })
} else {
  test_failed("MSigDB script", "File not found")
}

# Test example data script functions
if (file.exists("data-raw/example_GSE150392.R")) {
  tryCatch({
    source("data-raw/example_GSE150392.R", local = TRUE)
    
    # Test create_metadata function
    sample_names <- paste0("Sample", 1:10)
    metadata <- create_metadata(sample_names, treatment_pattern = "5|6|7|8|9|10")
    if (is.data.frame(metadata) && nrow(metadata) == 10) {
      test_passed("create_metadata() function")
    } else {
      test_failed("create_metadata() function", "Invalid return value")
    }
    
    # Test validate_dataset function with mock data
    mock_counts <- matrix(rpois(100, 10), nrow = 10, ncol = 10)
    colnames(mock_counts) <- sample_names
    rownames(mock_counts) <- paste0("Gene", 1:10)
    
    validation_result <- validate_dataset(mock_counts, metadata)
    if (isTRUE(validation_result)) {
      test_passed("validate_dataset() function")
    } else {
      test_failed("validate_dataset() function", "Validation failed")
    }
    
  }, error = function(e) {
    test_failed("Example data script loading", e$message)
  })
} else {
  test_failed("Example data script", "File not found")
}

# Test 4: Test master script
cat("\n=== Test 4: Master Script Testing ===\n")
if (file.exists("data-raw/generate_all_data.R")) {
  tryCatch({
    source("data-raw/generate_all_data.R", local = TRUE)
    
    # Test check_and_install_packages function
    check_and_install_packages()
    test_passed("check_and_install_packages() function")
    
    # Test validation function
    validation_result <- validate_generated_data()
    if (is.list(validation_result)) {
      test_passed("validate_generated_data() function")
    } else {
      test_failed("validate_generated_data() function", "Invalid return value")
    }
    
  }, error = function(e) {
    test_failed("Master script loading", e$message)
  })
} else {
  test_failed("Master script", "File not found")
}

# Test 5: Check existing data files
cat("\n=== Test 5: Existing Data Files ===\n")
data_files <- c(
  "data/geneSet_msig.rda",
  "data/example_counts.rda",
  "data/example_metadata.rda"
)

for (file in data_files) {
  if (file.exists(file)) {
    file_size <- round(file.size(file) / (1024^2), 2)
    test_passed(paste(basename(file), "(", file_size, "MB)"))
  } else {
    test_skipped(basename(file), "Not generated yet")
  }
}

# Test 6: Test data loading (if files exist)
cat("\n=== Test 6: Data Loading Test ===\n")
if (file.exists("data/example_counts.rda") && file.exists("data/example_metadata.rda")) {
  tryCatch({
    load("data/example_counts.rda")
    load("data/example_metadata.rda")
    
    if (exists("example_counts") && exists("example_metadata")) {
      # Check data structure
      if (is.matrix(example_counts) || is.data.frame(example_counts)) {
        test_passed("example_counts data structure")
      } else {
        test_failed("example_counts data structure", "Invalid data type")
      }
      
      if (is.data.frame(example_metadata)) {
        test_passed("example_metadata data structure")
      } else {
        test_failed("example_metadata data structure", "Invalid data type")
      }
      
      # Check dimensions match
      if (ncol(example_counts) == nrow(example_metadata)) {
        test_passed("Data dimensions consistency")
      } else {
        test_failed("Data dimensions consistency", "Dimension mismatch")
      }
    } else {
      test_failed("Data loading", "Objects not found after loading")
    }
  }, error = function(e) {
    test_failed("Data loading", e$message)
  })
} else {
  test_skipped("Data loading test", "Data files not available")
}

if (file.exists("data/geneSet_msig.rda")) {
  tryCatch({
    load("data/geneSet_msig.rda")
    
    if (exists("geneSet_msig")) {
      if (is.list(geneSet_msig) && length(geneSet_msig) > 0) {
        test_passed("geneSet_msig data structure")
        
        # Check first collection structure
        first_collection <- geneSet_msig[[1]]
        if (is.data.frame(first_collection) && 
            all(c("gs_name", "gene_symbol") %in% colnames(first_collection))) {
          test_passed("Gene set collection structure")
        } else {
          test_failed("Gene set collection structure", "Invalid column names")
        }
      } else {
        test_failed("geneSet_msig data structure", "Invalid data type")
      }
    } else {
      test_failed("Gene set loading", "Object not found after loading")
    }
  }, error = function(e) {
    test_failed("Gene set loading", e$message)
  })
} else {
  test_skipped("Gene set loading test", "Gene set file not available")
}

# Summary
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Data Generation Testing Summary\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("All tests completed!\n")
cat("\nTo generate data, run:\n")
cat("source('data-raw/generate_all_data.R')\n")
cat("main()\n")
cat("\nNote: Data generation requires internet connection and may take 10-20 minutes\n") 