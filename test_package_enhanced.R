# Enhanced test script for BulkS4 package
library(devtools)

# Load the package
cat("Loading BulkS4 package...\n")
load_all(".")

# Test helper function
test_passed <- function(test_name) {
  cat("✓", test_name, "PASSED\n")
}

test_failed <- function(test_name, error) {
  cat("✗", test_name, "FAILED:", error, "\n")
}

# Create test data
cat("\n=== Creating test data ===\n")
counts_mat <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
rownames(counts_mat) <- paste0("Gene", 1:100)
colnames(counts_mat) <- paste0("Sample", 1:10)

metadata <- data.frame(
  condition = rep(c("Control", "Treatment"), each = 5),
  batch = rep(c("Batch1", "Batch2"), times = 5),
  row.names = colnames(counts_mat)
)

# Test 1: Object creation
cat("\n=== Test 1: Object Creation ===\n")
tryCatch({
  bulk_obj <- BulkRNAseq(counts_mat, metadata, add_msig.geneSet = FALSE)
  test_passed("BulkRNAseq object creation")
}, error = function(e) {
  test_failed("BulkRNAseq object creation", e$message)
  stop("Cannot proceed without valid object")
})

# Test 2: Basic methods
cat("\n=== Test 2: Basic Methods ===\n")
tryCatch({
  print(bulk_obj)
  test_passed("show() method")
}, error = function(e) {
  test_failed("show() method", e$message)
})

tryCatch({
  dims <- dim(bulk_obj)
  if (identical(dims, c(100L, 10L))) {
    test_passed("dim() method")
  } else {
    test_failed("dim() method", "Incorrect dimensions")
  }
}, error = function(e) {
  test_failed("dim() method", e$message)
})

# Test 3: Accessor functions
cat("\n=== Test 3: Accessor Functions ===\n")
tryCatch({
  counts <- getCounts(bulk_obj)
  if (identical(dim(counts), dim(counts_mat))) {
    test_passed("getCounts()")
  } else {
    test_failed("getCounts()", "Dimension mismatch")
  }
}, error = function(e) {
  test_failed("getCounts()", e$message)
})

tryCatch({
  data <- getData(bulk_obj)
  if (identical(dim(data), dim(counts_mat))) {
    test_passed("getData()")
  } else {
    test_failed("getData()", "Dimension mismatch")
  }
}, error = function(e) {
  test_failed("getData()", e$message)
})

tryCatch({
  meta <- getMetadata(bulk_obj)
  if (identical(dim(meta), dim(metadata))) {
    test_passed("getMetadata()")
  } else {
    test_failed("getMetadata()", "Dimension mismatch")
  }
}, error = function(e) {
  test_failed("getMetadata()", e$message)
})

# Test 4: Subsetting
cat("\n=== Test 4: Subsetting ===\n")
tryCatch({
  subset_obj <- bulk_obj[1:50, 1:5]
  if (identical(dim(subset_obj), c(50L, 5L))) {
    test_passed("Subsetting operation")
  } else {
    test_failed("Subsetting operation", "Incorrect subset dimensions")
  }
}, error = function(e) {
  test_failed("Subsetting operation", e$message)
})

# Test 5: Metadata access
cat("\n=== Test 5: Metadata Access ===\n")
tryCatch({
  conditions <- bulk_obj$condition
  if (length(conditions) == 10 && all(conditions %in% c("Control", "Treatment"))) {
    test_passed("$ operator for metadata access")
  } else {
    test_failed("$ operator for metadata access", "Incorrect condition values")
  }
}, error = function(e) {
  test_failed("$ operator for metadata access", e$message)
})

# Test 6: Data manipulation
cat("\n=== Test 6: Data Manipulation ===\n")
tryCatch({
  # Test setData function
  normalized_data <- log2(counts_mat + 1)
  bulk_obj_norm <- setData(bulk_obj, normalized_data)
  test_passed("setData() function")
}, error = function(e) {
  test_failed("setData() function", e$message)
})

tryCatch({
  # Test addDiffResults function
  diff_result <- data.frame(
    gene = rownames(counts_mat)[1:10],
    logFC = rnorm(10),
    pvalue = runif(10),
    padj = runif(10)
  )
  bulk_obj_diff <- addDiffResults(bulk_obj, diff_result, "test_comparison")
  test_passed("addDiffResults() function")
}, error = function(e) {
  test_failed("addDiffResults() function", e$message)
})

# Test 7: Error handling
cat("\n=== Test 7: Error Handling ===\n")
tryCatch({
  # Test invalid metadata column access
  invalid_col <- bulk_obj$nonexistent_column
  test_failed("Invalid column access", "Should have thrown an error")
}, error = function(e) {
  test_passed("Invalid column access error handling")
})

# Test 8: List input
cat("\n=== Test 8: List Input ===\n")
tryCatch({
  data_list <- list(
    counts = counts_mat,
    metadata = metadata,
    exprSet = log2(counts_mat + 1),
    allDiff = data.frame(gene = "test", logFC = 1)
  )
  bulk_obj_list <- BulkRNAseq(data_list, add_msig.geneSet = FALSE)
  test_passed("List input creation")
}, error = function(e) {
  test_failed("List input creation", e$message)
})

cat("\n=== Test Summary ===\n")
cat("All core functionality tests completed!\n")
cat("If you see any FAILED tests above, please check the implementation.\n") 