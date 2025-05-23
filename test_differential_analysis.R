# Test script for differential expression analysis functions
library(devtools)

# Load the package
cat("Loading BulkS4 package...\n")
load_all(".")

# Create test data
cat("\n=== Creating test data ===\n")
set.seed(123)

# Create count matrix with some differential expression
n_genes <- 200
n_samples <- 12
counts_mat <- matrix(rpois(n_genes * n_samples, 50), nrow = n_genes, ncol = n_samples)

# Add some differential expression for first 50 genes
treat_samples <- 7:12
for (i in 1:50) {
  # Increase expression in treatment samples
  counts_mat[i, treat_samples] <- counts_mat[i, treat_samples] + rpois(6, 30)
}

rownames(counts_mat) <- paste0("Gene", 1:n_genes)
colnames(counts_mat) <- paste0("Sample", 1:n_samples)

# Create metadata
metadata <- data.frame(
  sample_id = colnames(counts_mat),
  condition = rep(c("Control", "Treatment"), each = 6),
  batch = rep(c("Batch1", "Batch2", "Batch3"), times = 4),
  row.names = colnames(counts_mat)
)

cat("Created count matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")
cat("Metadata conditions:", table(metadata$condition), "\n")

# Test 1: Create BulkRNAseq object
cat("\n=== Test 1: Creating BulkRNAseq object ===\n")
bulk_obj <- BulkRNAseq(counts_mat, metadata, add_msig.geneSet = FALSE)
print(bulk_obj)

# Test 2: Test individual analysis functions
cat("\n=== Test 2: Testing individual analysis functions ===\n")

# Prepare metadata for analysis functions
analysis_metadata <- metadata
analysis_metadata$group <- ifelse(analysis_metadata$condition == "Treatment", "treat", "con")

cat("Testing dedup_genes function...\n")
# Create test data with duplicates
test_df <- data.frame(
  gene_id = c("GeneA", "GeneB", "GeneA", "GeneC", "GeneB"),
  sample1 = c(5.2, 3.1, 4.8, 1.0, 2.5),
  sample2 = c(6.3, 2.9, 5.7, 3.2, 3.1),
  sample3 = c(2.0, 2.5, 3.9, 2.8, 2.2)
)

dedup_result <- dedup_genes(test_df, method = "first")
cat("Original rows:", nrow(test_df), "-> Deduplicated rows:", nrow(dedup_result), "\n")

# Test DESeq2 analysis (if available)
if (requireNamespace("DESeq2", quietly = TRUE)) {
  cat("\nTesting DESeq2 analysis...\n")
  tryCatch({
    deseq2_result <- run_DESeq2(counts_mat, analysis_metadata)
    cat("DESeq2 analysis completed successfully\n")
    cat("Total genes:", nrow(deseq2_result$allDiff), "\n")
    cat("Significant genes:", nrow(deseq2_result$fltDiff), "\n")
  }, error = function(e) {
    cat("DESeq2 analysis failed:", e$message, "\n")
  })
} else {
  cat("DESeq2 not available, skipping test\n")
}

# Test edgeR analysis (if available)
if (requireNamespace("edgeR", quietly = TRUE)) {
  cat("\nTesting edgeR analysis...\n")
  tryCatch({
    edger_result <- run_edgeR(counts_mat, analysis_metadata)
    cat("edgeR analysis completed successfully\n")
    cat("Total genes:", nrow(edger_result$allDiff), "\n")
    cat("Significant genes:", nrow(edger_result$fltDiff), "\n")
  }, error = function(e) {
    cat("edgeR analysis failed:", e$message, "\n")
  })
} else {
  cat("edgeR not available, skipping test\n")
}

# Test limma-voom analysis (if available)
if (requireNamespace("limma", quietly = TRUE) && requireNamespace("edgeR", quietly = TRUE)) {
  cat("\nTesting limma-voom analysis...\n")
  tryCatch({
    limma_result <- run_limma_voom(counts_mat, analysis_metadata)
    cat("limma-voom analysis completed successfully\n")
    cat("Total genes:", nrow(limma_result$allDiff), "\n")
    cat("Significant genes:", nrow(limma_result$fltDiff), "\n")
  }, error = function(e) {
    cat("limma-voom analysis failed:", e$message, "\n")
  })
} else {
  cat("limma/edgeR not available, skipping test\n")
}

# Test 3: Test integrated analysis methods
cat("\n=== Test 3: Testing integrated analysis methods ===\n")

# Test runDiffAnalysis function
available_methods <- c()
if (requireNamespace("DESeq2", quietly = TRUE)) available_methods <- c(available_methods, "DESeq2")
if (requireNamespace("edgeR", quietly = TRUE)) available_methods <- c(available_methods, "edgeR")
if (requireNamespace("limma", quietly = TRUE) && requireNamespace("edgeR", quietly = TRUE)) {
  available_methods <- c(available_methods, "limma_voom")
}

if (length(available_methods) > 0) {
  cat("Available methods:", paste(available_methods, collapse = ", "), "\n")
  
  for (method in available_methods) {
    cat("\nTesting", method, "with BulkRNAseq object...\n")
    tryCatch({
      bulk_obj <- runDiffAnalysis(
        bulk_obj, 
        method = method,
        group_col = "condition",
        treat_level = "Treatment",
        control_level = "Control"
      )
      cat(method, "analysis completed successfully\n")
    }, error = function(e) {
      cat(method, "analysis failed:", e$message, "\n")
    })
  }
  
  # Test comparison of methods
  if (length(available_methods) > 1) {
    cat("\n=== Test 4: Comparing multiple methods ===\n")
    tryCatch({
      bulk_obj_compare <- BulkRNAseq(counts_mat, metadata, add_msig.geneSet = FALSE)
      bulk_obj_compare <- compareDiffMethods(
        bulk_obj_compare,
        methods = available_methods[1:min(2, length(available_methods))],
        group_col = "condition",
        treat_level = "Treatment",
        control_level = "Control"
      )
      cat("Method comparison completed successfully\n")
    }, error = function(e) {
      cat("Method comparison failed:", e$message, "\n")
    })
  }
  
  # Test result extraction
  cat("\n=== Test 5: Testing result extraction ===\n")
  tryCatch({
    analysis_names <- names(bulk_obj@allDiff)
    if (length(analysis_names) > 0) {
      cat("Available analyses:", paste(analysis_names, collapse = ", "), "\n")
      
      # Extract results
      all_results <- extractDiffResults(bulk_obj, analysis_names[1])
      sig_results <- extractDiffResults(bulk_obj, analysis_names[1], significant_only = TRUE)
      
      cat("Total results:", nrow(all_results), "\n")
      cat("Significant results:", nrow(sig_results), "\n")
      
      # Get summary
      summary_df <- getAnalysisSummary(bulk_obj)
      print(summary_df)
    }
  }, error = function(e) {
    cat("Result extraction failed:", e$message, "\n")
  })
  
} else {
  cat("No analysis packages available for testing\n")
}

cat("\n=== All tests completed ===\n")
cat("Note: Some tests may be skipped if required packages are not installed\n")
cat("To install required packages, run:\n")
cat("BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'clusterProfiler'))\n") 