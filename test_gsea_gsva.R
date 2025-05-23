# Test script for GSEA and GSVA functionality in BulkS4 package
# This script tests the gene set enrichment and variation analysis methods

cat("Testing BulkS4 GSEA and GSVA Functions\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Load required packages
library(BulkS4)

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

# Test 1: Check if methods are available
cat("\n=== Test 1: Method Availability ===\n")
tryCatch({
  if (exists("gsea") && is.function(gsea)) {
    test_passed("gsea generic function exists")
  } else {
    test_failed("gsea generic function", "Function not found")
  }
}, error = function(e) {
  test_failed("gsea function check", e$message)
})

tryCatch({
  if (exists("gsva") && is.function(gsva)) {
    test_passed("gsva generic function exists")
  } else {
    test_failed("gsva generic function", "Function not found")
  }
}, error = function(e) {
  test_failed("gsva function check", e$message)
})

# Test 2: Check required packages
cat("\n=== Test 2: Package Dependencies ===\n")
required_packages <- c("clusterProfiler", "GSVA", "limma")

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    test_passed(paste("Package available:", pkg))
  } else {
    test_skipped(paste("Package missing:", pkg), "Install with BiocManager::install()")
  }
}

# Test 3: Create test data
cat("\n=== Test 3: Test Data Creation ===\n")
tryCatch({
  # Create mock expression data
  set.seed(123)
  n_genes <- 100
  n_samples <- 12
  
  # Generate mock count data
  counts <- matrix(
    rpois(n_genes * n_samples, lambda = 50),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Sample", 1:n_samples)
  
  # Generate mock normalized data (log2 transformed)
  normalized_data <- log2(counts + 1)
  
  # Create metadata
  metadata <- data.frame(
    sample_id = colnames(counts),
    group = factor(rep(c("control", "treatment"), each = 6)),
    batch = factor(rep(1:3, times = 4)),
    row.names = colnames(counts)
  )
  
  # Create mock differential expression results
  diff_results <- data.frame(
    gene_id = rownames(counts),
    logFC = rnorm(n_genes, mean = 0, sd = 1),
    AveExpr = rowMeans(normalized_data),
    t = rnorm(n_genes),
    P.Value = runif(n_genes),
    adj.P.Val = runif(n_genes),
    B = rnorm(n_genes),
    row.names = rownames(counts)
  )
  
  # Create mock gene sets
  gene_sets <- data.frame(
    gs_name = rep(paste0("GeneSet", 1:10), each = 10),
    gene_symbol = sample(rownames(counts), 100, replace = TRUE)
  )
  
  test_passed("Mock test data created")
  
}, error = function(e) {
  test_failed("Test data creation", e$message)
})

# Test 4: Create BulkRNAseq object
cat("\n=== Test 4: BulkRNAseq Object Creation ===\n")
tryCatch({
  # Create BulkRNAseq object
  bulk_obj <- BulkRNAseq(
    counts = counts,
    data = normalized_data,
    metadata = metadata
  )
  
  # Add differential results
  bulk_obj@allDiff[["test_analysis"]] <- diff_results
  
  # Add gene sets
  bulk_obj@geneSet[["test_geneset"]] <- gene_sets
  
  test_passed("BulkRNAseq object created with test data")
  
}, error = function(e) {
  test_failed("BulkRNAseq object creation", e$message)
})

# Test 5: Test GSEA functionality
cat("\n=== Test 5: GSEA Analysis Testing ===\n")
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  tryCatch({
    # Test basic GSEA
    bulk_obj_gsea <- gsea(bulk_obj, 
                         geneSet_name = "test_geneset",
                         analysis_name = "test_analysis",
                         verbose = FALSE)
    
    # Check if results were stored
    if ("test_geneset" %in% names(bulk_obj_gsea@gsea)) {
      test_passed("GSEA analysis completed and results stored")
    } else {
      test_failed("GSEA results storage", "Results not found in gsea slot")
    }
    
    # Test parameter validation
    tryCatch({
      gsea(bulk_obj, pvalueCutoff = 1.5)  # Should fail
      test_failed("GSEA parameter validation", "Invalid pvalue not caught")
    }, error = function(e) {
      test_passed("GSEA parameter validation (pvalueCutoff)")
    })
    
  }, error = function(e) {
    test_failed("GSEA analysis", e$message)
  })
} else {
  test_skipped("GSEA analysis", "clusterProfiler package not available")
}

# Test 6: Test GSVA functionality
cat("\n=== Test 6: GSVA Analysis Testing ===\n")
if (requireNamespace("GSVA", quietly = TRUE)) {
  tryCatch({
    # Test basic GSVA
    bulk_obj_gsva <- gsva(bulk_obj, 
                         geneSet_name = "test_geneset",
                         gsva_method = "gsva",
                         verbose = FALSE)
    
    # Check if results were stored
    if ("test_geneset" %in% names(bulk_obj_gsva@gsva)) {
      test_passed("GSVA analysis completed and results stored")
    } else {
      test_failed("GSVA results storage", "Results not found in gsva slot")
    }
    
    # Test ssGSEA method
    bulk_obj_ssgsea <- gsva(bulk_obj, 
                           geneSet_name = "test_geneset",
                           gsva_method = "ssgsea",
                           verbose = FALSE)
    
    if ("test_geneset" %in% names(bulk_obj_ssgsea@gsva)) {
      test_passed("ssGSEA analysis completed")
    } else {
      test_failed("ssGSEA analysis", "Results not stored")
    }
    
    # Test parameter validation
    tryCatch({
      gsva(bulk_obj, min_sz = -1)  # Should fail
      test_failed("GSVA parameter validation", "Invalid min_sz not caught")
    }, error = function(e) {
      test_passed("GSVA parameter validation (min_sz)")
    })
    
  }, error = function(e) {
    test_failed("GSVA analysis", e$message)
  })
} else {
  test_skipped("GSVA analysis", "GSVA package not available")
}

# Test 7: Test multiple gene sets
cat("\n=== Test 7: Multiple Gene Sets Testing ===\n")
tryCatch({
  # Add another gene set
  gene_sets2 <- data.frame(
    gs_name = rep(paste0("GeneSet2_", 1:5), each = 8),
    gene_symbol = sample(rownames(counts), 40, replace = TRUE)
  )
  bulk_obj@geneSet[["test_geneset2"]] <- gene_sets2
  
  # Test multiple gene sets with GSVA (if available)
  if (requireNamespace("GSVA", quietly = TRUE)) {
    bulk_obj_multi <- gsva(bulk_obj, 
                          geneSet_name = c("test_geneset", "test_geneset2"),
                          verbose = FALSE)
    
    if (all(c("test_geneset", "test_geneset2") %in% names(bulk_obj_multi@gsva))) {
      test_passed("Multiple gene sets GSVA analysis")
    } else {
      test_failed("Multiple gene sets GSVA", "Not all results stored")
    }
  } else {
    test_skipped("Multiple gene sets GSVA", "GSVA package not available")
  }
  
}, error = function(e) {
  test_failed("Multiple gene sets testing", e$message)
})

# Test 8: Error handling
cat("\n=== Test 8: Error Handling ===\n")
tryCatch({
  # Test with non-existent gene set
  tryCatch({
    gsea(bulk_obj, geneSet_name = "nonexistent")
    test_failed("Error handling", "Non-existent gene set not caught")
  }, error = function(e) {
    test_passed("Error handling for non-existent gene set")
  })
  
  # Test with invalid object
  tryCatch({
    gsea("not_an_object")
    test_failed("Error handling", "Invalid object not caught")
  }, error = function(e) {
    test_passed("Error handling for invalid object")
  })
  
}, error = function(e) {
  test_failed("Error handling tests", e$message)
})

# Summary
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("GSEA and GSVA Testing Summary\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("All tests completed!\n")
cat("\nTo run GSEA/GSVA analysis on real data:\n")
cat("1. Ensure clusterProfiler and GSVA packages are installed\n")
cat("2. Load gene sets using data generation scripts\n")
cat("3. Run differential expression analysis first\n")
cat("4. Use gsea() and gsva() methods on your BulkRNAseq object\n")

# Example usage
cat("\nExample usage:\n")
cat("# GSEA analysis\n")
cat("bulk_obj <- gsea(bulk_obj, geneSet_name = 'H')\n")
cat("\n# GSVA analysis\n")
cat("bulk_obj <- gsva(bulk_obj, geneSet_name = 'H', gsva_method = 'gsva')\n")
cat("\n# Multiple gene sets\n")
cat("bulk_obj <- gsea(bulk_obj, geneSet_name = c('H', 'C2_CP_KEGG'))\n") 