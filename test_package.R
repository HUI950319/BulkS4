# Test script for BulkS4 package
library(devtools)

# Load the package
load_all(".")

# Create test data
counts_mat <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
rownames(counts_mat) <- paste0("Gene", 1:100)
colnames(counts_mat) <- paste0("Sample", 1:10)

metadata <- data.frame(
  condition = rep(c("Control", "Treatment"), each = 5),
  batch = rep(c("Batch1", "Batch2"), times = 5),
  row.names = colnames(counts_mat)
)

# Test object creation
cat("Testing BulkRNAseq object creation...\n")
bulk_obj <- BulkRNAseq(counts_mat, metadata, add_msig.geneSet = FALSE)

# Test basic methods
cat("Testing basic methods...\n")
print(bulk_obj)
cat("Dimensions:", dim(bulk_obj), "\n")
cat("Number of genes:", length(rownames(bulk_obj)), "\n")
cat("Number of samples:", length(colnames(bulk_obj)), "\n")

# Test accessor functions
cat("Testing accessor functions...\n")
counts <- getCounts(bulk_obj)
data <- getData(bulk_obj)
meta <- getMetadata(bulk_obj)

cat("Counts matrix dimensions:", dim(counts), "\n")
cat("Data matrix dimensions:", dim(data), "\n")
cat("Metadata dimensions:", dim(meta), "\n")

# Test subsetting
cat("Testing subsetting...\n")
subset_obj <- bulk_obj[1:50, 1:5]
cat("Subset dimensions:", dim(subset_obj), "\n")

# Test metadata access
cat("Testing metadata access...\n")
conditions <- bulk_obj$condition
cat("Conditions:", unique(conditions), "\n")

cat("All tests completed successfully!\n") 
