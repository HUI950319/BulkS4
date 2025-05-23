# BulkS4: S4 Classes for Bulk RNA-seq Data Analysis

## Overview

BulkS4 is an R package that provides S4 classes and methods for bulk RNA-seq data analysis. The package is designed to store, analyze, and visualize bulk RNA-seq data, including differential expression analysis (using DESeq2, edgeR, limma-voom), gene set enrichment analysis (GSEA), and gene set variation analysis (GSVA).

## Key Features

- **BulkRNAseq S4 Class**: Unified data structure for storing raw counts, normalized data, sample metadata, and analysis results
- **Data Validation**: Built-in data consistency checks and validation mechanisms
- **Flexible Constructor**: Support for creating objects from matrices or preprocessed lists
- **Rich Methods**: Complete set of accessor and manipulation methods
- **MSigDB Integration**: Automatic loading of MSigDB gene sets (optional)
- **Differential Expression Analysis**: Integrated support for DESeq2, edgeR, and limma-voom
- **Gene Set Analysis**: GSEA and GSVA analysis capabilities

## Installation

```r
# Install from local source
# devtools::install_local("path/to/BulkS4")

# the package is on GitHub
devtools::install_github("HUI950319/BulkS4")
```

## Quick Start

### Load Example Data

```r
library(BulkS4)

# Load example datasets included in the package
data(counts)      # Example RNA-seq count matrix
data(metadata)    # Example sample metadata
data(geneSet_msig) # MSigDB gene sets

# Examine the example data
head(counts)
head(metadata)
length(geneSet_msig)
```

### Creating BulkRNAseq Objects

```r
# Create BulkRNAseq object from example data
bulk_obj <- BulkRNAseq(counts, metadata)

# Or create without MSigDB gene sets
bulk_obj <- BulkRNAseq(counts, metadata, add_msig.geneSet = FALSE)
```

### Basic Operations

```r
# View object summary
bulk_obj

# Get dimensions
dim(bulk_obj)

# Get gene and sample names
rownames(bulk_obj)
colnames(bulk_obj)

# Access metadata columns
bulk_obj$condition

# Subsetting
subset_obj <- bulk_obj[1:100, 1:5]  # First 100 genes, first 5 samples
```

### Data Access

```r
# Get raw count matrix
raw_counts <- getCounts(bulk_obj)

# Get normalized data matrix
norm_data <- getData(bulk_obj)

# Get metadata
sample_meta <- getMetadata(bulk_obj)

# Get differential analysis results (if available)
diff_results <- getDiffResults(bulk_obj)

# Get gene sets
gene_sets <- getGeneSets(bulk_obj)
```

### Data Manipulation

```r
# Set normalized data
normalized_data <- log2(getCounts(bulk_obj) + 1)  # Simple log2 transformation
bulk_obj <- setData(bulk_obj, normalized_data)

# Add differential analysis results
diff_result <- data.frame(
  gene = rownames(counts)[1:10],
  logFC = rnorm(10),
  pvalue = runif(10, 0, 0.05),
  padj = runif(10, 0, 0.05)
)
bulk_obj <- addDiffResults(bulk_obj, diff_result, "Control_vs_Treatment")
```

### Differential Expression Analysis

```r
# Run differential expression analysis using multiple methods
bulk_obj <- runDiffAnalysis(bulk_obj, 
                           group_col = "condition",
                           methods = c("DESeq2", "edgeR", "limma"))

# Compare methods
comparison <- compareDiffMethods(bulk_obj)

# Extract results
deseq2_results <- extractDiffResults(bulk_obj, method = "DESeq2")
```

### Gene Set Analysis

```r
# Run GSEA analysis
bulk_obj <- gsea(bulk_obj, 
                group_col = "condition",
                gene_sets = "hallmark")

# Run GSVA analysis  
bulk_obj <- gsva(bulk_obj,
                gene_sets = "hallmark",
                method = "gsva")

# Get analysis summary
summary(bulk_obj)
```

## S4 Class Structure

### BulkRNAseq Class Slots

- `counts`: Raw count matrix (genes × samples)
- `data`: Normalized data matrix (genes × samples)  
- `metadata`: Sample metadata (data.frame with sample names as rownames)
- `allDiff`: List of differential analysis results
- `geneSet`: List of gene sets for GSEA/GSVA analysis
- `gsea`: List of GSEA analysis results
- `gsva`: List of GSVA analysis results

### Available Methods

- `show()`: Display object summary
- `dim()`: Get matrix dimensions
- `rownames()`, `colnames()`: Get row and column names
- `$`: Access metadata columns
- `[`: Subset objects
- `summary()`: Detailed summary information

## Example Datasets

The package includes several example datasets:

- **counts**: Example RNA-seq count matrix with gene expression data
- **metadata**: Sample metadata with experimental conditions and batch information
- **geneSet_msig**: Curated gene sets from the Molecular Signatures Database (MSigDB)

```r
# Explore example data
data(counts)
dim(counts)          # Check dimensions
head(rownames(counts))  # View gene names

data(metadata) 
str(metadata)        # Check metadata structure
table(metadata$condition)  # View experimental groups

data(geneSet_msig)
length(geneSet_msig)     # Number of gene sets
names(geneSet_msig)[1:5] # First few gene set names
```

## Best Practices

1. **Data Consistency**: Ensure count matrix column names match metadata row names
2. **Gene Names**: Use standard gene symbols as row names
3. **Metadata**: Include all relevant experimental design information
4. **Normalization**: Properly normalize data before analysis
5. **Documentation**: Add meaningful names to analysis results
6. **Quality Control**: Check data quality before proceeding with analysis

## Dependencies

### Required
- R (>= 4.0.0)
- methods
- cli
- GSVA
- clusterProfiler
- limma

### Suggested
- DESeq2
- edgeR
- biomaRt
- SummarizedExperiment
- org.Hs.eg.db
- msigdbr

## License

MIT License

## Contributing

Issues and pull requests are welcome.

## Changelog

### v0.1.1
- Added comprehensive differential expression analysis
- Integrated GSEA and GSVA functionality
- Enhanced S4 class with additional slots
- Improved data validation and error handling
- Added example datasets

### v0.1.0
- Initial release
- Basic BulkRNAseq S4 class
- Core methods and utility functions
- MSigDB gene set integration 