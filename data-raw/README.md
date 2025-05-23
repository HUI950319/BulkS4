# Data Generation Scripts for BulkS4 Package

This directory contains scripts for generating package datasets used in the BulkS4 package.

## Overview

The BulkS4 package includes several datasets:

1. **MSigDB Gene Sets** (`geneSet_msig`) - Comprehensive collection of gene sets from the Molecular Signatures Database
2. **Example RNA-seq Data** (`example_counts`, `example_metadata`) - Subset of GSE150392 dataset for demonstrations

## Scripts

### 1. `generate_all_data.R` - Master Script

The main script that coordinates all data generation processes.

**Usage:**
```r
# Generate all datasets
source("data-raw/generate_all_data.R")
main()

# Custom options
main(
  skip_existing = FALSE,      # Regenerate existing data
  generate_msigdb = TRUE,     # Generate MSigDB gene sets
  generate_examples = TRUE,   # Generate example data
  subset_size = 1000,         # Number of genes in example
  create_docs = TRUE          # Create data documentation
)
```

**Features:**
- Automatic package dependency checking and installation
- Progress tracking and error handling
- Data validation and size reporting
- Automatic documentation generation

### 2. `msigdbr_geneset.R` - MSigDB Gene Sets

Downloads and processes gene sets from the Molecular Signatures Database using the `msigdbr` package.

**Features:**
- Downloads all available MSigDB collections
- Parallel processing for faster downloads
- Data validation and statistics
- Memory-efficient processing

**Output:**
- `geneSet_msig` - Named list of gene sets organized by collection

**Requirements:**
- `msigdbr` package
- `pbapply` package (optional, for progress bars)
- Internet connection

### 3. `example_GSE150392.R` - Example Dataset

Downloads and processes GSE150392 data from GEO to create example datasets.

**Features:**
- Automatic GEO data download
- Gene ID cleaning and deduplication
- Intelligent gene selection for examples
- Comprehensive data validation

**Output:**
- `example_counts` - Count matrix (1000 genes × samples)
- `example_metadata` - Sample metadata with group information

**Requirements:**
- `GEOquery` package
- Internet connection
- BulkS4 package functions (for gene deduplication)

## Data Generation Process

### Step 1: Install Dependencies

The scripts will automatically check and install required packages:

```r
# CRAN packages
install.packages(c("usethis", "devtools", "pbapply"))

# Bioconductor packages
BiocManager::install(c("GEOquery", "msigdbr"))
```

### Step 2: Run Data Generation

From the package root directory:

```r
# Option 1: Use master script (recommended)
source("data-raw/generate_all_data.R")
main()

# Option 2: Run individual scripts
source("data-raw/msigdbr_geneset.R")
source("data-raw/example_GSE150392.R")
```

### Step 3: Verify Generated Data

Check the `data/` directory for generated `.rda` files:

- `data/geneSet_msig.rda` (~80 MB)
- `data/example_counts.rda` (~1 MB)
- `data/example_metadata.rda` (~1 KB)

## Data Specifications

### MSigDB Gene Sets (`geneSet_msig`)

```r
# Structure
str(geneSet_msig, max.level = 2)
# List of 23 collections
# Each collection contains data.frame with:
#   - gs_name: Gene set name
#   - gene_symbol: Gene symbol

# Example usage
hallmark_sets <- geneSet_msig$H
pathway_sets <- geneSet_msig$C2_CP_KEGG
```

### Example Data (`example_counts`, `example_metadata`)

```r
# Count matrix
dim(example_counts)  # 1000 genes × 24 samples
class(example_counts)  # matrix

# Metadata
dim(example_metadata)  # 24 samples × 4 variables
colnames(example_metadata)  # sample_id, group, condition, batch
```

## Troubleshooting

### Common Issues

1. **Internet Connection Required**
   - Both MSigDB and GEO data require internet access
   - Consider using institutional networks for better stability

2. **Memory Limitations**
   - MSigDB download can use significant memory
   - Reduce parallel cores if experiencing memory issues
   - Consider generating datasets separately

3. **Package Dependencies**
   - Ensure Bioconductor is properly installed
   - Some packages may require system dependencies

### Error Solutions

```r
# If msigdbr fails
install.packages("msigdbr")

# If GEOquery fails
BiocManager::install("GEOquery")

# If parallel processing fails
# Edit scripts to use use_parallel = FALSE

# If memory issues occur
# Reduce n_cores in msigdbr_geneset.R
```

## Customization

### Modifying Example Dataset Size

Edit `example_GSE150392.R`:

```r
# Change subset size
create_example_subset(counts, metadata, n_genes = 500)  # Smaller dataset
create_example_subset(counts, metadata, n_genes = 2000) # Larger dataset
```

### Selecting Specific MSigDB Collections

Edit `msigdbr_geneset.R`:

```r
# Filter collections before processing
collections <- get_msigdb_collections()
collections <- collections[collections$gs_collection %in% c("H", "C2"), ]
```

### Adding Custom Datasets

1. Create new script in `data-raw/`
2. Add to `generate_all_data.R` main function
3. Update documentation in `R/data.R`

## Performance Notes

- **MSigDB Generation**: 5-15 minutes depending on internet speed
- **Example Data**: 2-5 minutes depending on internet speed
- **Total Size**: ~80-100 MB for all datasets
- **Parallel Processing**: Recommended for MSigDB generation

## Data Sources

- **MSigDB**: https://www.gsea-msigdb.org/gsea/msigdb/
- **GSE150392**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150392

## References

1. Liberzon, A., et al. (2015). The molecular signatures database hallmark gene set collection. Cell systems, 1(6), 417-425.
2. Overmyer, K. A., et al. (2021). Large-scale multi-omic analysis of COVID-19 severity. Cell systems, 12(1), 23-40. 