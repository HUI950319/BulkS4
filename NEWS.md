# BulkS4 0.1.1

## New Features

* Added comprehensive S4 classes for bulk RNA-seq data analysis
* Integrated support for differential expression analysis using DESeq2, edgeR, and limma-voom
* Added gene set enrichment analysis (GSEA) functionality
* Added gene set variation analysis (GSVA) functionality
* Comprehensive visualization functions including:
  - `bar()`: Gene expression bar/box plots
  - `hmap()`: Heatmaps with clustering
  - `volcano()`: Volcano plots for differential expression
  - `MA()`: MA plots
  - `pca_hc()`: PCA and hierarchical clustering plots

## Data Management

* BulkRNAseq S4 class for unified data storage
* Automatic MSigDB gene set integration
* Built-in data validation and consistency checks
* Example datasets included: `counts`, `metadata`, `geneSet_msig`

## Documentation

* Complete package documentation with roxygen2
* Three comprehensive vignettes:
  - Getting Started guide
  - Differential Expression Analysis tutorial
  - Visualization tutorial
* Professional pkgdown website with modern UI

## Bug Fixes

* Fixed LazyData compression warnings
* Resolved package check issues
* Improved error handling and input validation

---

# BulkS4 0.1.0

## Initial Release

* Initial development version
* Basic S4 class structure
* Core functionality for RNA-seq analysis 