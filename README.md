# BulkS4: S4 Classes for Bulk RNA-seq Data Analysis

<!-- badges: start -->
[![R build status](https://github.com/HUI950319/BulkS4/workflows/R-CMD-check/badge.svg)](https://github.com/HUI950319/BulkS4/actions)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-universe status badge](https://hui950319.r-universe.dev/badges/BulkS4)](https://hui950319.r-universe.dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![pkgdown](https://github.com/HUI950319/BulkS4/workflows/pkgdown/badge.svg)](https://hui950319.github.io/BulkS4/)
<!-- badges: end -->

## Overview

**BulkS4** is a comprehensive R package that provides S4 classes and methods for bulk RNA-seq data analysis. The package offers a unified framework for storing, analyzing, and visualizing bulk RNA-seq data with integrated support for:

- 🧬 **Differential Expression Analysis** using DESeq2, edgeR, and limma-voom
- 📊 **Gene Set Enrichment Analysis (GSEA)** 
- 🎯 **Gene Set Variation Analysis (GSVA)**
- 📈 **Rich Visualization Tools** for publication-ready plots
- 🗂️ **Integrated Data Management** with S4 classes

## ✨ Key Features

- **BulkRNAseq S4 Class**: Unified data structure for storing raw counts, normalized data, sample metadata, and analysis results
- **Data Validation**: Built-in data consistency checks and validation mechanisms
- **Flexible Constructor**: Support for creating objects from matrices or preprocessed lists
- **Rich Visualization**: Complete set of plotting functions for expression data, heatmaps, volcano plots, and more
- **MSigDB Integration**: Automatic loading of MSigDB gene sets (Hallmark, KEGG, Reactome, GO)
- **Multiple DE Methods**: Integrated support for DESeq2, edgeR, and limma-voom with method comparison
- **Gene Set Analysis**: GSEA and GSVA analysis capabilities with multiple gene set collections

## 📦 Installation

```r
# Install from GitHub (recommended)
if (!require("devtools")) install.packages("devtools")
devtools::install_github("HUI950319/BulkS4")

# Load the package
library(BulkS4)
```

## 🚀 Quick Start

### Load Example Data

```r
library(BulkS4)

# Load example datasets included in the package
data(counts)      # Example RNA-seq count matrix (500 genes × 6 samples)
data(metadata)    # Example sample metadata
data(geneSet_msig) # MSigDB gene sets

# Examine the example data
head(counts[, 1:4])
#>           Cov1 Cov2 Cov3 Mock1
#> CD36       156  142  123   89
#> GAPDH     1245 1198 1356 1289
#> ACTB      2134 2287 2156 2245
#> GUSB       234  198  267  198

head(metadata)
#>       sample group
#> Cov1    Cov1 treat
#> Cov2    Cov2 treat
#> Cov3    Cov3 treat
#> Mock1  Mock1   con
```

### Create BulkRNAseq Object

```r
# Create BulkRNAseq object from example data
bulk_obj <- BulkRNAseq(counts, metadata)

# View object summary
bulk_obj
#> BulkRNAseq object with:
#> - 500 genes and 6 samples
#> - Metadata: 2 columns (sample, group)
#> - Gene sets: 6 collections (18459 total sets)
#> - Differential analysis: Not performed
#> - GSEA results: Not available
#> - GSVA results: Not available

# Set normalized data (log2 transformation)
log_counts <- log2(getCounts(bulk_obj) + 1)
bulk_obj <- setData(bulk_obj, log_counts)
```

## 📊 Visualization Examples

BulkS4 provides comprehensive visualization functions for RNA-seq data analysis:

### 1. Gene Expression Bar/Box Plots

```r
# Single gene expression
bar(bulk_obj, genes = "CD36", group = "group")

# Multiple genes
bar(bulk_obj, genes = c("CD36", "GAPDH", "ACTB"), group = "group")
```

### 2. Heatmaps

```r
# Heatmap of specific genes
selected_genes <- rownames(bulk_obj)[1:20]
hmap(bulk_obj, genes = selected_genes, group = "group", 
     show_rownames = TRUE, scale = "row")

# Heatmap of top differential genes (after DE analysis)
hmap(bulk_obj, genes = 15, group = "group", scale = "row")
```

### 3. Principal Component Analysis

```r
# PCA with hierarchical clustering
pca_hc(bulk_obj, k = 2, add_ellipses = TRUE)
```

## 🧬 Differential Expression Analysis

Run comprehensive differential expression analysis with multiple methods:

```r
# Perform differential expression analysis
bulk_obj <- runDiffAnalysis(bulk_obj, 
                           group_col = "group",
                           methods = c("DESeq2", "edgeR", "limma"))

# Compare methods
comparison <- compareDiffMethods(bulk_obj)
print(comparison$summary)

# Extract results from specific method
deseq2_results <- extractDiffResults(bulk_obj, method = "DESeq2")
head(deseq2_results)
```

### Volcano Plot

```r
# Create volcano plot
volcano(bulk_obj, data = "allDiff", 
        fc_threshold = 1, 
        adj.P_value_threshold = 0.05,
        highlight_fc_threshold = 2)
```

### MA Plot

```r
# Create MA plot
MA(bulk_obj, data = "allDiff",
   fc_threshold = 1,
   adj.P_value_threshold = 0.05)
```

## 🎯 Gene Set Analysis

### Gene Set Enrichment Analysis (GSEA)

```r
# Run GSEA with Hallmark gene sets
bulk_obj <- gsea(bulk_obj, 
                group_col = "group",
                gene_sets = "hallmark",
                pvalueCutoff = 0.05)

# View GSEA results
gsea_results <- bulk_obj@gsea$hallmark
head(gsea_results)
```

### Gene Set Variation Analysis (GSVA)

```r
# Run GSVA analysis
bulk_obj <- gsva(bulk_obj,
                gene_sets = "hallmark",
                method = "gsva")

# Visualize pathway scores
pathway_names <- rownames(bulk_obj@gsva$hallmark)[1:10]
bar(bulk_obj, data = "hallmark", genes = pathway_names, group = "group")

# Heatmap of pathway activities
hmap(bulk_obj, data = "hallmark", genes = 15, group = "group", scale = "row")
```

## 🗂️ S4 Class Structure

### BulkRNAseq Class Slots

| Slot | Type | Description |
|------|------|-------------|
| `counts` | matrix | Raw count matrix (genes × samples) |
| `data` | matrix | Normalized data matrix (genes × samples) |
| `metadata` | data.frame | Sample metadata with sample names as rownames |
| `allDiff` | list | Differential analysis results from multiple methods |
| `geneSet` | list | Gene sets for GSEA/GSVA analysis |
| `gsea` | list | GSEA analysis results |
| `gsva` | list | GSVA analysis results |

### Available Methods

| Method | Description |
|--------|-------------|
| `show()` | Display object summary |
| `dim()` | Get matrix dimensions |
| `rownames()`, `colnames()` | Get row and column names |
| `$` | Access metadata columns |
| `[` | Subset objects |
| `summary()` | Detailed summary information |

### Data Access Functions

| Function | Purpose |
|----------|---------|
| `getCounts(obj)` | Get raw count matrix |
| `getData(obj)` | Get normalized data matrix |
| `getMetadata(obj)` | Get sample metadata |
| `getDiffResults(obj)` | Get differential analysis results |
| `getGeneSets(obj)` | Get gene sets |

## 📈 Visualization Gallery

The package includes comprehensive visualization capabilities:

### Expression Visualization
- **Bar/Box plots**: Compare gene expression between groups
- **Heatmaps**: Visualize expression patterns across samples and genes
- **Violin plots**: Show expression distributions

### Differential Expression Plots
- **Volcano plots**: Visualize fold change vs significance
- **MA plots**: Show average expression vs fold change
- **Heatmaps**: Display differential gene expression patterns

### Sample Relationship Plots
- **PCA plots**: Principal component analysis with group information
- **Hierarchical clustering**: Sample clustering dendrograms
- **Sample correlation**: Correlation heatmaps

### Gene Set Analysis Plots
- **Pathway scores**: Bar plots of GSVA pathway activities
- **Enrichment plots**: GSEA enrichment curves
- **Network plots**: Gene set interaction networks

## 📚 Example Datasets

The package includes carefully curated example datasets:

| Dataset | Description | Dimensions |
|---------|-------------|------------|
| `counts` | RNA-seq count matrix | 500 genes × 6 samples |
| `metadata` | Sample metadata | 6 samples × 2 variables |
| `geneSet_msig` | MSigDB gene sets | 6 collections, 18,459 sets |

```r
# Explore example data structure
data(counts)
dim(counts)                    # [1] 500   6
class(counts)                  # [1] "matrix" "array"

data(metadata)
str(metadata)
# 'data.frame': 6 obs. of 2 variables:
#  $ sample: chr [1:6] "Cov1" "Cov2" "Cov3" "Mock1" ...
#  $ group : Factor w/ 2 levels "con","treat": 2 2 2 1 1 1

data(geneSet_msig)
names(geneSet_msig)           # Gene set collections
# [1] "hallmark" "kegg"     "reactome" "go_bp"    "go_cc"    "go_mf"
```

## 💡 Best Practices

### Data Preparation
1. **Data Consistency**: Ensure count matrix column names match metadata row names
2. **Gene Names**: Use standard gene symbols (HGNC) as row names
3. **Metadata**: Include all relevant experimental design information
4. **Quality Control**: Remove low-count genes and check for batch effects

### Analysis Workflow
1. **Normalization**: Properly normalize data before analysis (log2, TPM, or method-specific)
2. **Multiple Methods**: Compare results from different DE methods
3. **Statistical Thresholds**: Use appropriate p-value and fold-change cutoffs
4. **Gene Set Analysis**: Validate findings with pathway analysis

### Visualization
1. **Consistent Colors**: Use consistent color schemes across plots
2. **Clear Labels**: Include meaningful titles and axis labels
3. **Statistical Annotation**: Add p-values and effect sizes where appropriate
4. **Publication Quality**: Use high-resolution formats for final figures

## 🔧 Advanced Usage

### Custom Gene Sets

```r
# Create custom gene set
custom_geneset <- list(
  "DNA_REPAIR" = c("BRCA1", "BRCA2", "ATM", "TP53"),
  "CELL_CYCLE" = c("CCND1", "CDK4", "RB1", "E2F1")
)

# Add to BulkRNAseq object
bulk_obj@geneSet$custom <- custom_geneset

# Run GSVA with custom gene sets
bulk_obj <- gsva(bulk_obj, gene_sets = "custom")
```

### Batch Processing

```r
# Process multiple contrasts
contrasts <- list(
  "TreatmentA_vs_Control" = c("TreatmentA", "Control"),
  "TreatmentB_vs_Control" = c("TreatmentB", "Control")
)

results_list <- lapply(contrasts, function(contrast) {
  # Subset data for specific contrast
  subset_meta <- metadata[metadata$group %in% contrast, ]
  subset_counts <- counts[, subset_meta$sample]
  
  # Create object and analyze
  obj <- BulkRNAseq(subset_counts, subset_meta)
  obj <- runDiffAnalysis(obj, group_col = "group", methods = "DESeq2")
  return(obj)
})
```

## 📋 Dependencies

### Required Packages
- R (≥ 4.0.0)
- methods, cli
- GSVA, clusterProfiler, limma
- dplyr, tidyr, ggplot2

### Suggested Packages
- DESeq2, edgeR, biomaRt
- SummarizedExperiment
- org.Hs.eg.db, msigdbr
- pheatmap, ggpubr, viridis
- factoextra, ggsci, patchwork

## 📖 Documentation

- **Package Website**: https://hui950319.github.io/BulkS4/
- **Getting Started**: [Vignette](https://hui950319.github.io/BulkS4/articles/getting-started.html)
- **Differential Analysis**: [Tutorial](https://hui950319.github.io/BulkS4/articles/differential-analysis.html)
- **Visualization Guide**: [Examples](https://hui950319.github.io/BulkS4/articles/visualization.html)

## 🤝 Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

- **Bug Reports**: [GitHub Issues](https://github.com/HUI950319/BulkS4/issues)
- **Feature Requests**: [GitHub Discussions](https://github.com/HUI950319/BulkS4/discussions)
- **Pull Requests**: [GitHub PRs](https://github.com/HUI950319/BulkS4/pulls)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## 🙏 Citation

If you use BulkS4 in your research, please cite:

```r
citation("BulkS4")
```

## 📝 Changelog

### v0.1.1 (Current)
- ✨ Added comprehensive visualization functions (`bar`, `hmap`, `volcano`, `MA`, `pca_hc`)
- 🧬 Enhanced differential expression analysis with method comparison
- 🎯 Integrated GSEA and GSVA functionality
- 📊 Improved S4 class with additional slots and validation
- 📚 Added comprehensive documentation and vignettes
- 🔧 Enhanced data validation and error handling

### v0.1.0
- 🎉 Initial release
- 📦 Basic BulkRNAseq S4 class
- 🛠️ Core methods and utility functions
- 🗂️ MSigDB gene set integration

---

<div align="center">
  <strong>Built with ❤️ for the R community</strong><br>
  <a href="https://hui950319.github.io/BulkS4/">📖 Documentation</a> •
  <a href="https://github.com/HUI950319/BulkS4/issues">🐛 Report Bug</a> •
  <a href="https://github.com/HUI950319/BulkS4/discussions">💬 Discussion</a>
</div> 