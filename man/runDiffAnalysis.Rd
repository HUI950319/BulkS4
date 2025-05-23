% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_methods.R
\name{runDiffAnalysis}
\alias{runDiffAnalysis}
\title{Run Differential Expression Analysis on BulkRNAseq Object}
\usage{
runDiffAnalysis(
  object,
  method = c("DESeq2", "edgeR", "limma_voom"),
  group_col = "group",
  treat_level = "treat",
  control_level = "con",
  fc_cutoff = 1,
  padj_cutoff = 0.05,
  min_counts = 1,
  normalize_method = "quantile",
  result_name = NULL
)
}
\arguments{
\item{object}{A BulkRNAseq object}

\item{method}{Character, analysis method: "DESeq2", "edgeR", or "limma_voom"}

\item{group_col}{Character, column name in metadata for grouping, default "group"}

\item{treat_level}{Character, treatment level name, default "treat"}

\item{control_level}{Character, control level name, default "con"}

\item{fc_cutoff}{Numeric, log2 fold change threshold, default 1}

\item{padj_cutoff}{Numeric, adjusted p-value threshold, default 0.05}

\item{min_counts}{Numeric, minimum count threshold (DESeq2 only), default 1}

\item{normalize_method}{Character, normalization method (limma-voom only), default "quantile"}

\item{result_name}{Character, name for storing results in object, default NULL (auto-generated)}
}
\value{
Updated BulkRNAseq object with analysis results added
}
\description{
Perform differential expression analysis on BulkRNAseq object using DESeq2, edgeR, or limma-voom
}
\examples{
\dontrun{
# Create example object
counts_mat <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
rownames(counts_mat) <- paste0("Gene", 1:100)
colnames(counts_mat) <- paste0("Sample", 1:10)

metadata <- data.frame(
  group = rep(c("con", "treat"), each = 5),
  row.names = colnames(counts_mat)
)

bulk_obj <- BulkRNAseq(counts_mat, metadata)

# Run different analyses
bulk_obj <- runDiffAnalysis(bulk_obj, method = "DESeq2")
bulk_obj <- runDiffAnalysis(bulk_obj, method = "edgeR")
bulk_obj <- runDiffAnalysis(bulk_obj, method = "limma_voom")
}
}
