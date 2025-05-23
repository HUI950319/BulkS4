% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{$,BulkRNAseq-method}
\alias{$,BulkRNAseq-method}
\title{Extract Metadata Columns}
\usage{
\S4method{$}{BulkRNAseq}(x, name)
}
\arguments{
\item{x}{A BulkRNAseq object}

\item{name}{Character. Column name in metadata}
}
\value{
Vector of values from the specified metadata column
}
\description{
Extract columns from metadata using $ operator
}
\examples{
# Assuming bulk_obj is a BulkRNAseq object with metadata
# bulk_obj$condition  # Extract condition column
}
