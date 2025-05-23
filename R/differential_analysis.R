#' @importFrom stats complete.cases aggregate var
#' @importFrom utils head

#' @title Remove Duplicate Genes from Expression Matrix
#' @description Remove duplicate gene IDs from expression matrix using various strategies
#'
#' @param df Data frame or tibble containing gene expression data
#' @param gene_col Character, column name containing gene IDs, default "gene_id"
#' @param method Character, method for handling duplicate genes:
#'   \itemize{
#'     \item "first": Keep first occurrence of each gene ID
#'     \item "max_mean": Keep row with highest mean expression
#'     \item "mean": Take average of all rows for each gene ID
#'     \item "max_variance": Keep row with highest expression variance
#'   }
#' @param gene_to_rownames Logical, whether to convert gene ID column to rownames
#'
#' @return Data frame with deduplicated gene expression data
#' @export
#'
#' @examples
#' # Create example data
#' df <- data.frame(
#'   gene_id = c("GeneA", "GeneB", "GeneA", "GeneC"),
#'   sample1 = c(5.2, 3.1, 4.8, 1.0),
#'   sample2 = c(6.3, 2.9, 5.7, 3.2),
#'   sample3 = c(2.0, 2.5, 3.9, 2.8)
#' )
#'
#' # Remove duplicates using different methods
#' dedup_genes(df, method = "first")
#' dedup_genes(df, method = "max_mean")
dedup_genes <- function(df,
                        gene_col = "gene_id",
                        method = c("first", "max_mean", "mean", "max_variance"),
                        gene_to_rownames = TRUE) {

  # Input validation
  if (!is.data.frame(df)) {
    stop("Input 'df' must be a data frame or tibble", call. = FALSE)
  }
  if (!gene_col %in% names(df)) {
    stop("Column '", gene_col, "' not found", call. = FALSE)
  }

  method <- match.arg(method)

  # Convert to data frame and handle missing values
  df <- as.data.frame(df)

  # Move gene column to first position
  gene_idx <- which(names(df) == gene_col)
  df <- df[, c(gene_idx, setdiff(seq_along(df), gene_idx))]

  # Handle empty strings and convert expression columns to numeric
  df[df == "" | df == " "] <- NA
  expr_cols <- setdiff(names(df), gene_col)
  df[expr_cols] <- lapply(df[expr_cols], as.numeric)

  # Remove rows with missing gene IDs
  na_genes <- sum(is.na(df[[gene_col]]))
  if (na_genes > 0) {
    message("Found ", na_genes, " rows with NA in '", gene_col, "' column")
    df <- df[!is.na(df[[gene_col]]), ]
    message("Removed ", na_genes, " rows with missing gene IDs")
  }

  # Remove rows with any missing expression values
  complete_rows <- complete.cases(df[expr_cols])
  if (sum(!complete_rows) > 0) {
    message("Removed ", sum(!complete_rows), " rows with missing expression values")
    df <- df[complete_rows, ]
  }

  if (nrow(df) == 0) {
    stop("No data remaining after filtering missing values", call. = FALSE)
  }

  # Handle duplicate genes based on method
  if (method == "first") {
    dedup_df <- df[!duplicated(df[[gene_col]]), ]
  } else if (method == "max_mean") {
    df$row_mean <- rowMeans(df[expr_cols], na.rm = TRUE)
    df <- df[order(df[[gene_col]], -df$row_mean), ]
    dedup_df <- df[!duplicated(df[[gene_col]]), ]
    dedup_df$row_mean <- NULL
  } else if (method == "mean") {
    dedup_df <- aggregate(df[expr_cols], by = list(df[[gene_col]]), FUN = mean, na.rm = TRUE)
    names(dedup_df)[1] <- gene_col
  } else if (method == "max_variance") {
    df$row_var <- apply(df[expr_cols], 1, var, na.rm = TRUE)
    df <- df[order(df[[gene_col]], -df$row_var), ]
    dedup_df <- df[!duplicated(df[[gene_col]]), ]
    dedup_df$row_var <- NULL
  }

  # Convert gene column to rownames if requested
  if (gene_to_rownames) {
    rownames(dedup_df) <- dedup_df[[gene_col]]
    dedup_df[[gene_col]] <- NULL
  }

  return(dedup_df)
}

#' @title Convert Gene IDs in Expression Matrix
#' @description Convert gene IDs in expression matrix using clusterProfiler or biomaRt
#'
#' @param expr_matrix Gene expression matrix with gene IDs as rownames
#' @param from_type Source ID type (e.g., "ENSEMBL", "ENTREZID", "SYMBOL")
#' @param to_type Target ID type (e.g., "SYMBOL", "ENTREZID", "ENSEMBL")
#' @param OrgDb Annotation database to use, default "org.Hs.eg.db"
#' @param use_biomart Whether to use biomaRt for conversion, default FALSE
#' @param biomart_dataset biomaRt dataset name, default "hsapiens_gene_ensembl"
#' @param gene_to_rownames Whether to convert gene ID column to rownames
#'
#' @return Expression matrix with converted gene IDs
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example expression matrix
#' expr_matrix <- matrix(rpois(60, 10), nrow = 6, ncol = 10)
#' rownames(expr_matrix) <- c("ENSG00000223972", "ENSG00000227232",
#'                           "ENSG00000278267", "ENSG00000243485",
#'                           "ENSG00000284332", "ENSG00000237613")
#' colnames(expr_matrix) <- paste0("Sample", 1:10)
#'
#' # Convert using clusterProfiler
#' convert_genes(expr_matrix)
#' }
convert_genes <- function(expr_matrix,
                          from_type = "ENSEMBL",
                          to_type = "SYMBOL",
                          OrgDb = "org.Hs.eg.db",
                          use_biomart = FALSE,
                          biomart_dataset = "hsapiens_gene_ensembl",
                          gene_to_rownames = TRUE) {

  # Input validation
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("expr_matrix must be a matrix or data frame", call. = FALSE)
  }

  # Convert to data frame with gene IDs as first column
  expr_df <- as.data.frame(expr_matrix)
  expr_df$ID <- rownames(expr_df)
  expr_df <- expr_df[, c("ID", setdiff(names(expr_df), "ID"))]

  if (use_biomart) {
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
      stop("Package 'biomaRt' is required for biomaRt conversion", call. = FALSE)
    }

    # biomaRt conversion
    ensembl <- biomaRt::useMart("ensembl", dataset = biomart_dataset)

    id_mapping <- list(
      "ENSEMBL" = "ensembl_gene_id",
      "ENTREZID" = "entrezgene_id",
      "SYMBOL" = "hgnc_symbol"
    )

    filter_name <- id_mapping[[from_type]]
    attr_name <- id_mapping[[to_type]]

    if (is.null(filter_name) || is.null(attr_name)) {
      stop("Unsupported ID type for biomaRt. Use ENSEMBL, ENTREZID, or SYMBOL", call. = FALSE)
    }

    converted_ids <- biomaRt::getBM(
      attributes = c(filter_name, attr_name),
      filters = filter_name,
      values = expr_df$ID,
      mart = ensembl
    )
    converted_ids <- converted_ids[!is.na(converted_ids[[attr_name]]), ]
    names(converted_ids) <- c("ID", "to_id")

  } else {
    # clusterProfiler conversion
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      stop("Package 'clusterProfiler' is required for ID conversion", call. = FALSE)
    }
    if (!requireNamespace(OrgDb, quietly = TRUE)) {
      stop("Package '", OrgDb, "' is not installed", call. = FALSE)
    }

    converted_ids <- clusterProfiler::bitr(
      expr_df$ID,
      fromType = from_type,
      toType = to_type,
      OrgDb = OrgDb
    )
    names(converted_ids) <- c("ID", "to_id")
  }

  # Print conversion statistics
  total_genes <- nrow(expr_df)
  converted_genes <- nrow(converted_ids)
  conversion_rate <- round(converted_genes / total_genes * 100, 2)

  message("\n=====================================")
  message(" Gene ID Conversion Summary")
  message("=====================================")
  message(" From: ", from_type, " To: ", to_type)
  message(" Total input genes: ", total_genes)
  message(" Successfully converted: ", converted_genes)
  message(" Conversion rate: ", conversion_rate, "%")
  message("=====================================\n")

  # Merge with expression data
  merged_df <- merge(converted_ids, expr_df, by = "ID", all.x = TRUE)

  # Remove original ID column and rename converted ID column
  merged_df$ID <- NULL
  names(merged_df)[1] <- "gene_id"

  # Handle duplicates using dedup_genes function
  final_df <- dedup_genes(merged_df, gene_col = "gene_id",
                         method = "max_mean", gene_to_rownames = gene_to_rownames)

  return(final_df)
}

#' @title Run DESeq2 Differential Expression Analysis
#' @description Perform differential expression analysis using DESeq2
#'
#' @param counts Gene expression count matrix (genes x samples)
#' @param metadata Sample information data frame with 'group' column ('treat' or 'con')
#' @param fc_cutoff Log2 fold change threshold for filtering, default 1
#' @param padj_cutoff Adjusted p-value threshold for filtering, default 0.05
#' @param min_counts Minimum count threshold for filtering low-expression genes, default 1
#'
#' @return List containing:
#'   \item{counts}{Original count matrix}
#'   \item{metadata}{Sample metadata}
#'   \item{exprSet}{Normalized expression matrix}
#'   \item{allDiff}{All genes differential analysis results}
#'   \item{fltDiff}{Filtered significant differential genes}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' counts <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
#' rownames(counts) <- paste0("Gene", 1:100)
#' colnames(counts) <- paste0("Sample", 1:10)
#'
#' metadata <- data.frame(
#'   sample = colnames(counts),
#'   group = rep(c("con", "treat"), each = 5)
#' )
#'
#' # Run DESeq2 analysis
#' result <- run_DESeq2(counts, metadata)
#' }
run_DESeq2 <- function(counts, metadata, fc_cutoff = 1, padj_cutoff = 0.05, min_counts = 1) {

  # Check required packages
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for this analysis", call. = FALSE)
  }

  # Input validation
  if (!("group" %in% colnames(metadata))) {
    stop("metadata must contain 'group' column", call. = FALSE)
  }
  if (!all(metadata$group %in% c("treat", "con"))) {
    stop("group column must only contain 'treat' and 'con' values", call. = FALSE)
  }

  # Ensure counts is numeric matrix
  if (is.character(counts[, 1])) {
    counts <- counts[, -1]
  }
  counts <- as.matrix(counts)
  mode(counts) <- "integer"

  # Create DESeq dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ group
  )

  # Filter low expression genes
  dds <- dds[rowSums(DESeq2::counts(dds)) > min_counts, ]

  # Variance stabilizing transformation
  vsd <- DESeq2::vst(dds, blind = FALSE)
  exprSet_vst <- as.data.frame(SummarizedExperiment::assay(vsd))

  # Differential expression analysis
  dds <- DESeq2::DESeq(dds)

  # Get results
  contrast <- c("group", "treat", "con")
  res <- DESeq2::results(dds, contrast = contrast, alpha = padj_cutoff)

  # Shrink log2 fold changes
  res_shrink <- DESeq2::lfcShrink(dds, contrast = contrast, res = res, type = "ashr")

  # Format results
  allDiff <- as.data.frame(res_shrink)
  allDiff$gene_id <- rownames(allDiff)
  allDiff <- allDiff[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
  names(allDiff) <- c("gene_id", "baseMean", "logFC", "lfcSE", "P.Value", "adj.P.Val")
  rownames(allDiff) <- allDiff$gene_id

  # Filter significant genes
  fltDiff <- allDiff[!is.na(allDiff$adj.P.Val) &
                     allDiff$adj.P.Val < padj_cutoff &
                     abs(allDiff$logFC) > fc_cutoff, ]

  return(list(
    counts = counts,
    metadata = metadata,
    exprSet = exprSet_vst,
    allDiff = allDiff,
    fltDiff = fltDiff
  ))
}

#' @title Run edgeR Differential Expression Analysis
#' @description Perform differential expression analysis using edgeR
#'
#' @param counts Gene expression count matrix (genes x samples)
#' @param metadata Sample information data frame with 'group' column ('treat' or 'con')
#' @param fc_cutoff Log2 fold change threshold for filtering, default 1
#' @param padj_cutoff Adjusted p-value threshold for filtering, default 0.05
#'
#' @return List containing analysis results (same structure as run_DESeq2)
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_edgeR(counts, metadata)
#' }
run_edgeR <- function(counts, metadata, fc_cutoff = 1, padj_cutoff = 0.05) {

  # Check required packages
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for this analysis", call. = FALSE)
  }

  # Input validation
  if (!("group" %in% colnames(metadata))) {
    stop("metadata must contain 'group' column", call. = FALSE)
  }
  if (!all(metadata$group %in% c("treat", "con"))) {
    stop("group column must only contain 'treat' and 'con' values", call. = FALSE)
  }

  # Create DGEList object
  y <- edgeR::DGEList(counts = counts, group = metadata$group)

  # Filter genes
  keep <- edgeR::filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]

  # Normalization
  y <- edgeR::calcNormFactors(y, method = "TMM")

  # Estimate dispersion
  design <- stats::model.matrix(~ group, data = metadata)
  rownames(design) <- colnames(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)

  # Fit model and test
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit)

  # Get results
  allDiff <- edgeR::topTags(qlf, n = nrow(y))$table
  allDiff$gene_id <- rownames(allDiff)
  allDiff <- allDiff[, c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")]
  names(allDiff) <- c("gene_id", "logFC", "logCPM", "F", "P.Value", "adj.P.Val")
  rownames(allDiff) <- allDiff$gene_id

  # Filter significant genes
  fltDiff <- allDiff[allDiff$adj.P.Val < padj_cutoff & abs(allDiff$logFC) > fc_cutoff, ]

  # Get normalized expression
  exprSet_norm <- edgeR::cpm(y, log = TRUE, prior.count = 1)

  return(list(
    counts = counts,
    metadata = metadata,
    exprSet = exprSet_norm,
    allDiff = allDiff,
    fltDiff = fltDiff
  ))
}

#' @title Run limma-voom Differential Expression Analysis
#' @description Perform differential expression analysis using limma-voom
#'
#' @param counts Gene expression count matrix (genes x samples)
#' @param metadata Sample information data frame with 'group' column ('treat' or 'con')
#' @param normalize_method Normalization method for voom, default "quantile"
#' @param fc_cutoff Log2 fold change threshold for filtering, default 1
#' @param padj_cutoff Adjusted p-value threshold for filtering, default 0.05
#'
#' @return List containing analysis results (same structure as run_DESeq2)
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_limma_voom(counts, metadata)
#' }
run_limma_voom <- function(counts, metadata,
                           normalize_method = "quantile",
                           fc_cutoff = 1,
                           padj_cutoff = 0.05) {

  # Check required packages
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("edgeR", quietly = TRUE)) {
    stop("Packages 'limma' and 'edgeR' are required for this analysis", call. = FALSE)
  }

  # Input validation
  if (!("group" %in% colnames(metadata))) {
    stop("metadata must contain 'group' column", call. = FALSE)
  }
  if (!all(metadata$group %in% c("treat", "con"))) {
    stop("group column must only contain 'treat' and 'con' values", call. = FALSE)
  }

  # Create design matrix
  group <- factor(metadata$group, levels = c("con", "treat"))
  design <- stats::model.matrix(~ group)
  rownames(design) <- colnames(counts)
  colnames(design) <- c("Intercept", "treat_vs_con")

  # Create DGEList and normalize
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)

  # voom transformation
  v <- limma::voom(dge, design = design, normalize.method = normalize_method)

  # Fit linear model
  fit <- limma::lmFit(v, design = design)
  fit <- limma::eBayes(fit)

  # Get results
  allDiff <- limma::topTable(fit, adjust.method = "BH", coef = 2, number = Inf)
  allDiff$gene_id <- rownames(allDiff)
  allDiff <- allDiff[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  rownames(allDiff) <- allDiff$gene_id

  # Filter significant genes
  fltDiff <- allDiff[allDiff$adj.P.Val < padj_cutoff & abs(allDiff$logFC) > fc_cutoff, ]

  return(list(
    counts = counts,
    metadata = metadata,
    exprSet = v$E,
    allDiff = allDiff,
    fltDiff = fltDiff
  ))
}
