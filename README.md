# BulkS4: S4 Classes for Bulk RNA-seq Data Analysis

## 概述

BulkS4是一个R包，提供了用于bulk RNA-seq数据分析的S4类和方法。该包设计用于存储和分析bulk RNA-seq数据，包括差异表达分析、基因集富集分析(GSEA)和基因集变异分析(GSVA)。

## 主要特性

- **BulkRNAseq S4类**: 统一的数据结构，用于存储原始计数、标准化数据、样本元数据和分析结果
- **数据验证**: 内置的数据一致性检查和验证机制
- **灵活的构造函数**: 支持从矩阵或预处理的列表创建对象
- **丰富的方法**: 提供完整的访问器和操作方法
- **MSigDB集成**: 自动加载MSigDB基因集（可选）

## 安装

```r
# 从本地安装
devtools::install_local("path/to/BulkS4")

# 或者如果包在GitHub上
# devtools::install_github("username/BulkS4")
```

## 快速开始

### 创建BulkRNAseq对象

```r
library(BulkS4)

# 创建示例数据
counts_mat <- matrix(rpois(1000, 10), nrow = 100, ncol = 10)
rownames(counts_mat) <- paste0("Gene", 1:100)
colnames(counts_mat) <- paste0("Sample", 1:10)

metadata <- data.frame(
  condition = rep(c("Control", "Treatment"), each = 5),
  batch = rep(c("Batch1", "Batch2"), times = 5),
  row.names = colnames(counts_mat)
)

# 创建BulkRNAseq对象
bulk_obj <- BulkRNAseq(counts_mat, metadata)
```

### 基本操作

```r
# 查看对象
bulk_obj

# 获取维度
dim(bulk_obj)

# 获取基因名和样本名
rownames(bulk_obj)
colnames(bulk_obj)

# 访问元数据列
bulk_obj$condition

# 子集化
subset_obj <- bulk_obj[1:50, 1:5]  # 前50个基因，前5个样本
```

### 数据访问

```r
# 获取原始计数矩阵
counts <- getCounts(bulk_obj)

# 获取标准化数据矩阵
data <- getData(bulk_obj)

# 获取元数据
metadata <- getMetadata(bulk_obj)

# 获取差异分析结果（如果有）
diff_results <- getDiffResults(bulk_obj)
```

### 数据设置

```r
# 设置标准化数据
normalized_data <- log2(counts_mat + 1)  # 简单的log2转换
bulk_obj <- setData(bulk_obj, normalized_data)

# 添加差异分析结果
diff_result <- data.frame(
  gene = rownames(counts_mat)[1:10],
  logFC = rnorm(10),
  pvalue = runif(10),
  padj = runif(10)
)
bulk_obj <- addDiffResults(bulk_obj, diff_result, "Control_vs_Treatment")
```

### 从预处理列表创建对象

```r
# 如果你有预处理的数据列表
data_list <- list(
  counts = counts_mat,
  metadata = metadata,
  exprSet = log2(counts_mat + 1),
  allDiff = diff_result
)

bulk_obj2 <- BulkRNAseq(data_list)
```

## S4类结构

### BulkRNAseq类的插槽(Slots)

- `counts`: 原始计数矩阵（基因 × 样本）
- `data`: 标准化数据矩阵（基因 × 样本）
- `metadata`: 样本元数据（data.frame，行名为样本名）
- `allDiff`: 差异分析结果列表
- `geneSet`: 基因集列表（用于GSEA/GSVA）
- `gsea`: GSEA分析结果列表
- `gsva`: GSVA分析结果列表

### 可用方法

- `show()`: 显示对象摘要
- `dim()`: 获取矩阵维度
- `rownames()`, `colnames()`: 获取行名和列名
- `$`: 访问元数据列
- `[`: 子集化对象
- `summary()`: 详细摘要信息

## 最佳实践

1. **数据一致性**: 确保计数矩阵的列名与元数据的行名一致
2. **基因名**: 使用标准的基因符号作为行名
3. **元数据**: 包含所有相关的实验设计信息
4. **标准化**: 在分析前适当标准化数据
5. **文档**: 为分析结果添加有意义的名称

## 依赖

- R (>= 4.0.0)
- methods
- cli

## 许可证

MIT License

## 贡献

欢迎提交问题和拉取请求。

## 更新日志

### v0.1.0
- 初始版本
- 基本的BulkRNAseq S4类
- 核心方法和实用函数
- MSigDB基因集集成 