# ===== PhosSight PCA Plot =====
# 展示主成分分析结果（使用归一化数据）

library(tidyverse)
library(data.table)
library(ggplot2)
library(readxl)
library(factoextra)

cat("Creating PhosSight PCA Plot\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ===== 配置路径 =====
phosSight_file_norm <- "../PhosSight_Merged_full_Normalized.tsv"
sample_info_file <- "../../mmc1.xlsx"
output_dir <- "../qcfenxi/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===== 1. 读取归一化数据 =====
cat("[1] Reading normalized data...\n")

PhosSight_norm <- read.delim(phosSight_file_norm, sep = "\t", check.names = FALSE, row.names = 1)

info_cols_norm <- c("quality", "AApos", "AA", "Modi", "GeneID", "Protein", "basename")
info_cols_norm <- info_cols_norm[info_cols_norm %in% colnames(PhosSight_norm)]
sample_cols_norm <- setdiff(colnames(PhosSight_norm), info_cols_norm)

# ===== 2. 读取样本信息 =====
cat("\n[2] Reading sample information...\n")

SampleInfo <- read_excel(sample_info_file)
SampleInfo <- as.data.frame(SampleInfo)
SampleInfo_unique <- SampleInfo[!duplicated(SampleInfo$Proteomics_Parent_Sample_IDs), ]
rownames(SampleInfo_unique) <- SampleInfo_unique$Proteomics_Parent_Sample_IDs

sample_cols_norm <- intersect(sample_cols_norm, rownames(SampleInfo_unique))
SampleInfo <- SampleInfo_unique[sample_cols_norm, , drop = FALSE]

SampleInfo$Group <- ifelse(SampleInfo$Proteomics_Tumor_Normal == "Tumor", "Tumor",
                           ifelse(SampleInfo$Proteomics_Tumor_Normal %in% 
                                  c("Adjacent_normal", "Myometrium_normal", "Enriched_normal"), 
                                  "Normal", NA))

cat(sprintf("  Matched samples: %d\n", length(sample_cols_norm)))
cat(sprintf("  Tumor samples: %d\n", sum(SampleInfo$Group == "Tumor", na.rm = TRUE)))
cat(sprintf("  Normal samples: %d\n", sum(SampleInfo$Group == "Normal", na.rm = TRUE)))

# ===== 3. 提取定量矩阵 =====
cat("\n[3] Extracting quantitative matrix...\n")

raw_matrix <- PhosSight_norm[, sample_cols_norm, drop = FALSE]
raw_matrix <- as.matrix(raw_matrix)

# 移除全NA的位点
raw_matrix <- raw_matrix[rowSums(!is.na(raw_matrix)) > 0, , drop = FALSE]

# 移除>50%缺失值的位点
missing_pct <- rowSums(is.na(raw_matrix)) / ncol(raw_matrix)
valid_sites <- which(missing_pct <= 0.5)
raw_matrix <- raw_matrix[valid_sites, , drop = FALSE]

# 用行中位数填充剩余NA
for (i in 1:nrow(raw_matrix)) {
  row_median <- median(raw_matrix[i, ], na.rm = TRUE)
  if (!is.na(row_median)) {
    raw_matrix[i, is.na(raw_matrix[i, ])] <- row_median
  }
}

# 移除仍有NA的位点
raw_matrix <- raw_matrix[complete.cases(raw_matrix), , drop = FALSE]

cat(sprintf("  After cleaning: %d sites × %d samples\n", nrow(raw_matrix), ncol(raw_matrix)))

# ===== 4. 执行PCA =====
cat("\n[4] Performing PCA...\n")

res.pca <- prcomp(t(raw_matrix), scale = TRUE, center = TRUE)

cat(sprintf("  PCA result: %d samples × %d PCs\n", nrow(res.pca$x), ncol(res.pca$x)))

# ===== 5. 准备样本标签 =====
cat("\n[5] Preparing sample labels...\n")

pca_sample_ids <- rownames(res.pca$x)
matched_indices <- match(pca_sample_ids, rownames(SampleInfo))
sample_labels <- SampleInfo$Group[matched_indices]
sample_labels <- as.factor(sample_labels)

# 设置因子水平顺序
if ("Normal" %in% levels(sample_labels) && "Tumor" %in% levels(sample_labels)) {
  levels(sample_labels) <- c("Normal", "Tumor")
}

# 移除NA标签的样本
valid_samples <- !is.na(sample_labels)
if (sum(!valid_samples) > 0) {
  res.pca$x <- res.pca$x[valid_samples, , drop = FALSE]
  sample_labels <- sample_labels[valid_samples]
}

cat(sprintf("  Groups for PCA plot: %s\n", paste(unique(sample_labels), collapse = ", ")))

# ===== 6. 绘制PCA图 =====
cat("\n[6] Drawing PCA plot...\n")

p <- fviz_pca_ind(res.pca, 
                  label = "none", 
                  habillage = sample_labels, 
                  addEllipses = TRUE, 
                  ellipse.level = 0.95,
                  palette = c("blue", "red"))

ggsave(file.path(output_dir, "PhosSight_PCA_Plot.pdf"), 
       plot = p, width = 6, height = 4, units = "in")
ggsave(file.path(output_dir, "PhosSight_PCA_Plot.png"), 
       plot = p, width = 1800, height = 1200, res = 300, units = "px")

cat(sprintf("  Saved: %s\n", file.path(output_dir, "PhosSight_PCA_Plot.png")))

cat("\nCompleted!\n")
