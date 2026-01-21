# ===== PhosSight Boxplot for All Samples =====
# 展示所有样本的箱线图（log2转换+中位数中心化后）

library(tidyverse)
library(data.table)
library(ggplot2)
library(readxl)
library(reshape2)

cat("Creating PhosSight Boxplot for All Samples\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ===== 配置路径 =====
phosSight_file <- "../PhosSight_Merged_full_WithSampleID.tsv"
sample_info_file <- "../../mmc1.xlsx"
output_dir <- "../qcfenxi/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===== 1. 读取数据 =====
cat("[1] Reading data...\n")

PhosSight <- read.delim(phosSight_file, sep = "\t", check.names = FALSE, row.names = 1)

info_cols <- c("quality", "AApos", "AA", "Modi", "GeneID", "Protein", "basename")
info_cols <- info_cols[info_cols %in% colnames(PhosSight)]
sample_cols <- setdiff(colnames(PhosSight), info_cols)
sample_cols <- sample_cols[!grepl("_Reference$", sample_cols)]

quant_matrix <- PhosSight[, sample_cols, drop = FALSE]
quant_matrix <- as.matrix(quant_matrix)

cat(sprintf("  Total sites: %d\n", nrow(quant_matrix)))
cat(sprintf("  Total samples: %d\n", length(sample_cols)))

# ===== 2. Log2转换 =====
cat("\n[2] Applying log2 transformation...\n")

quant_matrix_log2 <- quant_matrix
quant_matrix_log2[quant_matrix_log2 == 0] <- NA
quant_matrix_log2 <- log2(quant_matrix_log2)

# ===== 3. 中位数中心化 =====
cat("\n[3] Applying median centering...\n")

medians <- apply(quant_matrix_log2, 2, function(x) median(x, na.rm = TRUE))
mean_medians <- mean(medians, na.rm = TRUE)

quant_matrix_centered <- quant_matrix_log2
for (i in 1:ncol(quant_matrix_centered)) {
  quant_matrix_centered[, i] <- quant_matrix_log2[, i] + mean_medians - medians[i]
}

cat(sprintf("  After log2+median centering: range %.3f to %.3f\n", 
            min(quant_matrix_centered, na.rm = TRUE), 
            max(quant_matrix_centered, na.rm = TRUE)))

# ===== 4. 转换为长格式 =====
cat("\n[4] Converting to long format...\n")

data_long_box <- melt(quant_matrix_centered)
colnames(data_long_box) <- c("Site", "Sample", "Abundance")
data_long_box <- data_long_box[!is.na(data_long_box$Abundance), ]

cat(sprintf("  Boxplot data: %d rows\n", nrow(data_long_box)))

# ===== 5. 绘制箱线图 =====
cat("\n[5] Drawing boxplot...\n")

p <- ggplot(data_long_box, aes(x = Sample, y = Abundance)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  theme_minimal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "PhosSight", 
       x = "Sample ID", 
       y = "Log2(Intensity)")

ggsave(file.path(output_dir, "PhosSight_Boxplot_AllSamples.pdf"), 
       plot = p, width = 16, height = 6, units = "in")
ggsave(file.path(output_dir, "PhosSight_Boxplot_AllSamples.png"), 
       plot = p, width = 4800, height = 1800, res = 300, units = "px")

cat(sprintf("  Saved: %s\n", file.path(output_dir, "PhosSight_Boxplot_AllSamples.png")))

cat("\nCompleted!\n")
