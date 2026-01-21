# ===== PhosSight Density Plot =====
# 展示每个样本的密度分布图（log2转换+中位数中心化后）

library(tidyverse)
library(data.table)
library(ggplot2)
library(readxl)
library(reshape2)

cat("Creating PhosSight Density Plot\n")
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

quant_matrix_df <- as.data.frame(quant_matrix_centered)
data_long <- melt(quant_matrix_df, variable.name = "Sample", value.name = "Abundance")
data_long$Sample <- as.character(data_long$Sample)
data_long <- data_long[!is.na(data_long$Abundance), ]

cat(sprintf("  Density plot data: %d rows\n", nrow(data_long)))
cat(sprintf("  Number of samples: %d\n", length(unique(data_long$Sample))))

# ===== 5. 绘制密度图 =====
cat("\n[5] Drawing density plot...\n")

p <- ggplot(data_long, aes(x = Abundance, color = Sample)) +
  geom_density(show.legend = FALSE) +
  labs(title = "Density Plot of Each Sample", 
       x = "Normalized abundance", 
       y = "Density") +
  theme_minimal()

ggsave(file.path(output_dir, "PhosSight_Density_Plot.pdf"), 
       plot = p, width = 4, height = 4, units = "in")
ggsave(file.path(output_dir, "PhosSight_Density_Plot.png"), 
       plot = p, width = 2400, height = 1800, res = 300, units = "px")

cat(sprintf("  Saved: %s\n", file.path(output_dir, "PhosSight_Density_Plot.png")))

cat("\nCompleted!\n")
