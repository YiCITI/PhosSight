# ===== Missing Value Site Count Comparison =====
# 展示不同缺失值阈值下两种方法的位点数量对比

library(tidyverse)
library(data.table)
library(ggplot2)

cat("Creating Missing Value Site Count Comparison\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ===== 配置路径 =====
phosSight_file <- "../../PhosSight_Merged_full_WithSampleID.tsv"
phosphoRS_file <- "../../PhosphoRS_Merged_full_WithSampleID.tsv"
output_dir <- "../figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===== 1. 读取数据 =====
cat("[1] Reading data...\n")

PhosSight <- read.delim(phosSight_file, sep = '\t', row.names = 1, check.names = FALSE)
PhosphoRS <- read.delim(phosphoRS_file, sep = '\t', row.names = 1, check.names = FALSE)

PhosSight <- PhosSight[PhosSight$Modi == 'Phospho', ]
PhosphoRS <- PhosphoRS[PhosphoRS$Modi == 'Phospho', ]

# ===== 2. 提取样本列 =====
cat("\n[2] Extracting sample columns...\n")

info_cols <- c("Site", "quality", "AApos", "AA", "Modi", "GeneID", "Protein", "basename")
phosSight_sample_cols <- setdiff(colnames(PhosSight), info_cols)
phosSight_sample_cols <- phosSight_sample_cols[!grepl("_Reference$", phosSight_sample_cols)]

phosphoRS_sample_cols <- setdiff(colnames(PhosphoRS), info_cols)
phosphoRS_sample_cols <- phosphoRS_sample_cols[!grepl("_Reference$", phosphoRS_sample_cols)]

total_samples_ps <- length(phosSight_sample_cols)
total_samples_pr <- length(phosphoRS_sample_cols)

cat(sprintf("  PhosSight total samples: %d\n", total_samples_ps))
cat(sprintf("  PhosphoRS total samples: %d\n", total_samples_pr))

# ===== 3. 计算不同样本数量cutoff下的位点数量 =====
cat("\n[3] Calculating site counts at different sample count cutoffs...\n")

max_samples <- max(total_samples_ps, total_samples_pr)
sample_count_cutoffs <- seq(0, max_samples, by = 5)

phosSight_non_missing_count <- rowSums(!is.na(PhosSight[, phosSight_sample_cols]))
phosphoRS_non_missing_count <- rowSums(!is.na(PhosphoRS[, phosphoRS_sample_cols]))

results <- data.frame(
  SampleCountCutoff = sample_count_cutoffs,
  PhosSight_Count = 0,
  PhosphoRS_Count = 0
)

for (i in 1:length(sample_count_cutoffs)) {
  cutoff_samples <- sample_count_cutoffs[i]
  results$PhosSight_Count[i] <- sum(phosSight_non_missing_count >= cutoff_samples)
  results$PhosphoRS_Count[i] <- sum(phosphoRS_non_missing_count >= cutoff_samples)
}

results$Advantage <- results$PhosSight_Count - results$PhosphoRS_Count

# ===== 4. 绘制对比图 =====
cat("\n[4] Drawing comparison plot...\n")

plot_data <- data.frame(
  Count = c(results$PhosSight_Count, results$PhosphoRS_Count),
  SampleCountCutoff = c(results$SampleCountCutoff, results$SampleCountCutoff),
  Method = c(rep("PhosSight", nrow(results)), rep("PhosphoRS", nrow(results)))
)

plot_data$Method <- factor(plot_data$Method, levels = c("PhosSight", "PhosphoRS"))

p <- ggplot(plot_data, aes(x = Count, y = SampleCountCutoff, color = Method, group = Method)) +
  geom_line(linewidth = 1.5, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(
    values = c("PhosSight" = "#E74C3C", "PhosphoRS" = "#3498DB"),
    name = ""
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    text = element_text(size = 12)
  ) +
  labs(
    x = "Number of Phosphosites",
    y = "Non-Missing Value Cutoff (Number of Samples)"
  ) +
  scale_y_continuous(breaks = seq(0, max_samples, by = max(10, round(max_samples/10))), 
                     limits = c(0, max_samples))

ggsave(file.path(output_dir, "MissingValue_SiteCount_Comparison.png"),
       plot = p, width = 11, height = 7, dpi = 300, bg = "white")

ggsave(file.path(output_dir, "MissingValue_SiteCount_Comparison.svg"),
       plot = p, width = 11, height = 7, bg = "white")

cat(sprintf("  Saved: %s\n", file.path(output_dir, "MissingValue_SiteCount_Comparison.png")))
cat(sprintf("  Saved: %s\n", file.path(output_dir, "MissingValue_SiteCount_Comparison.svg")))

cat("\nCompleted!\n")
