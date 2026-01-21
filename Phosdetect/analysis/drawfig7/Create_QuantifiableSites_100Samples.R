# ===== Quantifiable Sites at 100 Sample Threshold =====
# 对比两种方法在100个样本阈值下的可定量位点数量

library(tidyverse)
library(data.table)
library(ggplot2)

cat("Creating Quantifiable Sites at 100 Sample Threshold\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ===== 配置路径 =====
phosSight_file <- "../../PhosSight_Merged_full_WithSampleID.tsv"
phosphoRS_file <- "../../PhosphoRS_Merged_full_WithSampleID.tsv"
output_dir <- "../figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

quantifiable_threshold <- 100

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

# ===== 3. 计算可定量位点 =====
cat("\n[3] Calculating quantifiable sites...\n")

phosSight_non_missing_count <- rowSums(!is.na(PhosSight[, phosSight_sample_cols]))
phosphoRS_non_missing_count <- rowSums(!is.na(PhosphoRS[, phosphoRS_sample_cols]))

phosSight_quantifiable <- sum(phosSight_non_missing_count >= quantifiable_threshold)
phosphoRS_quantifiable <- sum(phosphoRS_non_missing_count >= quantifiable_threshold)

cat(sprintf("  PhosSight quantifiable sites (>=%d samples): %d\n", quantifiable_threshold, phosSight_quantifiable))
cat(sprintf("  PhosphoRS quantifiable sites (>=%d samples): %d\n", quantifiable_threshold, phosphoRS_quantifiable))

# ===== 4. 分类可定量位点 =====
cat("\n[4] Classifying quantifiable sites...\n")

create_site_id <- function(df) {
  if ("Site" %in% colnames(df) && !all(is.na(df$Site))) {
    site_ids <- as.character(df$Site)
  } else if ("GeneID" %in% colnames(df) && "AApos" %in% colnames(df) && "AA" %in% colnames(df)) {
    site_ids <- paste(df$GeneID, df$AApos, df$AA, sep = "_")
  } else {
    site_ids <- rownames(df)
  }
  return(site_ids[!is.na(site_ids) & site_ids != ""])
}

ps_quantifiable_sites <- create_site_id(PhosSight[phosSight_non_missing_count >= quantifiable_threshold, ])
pr_quantifiable_sites <- create_site_id(PhosphoRS[phosphoRS_non_missing_count >= quantifiable_threshold, ])

ps_quantifiable_sites <- unique(ps_quantifiable_sites)
pr_quantifiable_sites <- unique(pr_quantifiable_sites)

not_identified_by_pr_count <- length(setdiff(ps_quantifiable_sites, pr_quantifiable_sites))
shared_quantifiable_count <- length(intersect(ps_quantifiable_sites, pr_quantifiable_sites))
loss_in_ps_count <- length(setdiff(pr_quantifiable_sites, ps_quantifiable_sites))

cat(sprintf("  Not identified by PhosphoRS: %d\n", not_identified_by_pr_count))
cat(sprintf("  Shared quantifiable: %d\n", shared_quantifiable_count))
cat(sprintf("  Loss in PhosSight: %d\n", loss_in_ps_count))

# ===== 5. 绘制堆叠柱状图 =====
cat("\n[5] Drawing stacked bar plot...\n")

plot_df_plot <- data.frame(
  Method = factor(c("PhosphoRS", "PhosSight"), levels = c("PhosphoRS", "PhosSight")),
  x_pos = c(0.3, 0.7),
  NotIdentified_value = c(0, not_identified_by_pr_count),
  PR_ge100_value = c(shared_quantifiable_count, shared_quantifiable_count),
  Loss_value = c(loss_in_ps_count, 0)
)

colors <- c("NotIdentified" = "#E58B00", "PR_ge100" = "#1F4E8C", "Loss" = "#7F7F7F")

p <- ggplot() +
  geom_rect(data = plot_df_plot[plot_df_plot$Method == "PhosphoRS", ],
            aes(xmin = x_pos - 0.1, xmax = x_pos + 0.1,
                ymin = -Loss_value, ymax = 0),
            fill = colors["Loss"], color = "black", linewidth = 0.5) +
  geom_rect(data = plot_df_plot[plot_df_plot$Method == "PhosphoRS", ],
            aes(xmin = x_pos - 0.1, xmax = x_pos + 0.1,
                ymin = 0, ymax = PR_ge100_value),
            fill = colors["PR_ge100"], color = "black", linewidth = 0.5) +
  geom_rect(data = plot_df_plot[plot_df_plot$Method == "PhosSight", ],
            aes(xmin = x_pos - 0.1, xmax = x_pos + 0.1,
                ymin = 0, ymax = PR_ge100_value),
            fill = colors["PR_ge100"], color = "black", linewidth = 0.5) +
  geom_rect(data = plot_df_plot[plot_df_plot$Method == "PhosSight", ],
            aes(xmin = x_pos - 0.1, xmax = x_pos + 0.1,
                ymin = PR_ge100_value, ymax = PR_ge100_value + NotIdentified_value),
            fill = colors["NotIdentified"], color = "black", linewidth = 0.5) +
  geom_text(data = plot_df_plot[plot_df_plot$Method == "PhosphoRS" & plot_df_plot$Loss_value > 0, ],
            aes(x = x_pos, y = -Loss_value - 500, label = Loss_value),
            size = 3, fontface = "bold", color = "black") +
  geom_text(data = plot_df_plot,
            aes(x = x_pos, y = PR_ge100_value / 2, label = PR_ge100_value),
            size = 3, fontface = "bold", color = "white") +
  geom_text(data = plot_df_plot[plot_df_plot$Method == "PhosSight" & plot_df_plot$NotIdentified_value > 0, ],
            aes(x = x_pos, y = PR_ge100_value + NotIdentified_value + 500, label = NotIdentified_value),
            size = 3, fontface = "bold", color = "black") +
  scale_x_continuous(
    breaks = c(0.3, 0.7),
    labels = c("PhosphoRS", "PhosSight"),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    breaks = seq(-5000, 30000, by = 5000),
    limits = c(-max(loss_in_ps_count + 1000, 5000), max(shared_quantifiable_count + not_identified_by_pr_count + 2000, 30000)),
    labels = scales::comma
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    text = element_text(size = 7),
    plot.title = element_text(size = 7, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "",
    y = "Quantifiable Sites (>=100 samples)",
    title = sprintf("Quantifiable Sites (>= %d Samples)", quantifiable_threshold)
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5)

ggsave(file.path(output_dir, "QuantifiableSites_100Samples.png"),
       plot = p, width = 4, height = 4, units = "in", dpi = 300, bg = "white")

ggsave(file.path(output_dir, "QuantifiableSites_100Samples.svg"),
       plot = p, width = 4, height = 4, units = "in", bg = "white")

cat(sprintf("  Saved: %s\n", file.path(output_dir, "QuantifiableSites_100Samples.png")))
cat(sprintf("  Saved: %s\n", file.path(output_dir, "QuantifiableSites_100Samples.svg")))

cat("\nCompleted!\n")
