# ===== PhosSight vs PhosphoRS Gain/Shared/Loss Comparison =====
# 展示两种方法的位点增益/共享/丢失对比

library(tidyverse)
library(data.table)
library(ggplot2)

cat("Creating PhosSight vs PhosphoRS Gain/Shared/Loss Comparison\n")
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

cat(sprintf("  PhosSight total sites: %d\n", nrow(PhosSight)))
cat(sprintf("  PhosphoRS total sites: %d\n", nrow(PhosphoRS)))

# ===== 2. 提取位点标识符 =====
cat("\n[2] Extracting site identifiers...\n")

create_site_id <- function(df) {
  if ("Site" %in% colnames(df) && !all(is.na(df$Site))) {
    site_ids <- as.character(df$Site)
  } else if ("GeneID" %in% colnames(df) && "AApos" %in% colnames(df) && "AA" %in% colnames(df)) {
    site_ids <- paste(df$GeneID, df$AApos, df$AA, sep = "_")
  } else {
    site_ids <- rownames(df)
  }
  return(unique(site_ids[!is.na(site_ids) & site_ids != ""]))
}

ps_sites <- create_site_id(PhosSight)
pr_sites <- create_site_id(PhosphoRS)

cat(sprintf("  PhosSight unique sites: %d\n", length(ps_sites)))
cat(sprintf("  PhosphoRS unique sites: %d\n", length(pr_sites)))

# ===== 3. 计算Gain/Shared/Loss =====
cat("\n[3] Calculating Gain/Shared/Loss...\n")

shared <- length(intersect(ps_sites, pr_sites))
gain <- length(setdiff(ps_sites, pr_sites))
loss <- length(setdiff(pr_sites, ps_sites))

cat(sprintf("  Shared: %d\n", shared))
cat(sprintf("  Gain: %d\n", gain))
cat(sprintf("  Loss: %d\n", loss))

# ===== 4. 准备绘图数据 =====
cat("\n[4] Preparing plot data...\n")

pr_total <- length(pr_sites)
scale_factor <- 100.0 / pr_total

plot_df_plot <- data.frame(
  Method = factor(c("PhosphoRS", "PhosSight"), levels = c("PhosphoRS", "PhosSight")),
  x_pos = c(0.3, 0.7),
  Shared_norm = shared * scale_factor,
  Loss_norm = c(loss * scale_factor, 0),
  Gain_norm = c(0, gain * scale_factor),
  Shared_value = shared,
  Loss_value = c(loss, 0),
  Gain_value = c(0, gain)
)

colors <- c("Loss" = "#7F7F7F", "Shared" = "#1F4E8C", "Gain" = "#E58B00")

# ===== 5. 绘制堆叠柱状图 =====
cat("\n[5] Drawing stacked bar plot...\n")

p <- ggplot() +
  geom_rect(data = plot_df_plot[plot_df_plot$Method == "PhosphoRS", ],
            aes(xmin = x_pos - 0.1, xmax = x_pos + 0.1,
                ymin = -Loss_norm, ymax = 0),
            fill = colors["Loss"], color = "black", linewidth = 0.5) +
  geom_rect(data = plot_df_plot,
            aes(xmin = x_pos - 0.1, xmax = x_pos + 0.1,
                ymin = 0, ymax = Shared_norm),
            fill = colors["Shared"], color = "black", linewidth = 0.5) +
  geom_rect(data = plot_df_plot[plot_df_plot$Method == "PhosSight", ],
            aes(xmin = x_pos - 0.1, xmax = x_pos + 0.1,
                ymin = Shared_norm, ymax = Shared_norm + Gain_norm),
            fill = colors["Gain"], color = "black", linewidth = 0.5) +
  geom_text(data = plot_df_plot[plot_df_plot$Method == "PhosphoRS" & plot_df_plot$Loss_value > 0, ],
            aes(x = x_pos, y = -Loss_norm - 6, label = Loss_value),
            size = 3, fontface = "bold", color = "black") +
  geom_text(data = plot_df_plot,
            aes(x = x_pos, y = Shared_norm / 2, label = Shared_value),
            size = 3, fontface = "bold", color = "white") +
  geom_text(data = plot_df_plot[plot_df_plot$Method == "PhosSight" & plot_df_plot$Gain_value > 0, ],
            aes(x = x_pos, y = Shared_norm + Gain_norm + 6, label = Gain_value),
            size = 3, fontface = "bold", color = "black") +
  scale_x_continuous(
    breaks = c(0.3, 0.7),
    labels = c("PhosphoRS", "PhosSight"),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    breaks = seq(-30, 160, by = 20),
    limits = c(-40, 175),
    labels = function(x) paste0(x, "%")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    text = element_text(size = 12)
  ) +
  labs(
    x = "",
    y = "Percentage of Sites (Normalized to PhosphoRS Total)",
    title = "PhosSight vs PhosphoRS: Gain/Shared/Loss"
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5)

ggsave(file.path(output_dir, "PhosSight_vs_PhosphoRS_GainSharedLoss.png"),
       plot = p, width = 6, height = 8, dpi = 300, bg = "white")

ggsave(file.path(output_dir, "PhosSight_vs_PhosphoRS_GainSharedLoss.svg"),
       plot = p, width = 6, height = 8, bg = "white")

cat(sprintf("  Saved: %s\n", file.path(output_dir, "PhosSight_vs_PhosphoRS_GainSharedLoss.png")))
cat(sprintf("  Saved: %s\n", file.path(output_dir, "PhosSight_vs_PhosphoRS_GainSharedLoss.svg")))

cat("\nCompleted!\n")
