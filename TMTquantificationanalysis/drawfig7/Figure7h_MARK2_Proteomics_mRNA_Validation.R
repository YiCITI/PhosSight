#Figure7h left
# ==============================================================================
# Script Name: MARK2 Proteomics Validation (Figure 7h Left)
# Description: Compares MARK2 protein abundance between Tumor and Normal tissues.
#              Calculates P-value, Log2 Fold Change (Diff), and 95% CI of the Difference.
#
# Usage:       Run from project root. Requires raw proteomics abundance files.
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ==============================================================================
# Configuration
# ==============================================================================
# Target Gene Info
TARGET_GENE   <- "MARK2"
TARGET_ENSEMBL <- "ENSG00000072518" 

# Input Files (Relative Paths)
base_raw  <- "./data/raw"
base_proc <- "./data/processed"

file_normal <- file.path(base_raw, "UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt")
file_tumor  <- file.path(base_raw, "UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt")

# Output Directory
output_dir <- file.path(base_proc, "Validation_Plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 1. Data Processing Function
# ==============================================================================
process_proteomics_data <- function(file_path, group_label) {
  message(paste0(">>> Reading: ", group_label, "..."))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # 1. Read Data
  dt <- fread(file_path, header = TRUE, sep = "\t")
  
  # 2. Process IDs (Remove version numbers like .15)
  # Assuming first column is 'idx'
  dt$EnsemblID_Base <- str_replace(dt$idx, "\\..*", "")
  
  # 3. Filter for Target Gene
  dt_sub <- dt %>% 
    filter(EnsemblID_Base == TARGET_ENSEMBL)
  
  if(nrow(dt_sub) == 0) {
    warning(paste0("Target gene ", TARGET_GENE, " not found in ", group_label))
    return(NULL)
  }
  
  # 4. Reshape to Long Format
  dt_long <- dt_sub %>%
    dplyr::select(-idx, -EnsemblID_Base) %>% 
    pivot_longer(
      cols = everything(),   # All remaining columns are samples
      names_to = "Sample", 
      values_to = "Intensity"
    ) %>%
    mutate(
      GeneSymbol = TARGET_GENE,
      Intensity = as.numeric(Intensity),
      Group = group_label
    ) %>%
    filter(!is.na(Intensity)) 
  
  return(dt_long)
}

# ==============================================================================
# 2. Statistical Analysis Function (Log2FC & Diff CI)
# ==============================================================================
calc_stats <- function(data) {
  # Check sample size
  if(length(unique(data$Group)) < 2 || min(table(data$Group)) < 2) {
    return(data.frame(Label = "Insufficient Data", MaxY = max(data$Intensity, na.rm=T)))
  }
  
  # Separate Vectors
  vec_tumor <- data$Intensity[data$Group == "Tumor"]
  vec_normal <- data$Intensity[data$Group == "Normal"]
  
  # Welch's T-test
  # By default, t.test calculates the Confidence Interval for the difference in means
  t_res <- t.test(vec_tumor, vec_normal)
  
  # 1. P-value
  p_val <- t_res$p.value
  
  # 2. Log2FC (Difference of Means)
  # t_res$estimate[1] is Mean of X (Tumor), [2] is Mean of Y (Normal)
  diff_mean <- t_res$estimate[1] - t_res$estimate[2] 
  
  # 3. 95% Confidence Interval of the Difference
  ci_lower <- t_res$conf.int[1]
  ci_upper <- t_res$conf.int[2]
  
  # Format Label
  label_text <- paste0(
    sprintf("p-value=%.1e", p_val), "\n",
    sprintf("Diff (Log2FC)=%.2f", diff_mean), "\n",
    sprintf("95%%C.I.=[%.2f, %.2f]", ci_lower, ci_upper)
  )
  
  return(data.frame(Label = label_text, MaxY = max(data$Intensity, na.rm = TRUE)))
}

# ==============================================================================
# 3. Main Execution
# ==============================================================================
data_normal <- process_proteomics_data(file_normal, "Normal")
data_tumor  <- process_proteomics_data(file_tumor,  "Tumor")

if(!is.null(data_normal) && !is.null(data_tumor)) {
  
  # Combine
  plot_df_all <- bind_rows(data_normal, data_tumor)
  
  # Set Factor Order: Tumor (Red) First, Normal (Blue) Second
  plot_df_all$Group <- factor(plot_df_all$Group, levels = c("Tumor", "Normal"))
  
  message(paste0(">>> Data ready for ", TARGET_GENE, ". Sample counts:"))
  print(table(plot_df_all$Group))
  
  # Calculate Statistics
  stats_df <- calc_stats(plot_df_all)
  
  # ==============================================================================
  # 4. Visualization
  # ==============================================================================
  message(">>> Generating Plot...")
  
  cols <- c("Tumor" = "#E41A1C", "Normal" = "#377EB8") 
  
  p <- ggplot(plot_df_all, aes(x = Group, y = Intensity)) +
    # Boxplot
    geom_boxplot(aes(fill = Group), width = 0.5, outlier.shape = NA, alpha = 0.8) +
    
    # Jitter Points
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
    
    # Stats Label
    geom_text(data = stats_df, 
              aes(x = 1.5, 
                  y = MaxY + (MaxY - min(plot_df_all$Intensity)) * 0.1, 
                  label = Label), 
              size = 4.5, lineheight = 1.1) +
    
    # Styling
    scale_fill_manual(values = cols) +
    theme_classic() +
    
    labs(
      title = TARGET_GENE, 
      y = "Log2 Protein Intensity", 
      x = ""
    ) +
    
    theme(
      # Box border style
      panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
      axis.line = element_blank(), # Remove default axis lines in favor of border
      
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "none",
      axis.text = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold")
    ) +
    
    # Expand Y axis for label space
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.4)))
  
  # Save
  output_file <- file.path(output_dir, paste0("Figure7h_", TARGET_GENE, "_Proteomics_Validation.pdf"))
  ggsave(output_file, p, width = 5, height = 6)
  message(paste("Plot saved to:", output_file))
  
} else {
  stop("Failed to merge Tumor and Normal data. Please check input files.")
}


#Figure7h right
# ==============================================================================
# Script Name: MARK2 mRNA Validation (Figure 7h Right)
# Description: Compares MARK2 mRNA expression (TPM) between Tumor and Normal tissues.
#              Reads from Excel, calculates P-value, Diff (Log2FC), and 95% CI.
#
# Usage:       Run from project root. Requires mRNA TPM Excel file.
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)

# ==============================================================================
# Configuration
# ==============================================================================
# Target Gene Info
TARGET_GENE     <- "MARK2"
TARGET_ENS_BASE <- "ENSG00000072518" 

# Input Files (Relative Paths)
base_raw  <- "./data/raw"
base_proc <- "./data/processed"

file_mrna <- file.path(base_raw, "UCEC-gene_TPM_removed_circRNA_tumor_normal_raw_log2(x+1)_BCM.xlsx")

# Output Directory
output_dir <- file.path(base_proc, "Validation_Plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 1. Data Loading & Processing
# ==============================================================================
message(paste0(">>> Reading mRNA Data from: ", basename(file_mrna)))

if (!file.exists(file_mrna)) stop("Input Excel file not found.")

# 1. Read Excel
dt <- read_excel(file_mrna)

# [CRITICAL]: Handle duplicate column names if they exist in Excel
colnames(dt) <- make.unique(colnames(dt))

message(paste("Data loaded. Dimensions:", paste(dim(dt), collapse = " x ")))

# 2. Process Gene IDs
# Assume first column is the ID column
id_col <- colnames(dt)[1] 
dt$EnsemblID_Base <- str_replace(dt[[id_col]], "\\..*", "")

# 3. Filter for Target Gene
dt_sub <- dt %>% filter(EnsemblID_Base == TARGET_ENS_BASE)

if(nrow(dt_sub) == 0) {
  stop(paste("Target gene", TARGET_GENE, "not found in dataset."))
}

# 4. Reshape to Long Format
# Exclude ID columns, keep all sample columns
cols_to_exclude <- c(id_col, "EnsemblID_Base")
sample_cols <- setdiff(colnames(dt), cols_to_exclude)

plot_df <- dt_sub %>%
  dplyr::select(all_of(sample_cols)) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Intensity") %>%
  mutate(Intensity = as.numeric(Intensity)) %>%
  filter(!is.na(Intensity))

# 5. Group Samples (Tumor vs Normal)
# Logic: Samples ending in "_T", "_T.1", etc. are Tumor. Others are Normal.
plot_df <- plot_df %>%
  mutate(
    Group = ifelse(grepl("_T($|\\.|\\d)", Sample), "Tumor", "Normal")
  )

# Set Factor Order: Tumor (Red) First
plot_df$Group <- factor(plot_df$Group, levels = c("Tumor", "Normal"))

message(">>> Sample Grouping Summary:")
print(table(plot_df$Group))

# ==============================================================================
# 2. Statistical Analysis
# ==============================================================================
calc_stats <- function(data) {
  # Check sample size
  if(length(unique(data$Group)) < 2 || min(table(data$Group)) < 2) {
    return(data.frame(Label = "Insufficient Data", MaxY = max(data$Intensity, na.rm=T)))
  }
  
  # Separate vectors
  vec_tumor <- data$Intensity[data$Group == "Tumor"]
  vec_normal <- data$Intensity[data$Group == "Normal"]
  
  # Welch's T-test
  t_res <- t.test(vec_tumor, vec_normal)
  
  # 1. P-value
  p_val <- t_res$p.value
  
  # 2. Diff (Log2FC equivalent for Log2 data)
  diff_mean <- t_res$estimate[1] - t_res$estimate[2] 
  
  # 3. 95% CI
  ci_lower <- t_res$conf.int[1]
  ci_upper <- t_res$conf.int[2]
  
  # Format Label
  label_text <- paste0(
    sprintf("p-value=%.1e", p_val), "\n",
    sprintf("Diff (Log2FC)=%.2f", diff_mean), "\n",
    sprintf("95%%C.I.=[%.2f, %.2f]", ci_lower, ci_upper)
  )
  
  return(data.frame(Label = label_text, MaxY = max(data$Intensity, na.rm = TRUE)))
}

stats_df <- calc_stats(plot_df)

# ==============================================================================
# 3. Visualization
# ==============================================================================
message(">>> Generating Plot...")

cols <- c("Tumor" = "#E41A1C", "Normal" = "#377EB8")

p <- ggplot(plot_df, aes(x = Group, y = Intensity)) +
  
  # Boxplot
  geom_boxplot(aes(fill = Group), width = 0.5, outlier.shape = NA, alpha = 0.8) +
  
  # Jitter Points
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
  
  # Stats Label
  geom_text(data = stats_df, 
            aes(x = 1.5, 
                y = MaxY + (MaxY - min(plot_df$Intensity)) * 0.1, 
                label = Label), 
            size = 4.5, lineheight = 1.1) +
  
  # Styling
  scale_fill_manual(values = cols) +
  theme_classic() +
  
  labs(
    title = TARGET_GENE, 
    y = "Log2 (TPM + 1)", 
    x = ""
  ) +
  
  theme(
    # Box border style
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    axis.line = element_blank(),
    
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold")
  ) +
  
  # Expand Y axis for label space
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.4)))

# Save
output_file <- file.path(output_dir, paste0("Figure7h_", TARGET_GENE, "_mRNA_Validation.pdf"))
ggsave(output_file, p, width = 5, height = 6)
message(paste("Plot saved to:", output_file))
