#Figure7f
# ==============================================================================
# Script Name: Visualization of Representative Differentially Phosphorylated Sites
# Description: Generates Figure 7f - Boxplots comparing intensity levels between 
#              Tumor and Normal tissues for specific high-confidence markers.
#              (Calculates P-value, Log2FC, and 95% C.I.)
#
# Usage:       Run from project root. Requires Metadata and Table X2 (Preprocessed).
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# ==============================================================================
# Configuration
# ==============================================================================
# Target Phosphosites to Plot
TARGET_IDS <- c("ESR1_167S", "GSK3B_9S", "AKT1_126S", "EZH2_367T")

# Input Files
base_raw  <- "./data/raw"
base_proc <- "./data/processed"

meta_file   <- file.path(base_raw, "mmc1.xlsx")
matrix_file <- file.path(base_proc, "Table_X2_PhosSight_50PercentMissingValueCutoff.txt")

# Output File
output_file <- file.path(base_proc, "Figure7f_Differential_Boxplots.pdf")

# ==============================================================================
# 1. Load and Clean Metadata
# ==============================================================================
message(">>> Loading Metadata...")

if (!file.exists(meta_file)) stop("Metadata file not found.")

meta_raw <- read_excel(meta_file)

meta_clean <- meta_raw %>%
  dplyr::select(
    Sample_List = Proteomics_Parent_Sample_IDs, 
    Type_Raw = Proteomics_Tumor_Normal
  ) %>%
  separate_rows(Sample_List, sep = ",") %>%
  mutate(Sample = trimws(Sample_List)) %>%
  mutate(Group = ifelse(Type_Raw == "Tumor", "Tumor", "Normal")) %>%
  dplyr::select(Sample, Group) %>%
  distinct() %>%
  na.omit()

message("Metadata Group Counts:")
print(table(meta_clean$Group))

# ==============================================================================
# 2. Load Matrix and Extract Targets
# ==============================================================================
message(">>> Loading Expression Matrix...")

if (!file.exists(matrix_file)) stop("Matrix file not found.")

raw_data <- fread(matrix_file, header = TRUE, sep = "\t")

# Ensure ID column exists
if (!"PhosphoSite" %in% colnames(raw_data)) stop("Column 'PhosphoSite' missing.")
raw_data$ID <- raw_data$PhosphoSite

# Extract Data for Targets
message(">>> Extracting Target Sites...")

# Identify sample columns
meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
sample_cols <- setdiff(colnames(raw_data), meta_cols)

# Filter and Reshape
matrix_long <- raw_data %>%
  filter(ID %in% TARGET_IDS) %>%
  dplyr::select(ID, all_of(sample_cols)) %>%
  pivot_longer(cols = -ID, names_to = "Sample", values_to = "Intensity") %>%
  mutate(Intensity = as.numeric(Intensity)) %>%
  filter(!is.na(Intensity)) # Keep negative values (Data is Median Centered)

# ==============================================================================
# 3. Merge Data
# ==============================================================================
plot_df <- matrix_long %>%
  inner_join(meta_clean, by = "Sample") 

# Set Factor Levels (Tumor first usually, or Normal first depending on preference)
# Here we set Tumor first to match common clinical plots
plot_df$Group <- factor(plot_df$Group, levels = c("Tumor", "Normal"))

message("Final Plotting Data Summary:")
print(table(plot_df$Group))

# ==============================================================================
# 4. Statistical Analysis
# ==============================================================================
# Note: Input data is already Log2 transformed.
# Log2FC = Mean(Tumor) - Mean(Normal)

calc_stats <- function(data) {
  # Safety check for sufficient data
  if(length(unique(data$Group)) < 2 || min(table(data$Group)) < 2) {
    return(data.frame(
      Label = "Insufficient Data", 
      MaxY = max(data$Intensity, na.rm = TRUE)
    ))
  }
  
  # Welch's T-test
  t_res <- t.test(Intensity ~ Group, data = data)
  p_val <- t_res$p.value
  
  # Calculate Means & Log2FC
  # t_res$estimate contains means of groups. 
  # Since Group factor is Tumor, Normal -> estimate[1] is Tumor
  mean_tumor <- t_res$estimate[1]
  mean_normal <- t_res$estimate[2]
  log2fc <- mean_tumor - mean_normal
  
  # 95% Confidence Interval of the DIFFERENCE
  ci_lower <- t_res$conf.int[1]
  ci_upper <- t_res$conf.int[2]
  
  # Format Label Text
  # Using scientific notation for small p-values
  p_text <- sprintf("p-value=%.1e", p_val)
  fc_text <- sprintf("Log2FC=%.2f", log2fc)
  ci_text <- sprintf("95%%C.I.=[%.2f, %.2f]", ci_lower, ci_upper)
  
  return(data.frame(
    Label = paste(p_text, fc_text, ci_text, sep = "\n"),
    MaxY = max(data$Intensity, na.rm = TRUE)
  ))
}

# Apply stats function to each ID
stats_df <- plot_df %>%
  group_by(ID) %>%
  do(calc_stats(.))

# ==============================================================================
# 5. Plotting
# ==============================================================================
message(">>> Generating Boxplots...")

# Define Colors (Red for Tumor, Blue for Normal)
cols <- c("Tumor" = "#E41A1C", "Normal" = "#377EB8")

p <- ggplot(plot_df, aes(x = Group, y = Intensity)) +
  # Boxplot
  geom_boxplot(aes(fill = Group), width = 0.6, outlier.shape = NA, alpha = 0.8) +
  
  # Jitter Points (Individual samples)
  geom_jitter(width = 0.2, size = 1, alpha = 0.4, color = "black") +
  
  # Faceting
  facet_wrap(~ID, scales = "free_y", nrow = 1) +
  
  # Statistical Annotations
  geom_text(data = stats_df, 
            aes(x = 1.5, 
                y = MaxY + (MaxY - min(plot_df$Intensity)) * 0.15, 
                label = Label), 
            size = 3.5, lineheight = 1.1) +
  
  # Aesthetics
  scale_fill_manual(values = cols) +
  theme_classic() +
  
  labs(
    y = "Log2 Intensity (Median Centered)", 
    x = ""
  ) +
  
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    # Remove X axis labels if redundant with legend, or keep them
    axis.text.x = element_text(angle = 0, hjust = 0.5) 
  ) +
  
  # Expand Y axis to fit text
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45)))

# ==============================================================================
# 6. Save Output
# ==============================================================================
# Ensure output directory exists
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

ggsave(output_file, plot = p, width = 10, height = 5) # Adjusted width for 1 row

message(paste("Plot saved to:", output_file))
