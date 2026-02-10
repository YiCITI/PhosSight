#Figure7b
# ==============================================================================
# Script Name: Data Completeness Analysis (Phosphosite Coverage)
# Description: Generates a completeness curve showing the number of identified 
#              phosphosites as a function of the non-missing value cutoff.
#              (Corresponds to Figure 7b)
#
# Usage:       Run from project root. Requires processed matrices (Table X1/X2).
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)

# ==============================================================================
# Define Processing Function
# ==============================================================================
process_matrix_completeness <- function(file_path, method_name, step_size = 10) {
  
  message(paste0("Processing: ", method_name, "..."))
  
  # 1. Read Data
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  df <- fread(file_path, header = TRUE, sep = "\t")
  
  # 2. Identify Sample Columns (Robust method)
  # Exclude standard metadata columns to identify sample intensity columns
  meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
  sample_cols <- setdiff(colnames(df), meta_cols)
  
  # Extract numerical matrix
  mat <- as.matrix(df[, ..sample_cols])
  
  total_samples <- ncol(mat)
  message(paste0("  - Total Samples: ", total_samples))
  
  # 3. Calculate Non-missing Values per Site
  # rowSums(!is.na) counts how many samples have a valid value for each site
  n_detected <- rowSums(!is.na(mat))
  
  # 4. Generate Cutoff Sequence
  # Sequence from 0 to Total Samples with specified step_size
  cutoffs <- seq(0, total_samples, by = step_size)
  
  # Optional: Ensure the exact maximum (total samples) is included if not covered by step
  if (tail(cutoffs, 1) != total_samples) {
    cutoffs <- c(cutoffs, total_samples)
  }
  
  # 5. Calculate Counts for each Cutoff
  counts <- numeric(length(cutoffs))
  
  for (i in seq_along(cutoffs)) {
    cutoff_val <- cutoffs[i]
    # Count sites present in at least 'cutoff_val' samples
    counts[i] <- sum(n_detected >= cutoff_val)
  }
  
  # Return formatted data frame
  return(data.frame(
    Cutoff = cutoffs,
    Phosphosites = counts,
    Method = method_name
  ))
}

# ==============================================================================
# Main Execution
# ==============================================================================

# Define relative paths
# Note: Input usually comes from Original Matrix (Table X1) to show full drop-off
base_dir <- "./data/raw" 
output_file <- "./data/processed/Figure7b_Phosphosite_Coverage.pdf"

# Input Files
file_phossight <- file.path(base_dir, "Table_X1_PhosSight_OriginalMatrix.txt")
file_phosphors <- file.path(base_dir, "Table_X1_PhosphoRS_OriginalMatrix.txt")

# Process Data
plot_data_1 <- process_matrix_completeness(file_phossight, "PhosSight", step_size = 10)
plot_data_2 <- process_matrix_completeness(file_phosphors, "PhosphoRS", step_size = 10)

# Merge
plot_data <- rbind(plot_data_1, plot_data_2)

# ==============================================================================
# Plotting
# ==============================================================================

if (!is.null(plot_data)) {
  
  # Define Colors (Using your specified hex codes)
  custom_colors <- c("PhosSight" = "#D68D8E", "PhosphoRS" = "#8EBFD4")
  
  p <- ggplot(plot_data, aes(x = Phosphosites, y = Cutoff, color = Method)) +
    
    # Line
    geom_path(size = 1.2) + # geom_path follows the order of data
    
    # Points (optional, can adjust size or remove if too crowded)
    geom_point(size = 2.5, alpha = 0.9) +
    
    # 50% Cutoff Line (92 samples)
    geom_hline(yintercept = 92, linetype = "dashed", color = "grey50", size = 0.6) +
    
    # Text Annotation for 50% line
    annotate("text", x = max(plot_data$Phosphosites)*0.85, y = 98, 
             label = "50% Cutoff (92 samples)", color = "grey40", size = 3.5, hjust = 0) +
    
    # Colors
    scale_color_manual(values = custom_colors) +
    
    # Labels
    labs(x = "# Identified Phosphosites", 
         y = "Non-missing Value Cutoff (# Samples)") +
    
    # Theme
    theme_classic() +
    theme(
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.line = element_line(size = 0.8, color = "black"),
      legend.position = c(0.80, 0.60), # Adjusted position
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.background = element_rect(fill = "transparent")
    )
  
  # Save
  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ggsave(output_file, plot = p, width = 6, height = 6)
  
  message(paste("Plot saved to:", output_file))
  
} else {
  message("Error: No data available for plotting.")
}











