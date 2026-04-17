# Step1 数据预处理
# ==============================================================================
# Script Name: Phosphoproteomics Data Preprocessing
# Description: This script performs standard preprocessing steps for TMT 
#              phosphoproteomics data, including:
#              1. Filtering phosphosites with >50% missing values.
#              2. Log2 transformation (if input is raw intensity).
#              3. Median Centering Normalization (sample-wise).
#
# Usage:       Run this script from the project root directory.
#              Input:  ./data/raw/
#              Output: ./data/processed/
# ==============================================================================
library(data.table)
library(dplyr)

# ==============================================================================
# Function Definition: Filter NA -> Log2 -> Median Center
# ==============================================================================
process_proteomics_data <- function(input_file, output_file) {
  
  message(paste0("\n>>> Processing file: ", basename(input_file)))
  
  # 1. Check Input File
  # ----------------------------------------------------------------------------
  if (!file.exists(input_file)) {
    stop(paste("Error: Input file not found:", input_file, 
               "\nPlease ensure the working directory is set to the project root."))
  }
  
  # Ensure the output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    message(paste("Creating output directory:", output_dir))
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 2. Read Data
  # ----------------------------------------------------------------------------
  dt <- fread(input_file, header = TRUE, sep = "\t")
  
  # Identify metadata columns
  # Standard metadata columns to exclude from numerical operations
  meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
  
  # Identify existing metadata columns in the current dataset
  existing_meta <- intersect(names(dt), meta_cols)
  
  # The remaining columns are treated as sample intensity columns
  sample_cols <- setdiff(names(dt), existing_meta)
  
  message(paste0("Detected sample count: ", length(sample_cols)))
  message(paste0("Original number of sites: ", nrow(dt)))
  
  # Extract the numerical matrix for processing
  mat <- as.matrix(dt[, ..sample_cols])
  
  # 3. 50% Missing Value Cutoff
  # ----------------------------------------------------------------------------
  # Calculate the number of NAs per row (phosphosite)
  na_counts <- rowSums(is.na(mat))
  
  # Threshold: Retain sites present in at least 50% of samples
  # (i.e., missing values must be <= 50% of total samples)
  threshold <- ncol(mat) * 0.5
  
  # Identify indices to keep
  keep_indices <- which(na_counts <= threshold)
  
  # Filter the matrix and metadata
  mat_filtered <- mat[keep_indices, ]
  dt_meta_filtered <- dt[keep_indices, ..existing_meta]
  
  removal_rate <- round((1 - nrow(mat_filtered) / nrow(dt)) * 100, 2)
  message(paste0("Sites retained after filtering: ", nrow(mat_filtered), 
                 " (Removal rate: ", removal_rate, "%)"))
  
  # 4. Log2 Transformation
  # ----------------------------------------------------------------------------
  # Check if data requires log transformation.
  # Assumption: If max value > 100, data is likely raw intensity.
  if (max(mat_filtered, na.rm = TRUE) > 100) {
    message("Performing Log2 transformation...")
    mat_log <- log2(mat_filtered)
    
    # Handle potential -Inf values (log2(0)) by converting them to NA
    mat_log[is.infinite(mat_log)] <- NA
  } else {
    message("Warning: Data values appear small. Assuming already log-transformed. Skipping Log2 step.")
    mat_log <- mat_filtered
  }
  
  # 5. Median Center Normalization
  # ----------------------------------------------------------------------------
  message("Performing Median Center Normalization...")
  
  # Calculate the median for each sample (column), ignoring NAs
  col_medians <- apply(mat_log, 2, median, na.rm = TRUE)
  
  # Print median range for quality check
  message(paste("Range of sample medians: ", round(min(col_medians), 2), 
                " - ", round(max(col_medians), 2)))
  
  # Center the data: Subtract the column median from each value
  mat_norm <- sweep(mat_log, 2, col_medians, FUN = "-")
  
  # 6. Merge and Save
  # ----------------------------------------------------------------------------
  # Combine metadata with the normalized matrix
  final_dt <- cbind(dt_meta_filtered, as.data.table(mat_norm))
  
  # Write to file
  fwrite(final_dt, output_file, sep = "\t", quote = FALSE, na = "NA")
  
  message(paste0("Processing complete! Output saved to: ", output_file))
}

# ==============================================================================
# Execution
# ==============================================================================

# Define relative paths
# NOTE: Ensure your R session working directory is the project root.

# Define directories
raw_dir <- "./data/raw"
proc_dir <- "./data/processed"

# Task 1: Process PhosphoRS Data
# ------------------------------
input_rs <- file.path(raw_dir, "Table_X1_PhosphoRS_OriginalMatrix.txt")
output_rs <- file.path(proc_dir, "Table_X2_PhosphoRS_50PercentMissingValueCutoff.txt")

process_proteomics_data(input_rs, output_rs)

# Task 2: Process PhosSight Data
# ------------------------------
input_ps <- file.path(raw_dir, "Table_X1_PhosSight_OriginalMatrix.txt")
output_ps <- file.path(proc_dir, "Table_X2_PhosSight_50PercentMissingValueCutoff.txt")

process_proteomics_data(input_ps, output_ps)