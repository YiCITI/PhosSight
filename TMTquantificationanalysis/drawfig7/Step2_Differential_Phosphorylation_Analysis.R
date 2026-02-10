# Step2 肿瘤组织与正常组织位点差异分析
# ==============================================================================
# Script Name: Differential Phosphorylation Analysis
# Description: This script performs differential expression analysis between 
#              Tumor and Normal samples using Welch's t-test.
#              1. Cleans metadata and aligns samples.
#              2. Calculates Log2 Fold Change (Log2FC).
#              3. Computes P-values (Welch's t-test) and FDR (Benjamini-Hochberg).
#
# Usage:       Run this script from the project root directory.
#              Requires: Metadata file (xlsx) and Preprocessed Data matrices (txt).
# ==============================================================================

library(data.table)
library(dplyr)
library(readxl)
library(tidyr)
library(stringr)

# ==============================================================================
# 1. Prepare Metadata (Clinical Grouping)
# ==============================================================================

# Define path to metadata (Adjust filename if necessary)
# Assumption: Metadata is stored in the 'data/raw' folder
meta_file <- "./data/raw/mmc1.xlsx" 

if (!file.exists(meta_file)) {
  stop(paste("Error: Metadata file not found at", meta_file))
}

message("Reading and cleaning metadata...")

meta_raw <- read_excel(meta_file)

# Cleaning and Grouping Logic
meta_clean <- meta_raw %>%
  dplyr::select(
    Sample_List = Proteomics_Parent_Sample_IDs, 
    Type_Raw = Proteomics_Tumor_Normal
  ) %>%
  # [CRITICAL]: Handle comma-separated sample IDs (e.g., "ID1, ID2")
  separate_rows(Sample_List, sep = ",") %>%
  mutate(Sample = trimws(Sample_List)) %>% # Remove potential whitespace
  # Define Groups: "Tumor" stays "Tumor", everything else is "Normal"
  mutate(Group = ifelse(Type_Raw == "Tumor", "Tumor", "Normal")) %>%
  dplyr::select(Sample, Group) %>%
  distinct()

# Check group statistics
message("Metadata Group Statistics:")
print(table(meta_clean$Group))


# ==============================================================================
# 2. Define Differential Analysis Function
# ==============================================================================
run_diff_analysis <- function(input_file, output_file, metadata) {
  
  message(paste0("\n------------------------------------------------"))
  message(paste0("Analyzing file: ", basename(input_file)))
  
  # 1. Read Expression Matrix
  # ----------------------------------------------------------------------------
  if (!file.exists(input_file)) {
    stop(paste("Error: Input file not found:", input_file))
  }
  
  data <- fread(input_file, header = TRUE, sep = "\t")
  
  # Identify sample columns (exclude non-numeric annotation columns)
  meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
  existing_meta_cols <- intersect(names(data), meta_cols)
  sample_cols <- setdiff(names(data), existing_meta_cols)
  
  # 2. Match Samples (Alignment)
  # ----------------------------------------------------------------------------
  # Find samples present in both the matrix and the metadata
  valid_samples <- intersect(sample_cols, metadata$Sample)
  
  if(length(valid_samples) < 4) {
    stop("Error: Too few valid matched samples (<4). Please check sample naming formats!")
  }
  
  # Extract sub-matrix for calculation
  mat <- as.matrix(data[, ..valid_samples])
  
  # Align group vector
  # Match ensures the order of group_vec corresponds exactly to the columns of mat
  group_vec <- metadata$Group[match(valid_samples, metadata$Sample)]
  
  message(paste0("Number of samples included in analysis: ", length(valid_samples)))
  message(paste0("Tumor: ", sum(group_vec == "Tumor"), 
                 ", Normal: ", sum(group_vec == "Normal")))
  
  # 3. Loop Calculation (T-test & Log2FC)
  # ----------------------------------------------------------------------------
  # Pre-allocate memory
  results_list <- vector("list", nrow(data))
  
  # Pre-calculate indices for speed
  idx_tumor <- which(group_vec == "Tumor")
  idx_normal <- which(group_vec == "Normal")
  
  message("Performing differential calculations (this may take a moment)...")
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)
  
  for (i in 1:nrow(data)) {
    # Extract one row of data
    vals <- mat[i, ]
    
    # Split data by group
    vals_tumor <- vals[idx_tumor]
    vals_normal <- vals[idx_normal]
    
    # Remove NAs
    vals_tumor <- vals_tumor[!is.na(vals_tumor)]
    vals_normal <- vals_normal[!is.na(vals_normal)]
    
    # Check if sample size is sufficient (at least 2 per group)
    if (length(vals_tumor) < 2 || length(vals_normal) < 2) {
      res <- c(NA, NA, NA, NA) # Mean_T, Mean_N, FC, Pval
    } else {
      # Calculate Means
      mean_t <- mean(vals_tumor)
      mean_n <- mean(vals_normal)
      
      # Calculate Log2 FoldChange
      # Since data is already Log2 transformed: Log2FC = Mean(Tumor) - Mean(Normal)
      log2fc <- mean_t - mean_n 
      
      # T-test (Two-sided, Welch's t-test for unequal variance)
      # tryCatch prevents the loop from breaking if data is constant
      p_val <- tryCatch({
        t.test(vals_tumor, vals_normal)$p.value
      }, error = function(e) NA)
      
      res <- c(mean_t, mean_n, log2fc, p_val)
    }
    
    results_list[[i]] <- res
    if(i %% 1000 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 4. Compile Results
  # ----------------------------------------------------------------------------
  results_df <- do.call(rbind, results_list)
  colnames(results_df) <- c("Mean_Tumor", "Mean_Normal", "log2FC", "pvalue")
  
  # Merge back with original metadata (GeneSymbol, PhosphoSite, etc.)
  final_df <- cbind(data[, ..existing_meta_cols], as.data.frame(results_df))
  
  # 5. Calculate FDR (BH adjustment)
  # ----------------------------------------------------------------------------
  final_df$adjust_pvalue <- NA
  valid_idx <- !is.na(final_df$pvalue)
  final_df$adjust_pvalue[valid_idx] <- p.adjust(final_df$pvalue[valid_idx], method = "BH")
  
  # 6. Generate Linear Fold Change (Optional)
  # ----------------------------------------------------------------------------
  # Calculated as 2^Log2FC for easier interpretation
  final_df$fold_change <- 2^final_df$log2FC
  
  # Sort by P-value
  final_df <- final_df %>% arrange(pvalue)
  
  # 7. Save Results
  # ----------------------------------------------------------------------------
  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  write.table(final_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message(paste0("Analysis complete! Results saved to: ", output_file))
  
  # Print summary of significant results (P < 0.05)
  n_sig <- sum(final_df$pvalue < 0.05, na.rm = TRUE)
  message(paste0("Number of significant sites (P < 0.05): ", n_sig))
}


# ==============================================================================
# Execution
# ==============================================================================

# Define directories using relative paths
proc_dir <- "./data/processed"

# 1. PhosSight Differential Analysis
# ----------------------------------
# Input comes from the output of the previous preprocessing script (Table X2)
input_ps <- file.path(proc_dir, "Table_X2_PhosSight_50PercentMissingValueCutoff.txt")
output_ps <- file.path(proc_dir, "Table_X3_PhosSight_DifferentialAnalysis.txt")

run_diff_analysis(input_ps, output_ps, meta_clean)

# 2. PhosphoRS Differential Analysis
# ----------------------------------
input_rs <- file.path(proc_dir, "Table_X2_PhosphoRS_50PercentMissingValueCutoff.txt")
output_rs <- file.path(proc_dir, "Table_X3_PhosphoRS_DifferentialAnalysis.txt")

run_diff_analysis(input_rs, output_rs, meta_clean)
