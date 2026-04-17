# Step4 ssKSEA分析
# ==============================================================================
# Script Name: Single-Sample Kinase-Substrate Enrichment Analysis (ssKSEA)
# Description: This script infers upstream kinase activities from phosphoproteomics data.
#              1. Loads and cleans the kinase-substrate database (KSData from KSEAapp).
#              2. Performs row-wise Z-score normalization on input data.
#              3. Calculates the mean Z-score of substrates for each kinase per sample.
#
# Usage:       Run this script from the project root directory.
#              Requires: Preprocessed Data (txt) from Step 1.
# ==============================================================================

library(data.table)
library(dplyr)
library(KSEAapp)

# ==============================================================================
# Define Main ssKSEA Function
# ==============================================================================
run_ssksea <- function(input_file, output_file) {
  
  message(paste0("\n----------------------------------------------------------"))
  message(paste0("Running ssKSEA on: ", basename(input_file)))
  
  # 1. Prepare Database
  # ----------------------------------------------------------------------------
  message("Loading and cleaning kinase-substrate database...")
  
  # Load KSData from KSEAapp package
  data(KSData) 
  ks_db_raw <- KSData
  
  # Clean and format database
  # Goal: Create Site_ID format (e.g., "SUBSTRATE_123S") to match input data
  ks_db <- ks_db_raw %>%
    filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human") %>%
    dplyr::select(
      Kinase = GENE,
      Substrate = SUB_GENE,
      Residue_Raw = SUB_MOD_RSD 
    ) %>%
    mutate(
      Residue_AA = substr(Residue_Raw, 1, 1), # Extract Amino Acid (e.g., S)
      Position = substring(Residue_Raw, 2),   # Extract Position (e.g., 315)
      # Construct ID: Substrate_PositionResidue (e.g., WRIP1_151S)
      Site_ID = paste0(Substrate, "_", Position, Residue_AA)
    ) %>%
    dplyr::select(Kinase, Site_ID) %>%
    distinct()
  
  # Filter kinases with fewer than 3 known substrates in the database
  valid_kinases <- ks_db %>% count(Kinase) %>% filter(n >= 3) %>% pull(Kinase)
  ks_db <- ks_db %>% filter(Kinase %in% valid_kinases)
  
  message(paste("Database ready. Contains", length(unique(ks_db$Kinase)), "kinases."))
  
  # 2. Read Data
  # ----------------------------------------------------------------------------
  if (!file.exists(input_file)) stop(paste("Error: Input file not found:", input_file))
  
  raw_data <- fread(input_file, header = TRUE, sep = "\t")
  
  # Define metadata columns to exclude
  meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
  sample_cols <- setdiff(colnames(raw_data), meta_cols)
  
  message(paste("Detected sample count:", length(sample_cols)))
  
  # Extract numerical matrix
  mat <- as.matrix(raw_data[, ..sample_cols])
  
  # Set row names using PhosphoSite column (e.g., WRIP1_151S)
  if(!"PhosphoSite" %in% colnames(raw_data)) {
    stop("Error: Input file must contain a 'PhosphoSite' column for ID matching.")
  }
  rownames(mat) <- raw_data$PhosphoSite
  
  # 3. Data Preprocessing (Z-score)
  # ----------------------------------------------------------------------------
  message("Preprocessing data...")
  
  # Logic check: If values are small, assume they are already Log2 transformed
  # (Input from Step 1 is typically already Log2 + Median Centered)
  if(max(mat, na.rm=TRUE) > 50) { 
    message("Large values detected. Performing Log2 transformation...")
    mat <- log2(mat + 1)
  } else {
    message("Data appears to be Log2 transformed. Skipping Log2 step.")
  }
  
  # [CRITICAL]: Row-wise Z-score Standardization
  # Even if data is median-centered by sample, ssKSEA requires that 
  # each phosphosite is standardized across samples so substrates contribute equally.
  message("Performing Row-wise Z-score Standardization...")
  mat_z <- t(scale(t(mat)))
  
  # 4. Calculate ssKSEA Scores
  # ----------------------------------------------------------------------------
  message("Calculating Single-Sample KSEA scores...")
  
  kinases <- unique(ks_db$Kinase)
  ksea_matrix <- matrix(NA, nrow = length(kinases), ncol = ncol(mat_z))
  rownames(ksea_matrix) <- kinases
  colnames(ksea_matrix) <- colnames(mat_z)
  
  count <- 0
  for (k in kinases) {
    # Find substrates for current kinase
    targets <- ks_db$Site_ID[ks_db$Kinase == k]
    
    # Match with available data
    matched_targets <- intersect(targets, rownames(mat_z))
    
    # Calculate Mean Z-score if enough substrates are found (Threshold >= 3)
    if (length(matched_targets) >= 3) {
      sub_mat <- mat_z[matched_targets, , drop = FALSE]
      # Mean of substrates for each sample, ignoring NAs
      ksea_matrix[k, ] <- colMeans(sub_mat, na.rm = TRUE)
    }
    
    count <- count + 1
    if(count %% 50 == 0) message(paste("Processed", count, "kinases..."))
  }
  
  # Remove rows that are entirely NA (Kinases not inferred)
  valid_rows <- rowSums(is.na(ksea_matrix)) != ncol(ksea_matrix)
  ksea_matrix <- ksea_matrix[valid_rows, ]
  
  # 5. Save Results
  # ----------------------------------------------------------------------------
  ksea_df <- as.data.frame(ksea_matrix)
  ksea_df$Kinase <- rownames(ksea_df)
  # Move Kinase column to front
  ksea_df <- ksea_df %>% dplyr::select(Kinase, everything())
  
  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  write.table(ksea_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("------------------------------------------------")
  message("ssKSEA Complete!")
  message(paste0("Number of quantifiable kinases: ", nrow(ksea_df)))
  message(paste0("Results saved to: ", output_file))
  message("------------------------------------------------")
}


# ==============================================================================
# Execution
# ==============================================================================

# Define directories using relative paths
proc_dir <- "./data/processed"

# 1. PhosSight ssKSEA Analysis
# ----------------------------
input_ps <- file.path(proc_dir, "Table_X2_PhosSight_50PercentMissingValueCutoff.txt")
output_ps <- file.path(proc_dir, "Table_X7_PhosSight_KSEA_Scores_50pct.txt")

run_ssksea(input_ps, output_ps)

# 2. PhosphoRS ssKSEA Analysis
# ----------------------------
input_rs <- file.path(proc_dir, "Table_X2_PhosphoRS_50PercentMissingValueCutoff.txt")
output_rs <- file.path(proc_dir, "Table_X7_PhosphoRS_KSEA_Scores_50pct.txt")

run_ssksea(input_rs, output_rs)
