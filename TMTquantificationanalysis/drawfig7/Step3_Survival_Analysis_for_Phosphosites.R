# Step3 位点表达量中位数分组的生存分析
# ==============================================================================
# Script Name: Survival Analysis for Phosphosites
# Description: This script performs univariate Cox proportional hazards regression 
#              to identify prognosis-associated phosphosites.
#              1. Cleans metadata and maps samples to clinical cases.
#              2. Averages replicates for each patient.
#              3. Stratifies patients into High/Low groups based on median expression.
#              4. Calculates Hazard Ratios (HR) and P-values (Cox model).
#
# Usage:       Run this script from the project root directory.
#              Requires: Metadata (xlsx), Survival Data (txt), and Preprocessed Data (txt).
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(readxl)

# ==============================================================================
# Define Survival Analysis Function
# ==============================================================================
run_site_survival <- function(input_file, survival_file, meta_file, output_file) {
  
  message(paste0("\n----------------------------------------------------------"))
  message(paste0("Processing file: ", basename(input_file)))
  
  # 1. Process Metadata (Build Sample -> Case Mapping)
  # ----------------------------------------------------------------------------
  if (!file.exists(meta_file)) stop(paste("Error: Metadata file not found at", meta_file))
  
  message("Reading and cleaning metadata...")
  meta_df <- read_excel(meta_file)
  
  # A. Filter for Tumor samples only
  # B. Extract Participant ID (Case) and Sample ID
  # C. Handle comma-separated sample lists (e.g., "ID1, ID2")
  sample_map <- meta_df %>%
    filter(Proteomics_Tumor_Normal == "Tumor") %>% 
    dplyr::select(
      case_id = Proteomics_Participant_ID, 
      sample_ids = Proteomics_Parent_Sample_IDs
    ) %>%
    separate_rows(sample_ids, sep = ",") %>%
    mutate(sample_ids = trimws(sample_ids)) %>% # Remove whitespace
    distinct()
  
  # Only samples in this list are valid tumor samples for analysis
  valid_tumor_samples <- sample_map$sample_ids
  
  message(paste0("Valid Tumor samples in metadata: ", length(valid_tumor_samples)))
  
  # 2. Read Expression Matrix
  # ----------------------------------------------------------------------------
  if (!file.exists(input_file)) stop(paste("Error: Input file not found at", input_file))
  
  data_mat <- fread(input_file, header = TRUE, sep = "\t")
  
  # 3. Construct Unique ID
  # ----------------------------------------------------------------------------
  if (!"ID" %in% colnames(data_mat)) {
    data_mat$ID <- paste(data_mat$GeneSymbol, data_mat$PhosphoSite, sep = "_")
  }
  
  # 4. Reshape and Filter Data
  # ----------------------------------------------------------------------------
  # Define standard metadata columns to exclude during pivoting
  meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
  existing_meta <- intersect(names(data_mat), meta_cols)
  
  # Convert to Long Format
  data_long <- data_mat %>%
    dplyr::select(ID, everything(), -any_of(setdiff(existing_meta, "ID"))) %>%
    pivot_longer(cols = -ID, names_to = "Sample_Full", values_to = "Value")
  
  # [CRITICAL]: Filter using inner_join with sample_map
  # This simultaneously filters out Normal samples and maps Sample IDs to Case IDs
  data_clean <- data_long %>%
    inner_join(sample_map, by = c("Sample_Full" = "sample_ids")) %>%
    filter(!is.na(Value))
  
  # Handle Replicates: Average multiple samples from the same patient
  data_final <- data_clean %>%
    group_by(ID, case_id) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  message(paste0("Number of unique cases (patients) after matching: ", length(unique(data_final$case_id))))
  
  # 5. Merge Survival Data
  # ----------------------------------------------------------------------------
  if (!file.exists(survival_file)) stop(paste("Error: Survival file not found at", survival_file))
  
  surv_info <- fread(survival_file, header = TRUE, sep = "\t")
  # Ensure case_id is character type for matching
  surv_info$case_id <- as.character(surv_info$case_id)
  
  analysis_df <- data_final %>%
    inner_join(surv_info, by = "case_id") %>%
    filter(!is.na(OS_days) & !is.na(OS_event))
  
  # 6. Loop Calculation (Cox Regression)
  # ----------------------------------------------------------------------------
  unique_ids <- unique(analysis_df$ID)
  results_list <- list()
  
  message(paste0("Number of phosphosites to analyze: ", length(unique_ids)))
  pb <- txtProgressBar(min = 0, max = length(unique_ids), style = 3)
  
  for (i in seq_along(unique_ids)) {
    cur_id <- unique_ids[i]
    sub_data <- analysis_df %>% filter(ID == cur_id)
    
    # Remove NAs
    sub_data <- na.omit(sub_data)
    
    # Filter: Minimum sample size required (e.g., detected in at least 20 patients)
    if(nrow(sub_data) < 20) next
    
    # Stratify patients by Median
    median_val <- median(sub_data$Value)
    sub_data <- sub_data %>%
      mutate(Group = ifelse(Value > median_val, "High", "Low"))
    
    n_high <- sum(sub_data$Group == "High")
    n_low  <- sum(sub_data$Group == "Low")
    
    # Filter: Minimum group size
    if(n_high < 5 || n_low < 5) next
    
    tryCatch({
      # Factorize group (Low as reference)
      sub_data$Group <- factor(sub_data$Group, levels = c("Low", "High"))
      
      # Cox Proportional Hazards Model
      fit <- coxph(Surv(OS_days, OS_event) ~ Group, data = sub_data)
      fit_summary <- summary(fit)
      
      res <- data.frame(
        ID = cur_id,
        median_value = median_val,
        n_high = n_high,
        n_low = n_low,
        pvalue = fit_summary$coefficients[1, "Pr(>|z|)"],
        hr = fit_summary$conf.int[1, "exp(coef)"],
        lower_95 = fit_summary$conf.int[1, "lower .95"],
        upper_95 = fit_summary$conf.int[1, "upper .95"]
      )
      results_list[[i]] <- res
    }, error = function(e) return(NULL))
    
    if(i %% 100 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 7. Compile and Save Results
  # ----------------------------------------------------------------------------
  final_res <- do.call(rbind, results_list)
  
  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if(!is.null(final_res) && nrow(final_res) > 0) {
    # Map back to GeneSymbol and PhosphoSite for readability
    meta_map <- data_mat %>% dplyr::select(ID, GeneSymbol, PhosphoSite) %>% distinct()
    final_res <- final_res %>% left_join(meta_map, by = "ID")
    
    # BH Adjustment for P-values
    final_res$adjust_pvalue <- p.adjust(final_res$pvalue, method = "BH")
    
    # Reorder columns
    final_res <- final_res %>%
      dplyr::select(GeneSymbol, PhosphoSite, median_value, n_high, n_low, 
                    pvalue, adjust_pvalue, hr, lower_95, upper_95) %>%
      arrange(pvalue)
    
    write.table(final_res, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste0("\nAnalysis complete! Saved to: ", output_file))
    message(paste0("Significant sites (P < 0.05): ", sum(final_res$pvalue < 0.05)))
  } else {
    message("\nWarning: No valid results generated.")
  }
}

# ==============================================================================
# Execution
# ==============================================================================

# Define directories using relative paths
# NOTE: Ensure you are in the project root directory
raw_dir <- "./data/raw"
proc_dir <- "./data/processed"

# Input Files
# ------------------------------
surv_file <- file.path(raw_dir, "UCEC_survival.txt")
meta_file <- file.path(raw_dir, "mmc1.xlsx")

# 1. PhosSight Survival Analysis
# ------------------------------
input_ps <- file.path(proc_dir, "Table_X2_PhosSight_50PercentMissingValueCutoff.txt")
output_ps <- file.path(proc_dir, "Table_X4_PhosSight_SurvivalAnalysis_OS.txt")

run_site_survival(input_ps, surv_file, meta_file, output_ps)

# 2. PhosphoRS Survival Analysis
# ------------------------------
input_rs <- file.path(proc_dir, "Table_X2_PhosphoRS_50PercentMissingValueCutoff.txt")
output_rs <- file.path(proc_dir, "Table_X4_PhosphoRS_SurvivalAnalysis_OS.txt")

run_site_survival(input_rs, surv_file, meta_file, output_rs)
