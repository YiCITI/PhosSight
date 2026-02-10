#Step5 ssKSEA推断激酶的生存分析
# ==============================================================================
# Script Name: Kinase Survival Analysis (Optimal Cutpoint)
# Description: This script performs univariate Cox proportional hazards regression 
#              to identify prognosis-associated kinases inferred by ssKSEA.
#              
#              Key Feature:
#              Instead of a median split, this script uses the 'survminer' package 
#              to determine the "Optimal Cutpoint" (maximally selected rank statistics)
#              for stratifying patients into High/Low kinase activity groups.
#
# Usage:       Run this script from the project root directory.
#              Requires: ssKSEA Scores (txt) from Step 4.
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(readxl)
library(survminer) # Required for surv_cutpoint

# ==============================================================================
# Define Survival Analysis Function (Optimal Cutpoint)
# ==============================================================================
run_ksea_survival_opt <- function(ksea_file, survival_file, meta_file, output_file) {
  
  message(paste0("\n----------------------------------------------------------"))
  message(paste0("Processing file: ", basename(ksea_file)))
  
  # 1. Process Metadata (Build Sample -> Case Mapping)
  # ----------------------------------------------------------------------------
  if (!file.exists(meta_file)) stop(paste("Error: Metadata file not found at", meta_file))
  
  message("Reading and cleaning metadata...")
  meta_df <- read_excel(meta_file)
  
  # Filter Tumor samples and map to Case IDs
  sample_map <- meta_df %>%
    filter(Proteomics_Tumor_Normal == "Tumor") %>% 
    dplyr::select(
      case_id = Proteomics_Participant_ID, 
      sample_ids = Proteomics_Parent_Sample_IDs
    ) %>%
    separate_rows(sample_ids, sep = ",") %>%
    mutate(sample_ids = trimws(sample_ids)) %>% 
    distinct()
  
  valid_tumor_samples <- sample_map$sample_ids
  message(paste0("Valid Tumor samples in metadata: ", length(valid_tumor_samples)))
  
  # 2. Read KSEA Matrix and Reshape
  # ----------------------------------------------------------------------------
  if (!file.exists(ksea_file)) stop(paste("Error: KSEA file not found at", ksea_file))
  if (!file.exists(survival_file)) stop(paste("Error: Survival file not found at", survival_file))
  
  ksea_data <- fread(ksea_file, header = TRUE, sep = "\t")
  surv_info <- fread(survival_file, header = TRUE, sep = "\t")
  
  # Convert wide matrix to long format
  ksea_long <- ksea_data %>%
    pivot_longer(cols = -Kinase, names_to = "Sample_Full", values_to = "Score")
  
  # 3. Sample Matching and Averaging
  # ----------------------------------------------------------------------------
  # Match samples to metadata
  ksea_clean <- ksea_long %>%
    inner_join(sample_map, by = c("Sample_Full" = "sample_ids")) %>%
    filter(!is.na(Score)) 
  
  # Average scores for replicates from the same patient (Case)
  ksea_final <- ksea_clean %>%
    group_by(Kinase, case_id) %>%
    summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  # 4. Merge Survival Data
  # ----------------------------------------------------------------------------
  surv_info$case_id <- as.character(surv_info$case_id)
  
  analysis_df <- ksea_final %>%
    inner_join(surv_info, by = "case_id") %>%
    filter(!is.na(OS_days) & !is.na(OS_event)) 
  
  total_patients <- length(unique(analysis_df$case_id))
  message(paste0("Number of valid patients (cases) for analysis: ", total_patients))
  
  if(total_patients < 20) {
    warning("Warning: Sample size is very small (<20). Results may be unreliable.")
  }
  
  # 5. Loop Cox Regression (Optimal Cutpoint)
  # ----------------------------------------------------------------------------
  all_kinases <- unique(analysis_df$Kinase)
  results_list <- list()
  
  message(paste0("Number of kinases to analyze: ", length(all_kinases)))
  pb <- txtProgressBar(min = 0, max = length(all_kinases), style = 3)
  
  for (i in seq_along(all_kinases)) {
    kinase <- all_kinases[i]
    
    sub_data <- analysis_df %>% filter(Kinase == kinase)
    sub_data <- na.omit(sub_data)
    
    # Basic filter for sample size
    if(nrow(sub_data) < 20) next
    
    # --- [KEY STEP] Determine Optimal Cutpoint ---
    tryCatch({
      # surv_cutpoint finds the cutoff that yields the most significant split
      # minprop = 0.20 ensures each group has at least 20% of samples (prevents extreme splits)
      res.cut <- surv_cutpoint(sub_data, time = "OS_days", event = "OS_event",
                               variables = "Score", minprop = 0.20)
      
      # Extract the specific cutoff value
      cut_point <- res.cut$cutpoint$cutpoint[1]
      
      # Stratify patients based on this optimal cutoff
      # > cut_point = High, <= cut_point = Low
      sub_data <- sub_data %>%
        mutate(Group = ifelse(Score > cut_point, "High", "Low"))
      
      # Count group sizes
      n_high <- sum(sub_data$Group == "High")
      n_low  <- sum(sub_data$Group == "Low")
      
      # Double check group sizes (redundant but safe)
      if(n_high < 5 || n_low < 5) next
      
      # --- Run Cox Regression ---
      sub_data$Group <- factor(sub_data$Group, levels = c("Low", "High"))
      
      fit <- coxph(Surv(OS_days, OS_event) ~ Group, data = sub_data)
      fit_summary <- summary(fit)
      
      results_list[[i]] <- data.frame(
        Kinase = kinase,
        cutpoint_value = cut_point, 
        n_high = n_high,
        n_low = n_low,
        pvalue = fit_summary$coefficients[1, "Pr(>|z|)"],
        hr = fit_summary$conf.int[1, "exp(coef)"],
        lower_95 = fit_summary$conf.int[1, "lower .95"],
        upper_95 = fit_summary$conf.int[1, "upper .95"]
      )
      
    }, error = function(e) { return(NULL) })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 6. Compile and Save Results
  # ----------------------------------------------------------------------------
  final_results <- do.call(rbind, results_list)
  
  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if(!is.null(final_results) && nrow(final_results) > 0) {
    
    # BH Adjustment for P-values
    final_results$adjust_pvalue <- p.adjust(final_results$pvalue, method = "BH")
    
    # Reorder columns
    final_results <- final_results %>% 
      dplyr::select(Kinase, cutpoint_value, n_high, n_low, pvalue, adjust_pvalue, hr, lower_95, upper_95) %>%
      arrange(pvalue)
    
    write.table(final_results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message(paste0("\nAnalysis complete! Results saved to: ", output_file))
    
    message("Top 5 Significant Kinases (Optimal Cutpoint):")
    print(head(final_results, 5))
    
  } else {
    message("\nWarning: No valid results generated.")
  }
}

# ==============================================================================
# Execution
# ==============================================================================

# Define directories using relative paths
raw_dir <- "./data/raw"
proc_dir <- "./data/processed"

# Input Files
surv_file <- file.path(raw_dir, "UCEC_survival.txt")
meta_file <- file.path(raw_dir, "mmc1.xlsx")

# 1. PhosSight Dataset
# --------------------
# Input comes from Step 4 (KSEA Scores)
ps_ksea_file <- file.path(proc_dir, "Table_X7_PhosSight_KSEA_Scores_50pct.txt")
ps_out_file  <- file.path(proc_dir, "Table_X8_PhosSight_KSEA_Survival_OS_OptCut.txt")

run_ksea_survival_opt(ps_ksea_file, surv_file, meta_file, ps_out_file)

# 2. PhosphoRS Dataset
# --------------------
rs_ksea_file <- file.path(proc_dir, "Table_X7_PhosphoRS_KSEA_Scores_50pct.txt")
rs_out_file  <- file.path(proc_dir, "Table_X8_PhosphoRS_KSEA_Survival_OS_OptCut.txt")

run_ksea_survival_opt(rs_ksea_file, surv_file, meta_file, rs_out_file)
