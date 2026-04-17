###Figure7i
# ==============================================================================
# Script Name: Kinase Survival Validation (Figure 7i)
# Description: Generates Kaplan-Meier survival curves for specific kinase activities
#              (e.g., MARK2) inferred by PhosSight (ssKSEA).
#              Allows for custom High/Low group splitting (e.g., Optimal Cutpoint counts).
#
# Usage:       Run from project root. Requires ssKSEA Scores (Table X7) and Survival Data.
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(readxl)
library(tibble)

# ==============================================================================
# Configuration
# ==============================================================================
# Target Kinase and Split Parameters (From Optimal Cutpoint Analysis in Step 5)
TARGET_KINASE <- "MARK2"
N_HIGH_TARGET <- 19  # Adjust based on your optimal cutpoint results
N_LOW_TARGET  <- 75  # Adjust based on your optimal cutpoint results

# Input Files (Relative Paths)
base_raw  <- "./data/raw"
base_proc <- "./data/processed"

surv_file <- file.path(base_raw, "UCEC_survival.txt")
meta_file <- file.path(base_raw, "mmc1.xlsx")
# Use PhosSight KSEA scores
ksea_file <- file.path(base_proc, "Table_X7_PhosSight_KSEA_Scores_50pct.txt")

# Output Directory
output_dir <- file.path(base_proc, "Figure7i_KM_Plots_Kinase")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_pdf <- file.path(output_dir, paste0("Figure7i_", TARGET_KINASE, "_OS_Survival.pdf"))

# ==============================================================================
# 1. Plotting Function (Strict Count Split)
# ==============================================================================
plot_ksea_survival_forced_os <- function(target_kinase, n_high_req, n_low_req, ksea_file, survival_file, meta_file, output_path) {
  
  message(paste0("\n----------------------------------------------------------"))
  message(paste0("Processing Kinase (OS): ", target_kinase))
  message(paste0("Target Split: High = ", n_high_req, ", Low = ", n_low_req))
  
  # 1. Clean Metadata
  if (!file.exists(meta_file)) stop("Metadata file missing.")
  meta_df <- read_excel(meta_file)
  
  sample_map <- meta_df %>%
    filter(Proteomics_Tumor_Normal == "Tumor") %>% 
    dplyr::select(case_id = Proteomics_Participant_ID, sample_ids = Proteomics_Parent_Sample_IDs) %>%
    separate_rows(sample_ids, sep = ",") %>%
    mutate(sample_ids = trimws(sample_ids)) %>% 
    distinct()
  
  # 2. Read KSEA Data
  if (!file.exists(ksea_file)) stop("KSEA file missing.")
  ksea_data <- fread(ksea_file, header = TRUE, sep = "\t")
  
  if(!target_kinase %in% ksea_data$Kinase) stop(paste("Kinase", target_kinase, "not found in KSEA data."))
  
  ksea_sub <- ksea_data %>% 
    filter(Kinase == target_kinase) %>%
    pivot_longer(cols = -Kinase, names_to = "Sample_Full", values_to = "Score")
  
  # 3. Match and Average (Sample -> Case)
  ksea_clean <- ksea_sub %>%
    inner_join(sample_map, by = c("Sample_Full" = "sample_ids")) %>%
    filter(!is.na(Score)) %>%
    group_by(Kinase, case_id) %>%
    summarise(Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  # 4. Merge Survival Data
  if (!file.exists(survival_file)) stop("Survival file missing.")
  surv_info <- fread(survival_file, header = TRUE, sep = "\t")
  surv_info$case_id <- as.character(surv_info$case_id)
  
  # Ensure numeric columns
  if("OS_days" %in% colnames(surv_info)) surv_info$OS_days <- as.numeric(surv_info$OS_days)
  if("OS_event" %in% colnames(surv_info)) surv_info$OS_event <- as.numeric(surv_info$OS_event)
  
  plot_data <- ksea_clean %>%
    inner_join(surv_info, by = "case_id") %>%
    filter(!is.na(OS_days) & !is.na(OS_event))
  
  # 5. Sample Size Validation
  n_real <- nrow(plot_data)
  n_target_total <- n_high_req + n_low_req
  
  message(paste0("Valid OS Samples available: ", n_real))
  
  # Warn if counts don't match (e.g., due to missing survival data for some patients)
  if (n_real != n_target_total) {
    warning(paste0("Mismatch: Real samples (", n_real, ") vs Target total (", n_target_total, ")."))
    warning("Adjusting split proportionally based on Score ranking.")
    
    # Proportional Adjustment
    prop_high <- n_high_req / n_target_total
    n_high_final <- round(n_real * prop_high)
  } else {
    n_high_final <- n_high_req
  }
  
  # 6. Apply Split
  plot_data <- plot_data %>% arrange(desc(Score)) # Sort descending
  plot_data$Group <- "Low"
  
  # Assign top N to High
  if (n_high_final > 0 && n_high_final < n_real) {
    plot_data$Group[1:n_high_final] <- "High"
  } else {
    stop("Invalid split size calculation.")
  }
  
  real_n_high <- sum(plot_data$Group == "High")
  real_n_low  <- sum(plot_data$Group == "Low")
  message(paste0("Final Split -> High: ", real_n_high, " / Low: ", real_n_low))
  
  plot_data$Group <- factor(plot_data$Group, levels = c("Low", "High"))
  
  # 7. Cox Regression Statistics
  cox_fit <- coxph(Surv(OS_days, OS_event) ~ Group, data = plot_data)
  cox_sum <- summary(cox_fit)
  
  hr_val   <- round(cox_sum$conf.int[1, "exp(coef)"], 2)
  hr_lower <- round(cox_sum$conf.int[1, "lower .95"], 2)
  hr_upper <- round(cox_sum$conf.int[1, "upper .95"], 2)
  
  p_val_raw <- cox_sum$coefficients[1, "Pr(>|z|)"]
  p_text <- ifelse(p_val_raw < 0.001, "p < 0.001", paste0("p = ", round(p_val_raw, 3)))
  
  # Label Format: HR = 1.52 (1.05-2.20)
  label_text <- paste0("PhosSight\nHR = ", hr_val, " (", hr_lower, "-", hr_upper, ")\n", p_text)
  
  # 8. Generate Plot
  fit <- survfit(Surv(OS_days, OS_event) ~ Group, data = plot_data)
  
  p_surv <- ggsurvplot(
    fit,
    data = plot_data,
    palette = c("#377EB8", "#E41A1C"), # Blue (Low), Red (High)
    size = 1.2,
    xlab = "Time (Days)", 
    ylab = "Survival Probability (OS)", 
    legend.title = "",
    legend.labs = c(paste(target_kinase, "low"), paste(target_kinase, "high")),
    pval = FALSE, # Custom annotation used instead
    risk.table = FALSE,
    ggtheme = theme_classic() + 
      theme(
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = c(0.75, 0.85),
        legend.background = element_blank(),
        panel.grid = element_blank()
      )
  )
  
  final_plot <- p_surv$plot +
    annotate("text", 
             x = max(plot_data$OS_days, na.rm=T) * 0.05, 
             y = 0.15, 
             label = label_text, 
             hjust = 0, 
             size = 4.5,
             fontface = "bold")
  
  # Save
  ggsave(output_path, final_plot, width = 5, height = 4.5)
  message(paste0("Plot saved to: ", output_path))
}

# ==============================================================================
# 2. Execution
# ==============================================================================
plot_ksea_survival_forced_os(
  target_kinase = TARGET_KINASE, 
  n_high_req = N_HIGH_TARGET, 
  n_low_req = N_LOW_TARGET, 
  ksea_file = ksea_file, 
  survival_file = surv_file, 
  meta_file = meta_file, 
  output_path = output_pdf
)

