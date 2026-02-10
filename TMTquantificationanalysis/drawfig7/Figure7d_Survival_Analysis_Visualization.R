#Figure7d
# ==============================================================================
# Script Name: Survival Analysis Visualization (Figure 7d)
# Description: Generates Kaplan-Meier survival curves for specific prognostic 
#              phosphosites (STMN1_46S, PARP1_368T).
#
# Usage:       Run from project root. Requires processed matrix (Table X2) and 
#              survival data (UCEC_survival.txt).
# ==============================================================================

library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(readxl)
library(ggplot2)

# ==============================================================================
# Configuration
# ==============================================================================
# Target Sites (Filtered as requested)
TARGET_SITES <- c("STMN1_46S", "PARP1_368T")

# Input Files
base_raw  <- "./data/raw"
base_proc <- "./data/processed"

surv_file  <- file.path(base_raw, "UCEC_survival.txt")
meta_file  <- file.path(base_raw, "mmc1.xlsx")
input_file <- file.path(base_proc, "Table_X2_PhosSight_50PercentMissingValueCutoff.txt")

# Output Directory
output_dir <- file.path(base_proc, "Figure7d_KM_Plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 1. Load and Preprocess Data
# ==============================================================================
message(">>> Loading and Preprocessing Data...")

# A. Process Metadata (Sample -> Case Mapping)
if (!file.exists(meta_file)) stop("Metadata file missing.")

meta_df <- read_excel(meta_file)
sample_map <- meta_df %>%
  filter(Proteomics_Tumor_Normal == "Tumor") %>%
  dplyr::select(case_id = Proteomics_Participant_ID, sample_ids = Proteomics_Parent_Sample_IDs) %>%
  separate_rows(sample_ids, sep = ",") %>%
  mutate(sample_ids = trimws(sample_ids)) %>%
  distinct()

# B. Read Expression Matrix
if (!file.exists(input_file)) stop("Expression Matrix file missing.")
data_mat <- fread(input_file, header = TRUE, sep = "\t")

# Check PhosphoSite column
if (!"PhosphoSite" %in% colnames(data_mat)) stop("Column 'PhosphoSite' missing in matrix.")

# C. Filter for Target Sites
message(paste0(">>> Filtering for targets: ", paste(TARGET_SITES, collapse = ", ")))
target_data <- data_mat %>% filter(PhosphoSite %in% TARGET_SITES)

if(nrow(target_data) == 0) stop("No target sites found in the matrix.")

# D. Reshape and Merge (Sample -> Case Averaging)
# 1. Identify valid sample columns
matrix_cols <- colnames(data_mat)
valid_samples <- intersect(matrix_cols, sample_map$sample_ids)

if(length(valid_samples) == 0) stop("No matching sample IDs found between matrix and metadata.")

# 2. Pivot and Average Replicates
clean_data <- target_data %>%
  dplyr::select(PhosphoSite, all_of(valid_samples)) %>%
  pivot_longer(cols = -PhosphoSite, names_to = "Sample_ID", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  inner_join(sample_map, by = c("Sample_ID" = "sample_ids")) %>%
  # [CRITICAL]: Average multiple samples per patient (Case)
  group_by(PhosphoSite, case_id) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# E. Merge Survival Info
if (!file.exists(surv_file)) stop("Survival file missing.")
surv_info <- fread(surv_file, header = TRUE, sep = "\t")
surv_info$case_id <- as.character(surv_info$case_id)
surv_info$OS_days <- as.numeric(surv_info$OS_days)
surv_info$OS_event <- as.numeric(surv_info$OS_event)

plot_ready_df <- clean_data %>%
  inner_join(surv_info, by = "case_id") %>%
  filter(!is.na(OS_days) & !is.na(OS_event))

message(paste0(">>> Data ready. Total patients (cases): ", length(unique(plot_ready_df$case_id))))

# ==============================================================================
# 2. Generate Survival Plots
# ==============================================================================

for (site in TARGET_SITES) {
  
  # Extract data for current site
  sub_data <- plot_ready_df %>% filter(PhosphoSite == site)
  
  # Basic Validation
  if (nrow(sub_data) < 20) {
    message(paste0("Skipping ", site, ": Insufficient data (<20 samples)."))
    next
  }
  
  # Stratify by Median
  med_val <- median(sub_data$Value)
  sub_data$Group <- ifelse(sub_data$Value > med_val, "High", "Low")
  sub_data$Group <- factor(sub_data$Group, levels = c("Low", "High")) 
  
  # Group Balance Check
  n_high <- sum(sub_data$Group == "High")
  n_low <- sum(sub_data$Group == "Low")
  
  if (n_high < 5 || n_low < 5) {
    message(paste0("Skipping ", site, ": Unbalanced groups (High=", n_high, ", Low=", n_low, ")."))
    next
  }
  
  message(paste0("Plotting: ", site, " (n=", nrow(sub_data), ")"))
  
  # ----------------------------------------------------------------------------
  # (1) Cox Regression Statistics
  # ----------------------------------------------------------------------------
  cox_fit <- coxph(Surv(OS_days, OS_event) ~ Group, data = sub_data)
  res_sum <- summary(cox_fit)
  
  p_val <- res_sum$coefficients[1, "Pr(>|z|)"]
  hr_val <- res_sum$conf.int[1, "exp(coef)"]
  lower_ci <- res_sum$conf.int[1, "lower .95"]
  upper_ci <- res_sum$conf.int[1, "upper .95"]
  
  # Format Label Text
  p_text <- ifelse(p_val < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", p_val)))
  
  label_text <- paste0(
    "HR = ", sprintf("%.2f", hr_val), "\n",
    "95% CI: ", sprintf("%.2f", lower_ci), "-", sprintf("%.2f", upper_ci), "\n",
    p_text
  )
  
  # ----------------------------------------------------------------------------
  # (2) Kaplan-Meier Plot
  # ----------------------------------------------------------------------------
  km_fit <- survfit(Surv(OS_days, OS_event) ~ Group, data = sub_data)
  
  p_surv <- ggsurvplot(
    km_fit,
    data = sub_data,
    palette = c("#377EB8", "#E41A1C"),  # Blue (Low), Red (High)
    size = 1.2,                         # Line thickness
    xlab = "Time (Days)", 
    ylab = "Survival Probability (OS)", 
    legend.title = "",
    legend.labs = c(paste(site, "low"), paste(site, "high")),
    pval = FALSE,                       # Custom annotation used instead
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
  
  # Add Statistical Annotation
  final_plot <- p_surv$plot +
    annotate("text", 
             x = max(sub_data$OS_days, na.rm=T) * 0.05, 
             y = 0.15, 
             label = label_text, 
             hjust = 0, 
             size = 4.5,
             fontface = "bold")
  
  # ----------------------------------------------------------------------------
  # (3) Save
  # ----------------------------------------------------------------------------
  safe_site_name <- gsub("[^A-Za-z0-9_]", "", site) 
  output_pdf <- file.path(output_dir, paste0("Figure7e_", safe_site_name, ".pdf"))
  
  ggsave(output_pdf, final_plot, width = 5, height = 4.5)
  message(paste0("  -> Saved: ", output_pdf))
}

message("\n>>> Figure 7e Generation Complete!")






