#Figure7c
# ==============================================================================
# Script Name: Comparative Analysis of Prognostic Phosphosites
# Description: Generates Figure 7c - A stacked bar chart classifying sites 
#              identified as prognostic by PhosSight based on their detection 
#              and significance status in the PhosphoRS workflow.
#
# Usage:       Run from project root. Requires processed matrices (Table X1) and 
#              survival analysis results (Table X4).
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)

# ==============================================================================
# Configuration
# ==============================================================================
# Thresholds
P_THRESHOLD  <- 0.05
HR_THRESHOLD <- 2.0  # Defines "Poor" (HR > 2) and "Good" (HR < 0.5)

# Input Files (Relative Paths)
base_raw  <- "./data/raw"
base_proc <- "./data/processed"

file_mat_rs  <- file.path(base_raw, "Table_X1_PhosphoRS_OriginalMatrix.txt")
file_surv_ps <- file.path(base_proc, "Table_X4_PhosSight_SurvivalAnalysis_OS.txt")
file_surv_rs <- file.path(base_proc, "Table_X4_PhosphoRS_SurvivalAnalysis_OS.txt")

# Output File
output_file <- file.path(base_proc, "Figure7c_Prognostic_Comparison.pdf")

# ==============================================================================
# 1. Load Data
# ==============================================================================
message(">>> Loading Data...")

if(!file.exists(file_mat_rs)) stop("Missing PhosphoRS Matrix file.")
if(!file.exists(file_surv_ps)) stop("Missing PhosSight Survival file.")

mat_rs  <- fread(file_mat_rs, header = TRUE, sep = "\t")
surv_ps <- fread(file_surv_ps, header = TRUE, sep = "\t")
surv_rs <- fread(file_surv_rs, header = TRUE, sep = "\t")

# Create Unique IDs
# Note: Assuming GeneSymbol and PhosphoSite columns exist
mat_rs$ID  <- paste(mat_rs$GeneSymbol, mat_rs$PhosphoSite, sep = "_")
surv_ps$ID <- paste(surv_ps$GeneSymbol, surv_ps$PhosphoSite, sep = "_")
surv_rs$ID <- paste(surv_rs$GeneSymbol, surv_rs$PhosphoSite, sep = "_")

# ==============================================================================
# 2. Build Lookup Lists
# ==============================================================================
message(">>> Building Lookup Lists...")

# List A: All IDs identified in PhosphoRS original matrix
rs_all_matrix_ids <- unique(mat_rs$ID)

# List B: IDs that were actually tested in PhosphoRS survival analysis
# (This implies they passed the missing value cutoff and were "Quantifiable")
rs_tested_ids <- unique(surv_rs$ID)

# List C: Significant prognostic sites in PhosphoRS
rs_significant_ids <- surv_rs %>% 
  filter(pvalue < P_THRESHOLD) %>% 
  pull(ID) %>% 
  unique()

# ==============================================================================
# 3. Process PhosSight Prognostic Sites
# ==============================================================================
message(">>> Processing PhosSight Prognostic Sites...")

# 3.1 Filter Significant PhosSight Sites
sig_ps_raw <- surv_ps %>% filter(pvalue < P_THRESHOLD)

# 3.2 Remove Duplicates
# (Keeps the first instance if a site appears multiple times)
sig_ps <- sig_ps_raw %>% distinct(ID, .keep_all = TRUE)
message(paste("  PhosSight Significant Sites (Unique):", nrow(sig_ps)))

# 3.3 Define Groups (Poor vs Good Prognosis)
# Poor: HR > Threshold
# Good: HR < (1 / Threshold)
sig_ps_poor <- sig_ps %>% filter(hr > HR_THRESHOLD)
sig_ps_good <- sig_ps %>% filter(hr < (1/HR_THRESHOLD))

# ==============================================================================
# 4. Classification Logic
# ==============================================================================
classify_status <- function(target_ids) {
  case_when(
    # 1. Significant in PhosphoRS
    target_ids %in% rs_significant_ids ~ "Significant",
    
    # 2. Not significant (but was tested/quantifiable)
    target_ids %in% rs_tested_ids ~ "Not significant",
    
    # 3. Not quantifiable (in matrix but filtered out)
    target_ids %in% rs_all_matrix_ids ~ "Not quantifiable",
    
    # 4. Not identified (not in matrix)
    TRUE ~ "Not identified"
  )
}

# Apply Classification
if(nrow(sig_ps_poor) > 0) {
  df_poor <- data.frame(ID = sig_ps_poor$ID)
  df_poor$Status <- classify_status(df_poor$ID)
} else {
  df_poor <- data.frame(ID = character(), Status = character())
}

if(nrow(sig_ps_good) > 0) {
  df_good <- data.frame(ID = sig_ps_good$ID)
  df_good$Status <- classify_status(df_good$ID)
} else {
  df_good <- data.frame(ID = character(), Status = character())
}

# ==============================================================================
# 5. Prepare Plot Data
# ==============================================================================
prepare_stats <- function(df, group_label) {
  if(nrow(df) == 0) return(NULL)
  
  df %>%
    group_by(Status) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    mutate(
      Total = sum(Count),
      Prop = Count / Total,
      Percentage = paste0(round(Prop * 100), "%"),
      Group = group_label
    )
}

label_poor <- paste0("Poor Survival\n(n=", nrow(df_poor), ")")
label_good <- paste0("Good Survival\n(n=", nrow(df_good), ")")

plot_df <- rbind(
  prepare_stats(df_poor, label_poor),
  prepare_stats(df_good, label_good)
)

# Set Factor Levels
plot_df$Status <- factor(plot_df$Status, 
                         levels = c("Not identified", "Not quantifiable", "Not significant", "Significant"))

# ==============================================================================
# 6. Plotting
# ==============================================================================
if (!is.null(plot_df)) {
  
  # Colors
  cols <- c(
    "Significant"      = "#2171B5", 
    "Not significant"  = "#6BAED6", 
    "Not quantifiable" = "#DEEBF7", 
    "Not identified"   = "#F7FBFF" 
  )
  
  # Calculate Label Positions
  plot_df <- plot_df %>%
    group_by(Group) %>%
    arrange(desc(Status)) %>%
    mutate(pos = cumsum(Count) - (0.5 * Count))
  
  p <- ggplot(plot_df, aes(x = Group, y = Count, fill = Status)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.3) +
    
    # Text Labels (Count + Percentage)
    # Filter out small segments (<3%) to prevent clutter
    geom_text(aes(y = pos, 
                  label = ifelse(Prop > 0.03, paste0(Count, "\n(", Percentage, ")"), "")),
              size = 3.5, color = "black") +
    
    scale_fill_manual(values = cols) +
    
    labs(
      title = "Status in PhosphoRS",
      x = "",
      y = "Number of Prognostic Phosphosites",
      fill = "Status in PhosphoRS"
    ) +
    
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(size = 11, face = "bold", color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      legend.position = "right",
      legend.title = element_blank()
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  # Save
  # Ensure output directory exists
  if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
  
  ggsave(output_file, p, width = 6, height = 6)
  message(paste("Plot saved to:", output_file))
  
} else {
  message("Warning: No prognostic sites found matching criteria. Plot not generated.")
}













