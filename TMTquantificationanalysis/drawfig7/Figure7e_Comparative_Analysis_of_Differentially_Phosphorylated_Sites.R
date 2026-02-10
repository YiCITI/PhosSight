#Figure7e
# ==============================================================================
# Script Name: Comparative Analysis of Differentially Phosphorylated Sites
# Description: Generates Figure 7e - A stacked bar chart classifying sites 
#              identified as differentially phosphorylated by PhosSight based on 
#              their status in the PhosphoRS workflow.
#
# Usage:       Run from project root. Requires processed matrices (Table X1) and 
#              differential analysis results (Table X3).
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)

# ==============================================================================
# Configuration
# ==============================================================================
# Define Threshold for "Significance"
ADJ_P_THRESHOLD <- 0.05 
# Define Threshold for "Quantifiable" (Number of samples)
QUANTIFIABLE_CUTOFF <- 92 # 50% of 183 samples

# Input Files (Relative Paths)
base_raw  <- "./data/raw"
base_proc <- "./data/processed"

file_mat_rs  <- file.path(base_raw, "Table_X1_PhosphoRS_OriginalMatrix.txt")
file_diff_ps <- file.path(base_proc, "Table_X3_PhosSight_DifferentialAnalysis.txt")
file_diff_rs <- file.path(base_proc, "Table_X3_PhosphoRS_DifferentialAnalysis.txt")

# Output File
output_file <- file.path(base_proc, "Figure7e_Differential_Comparison.pdf")

# ==============================================================================
# 1. Load Data
# ==============================================================================
message(">>> Loading Data...")

if(!file.exists(file_mat_rs)) stop("Missing PhosphoRS Matrix file.")
if(!file.exists(file_diff_ps)) stop("Missing PhosSight Diff Result file.")

# Read Matrix (For identification/quantification status)
mat_rs <- fread(file_mat_rs, header = TRUE, sep = "\t")

# Read Differential Results
diff_ps <- fread(file_diff_ps, header = TRUE, sep = "\t")
diff_rs <- fread(file_diff_rs, header = TRUE, sep = "\t")

# Create Unique IDs (Gene_Site)
mat_rs$ID  <- paste(mat_rs$GeneSymbol, mat_rs$PhosphoSite, sep = "_")
diff_ps$ID <- paste(diff_ps$GeneSymbol, diff_ps$PhosphoSite, sep = "_")
diff_rs$ID <- paste(diff_rs$GeneSymbol, diff_rs$PhosphoSite, sep = "_")

# ==============================================================================
# 2. Preprocess PhosphoRS Status
# ==============================================================================
message(">>> Preprocessing PhosphoRS Status...")

# 2.1 Identify Sample Columns in Matrix
meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
sample_cols_rs <- setdiff(colnames(mat_rs), meta_cols)

# 2.2 Calculate Observation Count per Site
# Count non-NA values for each row in PhosphoRS matrix
rs_counts <- rowSums(!is.na(as.matrix(mat_rs[, ..sample_cols_rs])))
names(rs_counts) <- mat_rs$ID

# 2.3 Define Sets
# Set A: Identified in PhosphoRS (Any detection)
set_rs_identified <- mat_rs$ID

# Set B: Quantifiable in PhosphoRS (>= 50% samples)
set_rs_quantifiable <- names(rs_counts)[rs_counts >= QUANTIFIABLE_CUTOFF]

# Set C: Significant in PhosphoRS (Diff Analysis, adj.p < 0.05)
set_rs_significant <- diff_rs %>% 
  filter(adjust_pvalue < ADJ_P_THRESHOLD) %>% 
  pull(ID)

# ==============================================================================
# 3. Classify PhosSight Significant Sites
# ==============================================================================
message(">>> Classifying PhosSight Hits...")

# 3.1 Filter PhosSight Significant Sites
sig_ps <- diff_ps %>% filter(adjust_pvalue < ADJ_P_THRESHOLD)

# 3.2 Define Classification Function
get_status <- function(ids) {
  sapply(ids, function(x) {
    if (!x %in% set_rs_identified) {
      return("Not identified")
    } else if (!x %in% set_rs_quantifiable) {
      return("Not quantifiable")
    } else if (x %in% set_rs_significant) {
      return("Significant")
    } else {
      return("Not significant")
    }
  })
}

# 3.3 Apply Classification
sig_ps$RS_Status <- get_status(sig_ps$ID)

# 3.4 Define Direction (Up/Down)
# Log2FC > 0 -> Higher in Tumor
# Log2FC < 0 -> Lower in Tumor
sig_ps$Direction <- ifelse(sig_ps$log2FC > 0, "Higher in Tumor", "Lower in Tumor")

# ==============================================================================
# 4. Prepare Plot Data
# ==============================================================================
plot_data <- sig_ps %>%
  group_by(Direction, RS_Status) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate Totals for Labels
plot_data <- plot_data %>%
  group_by(Direction) %>%
  mutate(Total_Group = sum(Count)) %>%
  ungroup() %>%
  mutate(
    # Create X-axis label with total count (e.g., "Higher in Tumor\n(n=4617)")
    Label = paste0(Direction, "\n(n=", Total_Group, ")"),
    # Calculate Percentage for text labels
    Prop = Count / Total_Group,
    Percentage = paste0(round(Prop * 100), "%")
  )

# Set Factor Levels for Stacking Order
# Order: Not identified (Top) -> Not quantifiable -> Not significant -> Significant (Bottom)
plot_data$RS_Status <- factor(plot_data$RS_Status, 
                              levels = c("Not identified", "Not quantifiable", "Not significant", "Significant"))

# Define Colors (Matching Figure 7d style)
# Significant: Dark Blue
# Not significant: Medium Blue
# Not quantifiable: Light Grey
# Not identified: White/Very Light Grey
custom_colors <- c(
  "Significant"      = "#2171B5", 
  "Not significant"  = "#9ECAE1", 
  "Not quantifiable" = "#D9D9D9", 
  "Not identified"   = "#F7F7F7"  
)

# Calculate label positions for text inside bars
plot_data <- plot_data %>%
  arrange(Direction, desc(RS_Status)) %>%
  group_by(Direction) %>%
  mutate(pos = cumsum(Count) - (0.5 * Count))

# ==============================================================================
# 5. Plotting (Stacked Bar Chart)
# ==============================================================================
message(">>> Generating Plot...")

p <- ggplot(plot_data, aes(x = Label, y = Count, fill = RS_Status)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, size = 0.3) +
  
  # Add Text Labels inside bars (Count + Percentage)
  # Only show label if the segment is large enough (> 3%) to avoid clutter
  geom_text(aes(y = pos, 
                label = ifelse(Prop > 0.03, paste0(Count, "\n(", Percentage, ")"), "")),
            size = 3.5, color = "black") +
  
  scale_fill_manual(values = custom_colors) +
  
  labs(
    title = "",
    x = "",
    y = "Number of Phosphosites",
    fill = "Status in PhosphoRS"
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  
  # Add small buffer to top of Y axis
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Save
# Ensure output directory exists
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

ggsave(output_file, p, width = 6, height = 6)
message(paste("Plot saved to:", output_file))





