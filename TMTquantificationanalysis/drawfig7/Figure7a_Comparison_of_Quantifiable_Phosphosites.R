#Figure7a 
# ==============================================================================
# Script Name: Comparison of Quantifiable Phosphosites (Figure 7a)
# Description: Generates a vertical stacked bar chart comparing the number of 
#              quantifiable sites (<50% missing values) between PhosphoRS and 
#              PhosSight.
#
# Usage:       Run from project root. Requires Original Matrices (Table X1).
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

# ==============================================================================
# Configuration
# ==============================================================================
# Cutoff for "Quantifiable" (50% of 183 samples = 91.5 -> 92)
CUTOFF <- 92 

# Input Files (Relative Paths)
base_raw <- "./data/raw"
base_proc <- "./data/processed"

file_phossight <- file.path(base_raw, "Table_X1_PhosSight_OriginalMatrix.txt")
file_phosphors <- file.path(base_raw, "Table_X1_PhosphoRS_OriginalMatrix.txt")

# Output File
output_file <- file.path(base_proc, "Figure7a_Quantifiable_Sites_Comparison.pdf")

# ==============================================================================
# 1. Data Loading & Preprocessing
# ==============================================================================
process_counts <- function(file_path) {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  df <- fread(file_path, header = TRUE, sep = "\t")
  
  # Create Unique ID
  # Assuming GeneSymbol and PhosphoSite columns exist
  df$ID <- paste(df$GeneSymbol, df$PhosphoSite, sep = "_")
  
  # Identify Sample Columns dynamically
  meta_cols <- c("GeneSymbol", "PhosphoSite", "Protein", "Gene", "Description", "ID")
  sample_cols <- setdiff(colnames(df), meta_cols)
  
  # Extract Matrix and Count Non-NA values
  data_mat <- as.matrix(df[, ..sample_cols])
  counts <- rowSums(!is.na(data_mat))
  
  return(data.frame(ID = df$ID, Count = counts, stringsAsFactors = FALSE))
}

message(">>> Loading PhosSight Data...")
df_ps <- process_counts(file_phossight)

message(">>> Loading PhosphoRS Data...")
df_rs <- process_counts(file_phosphors)

# ==============================================================================
# 2. Set Calculations
# ==============================================================================
message(">>> Calculating Sets...")

# Define Sets based on Cutoff
ids_rs_quant <- df_rs$ID[df_rs$Count >= CUTOFF]
ids_ps_quant <- df_ps$ID[df_ps$Count >= CUTOFF]

# All IDs identified in RS (regardless of quantification)
ids_rs_all <- df_rs$ID 

# --- Categorize PhosSight Quantifiable Sites ---

# 1. Shared: Quantifiable in BOTH
green_ids <- intersect(ids_ps_quant, ids_rs_quant)
n_shared <- length(green_ids)

# 2. Newly Quantifiable (Low in RS): 
# Quantifiable in PS, Present in RS matrix, but NOT quantifiable in RS (< 92)
light_pink_ids <- setdiff(intersect(ids_ps_quant, ids_rs_all), ids_rs_quant)
n_newly_quant <- length(light_pink_ids) 

# 3. Newly Identified (Not in RS):
# Quantifiable in PS, Completely missing from RS matrix
dark_pink_ids <- setdiff(ids_ps_quant, ids_rs_all)
n_newly_id <- length(dark_pink_ids) 

# --- Categorize PhosphoRS Loss ---

# 4. Loss: Quantifiable in RS but NOT in PS
loss_ids <- setdiff(ids_rs_quant, ids_ps_quant)
n_loss <- length(loss_ids)

# --- Totals ---
n_rs_total <- length(ids_rs_quant)
n_ps_total <- length(ids_ps_quant)

# Percentage Increase
pct_increase <- round((n_ps_total - n_rs_total) / n_rs_total * 100, 1)

message(paste("PhosSight Total:", n_ps_total))
message(paste("PhosphoRS Total:", n_rs_total))
message(paste("Increase:", pct_increase, "%"))

# ==============================================================================
# 3. Prepare Plot Data
# ==============================================================================

# Legend Labels
label_not_id  <- paste0("Newly identified: ", comma(n_newly_id))
label_low_rs  <- paste0("Newly quantifiable (Low in RS): ", comma(n_newly_quant))
label_shared  <- paste0("Shared quantifiable sites: ", comma(n_shared))
label_loss    <- paste0("Loss in PhosSight: ", comma(n_loss))

# Build Data Frame
plot_df <- rbind(
  # Reference Column (PhosphoRS)
  # Uses 'Reference' type for the main bar and 'Loss' for the negative part
  data.frame(Method="PhosphoRS", Type="Reference", Count=n_shared), 
  data.frame(Method="PhosphoRS", Type="Loss",      Count=-n_loss),
  
  # PhosSight Column
  data.frame(Method="PhosSight", Type="Shared",        Count=n_shared),
  data.frame(Method="PhosSight", Type="Newly_Quant",   Count=n_newly_quant),
  data.frame(Method="PhosSight", Type="Newly_ID",      Count=n_newly_id)
)

# Set Factor Levels for Ordering
plot_df$Method <- factor(plot_df$Method, levels = c("PhosphoRS", "PhosSight"))

# Stacking Order: Newly_Quant (Top) -> Newly_ID -> Shared (Bottom) -> Reference -> Loss
# Note: For the positive stack, ggplot stacks based on level order.
plot_df$Type <- factor(plot_df$Type, 
                       levels = c("Newly_Quant", "Newly_ID", "Shared", "Reference", "Loss"))

# Calculate Label Positions (Midpoints of bars)
plot_df <- plot_df %>%
  group_by(Method) %>%
  mutate(
    pos_label = case_when(
      Type == "Loss"        ~ Count / 2, # Negative
      Type == "Reference"   ~ Count / 2,
      Type == "Shared"      ~ Count / 2,
      Type == "Newly_ID"    ~ n_shared + (Count / 2),
      Type == "Newly_Quant" ~ n_shared + n_newly_id + (Count / 2)
    ),
    label_text = comma(abs(Count))
  ) %>%
  ungroup()

# ==============================================================================
# 4. Visualization
# ==============================================================================

# Define Colors
cols <- c(
  "Reference"   = "#205595",  # Blue (Same as Shared)
  "Shared"      = "#205595",  # Blue
  "Newly_ID"    = "#C93838",  # Brick Red
  "Newly_Quant" = "#EFA93A",  # Orange
  "Loss"        = "#7F7F7F"   # Dark Grey
)

# Map Types to Legend Labels
labels_map <- c(
  "Newly_Quant" = label_low_rs,
  "Newly_ID"    = label_not_id,
  "Shared"      = label_shared,
  "Loss"        = label_loss
)

# Define Y-axis limits buffer
limit_left_buffer <- if(n_loss > 0) -n_loss * 1.5 else -2000

message(">>> Generating Plot...")

p <- ggplot(plot_df, aes(x = Method, y = Count, fill = Type)) +
  
  # Bars
  geom_col(width = 0.55, color = "black", size = 0.3) +
  
  # Zero Line
  geom_hline(yintercept = 0, color = "black", size = 0.8) +
  
  # Text Labels inside bars
  geom_text(aes(y = pos_label, label = label_text), 
            color = "white", size = 3.8, fontface = "bold") +
  
  # Colors and Legend
  scale_fill_manual(values = cols, 
                    breaks = c("Newly_Quant", "Newly_ID", "Shared", "Loss"),
                    labels = labels_map) +
  
  # Y Axis formatting
  scale_y_continuous(
    labels = function(x) comma(abs(x)), 
    limits = c(limit_left_buffer, NA),
    expand = expansion(mult = c(0.1, 0.2)) 
  ) +
  
  # Total Annotation (PhosphoRS)
  annotate("text", x = 1, y = n_rs_total + (n_ps_total * 0.05), 
           label = comma(n_rs_total), 
           size = 5, fontface = "bold", hjust = 0.5) +
  
  # Total Annotation (PhosSight)
  annotate("text", x = 2, y = n_ps_total + (n_ps_total * 0.08),
           label = paste0(comma(n_ps_total), "\n(+", pct_increase, "%)"),
           size = 5, fontface = "bold", hjust = 0.5) +
  
  # Theme
  theme_minimal() +
  labs(x = "", y = "# Quantifiable sites (< 50% missing values)", fill = "") +
  
  theme(
    panel.grid = element_blank(), 
    # Add a border around the plot panel
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(size = 14, color = "black", face = "bold", margin = margin(t=5)),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.y = element_line(color = "black", size = 0.5), 
    
    # Legend Styling
    legend.position = "top",
    legend.justification = "left", # Align legend to left
    legend.direction = "vertical", # Stack legend items
    legend.text = element_text(size = 10),
    legend.box.margin = margin(b = 10, l = 10)
  )

# ==============================================================================
# 5. Save Output
# ==============================================================================
# Ensure output directory exists
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

ggsave(output_file, p, width = 6, height = 8)
message(paste("Plot saved to:", output_file))

