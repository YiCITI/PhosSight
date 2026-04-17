#Figure7g
# ==============================================================================
# Script Name: Comparison of Identified Kinases (Figure 7g)
# Description: Generates a stacked bar chart comparing the number of kinases 
#              inferred by PhosphoRS vs PhosSight (via ssKSEA).
#
# Usage:       Run from project root. Requires KSEA Score Tables (Table X7).
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

# ==============================================================================
# Configuration
# ==============================================================================
# Input Files (Relative Paths)
base_proc <- "./data/processed"
file_ps <- file.path(base_proc, "Table_X7_PhosSight_KSEA_Scores_50pct.txt")
file_rs <- file.path(base_proc, "Table_X7_PhosphoRS_KSEA_Scores_50pct.txt")

# Output File
output_file <- file.path(base_proc, "Figure7g_Kinase_Comparison.pdf")

# ==============================================================================
# 1. Data Loading & Set Analysis
# ==============================================================================
# Function to extract unique kinases from file
get_kinases <- function(fpath) {
  if(!file.exists(fpath)) stop(paste("File not found:", fpath))
  
  dt <- fread(fpath, header = TRUE)
  
  # Robust column detection
  if("Kinase" %in% colnames(dt)) {
    return(unique(dt$Kinase))
  } else {
    warning("Column 'Kinase' not found, using first column.")
    return(unique(dt[[1]])) 
  }
}

message(">>> Loading Kinase Lists...")
vec_ps <- get_kinases(file_ps)
vec_rs <- get_kinases(file_rs)

message(paste("PhosSight Total:", length(vec_ps)))
message(paste("PhosphoRS Total:", length(vec_rs)))

# --- Set Operations ---
# Shared: Intersection
shared_ids <- intersect(vec_rs, vec_ps)
n_shared <- length(shared_ids)

# Gain: PhosSight Only
gain_ids <- setdiff(vec_ps, vec_rs)
n_gain <- length(gain_ids)

# Loss: PhosphoRS Only (Optional check)
loss_ids <- setdiff(vec_rs, vec_ps)
n_loss <- length(loss_ids)

# Totals
total_rs <- length(vec_rs)
total_ps <- length(vec_ps)
pct_increase <- round((total_ps - total_rs) / total_rs * 100, 1)

if (n_loss > 0) {
  message(paste("Note: There are", n_loss, "kinases unique to PhosphoRS (Loss)."))
  message("Current plot logic focuses on Shared + Gain.")
}

# ==============================================================================
# 2. Prepare Plot Data
# ==============================================================================
# Label Text
label_shared <- paste0("Shared Kinases: ", comma(n_shared))
label_gain   <- paste0("Gain in PhosSight: ", comma(n_gain))

plot_df <- rbind(
  # PhosphoRS Bar (Reference)
  data.frame(Method="PhosphoRS", Type="Shared", Count=n_shared),
  
  # PhosSight Bar
  data.frame(Method="PhosSight", Type="Shared", Count=n_shared),
  data.frame(Method="PhosSight", Type="Gain",   Count=n_gain)
)

# Factor Ordering
plot_df$Method <- factor(plot_df$Method, levels = c("PhosphoRS", "PhosSight"))
# Stack Order: Shared (Bottom), Gain (Top)
plot_df$Type <- factor(plot_df$Type, levels = c("Gain", "Shared"))

# Calculate Label Positions (Midpoints)
plot_df <- plot_df %>%
  group_by(Method) %>%
  mutate(
    pos_label = case_when(
      Type == "Shared" ~ Count / 2,
      Type == "Gain"   ~ n_shared + (Count / 2)
    ),
    label_text = comma(abs(Count))
  ) %>%
  ungroup()

# ==============================================================================
# 3. Visualization
# ==============================================================================
# Colors
cols <- c(
  "Gain"   = "#EFA93A", # Orange
  "Shared" = "#205595"  # Blue
)

# Legend Map
labels_map <- c(
  "Gain"   = label_gain,
  "Shared" = label_shared
)

message(">>> Generating Plot...")

p <- ggplot(plot_df, aes(x = Method, y = Count, fill = Type)) +
  
  # Bars
  geom_col(width = 0.6, color = "black", size = 0.5) +
  
  # Text Labels inside bars
  geom_text(aes(y = pos_label, label = ifelse(Count != 0, label_text, "")), 
            color = "white", size = 4.5, fontface = "bold") +
  
  # Colors
  scale_fill_manual(values = cols, labels = labels_map,
                    breaks = c("Gain", "Shared")) +
  
  # Y Axis
  scale_y_continuous(
    labels = comma, 
    limits = c(0, NA), 
    expand = expansion(mult = c(0, 0.2)) # Add space at top for annotation
  ) +
  
  # Annotation: PhosphoRS Total
  annotate("text", x = 1, y = n_shared + (total_ps * 0.05), 
           label = comma(n_shared), 
           size = 5, fontface = "bold", hjust = 0.5) +
  
  # Annotation: PhosSight Total + Increase
  annotate("text", x = 2, y = n_shared + n_gain + (total_ps * 0.05),
           label = paste0("Total: ", comma(total_ps), "\n(+", pct_increase, "%)"),
           size = 5, fontface = "bold", hjust = 0.5) +
  
  # Theme: Classic provides L-shaped axis lines by default
  theme_classic() + 
  
  labs(x = "", y = "# Identified Kinases", fill = "") + 
  
  theme(
    # Axes
    axis.line = element_line(color = "black", size = 0.8), # Thicker L-shape
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    
    # Legend
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text = element_text(size = 11),
    legend.box.margin = margin(b = 10)
  )

# ==============================================================================
# 4. Save Output
# ==============================================================================
# Ensure output directory exists
if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

ggsave(output_file, p, width = 5, height = 7)
message(paste("Plot saved to:", output_file))
