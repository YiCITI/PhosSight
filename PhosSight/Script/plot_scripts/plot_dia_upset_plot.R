# Create UpSet plot using CSV data with command line arguments
# Method 1: Try svglite package
if (!require(svglite)) {
  install.packages("svglite", repos="https://cran.rstudio.com/")
  library(svglite)
} else {
  library(svglite)
}

# Method 2: Try Cairo package
if (!require(Cairo)) {
  install.packages("Cairo", repos="https://cran.rstudio.com/")
  library(Cairo)
} else {
  library(Cairo)
}

library(tidyverse)
library(data.table)

# Use UpSetR instead of ComplexUpset for better compatibility
if (!require(UpSetR)) {
  install.packages("UpSetR", repos="https://cran.rstudio.com/")
  library(UpSetR)
} else {
  library(UpSetR)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if CSV file path is provided
if (length(args) == 0) {
  cat("Usage: Rscript draw_upset_plot.R <csv_file_path>\n")
  cat("Example: Rscript draw_upset_plot.R csv4upsetplot/psm_upset_plot_202503_lvjiawei_V12_0h_24h_finetuned_cache_all.csv\n")
  stop("Please provide the CSV file path as a command line argument")
}

# Get CSV file path from command line argument
csv_data_path <- args[1]

# Check if file exists
if (!file.exists(csv_data_path)) {
  stop(paste("CSV file not found:", csv_data_path))
}

# Generate output SVG path by replacing .csv with .svg
output_file <- gsub("\\.csv$", ".svg", csv_data_path)

cat("Input CSV file:", csv_data_path, "\n")
cat("Output SVG file:", output_file, "\n")

# Read the CSV data
cat("Reading CSV data...\n")
upset_data <- fread(csv_data_path)

# Remove the first column (PSM identifiers) and keep only the binary data
# The first column contains the PSM identifiers, columns 2-11 contain the binary data
upset_matrix <- upset_data[, -1]  # Remove first column (PSM identifiers)

# Set row names to the PSM identifiers for reference
rownames(upset_matrix) <- upset_data[[1]]

# Convert to data frame and ensure binary format
upset_matrix <- as.data.frame(upset_matrix)
upset_matrix[upset_matrix == 1] <- 1
upset_matrix[upset_matrix == 0] <- 0

cat("Data loaded successfully. Dimensions:", nrow(upset_matrix), "x", ncol(upset_matrix), "\n")
cat("Column names:", colnames(upset_matrix), "\n")

# Print summary statistics for each method
cat("Summary of data:\n")
for (col_name in colnames(upset_matrix)) {
  count <- sum(upset_matrix[[col_name]] == 1)
  cat(paste(col_name, ":", count, "PSMs\n"))
}


#======================= Draw upset plot ============================#
# The data is already in the correct format for UpSetR
# Column order: Original, Top 10%, Top 20%, Top 30%, Top 40%, Top 50%, Top 60%, Top 70%, Top 80%, Top 90%
cat("Creating UpSet plot...\n")

# Try multiple SVG methods with error handling

# Method 1: Try svglite (most compatible)
tryCatch({
  svglite(output_file, width = 7, height = 4)
  cat("Using svglite for SVG output\n")
}, error = function(e) {
  cat("svglite failed, trying alternative methods...\n")
})

# Use UpSetR with the CSV data
# Column order: Original, Top 10%, Top 20%, Top 30%, Top 40%, Top 50%, Top 60%, Top 70%, Top 80%, Top 90%
upset(upset_matrix, 
      sets = c("Top 90%", "Top 80%", "Top 70%", "Top 60%", "Top 50%", "Top 40%", "Top 30%", "Top 20%", "Top 10%", "Original"),
      sets.bar.color = c("#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080"),  # All gray for left bars
      main.bar.color = "#000000",  # Default black for top bars
      matrix.color = "#000000",    # Default black for dots and lines
      order.by = "freq",
      keep.order = TRUE,
      mainbar.y.label = "Intersection size",
      sets.x.label = "#PSMs",
      point.size = 2.0,           # Increased point size for thicker bars
      line.size = 1.0,            # Increased line thickness
      text.scale = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.0),  # Increased text scale for all text elements
      show.numbers = "yes",       # Show numbers above bars (UpSetR expects "yes"/"no")
      mb.ratio = c(0.3, 0.7),  # Increase top bar height: 40% for main bars (with numbers), 60% for set bars
      # Reduce y-axis range by limiting intersections to decrease bar height
      nintersects = 20,           # Show top 20 intersections to reduce y-axis range
      # Use queries to set specific colors: black for all intersections
      queries = list(
        list(query = elements, params = list("Original"), color = "#000000", active = T),
        list(query = elements, params = list("Top 10%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 20%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 30%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 40%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 50%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 60%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 70%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 80%"), color = "#000000", active = T),
        list(query = elements, params = list("Top 90%"), color = "#000000", active = T)
      ))

dev.off()

# Print final output information
cat("Output saved to:", output_file, "\n")
if (file.exists(output_file)) {
  cat("✓ File successfully created\n")
} else {
  cat("✗ File creation failed\n")
}