# Create UpSet plot for A549 dataset using CSV data
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

# Check if required arguments are provided
if (length(args) < 2) {
  cat("Usage: Rscript generate_upset_plot_A549.R --peptides_input <path> --psms_input <path> --output_dir <path> --dataset_name <name>\n")
  cat("Example: Rscript generate_upset_plot_A549.R --peptides_input upset_plot_202503_A549_finetuned_phos_peptides.csv --psms_input upset_plot_202503_A549_finetuned_phos_PSMs.csv --output_dir ./output/ --dataset_name A549\n")
  stop("Please provide required arguments")
}

# Parse arguments
i <- 1
while (i <= length(args)) {
  if (args[i] == "--peptides_input") {
    peptides_input <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--psms_input") {
    psms_input <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--dataset_name") {
    dataset_name <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Dataset:", dataset_name, "\n")
cat("Peptides input:", peptides_input, "\n")
cat("PSMs input:", psms_input, "\n")
cat("Output directory:", output_dir, "\n")

# Function to generate upset plot from CSV
generate_upset_plot <- function(csv_data_path, output_suffix, plot_title) {
  # Check if file exists
  if (!file.exists(csv_data_path)) {
    cat(paste("Warning: CSV file not found:", csv_data_path, "\n"))
    return(NULL)
  }

  # Generate output SVG path
  output_file <- file.path(output_dir, paste0(dataset_name, "_", output_suffix, ".svg"))

  cat("Processing:", csv_data_path, "\n")
  cat("Output:", output_file, "\n")

  # Read the CSV data
  cat("Reading CSV data...\n")
  upset_data <- fread(csv_data_path)

  # Remove the first column (PSM/peptide identifiers) and keep only the binary data
  upset_matrix <- upset_data[, -1]

  # Set row names to the PSM/peptide identifiers for reference
  rownames(upset_matrix) <- upset_data[[1]]

  # Convert to data frame and ensure binary format
  upset_matrix <- as.data.frame(upset_matrix)
  upset_matrix[upset_matrix == 1] <- 1
  upset_matrix[upset_matrix == 0] <- 0

  cat("Data loaded successfully. Dimensions:", nrow(upset_matrix), "x", ncol(upset_matrix), "\n")
  cat("Column names:", colnames(upset_matrix), "\n")

  # Print summary statistics
  cat("Summary of data:\n")
  for (col_name in colnames(upset_matrix)) {
    count <- sum(upset_matrix[[col_name]] == 1)
    cat(paste(col_name, ":", count, "items\n"))
  }

  #======================= Draw upset plot ============================#
  # The data is already in the correct format for UpSetR
  # Column order: Original, Top 10%, Top 20%, Top 30%, Top 40%, Top 50%, Top 60%, Top 70%, Top 80%, Top 90%
  cat("Creating UpSet plot...\n")

  # Font size 7 is set via pointsize = 7 in graphics device
  # text.scale in UpSetR: c(intersection size title, intersection size tick labels,
  #                         set size title, set size tick labels, set names, numbers above bars)
  # text.scale = 1 means no additional scaling, so all text uses pointsize = 7

  # Set global graphics parameters to match matplotlib
  # Try to use Arial first, fallback to sans-serif
  if (Sys.info()["sysname"] == "Windows") {
    # On Windows, Arial is usually available
    par(family = "sans")  # Use sans-serif font family (Arial on Windows)
  } else {
    # On Linux/Mac, try to use Arial or DejaVu Sans
    par(family = "sans")  # Use sans-serif font family
  }

  # Explicitly set cex = 1 to ensure no global scaling
  # text.scale in UpSetR handles font scaling independently
  par(cex = 1)  # No global character expansion (ensure font size = 7 from pointsize)
  par(lwd = 1)  # Line width = 1 (axes.linewidth, lines.linewidth, patch.linewidth)
  par(tcl = -0.2)  # Tick mark length
  par(mgp = c(2, 0.5, 0))  # Margin line positions (axis title, axis labels, axis line)
  par(las = 1)  # Labels parallel to axis

  # Try multiple SVG methods with error handling
  # Apply matplotlib-style settings after opening device

  # Method 1: Try svglite (most compatible)
  # Reduced figure size to 3.5x3.5 to make fonts appear larger relative to figure size
  fig_width = 6.5
  fig_height = 4.5
  pointsize = 7
  tryCatch({
    svglite(output_file, width = fig_width, height = fig_height, pointsize = pointsize)  # Set point size to 7
    cat("Using svglite for SVG output\n")
    # Re-apply graphics parameters after opening device (cex = 1 ensures no scaling, font size = 7)
    par(family = "sans", cex = 1, lwd = 1, tcl = -0.2, mgp = c(2, 0.5, 0), las = 1)
  }, error = function(e) {
    cat("svglite failed, trying alternative methods...\n")

    # Method 2: Try Cairo SVG
    tryCatch({
      Cairo(file = output_file, type = "svg", width = fig_width, height = fig_height, pointsize = pointsize)
        cat("Using Cairo for SVG output\n")
      # Re-apply graphics parameters after opening device (cex = 1 ensures no scaling, font size = 7)
      par(family = "sans", cex = 1, lwd = 1, tcl = -0.2, mgp = c(2, 0.5, 0), las = 1)
    }, error = function(e2) {
      cat("Cairo SVG failed, trying basic svg()...\n")

      # Method 3: Try basic svg() function
      tryCatch({
        svg(output_file, width = fig_width, height = fig_height, pointsize = pointsize)
          cat("Using basic svg() function\n")
        # Re-apply graphics parameters after opening device (cex = 1 ensures no scaling, font size = 7)
        par(family = "sans", cex = 1, lwd = 1, tcl = -0.2, mgp = c(2, 0.5, 0), las = 1)
      }, error = function(e3) {
        cat("All SVG methods failed, falling back to PDF\n")
        # Fallback to PDF
        output_file <<- gsub("\\.svg$", ".pdf", output_file)
        pdf(output_file, width = fig_width, height = fig_height, pointsize = pointsize)
        # Re-apply graphics parameters after opening device (cex = 1 ensures no scaling, font size = 7)
        par(family = "sans", cex = 1, lwd = 1, tcl = -0.2, mgp = c(2, 0.5, 0), las = 1)
      })
    })
  })

  # Ensure font size parameters are set before calling upset()
  # This ensures all text elements, including numbers above bars, use font size 7
  par(cex = 1)  # No character expansion
  par(ps = 7)   # Explicitly set point size to 7

  upset(upset_matrix,
        sets = c("Top 90%", "Top 80%", "Top 70%", "Top 60%", "Top 50%", "Top 40%", "Top 30%", "Top 20%", "Top 10%", "Original"),
        sets.bar.color = c("#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080", "#808080"),  # All gray for left bars
        main.bar.color = "#000000",  # Default black for top bars
        matrix.color = "#000000",    # Default black for dots and lines
        order.by = "freq",
        keep.order = TRUE,
        mainbar.y.label = "Intersection size",
        sets.x.label = "#Items",
        point.size = 1.0,           # patch.linewidth = 1
        line.size = 1.0,            # lines.linewidth = 1
        # text.scale: c(intersection size title, intersection size tick labels,
        #               set size title, set size tick labels, set names, numbers above bars)
        # All elements set to 1 to use pointsize = 7 without additional scaling
        # The 6th element specifically controls numbers above bars
        text.scale = c(1, 1, 1, 1, 1, 1),  # All text elements (including numbers above bars) use pointsize = 7
        show.numbers = "yes",       # Show numbers above bars (UpSetR expects "yes"/"no")
        number.angles = 0,          # Set number angle to 0 (horizontal) for better readability
        mb.ratio = c(0.5, 0.5),  # (0.3, 0.7)
        # Reduce y-axis range by limiting intersections to decrease bar height
        nintersects = 20)

  dev.off()

  # Print final output information
  cat("Output saved to:", output_file, "\n")
  if (file.exists(output_file)) {
    cat("✓ File successfully created\n")
  } else {
    cat("✗ File creation failed\n")
  }
  return(output_file)
}

# Generate plots for both peptides and PSMs
cat("\n=================== Generating UpSet Plots ===================\n\n")

# Peptides upset plot
cat("Generating peptides upset plot...\n")
peptides_output <- generate_upset_plot(peptides_input, "peptides_upset", "Phosphopeptides")

# PSMs upset plot
cat("\nGenerating PSMs upset plot...\n")
psms_output <- generate_upset_plot(psms_input, "psms_upset", "Phospho-PSMs")

cat("\n=================== Summary ===================\n")
cat("Peptides plot:", ifelse(is.null(peptides_output), "Failed", peptides_output), "\n")
cat("PSMs plot:", ifelse(is.null(psms_output), "Failed", psms_output), "\n")
cat("\nAll plots generated successfully!\n")
