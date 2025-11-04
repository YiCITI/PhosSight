# Merge four methods PSM_Site level results together
# Load required packages
load_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

load_package("svglite")
load_package("Cairo")
load_package("UpSetR")
library(tidyverse)
library(data.table)

# File paths
all_features_path <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Features.Localization.entropy.Label.txt'
psm_pga_results_path <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/PGA/psm_level/pga-peptideSummary.txt'
percolator_results_path <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/DeepRescore2/DeepRescore2.psms.txt'
percolator_without_dlfeatures_path <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/DeepRescore2/DeepRescore2WithoutDLFeatures.psms.txt'
deeprescore2_with_phossight_results_path <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/DeepRescore2Results.Label.txt'
deeprescore2_results_path <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/results_v1/DeepRescore2Results.Label.txt'

# Load and prepare base data (only once)
all_features <- fread(all_features_path)
psm_pga_results <- fread(psm_pga_results_path) %>% select(index, peptide, evalue)
colnames(psm_pga_results)[1] <- "Title"
base_data <- left_join(psm_pga_results, all_features, by = "Title")
base_data <- base_data[base_data$PhosphoLabel == 1,]

# Function to create method from base data
create_method <- function(data, seq_label_col, order_cols) {
  method <- data[data[[seq_label_col]] == 'TRUE',]
  method <- method[order(method[[order_cols[1]]], method[[order_cols[2]]], decreasing = TRUE),]
  return(method)
}

# Create Method1-6
Method1 <- create_method(base_data, "phosphors_Sequence_Label", c("PhosphoRS_IsoformProbability", "PhosphoRS_IsoformScore"))
Method2_entropy <- create_method(base_data, "pdeep_Sequence_Label", c("pDeepProb", "PhosphoRS_IsoformScore"))
Method3_RTRatio <- create_method(base_data, "autort_Sequence_Label", c("AutoRTProb", "PhosphoRS_IsoformScore"))
Method4_entropy_RTRatio <- create_method(base_data, "autort_pdeep_Sequence_Label", c("autort_pDeep_Prob", "PhosphoRS_IsoformScore"))
Method5_PhosSight <- create_method(base_data, "phosSight_Sequence_Label", c("PhosSightProb", "PhosphoRS_IsoformScore"))
Method6_entropy_RTRatio_PhosSight <- create_method(base_data, "autort_pdeep_phosSight_Sequence_Label", c("autort_pDeep_phosSight_Prob", "PhosphoRS_IsoformScore"))

# Load DeepRescore2 results
deeprescore2_percolator_results <- fread(percolator_results_path)
rescore_percolator_results <- fread(percolator_without_dlfeatures_path)
deeprescore2_results <- fread(deeprescore2_results_path)
deeprescore2_with_phossight_results <- fread(deeprescore2_with_phossight_results_path)

# Filter percolator results
rescore_percolator_results <- rescore_percolator_results[rescore_percolator_results$`q-value` <= 0.01,]
deeprescore2_percolator_results <- deeprescore2_percolator_results[deeprescore2_percolator_results$`q-value` <= 0.01,]

# Create Method7-10
Method7 <- deeprescore2_results[deeprescore2_results$Title %in% rescore_percolator_results$PSMId,]
Method7 <- create_method(Method7, "phosphors_Sequence_Label", c("PhosphoRS_IsoformProbability", "PhosphoRS_IsoformScore"))

Method8 <- deeprescore2_results[deeprescore2_results$Title %in% deeprescore2_percolator_results$PSMId,]
Method8 <- create_method(Method8, "phosphors_Sequence_Label", c("PhosphoRS_IsoformProbability", "PhosphoRS_IsoformScore"))

Method9 <- deeprescore2_results[deeprescore2_results$Title %in% deeprescore2_percolator_results$PSMId,]
Method9 <- create_method(Method9, "autort_pdeep_Sequence_Label", c("autort_pDeep_Prob", "PhosphoRS_IsoformScore"))

Method10 <- deeprescore2_with_phossight_results[deeprescore2_with_phossight_results$Title %in% deeprescore2_percolator_results$PSMId,]
Method10 <- create_method(Method10, "autort_pdeep_phosSight_Sequence_Label", c("autort_pDeep_phosSight_Prob", "PhosphoRS_IsoformScore"))

# Function to calculate FLR and filter results
calculate_flr <- function(method_data, seq_label_col, site_label_col, method_name) {
  tl <- 0
  fl <- 0
  TL_vec <- c()
  FL_vec <- c()
  cat(method_name, ':', nrow(method_data), '\n')
  
  for (i in seq_len(nrow(method_data))) {
    if (method_data[[seq_label_col]][i] == 'TRUE' && method_data[[site_label_col]][i] == 'TRUE') {
      tl <- tl + 1
    } else {
      fl <- fl + 1
    }
    TL_vec <- c(TL_vec, tl)
    FL_vec <- c(FL_vec, fl)
  }
  
  FLR_vec <- FL_vec / (FL_vec + TL_vec)
  count <- max(TL_vec[FLR_vec <= 0.01], na.rm = TRUE)
  method_data$FLR <- FLR_vec
  cat('count:', count, '\n')
  
  results <- method_data[method_data[[seq_label_col]] == 'TRUE' & 
                         method_data[[site_label_col]] == 'TRUE' & 
                         method_data$FLR <= 0.01,]
  return(list(data = method_data, results = results))
}

# Calculate FLR for all methods
calc1 <- calculate_flr(Method1, "phosphors_Sequence_Label", "phosphors_Site_Label", "Method1 (PhosphoRS)")
Method1 <- calc1$data
results1 <- calc1$results

calc2 <- calculate_flr(Method2_entropy, "pdeep_Sequence_Label", "pdeep_Site_Label", "Method2 (PhosphoRS + pDeep3)")
Method2_entropy <- calc2$data
results2 <- calc2$results

calc3 <- calculate_flr(Method3_RTRatio, "autort_Sequence_Label", "autort_Site_Label", "Method3 (PhosphoRS + AutoRT)")
Method3_RTRatio <- calc3$data
results3 <- calc3$results

calc4 <- calculate_flr(Method4_entropy_RTRatio, "autort_pdeep_Sequence_Label", "autort_pdeep_Site_Label", "Method4 (PhosphoRS + AutoRT + pDeep3)")
Method4_entropy_RTRatio <- calc4$data
results4 <- calc4$results

calc5 <- calculate_flr(Method5_PhosSight, "phosSight_Sequence_Label", "phosSight_Site_Label", "Method5 (PhosphoRS + PhosSight)")
Method5_PhosSight <- calc5$data
results5 <- calc5$results

calc6 <- calculate_flr(Method6_entropy_RTRatio_PhosSight, "autort_pdeep_phosSight_Sequence_Label", "autort_pdeep_phosSight_Site_Label", "Method6 (PhosphoRS + AutoRT + pDeep3 + PhosSight)")
Method6_entropy_RTRatio_PhosSight <- calc6$data
results6 <- calc6$results

calc7 <- calculate_flr(Method7, "phosphors_Sequence_Label", "phosphors_Site_Label", "Method7 (PhosphoRS + Rescore)")
Method7 <- calc7$data
results7 <- calc7$results

calc8 <- calculate_flr(Method8, "phosphors_Sequence_Label", "phosphors_Site_Label", "Method8 (PhosphoRS + DeepRescore)")
Method8 <- calc8$data
results8 <- calc8$results

calc9 <- calculate_flr(Method9, "autort_pdeep_Sequence_Label", "autort_pdeep_Site_Label", "Method9 (Method4 + DeepRescore)")
Method9 <- calc9$data
results9 <- calc9$results

calc10 <- calculate_flr(Method10, "autort_pdeep_phosSight_Sequence_Label", "autort_pdeep_phosSight_Site_Label", "Method10 (Method6 + DeepRescore)")
Method10 <- calc10$data
results10 <- calc10$results


#=======================Draw upset plot============================#
# Only keep 5 methods: Method1, Method2 (was Method9), Method3 (was Method7), Method4 (was Method5), Method5 (was Method10)
key1 <- paste0(results1$Title,'_',results1$IsoformSequence_PhosphoRS)
key2 <- paste0(results9$Title,'_',results9$IsoformSequence_autort_pDeep)  # Method9 -> Method2
key3 <- paste0(results7$Title,'_',results7$IsoformSequence_PhosphoRS)     # Method7 -> Method3
key4 <- paste0(results5$Title,'_',results5$IsoformSequence_PhosSight)   # Method5 -> Method4
key5 <- paste0(results10$Title,'_',results10$IsoformSequence_autort_pDeep_phosSight)  # Method10 -> Method5

totalKeys <- Reduce(union, list(key1, key2, key3, key4, key5))

results <- data.frame(key = totalKeys)
results$Method1 <- 0
results$Method2 <- 0
results$Method3 <- 0
results$Method4 <- 0
results$Method5 <- 0
rownames(results) <- results$key
results$key <- NULL

results[key1,]$Method1 <- 1
results[key2,]$Method2 <- 1
results[key3,]$Method3 <- 1
results[key4,]$Method4 <- 1
results[key5,]$Method5 <- 1

# Generate SVG output
output_file <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Figure3C_v3.svg'

tryCatch({
  svglite(output_file, width = 5, height = 5)
}, error = function(e) {
  tryCatch({
    Cairo(file = output_file, type = "svg", width = 5, height = 5)
  }, error = function(e2) {
    svg(output_file, width = 5, height = 5)
  })
})

# Draw UpSet plot
upset(results, 
      sets = c("Method5", "Method4", "Method3", "Method2", "Method1"),
      sets.bar.color = "#808080",
      main.bar.color = "#000000",
      matrix.color = "#000000",
      order.by = "freq",
      keep.order = TRUE,
      mainbar.y.label = "Results of Methods",
      sets.x.label = "#PSMs under 1% FLR",
      point.size = 1.5,
      line.size = 0.6,
      text.scale = 1)

dev.off()
cat("Output saved to:", output_file, "\n")