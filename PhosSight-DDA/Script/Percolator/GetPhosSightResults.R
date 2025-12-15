# Get DeepRescore2 results
library(data.table)
library(tidyverse)

if (FALSE){
featurePath <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Features'
DeepRescore2ResultsPath <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/DeepRescore2'
outputPath <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138'
} else {
args <- commandArgs(T)
featurePath <- args[1]
DeepRescore2ResultsPath <- args[2]
outputPath <- args[3]
}

all_features <- fread(paste0(featurePath,'/Features.Localization.entropy.txt'))
DeepRescore2 <- fread(paste0(DeepRescore2ResultsPath,'/PhosSight.psms.txt'))
DeepRescore2 <- DeepRescore2[DeepRescore2$`q-value`<=0.01,]
DeepRescore2 <- all_features %>% filter(Title %in% DeepRescore2$PSMId)
DeepRescore2_Results <- DeepRescore2[DeepRescore2$autort_pDeep_phosSight_Prob>=0.75,]
write.table(DeepRescore2_Results,
            paste0(outputPath,'/PhosSightResults.txt'),
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)
cat('PhosSightResults.txt has been saved to', paste0(outputPath,'/PhosSightResults.txt'), '\n')