library(tidyverse)
library(data.table)

# 生成传统Percolator输入文件（不包含深度学习特征）
# 用于比较传统重评分和深度学习重评分的差异

if (FALSE){
features <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Features/Features.Localization.entropy.txt')
raw_psms <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/PGA/pga-rawPSMs.txt') %>% select(index)
output <- 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/Traditional.pin'
software='maxquant'
phosphors <- fread('E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/PhosphoRS/PhosphoRS.txt')
VariableMods = '1,Oxidation,M,15.994919,1;2,Phospho,S,79.966331,2;3,Phospho,T,79.966331,2;4,Phospho,Y,79.966331,2'
FixedMods = 'null'
} else {
args <- commandArgs(T)
features <- fread(args[1])
raw_psms <- fread(args[2]) %>% select(index)
output <- args[3]
software <- args[4]
phosphors <- fread(args[5])
VariableMods <- args[6]
FixedMods <- args[7]
}

VariableModsInfo = unlist(strsplit(VariableMods,';',fixed = TRUE))
if (FixedMods!='null'){
  FixedModsInfo = unlist(strsplit(FixedMods,';',fixed = TRUE))
}

# 合并特征和PSM数据
data <- left_join(raw_psms, features, by="index")
colnames(data)[1] <- "PSMId"

# 添加PhosphoRS结果
data <- left_join(data, phosphors, by = c("PSMId" = "Spectrum.Name"))

# 根据搜索引擎添加相应特征
if (software == "maxquant"){
  # MaxQuant特征
  data$score <- data$Score
  data$deltascore <- data$DeltaScore
  data$peplen <- data$Length
  data$charge <- data$Charge
  data$mass <- data$Mass
  data$deltamass <- data$MassError
  data$absdeltamass <- abs(data$MassError)
  data$isdecoy <- ifelse(grepl("^REV_", data$Proteins), 1, 0)
  
} else if (software == "comet"){
  # Comet特征
  data$score <- data$XCorr
  data$deltascore <- data$DeltaCn
  data$peplen <- data$Length
  data$charge <- data$Charge
  data$mass <- data$Mass
  data$deltamass <- data$MassError
  data$absdeltamass <- abs(data$MassError)
  data$isdecoy <- ifelse(grepl("^REV_", data$Proteins), 1, 0)
  
} else if (software == "msgf"){
  # MS-GF+特征
  data$score <- data$SpecEValue
  data$deltascore <- data$DeltaScore
  data$peplen <- data$Length
  data$charge <- data$Charge
  data$mass <- data$Mass
  data$deltamass <- data$MassError
  data$absdeltamass <- abs(data$MassError)
  data$isdecoy <- ifelse(grepl("^REV_", data$Proteins), 1, 0)
  
} else if (software == "xtandem"){
  # X!Tandem特征
  data$score <- data$Hyperscore
  data$deltascore <- data$DeltaScore
  data$peplen <- data$Length
  data$charge <- data$Charge
  data$mass <- data$Mass
  data$deltamass <- data$MassError
  data$absdeltamass <- abs(data$MassError)
  data$isdecoy <- ifelse(grepl("^REV_", data$Proteins), 1, 0)
}

# 添加磷酸化相关特征
data$phospho <- ifelse(data$PhosphoLabel == 1, 1, 0)
data$phosphors_score <- data$PhosphoRS_IsoformScore
data$phosphors_prob <- data$PhosphoRS_IsoformProbability

# 选择Percolator需要的列
percolator_data <- data %>%
  select(PSMId, score, deltascore, peplen, charge, mass, deltamass, absdeltamass, 
         isdecoy, phospho, phosphors_score, phosphors_prob)

# 处理缺失值
percolator_data[is.na(percolator_data)] <- 0

# 写入Percolator输入文件
write.table(percolator_data, output, sep = "\t", row.names = FALSE, quote = FALSE)

cat('Traditional Percolator input file saved to:', output, '\n')
