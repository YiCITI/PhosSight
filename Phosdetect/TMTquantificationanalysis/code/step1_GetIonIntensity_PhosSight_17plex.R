# R脚本：将MASCI定量结果与PhosSight鉴定结果匹配（17 plex完整版）
# TMT01-TMT17共17个批次，170个TMT通道
# 基于成功的 4plex 脚本修改

library(tidyverse)
library(data.table)

# ===== 配置路径 =====
# MASCI定量结果目录（TMT01-TMT17）
pathin_list <- c(
  'D:/MASCIresult/QuantificationResults01',
  'D:/MASCIresult/QuantificationResults02',
  'D:/MASCIresult/QuantificationResults03',
  'D:/MASCIresult/QuantificationResults04',
  'D:/MASCIresult/QuantificationResults05',
  'D:/MASCIresult/QuantificationResults06',
  'D:/MASCIresult/QuantificationResults07',
  'D:/MASCIresult/QuantificationResults08',
  'D:/MASCIresult/QuantificationResults09',
  'D:/MASCIresult/QuantificationResults10',
  'D:/MASCIresult/QuantificationResults11',
  'D:/MASCIresult/QuantificationResults12',
  'D:/MASCIresult/QuantificationResults13',
  'D:/MASCIresult/QuantificationResults14',
  'D:/MASCIresult/QuantificationResults15',
  'D:/MASCIresult/QuantificationResults16',
  'D:/MASCIresult/QuantificationResults17'
)

# plex列表
plex_list <- c('TMT01', 'TMT02', 'TMT03', 'TMT04', 'TMT05', 'TMT06', 'TMT07', 'TMT08', 
               'TMT09', 'TMT10', 'TMT11', 'TMT12', 'TMT13', 'TMT14', 'TMT15', 'TMT16', 'TMT17')

# PhosSight鉴定结果文件（需要先合并所有批次的鉴定结果）
pathin2 = 'PhosSightResults_TMT01_17.txt'

# 输出文件
pathout = 'PhosSight_Intensity_full.txt'

# ===== 读取PhosSight鉴定结果 =====
cat("========================================\n")
cat("Step 1: 读取PhosSight鉴定结果（TMT01-17）\n")
cat("========================================\n\n")

if (!file.exists(pathin2)) {
  stop(paste0("错误: 鉴定结果文件不存在: ", pathin2, "\n请先运行 merge_identification_results.R 合并鉴定结果"))
}

PhosSight <- fread(pathin2)

# 选择需要的列
cols_to_select <- c("Title", "Charge", "Peptide", "Proteins")

if ("IsoformSequence_PhosSight" %in% colnames(PhosSight)) {
  cols_to_select <- c(cols_to_select, "IsoformSequence_PhosSight")
} else if ("Mod_Sequence_for_phosphoRS" %in% colnames(PhosSight)) {
  cols_to_select <- c(cols_to_select, "Mod_Sequence_for_phosphoRS")
}

if ("Modifications_abbrev" %in% colnames(PhosSight)) {
  cols_to_select <- c(cols_to_select, "Modifications_abbrev")
} else if ("Modification" %in% colnames(PhosSight)) {
  cols_to_select <- c(cols_to_select, "Modification")
}

PhosSight <- PhosSight %>% select(all_of(cols_to_select))

# 构建Title2用于匹配
cat("构建Title2用于匹配...\n")
PhosSight$Title2 <- as.character(sapply(PhosSight$Title, function(x) {
  parts <- strsplit(x, '[.]')[[1]]
  if (length(parts) >= 3) {
    paste(parts[1], parts[2], parts[3], sep = '.')
  } else {
    NA_character_
  }
}, USE.NAMES = FALSE))

# 提取TMT plex信息（从Title前缀判断）
PhosSight$Plex <- sapply(PhosSight$Title, function(x) {
  prefix <- substr(x, 1, 2)
  plex_num <- as.numeric(prefix)
  if (!is.na(plex_num) && plex_num >= 1 && plex_num <= 17) {
    return(paste0("TMT", sprintf("%02d", plex_num)))
  } else {
    return(NA)
  }
})

cat("PhosSight总记录数:", nrow(PhosSight), "\n\n")
cat("各批次记录数:\n")
for (plex in plex_list) {
  count <- sum(PhosSight$Plex == plex, na.rm = TRUE)
  cat(paste0("  ", plex, ": ", count, "\n"))
}

# ===== 读取并处理每个MASCI批次的数据 =====
cat("\n========================================\n")
cat("Step 2: 读取MASCI定量文件\n")
cat("========================================\n\n")

ion_cols <- c('Ion_126.128', 'Ion_127.125', 'Ion_127.131', 'Ion_128.128', 'Ion_128.134',
              'Ion_129.131', 'Ion_129.138', 'Ion_130.135', 'Ion_130.141', 'Ion_131.138')

Results <- PhosSight

# 为每个plex初始化空列
for (plex in plex_list) {
  for (col in ion_cols) {
    new_col <- paste0(plex, "_", col)
    Results[[new_col]] <- NA_real_
  }
  Results[[paste0(plex, "_PeakArea")]] <- NA_real_
  Results[[paste0(plex, "_ParentIonIntensity")]] <- NA_real_
}

# 处理每个MASCI目录
total_matched <- 0
for (idx in seq_along(pathin_list)) {
  pathin1 <- pathin_list[idx]
  plex <- plex_list[idx]
  
  if (!dir.exists(pathin1)) {
    cat("  [跳过]", plex, "- 目录不存在:", pathin1, "\n")
    next
  }
  
  cat("  处理", plex, "- 从目录读取:", pathin1, "\n")
  
  # 读取MASCI数据（只读取 ReporterIons 和 SICstats 文件）
  totalData <- data.table()
  FileNames <- list.files(path = pathin1)
  
  for (i in seq_along(FileNames)) {
    FileName <- FileNames[i]
    
    if (grepl("_ReporterIons.txt", FileName, fixed = TRUE)) {
      name <- gsub("_ReporterIons.txt", "", FileName)
      FileName2 <- paste0(name, "_SICstats.txt")
      
      file_path1 <- file.path(pathin1, FileName)
      if (!file.exists(file_path1)) next
      data1 <- fread(file_path1)
      
      file_path2 <- file.path(pathin1, FileName2)
      if (!file.exists(file_path2)) next
      data2 <- fread(file_path2)
      
      data <- left_join(data1, data2, by = c('ScanNumber' = 'FragScanNumber')) %>%
        select(ScanNumber, all_of(ion_cols), PeakArea, ParentIonIntensity)
      
      # 构建Title2 (ScanNumber)
      data$ScanNumber <- paste0(name, '.', data$ScanNumber, '.', data$ScanNumber)
      
      if (nrow(totalData) == 0) {
        totalData <- data
      } else {
        totalData <- rbind(totalData, data, fill = TRUE)
      }
    }
  }
  
  cat("    MASCI数据记录数:", nrow(totalData), "\n")
  
  if (nrow(totalData) == 0) {
    cat("    [警告] 没有找到MASCI数据，跳过\n")
    next
  }
  
  # 计算归一化的TMT强度
  total_intensity <- rowSums(totalData[, ..ion_cols], na.rm = TRUE)
  total_intensity[total_intensity == 0] <- NA
  
  # 归一化
  for (col in ion_cols) {
    totalData[[col]] <- (totalData[[col]] / total_intensity) * totalData$PeakArea
  }
  
  # 匹配到对应plex的PSM
  # 优先使用Title前缀判断（更可靠），如果Title前缀匹配不到，再尝试Plex列
  plex_num <- as.numeric(substr(plex, 4, 5))
  prefix <- sprintf("%02d", plex_num)
  plex_mask <- substr(Results$Title, 1, 2) == prefix
  
  # 如果Title前缀匹配不到数据，尝试使用Plex列（作为fallback）
  if (sum(plex_mask, na.rm = TRUE) == 0 && 'Plex' %in% colnames(Results)) {
    plex_mask <- Results$Plex == plex & !is.na(Results$Plex)
    if (sum(plex_mask, na.rm = TRUE) > 0) {
      cat("    [注意] Title前缀", prefix, "匹配不到数据，使用Plex列匹配\n")
    }
  }
  
  plex_results <- Results[plex_mask, ]
  
  if (nrow(plex_results) == 0) {
    cat("    [警告] 没有找到", plex, "的PSM数据，跳过匹配\n\n")
    next
  }
  
  cat("    找到", plex, "的PSM数:", nrow(plex_results), "\n")
  
  # 匹配
  matched <- left_join(
    plex_results %>% select(Title2),
    totalData %>% rename(Title2 = ScanNumber),
    by = "Title2"
  )
  
  # 填充对应plex的列
  for (col in ion_cols) {
    new_col <- paste0(plex, "_", col)
    Results[[new_col]][plex_mask] <- matched[[col]]
  }
  Results[[paste0(plex, "_PeakArea")]][plex_mask] <- matched$PeakArea
  Results[[paste0(plex, "_ParentIonIntensity")]][plex_mask] <- matched$ParentIonIntensity
  
  matched_count <- sum(!is.na(matched$PeakArea))
  total_matched <- total_matched + matched_count
  cat("    成功匹配记录数:", matched_count, "\n\n")
}

# ===== 输出结果 =====
cat("========================================\n")
cat("Step 3: 保存结果\n")
cat("========================================\n\n")

output_cols <- c('Title', 'Charge', 'Peptide', 'Proteins', 'Plex')

if ("IsoformSequence_PhosSight" %in% colnames(Results)) {
  output_cols <- c(output_cols, 'IsoformSequence_PhosSight')
} else if ("Mod_Sequence_for_phosphoRS" %in% colnames(Results)) {
  output_cols <- c(output_cols, 'Mod_Sequence_for_phosphoRS')
}

if ("Modifications_abbrev" %in% colnames(Results)) {
  output_cols <- c(output_cols, 'Modifications_abbrev')
} else if ("Modification" %in% colnames(Results)) {
  output_cols <- c(output_cols, 'Modification')
}

# 添加所有TMT通道列
for (plex in plex_list) {
  output_cols <- c(output_cols, paste0(plex, "_PeakArea"), paste0(plex, "_ParentIonIntensity"))
  for (col in ion_cols) {
    output_cols <- c(output_cols, paste0(plex, "_", col))
  }
}

Results <- Results %>% select(all_of(output_cols))

cat("匹配后的总记录数:", nrow(Results), "\n")
cat("总匹配数:", total_matched, "\n\n")

cat("各批次匹配统计:\n")
for (plex in plex_list) {
  matched_count <- sum(!is.na(Results[[paste0(plex, "_PeakArea")]]))
  cat(paste0("  ", plex, ": ", matched_count, "\n"))
}

# 保存结果
cat("\n保存结果到:", pathout, "\n")
write.table(Results,
            pathout,
            sep = '\t',
            quote = FALSE,
            row.names = FALSE)

cat("\n========================================\n")
cat("完成！\n")
cat("========================================\n")
cat("输出文件:", pathout, "\n")
cat("总列数:", ncol(Results), "\n")
cat("TMT通道列:", length(ion_cols) * length(plex_list), "个（", length(plex_list), "个plex各", length(ion_cols), "个）\n")

