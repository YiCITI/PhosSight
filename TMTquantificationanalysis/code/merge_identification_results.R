# R脚本：合并TMT01-TMT17的PhosSight和PhosphoRS鉴定结果
# 运行此脚本前，需要确保各批次的鉴定结果文件存在

library(data.table)

cat("====================================\n")
cat("合并 TMT01-TMT17 鉴定结果\n")
cat("====================================\n\n")

# ===== 配置路径 =====
# 所有TMT01-17的鉴定结果文件都在这个目录下
result_dir <- "D:/phossight&phosphors/TMT_Results"

# 输出文件
phosSight_output <- "PhosSightResults_TMT01_17.txt"
phosphoRS_output <- "PhosphoRSResults_TMT01_17.txt"

# ===== 合并 PhosSight 结果 =====
cat("[1] 合并 PhosSight 结果...\n")

phosSight_all <- data.table()

# 读取TMT01-17的所有批次文件
for (plex_num in 1:17) {
  plex_str <- sprintf("%02d", plex_num)
  plex_name <- paste0("TMT", plex_str)
  
  # 文件名格式: TMTXX.PhosSightResults.txt
  file_path <- file.path(result_dir, paste0(plex_name, ".PhosSightResults.txt"))
  
  if (file.exists(file_path)) {
    cat("  读取 ", plex_name, ":", file_path, "\n")
    df <- fread(file_path)
    cat("    ", plex_name, " 记录数:", nrow(df), "\n")
    phosSight_all <- rbind(phosSight_all, df, fill = TRUE)
  } else {
    cat("  [警告] ", plex_name, " 文件不存在:", file_path, "\n")
  }
}

if (nrow(phosSight_all) > 0) {
  cat("  合并后总记录数:", nrow(phosSight_all), "\n")
  
  # 验证 Title 前缀分布
  cat("  验证 Title 前缀分布:\n")
  for (prefix in sprintf("%02d", 1:17)) {
    count <- sum(grepl(paste0("^", prefix), phosSight_all$Title))
    if (count > 0) {
      cat("    ", prefix, ":", count, "\n")
    }
  }
  
  # 保存
  cat("  保存合并后的文件:", phosSight_output, "\n")
  fwrite(phosSight_all, phosSight_output, sep = "\t")
  cat("  ✓ PhosSight 合并完成！\n")
} else {
  cat("  [错误] 没有找到任何 PhosSight 结果文件！\n")
}

# ===== 合并 PhosphoRS 结果 =====
cat("\n[2] 合并 PhosphoRS 结果...\n")

phosphoRS_all <- data.table()

# 读取TMT01-17的所有批次文件
for (plex_num in 1:17) {
  plex_str <- sprintf("%02d", plex_num)
  plex_name <- paste0("TMT", plex_str)
  
  # 文件名格式: TMTXX.Method1Results.txt
  file_path <- file.path(result_dir, paste0(plex_name, ".Method1Results.txt"))
  
  if (file.exists(file_path)) {
    cat("  读取 ", plex_name, ":", file_path, "\n")
    df <- fread(file_path)
    cat("    ", plex_name, " 记录数:", nrow(df), "\n")
    phosphoRS_all <- rbind(phosphoRS_all, df, fill = TRUE)
  } else {
    cat("  [警告] ", plex_name, " 文件不存在:", file_path, "\n")
  }
}

if (nrow(phosphoRS_all) > 0) {
  cat("  合并后总记录数:", nrow(phosphoRS_all), "\n")
  
  # 验证 Title 前缀分布
  cat("  验证 Title 前缀分布:\n")
  for (prefix in sprintf("%02d", 1:17)) {
    count <- sum(grepl(paste0("^", prefix), phosphoRS_all$Title))
    if (count > 0) {
      cat("    ", prefix, ":", count, "\n")
    }
  }
  
  # 保存
  cat("  保存合并后的文件:", phosphoRS_output, "\n")
  fwrite(phosphoRS_all, phosphoRS_output, sep = "\t")
  cat("  ✓ PhosphoRS 合并完成！\n")
} else {
  cat("  [错误] 没有找到任何 PhosphoRS 结果文件！\n")
}

# ===== 合并 DeepRescore2 结果 =====
cat("\n[3] 合并 DeepRescore2 结果...\n")

deepRescore2_all <- data.table()

# 读取TMT01-17的所有批次文件
for (plex_num in 1:17) {
  plex_str <- sprintf("%02d", plex_num)
  plex_name <- paste0("TMT", plex_str)
  
  # 文件名格式: TMTXX.DeepRescore2Results.txt
  file_path <- file.path(result_dir, paste0(plex_name, ".DeepRescore2Results.txt"))
  
  if (file.exists(file_path)) {
    cat("  读取 ", plex_name, ":", file_path, "\n")
    df <- fread(file_path)
    cat("    ", plex_name, " 记录数:", nrow(df), "\n")
    deepRescore2_all <- rbind(deepRescore2_all, df, fill = TRUE)
  } else {
    cat("  [警告] ", plex_name, " 文件不存在:", file_path, "\n")
  }
}

if (nrow(deepRescore2_all) > 0) {
  cat("  合并后总记录数:", nrow(deepRescore2_all), "\n")
  
  # 验证 Title 前缀分布
  cat("  验证 Title 前缀分布:\n")
  for (prefix in sprintf("%02d", 1:17)) {
    count <- sum(grepl(paste0("^", prefix), deepRescore2_all$Title))
    if (count > 0) {
      cat("    ", prefix, ":", count, "\n")
    }
  }
  
  # 保存
  cat("  保存合并后的文件:", deepRescore2_output, "\n")
  fwrite(deepRescore2_all, deepRescore2_output, sep = "\t")
  cat("  ✓ DeepRescore2 合并完成！\n")
} else {
  cat("  [错误] 没有找到任何 DeepRescore2 结果文件！\n")
}

cat("\n====================================\n")
cat("合并完成！\n")
cat("====================================\n\n")

cat("输出文件:\n")
cat("  - ", phosSight_output, "\n")
cat("  - ", phosphoRS_output, "\n")
cat("  - ", deepRescore2_output, "\n")
cat("\n下一步: 运行 step1_GetIonIntensity 脚本\n")
