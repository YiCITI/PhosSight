# drawsuppfigs8 - Data Files

## 数据文件说明

本文件夹包含用于生成补充图表的R脚本和所需的数据文件。由于数据文件较大（>100MB），这些文件已被添加到 `.gitignore` 中，不会上传到GitHub。

## 所需数据文件

运行脚本前，请确保以下数据文件存在于本文件夹中：

- `PhosSight_Merged_full_WithSampleID.tsv` - PhosSight合并数据（带样本ID）
- `PhosSight_Merged_full_Normalized.tsv` - PhosSight归一化数据
- `mmc1.xlsx` - 样本信息Excel文件

## 脚本说明

1. **Create_PhosSight_Boxplot_AllSamples.R** - 创建所有样本的箱线图（log2转换+中位数中心化后）
2. **Create_PhosSight_Density_Plot.R** - 创建每个样本的密度分布图
3. **Create_PhosSight_PCA_Plot.R** - 创建主成分分析（PCA）图

## 使用方法

1. 将所需的数据文件放置在本文件夹中
2. 在R中运行相应的脚本
3. 生成的图表将保存在脚本所在目录

## 注意事项

- 所有脚本已配置为自动查找脚本所在目录中的数据文件
- 如果数据文件不在脚本目录，脚本会尝试从当前工作目录查找
- 确保数据文件格式正确，列名与脚本中使用的列名一致
- PCA脚本需要 `mmc1.xlsx` 文件来获取样本分组信息（Tumor/Normal）
