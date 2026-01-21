# drawfig7 - Data Files

## 数据文件说明

本文件夹包含用于生成图表的R脚本和所需的数据文件。由于数据文件较大（>100MB），这些文件已被添加到 `.gitignore` 中，不会上传到GitHub。

## 所需数据文件

运行脚本前，请确保以下数据文件存在于本文件夹中：

- `PhosSight_Merged_full_WithSampleID.tsv` - PhosSight合并数据（带样本ID）
- `PhosphoRS_Merged_full_WithSampleID.tsv` - PhosphoRS合并数据（带样本ID）

## 脚本说明

1. **Create_QuantifiableSites_100Samples.R** - 创建100样本阈值下的可定量位点对比图
2. **Create_PhosSight_vs_PhosphoRS_GainSharedLoss.R** - 创建PhosSight vs PhosphoRS的增益/共享/丢失对比图
3. **Create_MissingValue_SiteCount_Comparison.R** - 创建不同缺失值阈值下的位点数量对比图

## 使用方法

1. 将所需的数据文件放置在本文件夹中
2. 在R中运行相应的脚本
3. 生成的图表将保存在脚本所在目录

## 注意事项

- 所有脚本已配置为自动查找脚本所在目录中的数据文件
- 如果数据文件不在脚本目录，脚本会尝试从当前工作目录查找
- 确保数据文件格式正确，列名与脚本中使用的列名一致
