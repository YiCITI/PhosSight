# drawfig7 - Data Files

## 数据文件说明

本文件夹包含用于生成图表的R脚本和所需的数据文件。由于数据文件较大（>100MB），这些文件已被添加到 `.gitignore` 中，不会上传到GitHub。

## 所需数据文件

运行脚本前，请确保以下数据文件存在于本文件夹中：

- `Table_X1_PhosSight_OriginalMatrix.txt` - PhosSight合并数据（带样本ID）
- `Table_X1_PhosphoRS_OriginalMatrix.txt` - PhosphoRS合并数据（带样本ID）
- `UCEC_survival.txt` - UCEC数据集的生存信息
- `mmc1.xlsx` - UCEC数据集的样本信息
- `UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt` - UCEC数据集正常组织的蛋白质组信息
- `UCEC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` - UCEC数据集肿瘤组织的蛋白质组信息
- `UCEC-gene_TPM_removed_circRNA_tumor_normal_raw_log2(x+1)_BCM.xlsx` - UCEC数据集的转录组信息

## 脚本说明

1. **Step1_Dat_Preprocessing** - 数据预处理
2. **Step2_Differential_Phosphorylation_Analysis** - 肿瘤组织与正常组织的差异磷酸化位点分析
3. **Step3_Survival_Analysis_for_Phosphosites** - 磷酸化位点中位数分组的生存分析
4. **Step4_ssKSEA_Analysis** - ssKSEA激酶推断
5. **Step5_Kinase_Survival_Analysis** - 激酶的生存分析
6. **Figure7a_Comparison_of_Quantifiable_Phosphosites** - PhosphoRS与PhosSight可定量位点在50%missing value截断值下的区别
7. **Figure7b_Data_Completeness_Analysis** - PhosphoRS与PhosSight可定量位点在不同截断值下的区别
8. **Figure7c_Comparative_Analysis_of Prognostic_Phosphosites** - PhosphoRS与PhosSight磷酸化位点的生存分析区别
9. **Figure7d_Survival_Analysis_Visualization** - PhosSight的特有生存分析显著位点的KM图(以"STMN1_46S", "PARP1_368T"为例)
10. **Figure7e_Comparative_Analysis_of_Differentially_Phosphorylated_Sites** - PhosphoRS与PhosSight磷酸化位点肿瘤与正常组织的差异分析区别
11. **Figure7f_Visualization_of_Representative_Differentially_Phosphorylated_Sites** - PhosSight的特有肿瘤与正常组织的差异显著位点(以"ESR1_167S", "GSK3B_9S", "AKT1_126S", "EZH2_367T"为例)
12. **Figure7g_Comparison_of_Identified_Kinases** - PhosphoRS与PhosSight的激酶推断数量区别
13. **Figure7h_MARK2_Proteomics_mRNA_Validation** - PhosSight特有生存分析显著激酶(以"MARK2"为例)的肿瘤组织与正常组织蛋白质组/转录组差异
14. **Figure7i_Kinase_Survival_Validation** - PhosSight特有生存分析显著激酶的KM图(以"MARK2"为例)
 
## 使用方法

1. 将所需的数据文件放置在本文件夹中
2. 在R中运行相应的脚本
3. 生成的图表将保存在脚本所在目录

## 注意事项

- 所有脚本已配置为自动查找脚本所在目录中的数据文件
- 如果数据文件不在脚本目录，脚本会尝试从当前工作目录查找
- 确保数据文件格式正确，列名与脚本中使用的列名一致
