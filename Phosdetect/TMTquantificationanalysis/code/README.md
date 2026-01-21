# TMT01-17 完整版流程（Full SiteQuant）

## 概述

这套流程处理 TMT01-TMT17 共 17 个批次的磷酸化蛋白组数据，使用**完整版 sitequant** 进行位点定量。

### 完整版 vs 简化版 SiteQuant

| 特性 | 完整版 | 简化版 |
|------|--------|--------|
| 位点标识 | 蛋白名_氨基酸_位置 (如 `sp\|Q13523\|PRP4B_HUMAN_S45`) | 修饰序列字符串 |
| FASTA依赖 | 需要 | 不需要 |
| 位点定位 | 精确（基于蛋白序列） | 近似（基于修饰序列） |
| 适用场景 | motif分析、位点注释、数据库比对 | 快速统计、缺失值分析 |

## 数据配置

- **TMT批次**: TMT01-TMT17（共17个）
- **TMT通道**: 170个（17 × 10）
- **样本通道**: 153个（排除17个126参考通道）
- **MASCI结果**: `E:\MASCIresult\QuantificationResults01` - `QuantificationResults17`
- **FASTA文件**: `swiss_prot_human_20190214_target_conts_96Libraries.fasta`

## 流程步骤

### Step 0: 合并鉴定结果
```bash
Rscript merge_identification_results.R
```
- 输入: 各批次的 PhosSight/PhosphoRS 鉴定结果
- 输出: `PhosSightResults_TMT01_17.txt`, `PhosphoRSResults_TMT01_17.txt`

### Step 1: GetIonIntensity
```bash
Rscript step1_GetIonIntensity_PhosSight_17plex.R
Rscript step1_GetIonIntensity_PhosphoRS_17plex.R
```
- 读取 MASCI 定量结果（ReporterIons + SICstats）
- 与鉴定结果匹配
- 输出: `*_Intensity_17plex.txt`（170个TMT通道）

### Step 2: SiteQuant（完整版）
```bash
python step2_sitequant_full_17plex.py -i PhosSight_Intensity_17plex.txt -f <FASTA> -o PhosSight_17plex_full_site_table.tsv
```
- 按批次拆分数据
- 对每个批次调用原版 sitequant（preprocess + runner + modisite_quant）
- 使用 FASTA 进行精确位点定位
- 输出: `*_17plex_full_site_table.tsv`

### Step 3-5: 后处理
```bash
python step3_AdjustSiteTable_17plex.py -i ... -o ...
python step4_MergeSiteTables_17plex.py -i ... -o ...
python step5_GetUniprotIDGeneName_17plex.py -i ... -o ...
```

### Step 6: 缺失值分析
```bash
Rscript analysis/Compare_PhosSight_vs_PhosphoRS_MissingValue_17plex.R
Rscript analysis/Plot_PhosSight_vs_PhosphoRS_17plex.R
```

## 一键运行

```bash
cd C:\Users\wk\Desktop\TMTQuantification\full
.\run_complete_workflow_17plex.bat
```

## 输出文件

| 步骤 | 文件 | 说明 |
|------|------|------|
| Step 0 | `*Results_TMT01_17.txt` | 合并后的鉴定结果 |
| Step 1 | `*_Intensity_17plex.txt` | PSM级定量数据（170通道） |
| Step 2 | `*_17plex_full_site_table.tsv` | 完整版位点表（精确位点定位） |
| Step 3 | `*_site_table_Adjust_17plex.tsv` | 调整后的位点表 |
| Step 4 | `*_Merged_17plex.tsv` | 合并后的位点表 |
| Step 5 | `*_UniprotID_GeneName_17plex.tsv` | 带基因名的位点表 |
| Step 6 | `analysis/MissingValueCutoff/*.png` | 缺失值分析图 |

## 注意事项

1. **运行前确保**:
   - MASCI 定量结果已生成（`E:\MASCIresult\QuantificationResults01-17`）
   - 各批次鉴定结果文件存在
   - FASTA 文件路径正确

2. **环境要求**:
   - R 环境: `R_env`（包含 tidyverse, data.table, ggplot2）
   - Python 环境: `base`（包含 pandas, numpy）

3. **运行时间**:
   - 完整版 sitequant 比简化版慢很多（需要 FASTA 映射）
   - 17 个批次预计需要数小时

4. **如果某个批次失败**:
   - 可以单独运行该批次
   - 使用 `--plexes` 参数指定批次范围，如 `--plexes 1-5` 或 `--plexes 1,2,6,7`












