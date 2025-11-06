# PhosDetect 肽段可检测性预测脚本使用说明

## 脚本功能

`predict_detectability.sh` 是一个用于从FASTA文件读取肽段序列，使用PhosSight模型预测detectability scores的Shell脚本。

## 使用方法

### 基本用法

```bash
# 使用默认输入文件
./predict_detectability.sh

# 指定输入FASTA文件
./predict_detectability.sh /path/to/peptides.fasta

# 指定输入和输出文件
./predict_detectability.sh /path/to/peptides.fasta /path/to/output.txt
```

### 查看帮助

```bash
./predict_detectability.sh --help
```

## 工作流程

1. **读取FASTA文件** - 从指定的FASTA文件中读取肽段序列
2. **加载模型** - 加载预训练的PhosSight模型
3. **预测** - 对每个肽段序列进行可检测性预测
4. **保存结果** - 将预测结果保存为TXT文件（格式：sequence\tdetectability_score）

## 输出格式

输出文件为制表符分隔的文本文件，包含以下列：

```
sequence    detectability_score
AAAKKK      0.823456
SSSTTT      0.745321
...
```

## 配置参数

在使用脚本前，请根据实际情况修改脚本中的以下配置参数：

- `PHOSSIGHT_MODEL_PATH`: PhosSight模型文件路径
- `PHOSSIGHT_CODE_PATH`: PhosSight模型代码路径
- `PYTHON_SCRIPT`: Python预测脚本路径
- `DEFAULT_INPUT_FASTA`: 默认输入FASTA文件路径

## 依赖要求

- Python 3.x
- PyTorch
- NumPy
- PhosSight模型文件

## 示例

```bash
# 预测示例FASTA文件
./predict_detectability.sh \
    /data0/wangb/wbscy/czy/0925/human_all_peptides_2_7_46.fasta \
    /data0/wangb/wbscy/czy/0925/human_peptides_detectability_scores.txt
```

## 注意事项

1. 确保模型文件路径正确
2. 确保Python环境已安装所需依赖
3. 输入FASTA文件格式应符合标准FASTA格式
4. 脚本会自动检查CUDA可用性，如果可用将使用GPU加速

## 错误处理

脚本会检查：
- 输入文件是否存在
- 模型文件是否存在
- Python和PyTorch是否已安装
- 预测过程是否成功

如果出现错误，脚本会显示详细的错误信息并退出。

