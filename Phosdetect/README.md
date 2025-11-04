

├── model.py              # 模型定义（包含V2改进）
├── train.py              # 训练脚本
├── test.py               # 测试脚本
├── run_training.py       # 便捷训练脚本
├── run_testing.py        # 便捷测试脚本
├── README.md             # 本文档
├── logs/                 # 训练日志目录
├── models/               # 模型保存目录
```


训练模型
```bash
# 使用便捷脚本
python run_training.py --data_path ../data/train/balanced_dataset_1.csv

# 或直接使用训练脚本
python train.py -p ../data/train/balanced_dataset_1.csv \
                -m ./models/best_model.pth \
                --model_type bigru_improved_v2 \
                -b 128 -e 100 -lr 0.0005 -d 0.3
```

测试模型
```bash
# 使用便捷脚本
python run_testing.py --model_path ./models/best_model.pth

# 或直接使用测试脚本
python test.py -p ../data/test/DeepDetect_human.csv \
               -m ./models/best_model.pth \
               --model_type bigru_improved_v2
```

## 参数说明

### 训练参数
- `--data_path`: 训练数据路径
- `--batch_size`: 批次大小（默认128）
- `--epochs`: 训练轮数（默认100）
- `--learning_rate`: 学习率（默认0.0005）
- `--dropout`: Dropout率（默认0.3）
- `--patience`: 早停耐心值（默认15）
- `--optimizer`: 优化器类型（adam/adamw/sgd）
- `--scheduler`: 学习率调度器（plateau/cosine/step）

### 模型参数
- `--in_features`: 输入特征维度（默认10）
- `--out_features`: 输出特征维度（默认20）
- `--num_layers`: GRU层数（默认2）

## 物理化学性质特征

### 氨基酸性质表
| 氨基酸 | 疏水性 | 电荷 | 极性 |
|--------|--------|------|------|
| A (Ala) | 1.8 | 0.0 | 0.0 |
| C (Cys) | 2.5 | 0.0 | 0.0 |
| D (Asp) | -3.5 | -1.0 | 1.0 |
| E (Glu) | -3.5 | -1.0 | 1.0 |
| F (Phe) | 2.8 | 0.0 | 0.0 |
| G (Gly) | -0.4 | 0.0 | 0.0 |
| H (His) | -3.2 | 0.1 | 1.0 |
| I (Ile) | 4.5 | 0.0 | 0.0 |
| K (Lys) | -3.9 | 1.0 | 1.0 |
| L (Leu) | 3.8 | 0.0 | 0.0 |
| M (Met) | 1.9 | 0.0 | 0.0 |
| N (Asn) | -3.5 | 0.0 | 1.0 |
| P (Pro) | -1.6 | 0.0 | 0.0 |
| Q (Gln) | -3.5 | 0.0 | 1.0 |
| R (Arg) | -4.5 | 1.0 | 1.0 |
| S (Ser) | -0.8 | 0.0 | 1.0 |
| T (Thr) | -0.7 | 0.0 | 1.0 |
| V (Val) | 4.2 | 0.0 | 0.0 |
| W (Trp) | -0.9 | 0.0 | 0.0 |
| Y (Tyr) | -1.3 | 0.0 | 1.0 |
| s (pSer) | -0.8 | -1.0 | 1.0 |
| t (pThr) | -0.7 | -1.0 | 1.0 |
| y (pTyr) | -1.3 | -1.0 | 1.0 |
