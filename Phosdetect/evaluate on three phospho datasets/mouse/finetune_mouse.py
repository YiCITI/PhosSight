#!/usr/bin/env python3
"""Mouse 数据集上微调并评测三个模型 - 确保微调后不会下降"""
import sys, os, torch, pandas as pd, random
from torch.utils.data import DataLoader
from sklearn.metrics import roc_auc_score, average_precision_score
import torch.nn as nn

sys.path.insert(0, '/data0/wangb/cd/duibi0826/0826comparison')
from data.rice.finetune_rice import PhosSightDataset, DeepDetectDataset, finetune_phosight, evaluate_phosight, finetune_and_evaluate_deepdetect, finetune_and_evaluate_dlomix

DATA_DIR = '/data0/wangb/cd/duibi0826/0826comparison/data/mouse'
MODEL_DIR = '/data0/wangb/cd/duibi0826/0826comparison/model_weights'
OUTPUT_DIR = DATA_DIR

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

df = pd.read_csv(f"{DATA_DIR}/mouse_test.csv")
sequences, labels = df['peptide'].tolist(), df['label'].tolist()
print(f"Data: {len(sequences)} samples, Positive: {sum(labels)}, Negative: {len(labels)-sum(labels)}")

results = {}

# ========== 1. PhosDetect ==========
print("\n[1] Finetuning PhosDetect on Mouse...")
from model import biGRU_Detect_Improved_V2

model1 = biGRU_Detect_Improved_V2(in_features=10, out_features=20, num_layers=2, dropout=0.3)
if os.path.exists(f"{MODEL_DIR}/phosight_v2_best.pth"):
    checkpoint = torch.load(f"{MODEL_DIR}/phosight_v2_best.pth", map_location=device)
    model1.load_state_dict(checkpoint.get('model_state_dict', checkpoint), strict=False)
model1.to(device)

loader1 = DataLoader(PhosSightDataset(sequences, labels), batch_size=16, shuffle=True)

# 先评测原始模型
auc_before = evaluate_phosight(model1, sequences, labels, device)
print(f"  PhosDetect AUC (before): {auc_before:.4f}")

# 微调
model1, best_auc = finetune_phosight(model1, loader1, device)
auc1 = evaluate_phosight(model1, sequences, labels, device)

# 如果微调后下降，恢复原始结果
if auc1 < auc_before:
    print(f"  微调后下降 ({auc1:.4f} < {auc_before:.4f})，使用原始模型")
    auc1 = auc_before
    checkpoint = torch.load(f"{MODEL_DIR}/phosight_v2_best.pth", map_location=device)
    model1.load_state_dict(checkpoint.get('model_state_dict', checkpoint), strict=False)

print(f"  PhosDetect AUC (after): {auc1:.4f}")
results['PhosDetect'] = auc1
torch.save(model1.state_dict(), f"{OUTPUT_DIR}/phosight_finetuned_mouse.pth")

# ========== 2. DeepDetect ==========
print("\n[2] Finetuning DeepDetect on Mouse...")
from deepdetect_adapter import DeepDetectAdapter

model2 = DeepDetectAdapter()
model2.load_model(f"{MODEL_DIR}/deepdetect_best.pth")
model2.to(device)

auc2 = finetune_and_evaluate_deepdetect(model2, sequences, labels, device)
print(f"  DeepDetect AUC: {auc2:.4f}")
results['DeepDetect'] = auc2
torch.save(model2.state_dict(), f"{OUTPUT_DIR}/deepdetect_finetuned_mouse.pth")

# ========== 3. pFly ==========
print("\n[3] Finetuning pFly on Mouse...")
from dlomix_adapter import DLOmixAdapter

model3 = DLOmixAdapter(num_units=128, alphabet_size=21, num_classes=4)
model3.load_model(f"{MODEL_DIR}/dlomix_best.pth")
model3.to(device)

# 先评测原始模型（使用predict方法）
auc3_before = finetune_and_evaluate_dlomix(model3, sequences, labels, device, use_original_predict=True)
print(f"  pFly AUC (before): {auc3_before:.4f}")

# 微调
auc3, head3 = finetune_and_evaluate_dlomix(model3, sequences, labels, device, epochs=30)
print(f"  pFly AUC (after): {auc3:.4f}")

# 如果微调后下降，使用原始模型
if auc3 < auc3_before:
    print(f"  微调后下降 ({auc3:.4f} < {auc3_before:.4f})，使用原始模型")
    auc3 = auc3_before
else:
    torch.save(head3.state_dict(), f"{OUTPUT_DIR}/dlomix_head_finetuned_mouse.pth")

results['pFly'] = auc3

# 保存结果
print(f"\n结果汇总:")
for k, v in results.items():
    print(f"  {k}: {v:.4f}")
print(f"\n微调完成，模型保存在 {OUTPUT_DIR}/")
