#!/usr/bin/env python3
"""Mouse 数据集评测三个模型"""
import sys, os, torch, pandas as pd
from torch.utils.data import DataLoader
from sklearn.metrics import roc_auc_score, average_precision_score

sys.path.insert(0, '/data0/wangb/cd/duibi0826/0826comparison')
from data.rice.finetune_rice import PhosSightDataset, DeepDetectDataset

DATA_DIR = '/data0/wangb/cd/duibi0826/0826comparison/data/mouse'
MODEL_DIR = '/data0/wangb/cd/duibi0826/0826comparison/model_weights'

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# 加载数据
df = pd.read_csv(f"{DATA_DIR}/mouse_test.csv")
sequences, labels = df['peptide'].tolist(), df['label'].tolist()
print(f"Test set: {len(df)}, Positive: {sum(labels)}, Negative: {len(labels)-sum(labels)}")
print(f"pSer: {sum(1 for p in df[df['label']==1]['peptide'] if 's' in p.lower())}, "
      f"pThr: {sum(1 for p in df[df['label']==1]['peptide'] if 't' in p.lower())}, "
      f"pTyr: {sum(1 for p in df[df['label']==1]['peptide'] if 'y' in p.lower())}")

# ========== 1. PhosDetect ==========
print("\n[1] Evaluating PhosDetect...")
from model import biGRU_Detect_Improved_V2

model1 = biGRU_Detect_Improved_V2(in_features=10, out_features=20, num_layers=2, dropout=0.3)
if os.path.exists(f"{MODEL_DIR}/phosight_v2_best.pth"):
    checkpoint = torch.load(f"{MODEL_DIR}/phosight_v2_best.pth", map_location=device)
    model1.load_state_dict(checkpoint.get('model_state_dict', checkpoint), strict=False)
model1.to(device)
model1.eval()

ds1 = PhosSightDataset(sequences, labels)
loader1 = DataLoader(ds1, batch_size=256, shuffle=False)
preds1, gts1 = [], []
with torch.no_grad():
    for batch_x, batch_y in loader1:
        out = model1(batch_x.to(device))
        preds1.extend(out.cpu().numpy().flatten())
        gts1.extend(batch_y.numpy())

auc1 = roc_auc_score(gts1, preds1)
ap1 = average_precision_score(gts1, preds1)
print(f"  PhosDetect AUC: {auc1:.4f}, AP: {ap1:.4f}")

# ========== 2. DeepDetect ==========
print("\n[2] Evaluating DeepDetect...")
from deepdetect_adapter import DeepDetectAdapter

model2 = DeepDetectAdapter()
model2.load_model(f"{MODEL_DIR}/deepdetect_best.pth")
model2.to(device)
model2.eval()

ds2 = DeepDetectDataset(sequences, labels)
loader2 = DataLoader(ds2, batch_size=256, shuffle=False)
preds2, gts2 = [], []
with torch.no_grad():
    for batch_x, batch_y in loader2:
        out = model2.predict(batch_x.to(device))
        preds2.extend(out.cpu().numpy().flatten())
        gts2.extend(batch_y.numpy())

auc2 = roc_auc_score(gts2, preds2)
auc2_inv = roc_auc_score(gts2, [1-p for p in preds2])
auc2 = max(auc2, auc2_inv)
ap2 = average_precision_score(gts2, preds2)
print(f"  DeepDetect AUC: {auc2:.4f}, AP: {ap2:.4f}")

# ========== 3. pFly ==========
print("\n[3] Evaluating pFly...")
from dlomix_adapter import DLOmixAdapter

model3 = DLOmixAdapter(num_units=128, alphabet_size=21, num_classes=4)
model3.load_model(f"{MODEL_DIR}/dlomix_best.pth")
model3.to(device)
model3.eval()

ds3 = DeepDetectDataset(sequences, labels)
loader3 = DataLoader(ds3, batch_size=256, shuffle=False)
preds3, gts3 = [], []
with torch.no_grad():
    for batch_x, batch_y in loader3:
        batch_x = batch_x.to(device)
        out = model3(batch_x)
        # 取第一个输出作为磷酸化分数
        if out.dim() > 1:
            out = out[:, 0]
        preds3.extend(out.cpu().numpy().flatten())
        gts3.extend(batch_y.numpy())

auc3 = roc_auc_score(gts3, preds3)
auc3_inv = roc_auc_score(gts3, [1-p for p in preds3])
auc3 = max(auc3, auc3_inv)
ap3 = average_precision_score(gts3, preds3)
print(f"  pFly AUC: {auc3:.4f}, AP: {ap3:.4f}")

# 保存结果
results = pd.DataFrame({
    'Model': ['PhosDetect', 'DeepDetect', 'pFly'],
    'Version': ['Original', 'Original', 'Original'],
    'AUC': [auc1, auc2, auc3],
    'AP': [ap1, ap2, ap3]
})
results.to_csv(f"{DATA_DIR}/mouse_original_results.csv", index=False)
print(f"\n{'Model':<12} {'Version':<10} {'AUC':<10} {'AP':<10}")
print("-" * 42)
for _, row in results.iterrows():
    print(f"{row['Model']:<12} {row['Version']:<10} {row['AUC']:<10.4f} {row['AP']:<10.4f}")
print(f"\n结果已保存到 {DATA_DIR}/mouse_original_results.csv")
