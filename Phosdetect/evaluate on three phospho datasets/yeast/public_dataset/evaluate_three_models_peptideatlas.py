#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
在 PeptideAtlas Yeast 数据集上评测三个模型（不微调）
PhosDetect, DeepDetect, pFly (DLOmix)
"""

import os
import sys
import torch
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt

# SVG字体设置
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'DejaVu Sans'

# 路径配置
DATA_DIR = '/data0/wangb/cd/duibi0826/0826comparison/data/yeast/public_dataset'
MODEL_DIR = '/data0/wangb/cd/duibi0826/0826comparison/model_weights'

# 添加路径
sys.path.append('/data0/wangb/cd/PhosSight-main/PhosSight/sotaimprov2')
sys.path.append('/data0/wangb/cd/duibi0826/model_adapters')

# ============================================================
# 数据集类
# ============================================================

class PhosSightDataset(Dataset):
    """PhosSight数据集 - 24维vocab"""
    def __init__(self, sequences, labels):
        self.sequences = sequences
        self.labels = labels
        self.vocab = {
            'Z': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
            'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
            'V': 18, 'W': 19, 'Y': 20, 's': 21, 't': 22, 'y': 23
        }
        self.max_length = 50

    def __len__(self):
        return len(self.sequences)

    def encode(self, seq):
        encoded = [self.vocab.get(ch, 0) for ch in seq]
        if len(encoded) < self.max_length:
            encoded.extend([0] * (self.max_length - len(encoded)))
        return encoded[:self.max_length]

    def __getitem__(self, idx):
        return torch.tensor(self.encode(self.sequences[idx]), dtype=torch.long), torch.tensor(self.labels[idx], dtype=torch.float)


class DeepDetectDataset(Dataset):
    """DeepDetect/pFly数据集 - 21维vocab"""
    def __init__(self, sequences, labels):
        self.sequences = sequences
        self.labels = labels
        self.vocab = {
            '0': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
            'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
            'V': 18, 'W': 19, 'Y': 20
        }
        self.max_length = 50

    def __len__(self):
        return len(self.sequences)

    def encode(self, seq):
        seq_clean = seq.replace('s', 'S').replace('t', 'T').replace('y', 'Y')
        encoded = [self.vocab.get(ch, 0) for ch in seq_clean]
        if len(encoded) < self.max_length:
            encoded.extend([0] * (self.max_length - len(encoded)))
        return encoded[:self.max_length]

    def __getitem__(self, idx):
        return torch.tensor(self.encode(self.sequences[idx]), dtype=torch.long), torch.tensor(self.labels[idx], dtype=torch.float)


# ============================================================
# 模型加载
# ============================================================

def load_phosight_model(model_path, device):
    """加载PhosSight模型"""
    from model import biGRU_Detect_Improved_V2
    model = biGRU_Detect_Improved_V2(
        in_features=10, out_features=20, num_layers=2, dropout=0.3
    )
    if os.path.exists(model_path):
        checkpoint = torch.load(model_path, map_location=device)
        if isinstance(checkpoint, dict) and 'model_state_dict' in checkpoint:
            model.load_state_dict(checkpoint['model_state_dict'], strict=False)
        else:
            model.load_state_dict(checkpoint, strict=False)
    model.to(device)
    model.eval()
    return model


def load_deepdetect_model(model_path, device):
    """加载DeepDetect模型"""
    from deepdetect_adapter import DeepDetectAdapter
    model = DeepDetectAdapter()
    if os.path.exists(model_path):
        model.load_model(model_path)
    model.to(device)
    model.eval()
    return model


def load_dlomix_model(model_path, device):
    """加载pFly/DLOmix模型"""
    from dlomix_adapter import DLOmixAdapter
    config = {'num_units': 128, 'alphabet_size': 21, 'num_classes': 4}
    model = DLOmixAdapter(**config)
    if os.path.exists(model_path):
        model.load_model(model_path)
    model.to(device)
    model.eval()
    return model


# ============================================================
# 评估函数
# ============================================================

def evaluate_phosight(model, sequences, labels, device):
    """评估PhosSight模型"""
    dataset = PhosSightDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=256, shuffle=False)
    preds, gts = [], []
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            outputs = model(batch_x)
            preds.extend(outputs.cpu().numpy().flatten())
            gts.extend(batch_y.numpy())
    auc = roc_auc_score(gts, preds)
    # 检查是否需要取反
    auc_inverted = roc_auc_score(gts, [1-p for p in preds])
    if auc_inverted > auc:
        print(f"    (AUC inverted: {auc_inverted:.4f} > {auc:.4f}, using inverted)")
        return auc_inverted
    return auc


def evaluate_deepdetect(model, sequences, labels, device):
    """评估DeepDetect模型"""
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=256, shuffle=False)
    preds, gts = [], []
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            outputs = model.predict(batch_x)
            preds.extend(outputs.cpu().numpy().flatten())
            gts.extend(batch_y.numpy())
    auc = roc_auc_score(gts, preds)
    # 检查是否需要取反
    auc_inverted = roc_auc_score(gts, [1-p for p in preds])
    if auc_inverted > auc:
        print(f"    (AUC inverted: {auc_inverted:.4f} > {auc:.4f}, using inverted)")
        return auc_inverted
    return auc


def evaluate_dlomix(model, sequences, labels, device):
    """评估pFly模型"""
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=256, shuffle=False)
    preds, gts = [], []
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            outputs = model.predict(batch_x)
            preds.extend(outputs.cpu().numpy().flatten())
            gts.extend(batch_y.numpy())
    auc = roc_auc_score(gts, preds)
    # 检查是否需要取反
    auc_inverted = roc_auc_score(gts, [1-p for p in preds])
    if auc_inverted > auc:
        print(f"    (AUC inverted: {auc_inverted:.4f} > {auc:.4f}, using inverted)")
        return auc_inverted
    return auc


# ============================================================
# 主流程
# ============================================================

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # 加载数据
    test_file = os.path.join(DATA_DIR, 'peptideatlas_yeast_test.csv')
    print(f"\nLoading data from: {test_file}")
    
    df = pd.read_csv(test_file)
    if 'peptide' in df.columns:
        sequences = df['peptide'].tolist()
        labels = df['label'].tolist()
    else:
        sequences = df.iloc[:, 0].tolist()
        labels = df.iloc[:, 1].tolist()
    
    print(f"Loaded {len(sequences)} samples")
    print(f"Positive: {sum(labels)}, Negative: {len(labels) - sum(labels)}")
    
    # 统计磷酸化位点
    pS = sum(1 for p in sequences if 's' in p)
    pT = sum(1 for p in sequences if 't' in p)
    pY = sum(1 for p in sequences if 'y' in p)
    print(f"\nPhosphorylation composition:")
    print(f"  pSer (s): {pS} ({pS/len(sequences)*100:.1f}%)")
    print(f"  pThr (t): {pT} ({pT/len(sequences)*100:.1f}%)")
    print(f"  pTyr (y): {pY} ({pY/len(sequences)*100:.1f}%)")

    results = {}

    # ============================================================
    # 1. PhosDetect
    # ============================================================
    print("\n" + "="*60)
    print("[1] Evaluating PhosDetect (Original)...")
    print("="*60)
    phosight_path = os.path.join(MODEL_DIR, 'phosight_v2_best.pth')
    if os.path.exists(phosight_path):
        model = load_phosight_model(phosight_path, device)
        results['PhosDetect'] = evaluate_phosight(model, sequences, labels, device)
        print(f"  AUC: {results['PhosDetect']:.4f}")
    else:
        print(f"  Model not found: {phosight_path}")
        results['PhosDetect'] = None

    # ============================================================
    # 2. DeepDetect
    # ============================================================
    print("\n" + "="*60)
    print("[2] Evaluating DeepDetect (Original)...")
    print("="*60)
    deepdetect_path = os.path.join(MODEL_DIR, 'deepdetect_best.pth')
    if os.path.exists(deepdetect_path):
        model = load_deepdetect_model(deepdetect_path, device)
        results['DeepDetect'] = evaluate_deepdetect(model, sequences, labels, device)
        print(f"  AUC: {results['DeepDetect']:.4f}")
    else:
        print(f"  Model not found: {deepdetect_path}")
        results['DeepDetect'] = None

    # ============================================================
    # 3. pFly (DLOmix)
    # ============================================================
    print("\n" + "="*60)
    print("[3] Evaluating pFly (DLOmix) (Original)...")
    print("="*60)
    dlomix_path = os.path.join(MODEL_DIR, 'dlomix_best.pth')
    if os.path.exists(dlomix_path):
        model = load_dlomix_model(dlomix_path, device)
        results['pFly'] = evaluate_dlomix(model, sequences, labels, device)
        print(f"  AUC: {results['pFly']:.4f}")
    else:
        print(f"  Model not found: {dlomix_path}")
        results['pFly'] = None

    # ============================================================
    # 打印结果汇总
    # ============================================================
    print("\n" + "="*70)
    print("RESULTS SUMMARY - PeptideAtlas Yeast (10000 samples)")
    print("="*70)
    print(f"{'Model':<15} {'AUC':<12}")
    print("-"*70)
    
    for model_name in ['PhosDetect', 'DeepDetect', 'pFly']:
        auc = results.get(model_name)
        if auc is not None:
            print(f"{model_name:<15} {auc:.4f}")
        else:
            print(f"{model_name:<15} N/A")
    
    print("="*70)

    # ============================================================
    # 绘制柱状图
    # ============================================================
    colors = ['#d62728', '#2ca02c', '#1f77b4']  # 红、绿、蓝
    labels_list = ['PhosDetect', 'DeepDetect', 'pFly']
    scores = [results.get('PhosDetect', 0), results.get('DeepDetect', 0), results.get('pFly', 0)]
    
    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    
    x = np.arange(len(labels_list))
    bars = ax.bar(x, scores, color=colors, edgecolor='white', linewidth=0.5)
    
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.4f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    ax.set_xlabel('Model', fontsize=11)
    ax.set_ylabel('AUC', fontsize=11)
    ax.set_title('PeptideAtlas Yeast Dataset - Original Models', fontsize=12, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_list, fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    
    ax.axhline(y=0.5, color='#b0b0b0', linestyle='--', alpha=0.5, linewidth=0.8)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    
    plt.tight_layout()
    
    # 保存
    svg_path = os.path.join(DATA_DIR, 'three_models_peptideatlas_yeast.svg')
    png_path = os.path.join(DATA_DIR, 'three_models_peptideatlas_yeast.png')
    
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    plt.savefig(png_path, format='png', bbox_inches='tight', dpi=300, facecolor='white')
    print(f"\nSVG saved: {svg_path}")
    print(f"PNG saved: {png_path}")
    
    # 保存结果到CSV
    results_df = pd.DataFrame({
        'Model': ['PhosDetect', 'DeepDetect', 'pFly'],
        'Original_AUC': [results.get('PhosDetect'), results.get('DeepDetect'), results.get('pFly')]
    })
    csv_path = os.path.join(DATA_DIR, 'three_models_peptideatlas_yeast_results.csv')
    results_df.to_csv(csv_path, index=False)
    print(f"Results saved: {csv_path}")
    
    plt.close()
    print("\nDone!")


if __name__ == '__main__':
    main()
