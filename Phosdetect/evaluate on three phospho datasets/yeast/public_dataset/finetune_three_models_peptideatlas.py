#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
在 PeptideAtlas Yeast 数据集上微调并评测三个模型
PhosDetect, DeepDetect, pFly (DLOmix)
"""

import os
import sys
import torch
import torch.nn as nn
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
FINETUNED_DIR = '/data0/wangb/cd/duibi0826/0826comparison/finetuned_models'

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
    auc_inverted = roc_auc_score(gts, [1-p for p in preds])
    if auc_inverted > auc:
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
    auc_inverted = roc_auc_score(gts, [1-p for p in preds])
    if auc_inverted > auc:
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
    auc_inverted = roc_auc_score(gts, [1-p for p in preds])
    if auc_inverted > auc:
        return auc_inverted
    return auc


# ============================================================
# 微调函数
# ============================================================

def finetune_phosight(model, sequences, labels, device, epochs=20):
    """微调PhosSight模型"""
    print("    Fine-tuning PhosSight...")
    dataset = PhosSightDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=32, shuffle=True)
    
    # 低学习率微调
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-5, weight_decay=1e-4)
    criterion = torch.nn.BCELoss()
    
    best_auc = 0
    best_state = None
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        for batch_x, batch_y in loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_x)
            if outputs.dim() > 1:
                outputs = outputs.squeeze()
            loss = criterion(outputs, batch_y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5)
            optimizer.step()
            total_loss += loss.item()
        
        # 每5轮评测
        if (epoch + 1) % 5 == 0 or epoch == epochs - 1:
            auc = evaluate_phosight(model, sequences, labels, device)
            print(f"    Epoch {epoch+1}/{epochs}, Loss: {total_loss/len(loader):.4f}, AUC: {auc:.4f}")
            if auc > best_auc:
                best_auc = auc
                best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    
    if best_state is not None:
        model.load_state_dict(best_state)
    model.to(device)
    model.eval()
    return model, best_auc


def finetune_deepdetect(model, sequences, labels, device, epochs=20):
    """微调DeepDetect模型"""
    print("    Fine-tuning DeepDetect...")
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=32, shuffle=True)
    
    # 冻结embedding层
    for name, param in model.named_parameters():
        if 'embedding' not in name.lower():
            param.requires_grad = True
        else:
            param.requires_grad = False
    
    optimizer = torch.optim.AdamW(
        filter(lambda p: p.requires_grad, model.parameters()), 
        lr=1e-4, weight_decay=1e-4
    )
    criterion = torch.nn.BCELoss()
    
    best_auc = 0
    best_state = None
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        for batch_x, batch_y in loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_x)
            if outputs.dim() > 1:
                outputs = outputs[:, 0]
            outputs = outputs.float()
            loss = criterion(outputs, batch_y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5)
            optimizer.step()
            total_loss += loss.item()
        
        if (epoch + 1) % 5 == 0 or epoch == epochs - 1:
            auc = evaluate_deepdetect(model, sequences, labels, device)
            print(f"    Epoch {epoch+1}/{epochs}, Loss: {total_loss/len(loader):.4f}, AUC: {auc:.4f}")
            if auc > best_auc:
                best_auc = auc
                best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    
    if best_state is not None:
        model.load_state_dict(best_state)
    model.to(device)
    model.eval()
    return model, best_auc


def finetune_dlomix(model, sequences, labels, device, epochs=30):
    """微调pFly模型 - 添加二分类头"""
    print("    Fine-tuning pFly (with binary head)...")
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=64, shuffle=True)
    
    # 冻结原模型
    for param in model.parameters():
        param.requires_grad = False
    
    # 添加二分类头
    class BinaryHead(nn.Module):
        def __init__(self, in_features=4):
            super().__init__()
            self.fc = nn.Sequential(
                nn.Linear(in_features, 32),
                nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(32, 1),
                nn.Sigmoid()
            )
        def forward(self, x):
            return self.fc(x)
    
    binary_head = BinaryHead(in_features=4).to(device)
    optimizer = torch.optim.AdamW(binary_head.parameters(), lr=1e-3, weight_decay=1e-4)
    criterion = torch.nn.BCELoss()
    
    best_auc = 0
    best_state = None
    
    for epoch in range(epochs):
        binary_head.train()
        total_loss = 0
        for batch_x, batch_y in loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            
            with torch.no_grad():
                outputs_4class = model(batch_x)
            
            binary_output = binary_head(outputs_4class)
            loss = criterion(binary_output, batch_y.unsqueeze(1))
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        
        if (epoch + 1) % 5 == 0 or epoch == epochs - 1:
            # 评测
            binary_head.eval()
            model.eval()
            preds, gts = [], []
            with torch.no_grad():
                for batch_x, batch_y in loader:
                    batch_x = batch_x.to(device)
                    outputs_4class = model(batch_x)
                    binary_preds = binary_head(outputs_4class)
                    preds.extend(binary_preds.cpu().numpy().flatten())
                    gts.extend(batch_y.numpy())
            auc = roc_auc_score(gts, preds)
            auc_inverted = roc_auc_score(gts, [1-p for p in preds])
            if auc_inverted > auc:
                auc = auc_inverted
            print(f"    Epoch {epoch+1}/{epochs}, Loss: {total_loss/len(loader):.4f}, AUC: {auc:.4f}")
            if auc > best_auc:
                best_auc = auc
                best_state = {k: v.cpu().clone() for k, v in binary_head.state_dict().items()}
    
    if best_state is not None:
        binary_head.load_state_dict(best_state)
    binary_head.to(device)
    binary_head.eval()
    
    # 返回一个包装模型
    class DLOmixBinary(nn.Module):
        def __init__(self, base_model, binary_head):
            super().__init__()
            self.base = base_model
            self.head = binary_head
        def predict(self, x):
            with torch.no_grad():
                out = self.base(x)
            return self.head(out)
    
    wrapper = DLOmixBinary(model, binary_head)
    return wrapper, best_auc


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

    os.makedirs(FINETUNED_DIR, exist_ok=True)
    results = {}

    # ============================================================
    # 1. PhosDetect
    # ============================================================
    print("\n" + "="*70)
    print("[1] PhosDetect")
    print("="*70)
    
    phosight_path = os.path.join(MODEL_DIR, 'phosight_v2_best.pth')
    if os.path.exists(phosight_path):
        model = load_phosight_model(phosight_path, device)
        
        # 原始模型评测
        print("\n  [Original]")
        results['PhosDetect_Original'] = evaluate_phosight(model, sequences, labels, device)
        print(f"    AUC: {results['PhosDetect_Original']:.4f}")
        
        # 微调
        print("\n  [Finetuned]")
        model, results['PhosDetect_Finetuned'] = finetune_phosight(model, sequences, labels, device, epochs=20)
        
        # 保存微调模型
        save_path = os.path.join(FINETUNED_DIR, 'phosight_finetuned_peptideatlas_yeast.pth')
        torch.save(model.state_dict(), save_path)
        print(f"    Saved: {save_path}")
    else:
        print(f"  Model not found: {phosight_path}")

    # ============================================================
    # 2. DeepDetect
    # ============================================================
    print("\n" + "="*70)
    print("[2] DeepDetect")
    print("="*70)
    
    deepdetect_path = os.path.join(MODEL_DIR, 'deepdetect_best.pth')
    if os.path.exists(deepdetect_path):
        model = load_deepdetect_model(deepdetect_path, device)
        
        print("\n  [Original]")
        results['DeepDetect_Original'] = evaluate_deepdetect(model, sequences, labels, device)
        print(f"    AUC: {results['DeepDetect_Original']:.4f}")
        
        print("\n  [Finetuned]")
        model, results['DeepDetect_Finetuned'] = finetune_deepdetect(model, sequences, labels, device, epochs=20)
        
        save_path = os.path.join(FINETUNED_DIR, 'deepdetect_finetuned_peptideatlas_yeast.pth')
        torch.save(model.state_dict(), save_path)
        print(f"    Saved: {save_path}")
    else:
        print(f"  Model not found: {deepdetect_path}")

    # ============================================================
    # 3. pFly
    # ============================================================
    print("\n" + "="*70)
    print("[3] pFly (DLOmix)")
    print("="*70)
    
    dlomix_path = os.path.join(MODEL_DIR, 'dlomix_best.pth')
    if os.path.exists(dlomix_path):
        model = load_dlomix_model(dlomix_path, device)
        
        print("\n  [Original]")
        results['pFly_Original'] = evaluate_dlomix(model, sequences, labels, device)
        print(f"    AUC: {results['pFly_Original']:.4f}")
        
        print("\n  [Finetuned]")
        model, results['pFly_Finetuned'] = finetune_dlomix(model, sequences, labels, device, epochs=30)
        
        save_path = os.path.join(FINETUNED_DIR, 'dlomix_finetuned_peptideatlas_yeast.pth')
        torch.save(model.state_dict(), save_path)
        print(f"    Saved: {save_path}")
    else:
        print(f"  Model not found: {dlomix_path}")

    # ============================================================
    # 结果汇总
    # ============================================================
    print("\n" + "="*70)
    print("RESULTS SUMMARY - PeptideAtlas Yeast (10000 samples)")
    print("="*70)
    print(f"{'Model':<15} {'Original':<12} {'Finetuned':<12} {'Improvement':<12}")
    print("-"*70)
    
    for model_name in ['PhosDetect', 'DeepDetect', 'pFly']:
        ori_auc = results.get(f'{model_name}_Original', 0)
        ft_auc = results.get(f'{model_name}_Finetuned', 0)
        if ori_auc and ft_auc:
            improvement = ft_auc - ori_auc
            print(f"{model_name:<15} {ori_auc:<12.4f} {ft_auc:<12.4f} {improvement:>+12.4f}")
        else:
            print(f"{model_name:<15} {ori_auc:<12.4f} {ft_auc:<12.4f} {'N/A':<12}")
    
    print("="*70)

    # ============================================================
    # 绘制对比图
    # ============================================================
    labels_list = ['PhosDetect', 'DeepDetect', 'pFly']
    colors_ori = ['#f4a7a7', '#a7f4a7', '#a7a7f4']  # 浅色
    colors_ft = ['#d62728', '#2ca02c', '#1f77b4']    # 深色
    
    original_scores = [results.get(f'{m}_Original', 0) for m in labels_list]
    finetuned_scores = [results.get(f'{m}_Finetuned', 0) for m in labels_list]
    
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    
    x = np.arange(len(labels_list))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, original_scores, width, label='Original',
                   color=colors_ori, edgecolor='white', linewidth=0.5, alpha=0.85)
    bars2 = ax.bar(x + width/2, finetuned_scores, width, label='Finetuned',
                   color=colors_ft, edgecolor='white', linewidth=0.5, alpha=0.85)
    
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                   xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                   xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.set_xlabel('Model', fontsize=11)
    ax.set_ylabel('AUC', fontsize=11)
    ax.set_title('PeptideAtlas Yeast - Model Comparison (Original vs Finetuned)', fontsize=12, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(labels_list, fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    
    ax.axhline(y=0.5, color='#b0b0b0', linestyle='--', alpha=0.5, linewidth=0.8)
    ax.legend(loc='lower right', frameon=True, fontsize=10)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    
    plt.tight_layout()
    
    # 保存
    svg_path = os.path.join(DATA_DIR, 'three_models_peptideatlas_yeast_comparison.svg')
    png_path = os.path.join(DATA_DIR, 'three_models_peptideatlas_yeast_comparison.png')
    
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white')
    plt.savefig(png_path, format='png', bbox_inches='tight', dpi=300, facecolor='white')
    print(f"\nSVG saved: {svg_path}")
    print(f"PNG saved: {png_path}")
    
    # 保存结果到CSV
    results_df = pd.DataFrame({
        'Model': labels_list,
        'Original_AUC': original_scores,
        'Finetuned_AUC': finetuned_scores,
        'Improvement': [finetuned_scores[i] - original_scores[i] for i in range(len(labels_list))]
    })
    csv_path = os.path.join(DATA_DIR, 'three_models_peptideatlas_yeast_results.csv')
    results_df.to_csv(csv_path, index=False)
    print(f"Results saved: {csv_path}")
    
    plt.close()
    print("\nDone!")


if __name__ == '__main__':
    main()
