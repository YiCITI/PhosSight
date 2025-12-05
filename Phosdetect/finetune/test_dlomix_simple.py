#!/usr/bin/env python3
"""
简洁的DLOmix模型测试脚本
"""

import os
import sys
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from tqdm import tqdm

# 添加项目路径
sys.path.append('/data0/wangb/cd/duibi0826')
from model_adapters.dlomix_adapter import DLOmixAdapter

class PeptideDataset(Dataset):
    """肽段数据集 - 使用与训练时相同的预处理"""
    
    def __init__(self, sequences, labels, max_length=30):
        self.sequences = sequences
        self.labels = labels
        self.max_length = max_length
        
        # DLOmix词汇表：21维
        self.vocab = {
            '0': 0,  # padding
            'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
            'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
            'V': 18, 'W': 19, 'Y': 20
        }
    
    def __len__(self):
        return len(self.sequences)
    
    def encode_sequence(self, sequence):
        """编码氨基酸序列"""
        # 处理磷酸化氨基酸
        sequence = sequence.replace('S[UNIMOD:21]', 'S')
        sequence = sequence.replace('T[UNIMOD:21]', 'T')
        sequence = sequence.replace('Y[UNIMOD:21]', 'Y')
        
        # 编码序列
        encoded = []
        for aa in sequence:
            if aa in self.vocab:
                encoded.append(self.vocab[aa])
            else:
                encoded.append(self.vocab['0'])  # 未知氨基酸用padding代替
        
        # 填充到固定长度
        if len(encoded) < self.max_length:
            encoded.extend([self.vocab['0']] * (self.max_length - len(encoded)))
        else:
            encoded = encoded[:self.max_length]
        
        return torch.tensor(encoded, dtype=torch.long)
    
    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        label = self.labels[idx]
        
        encoded_seq = self.encode_sequence(sequence)
        return encoded_seq, torch.tensor(label, dtype=torch.float32)

def test_dlomix_on_datasets():
    """在七个测试集上测试DLOmix模型"""
    
    # 设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # 模型配置（与训练时相同）
    config = {
        'num_units': 128,
        'alphabet_size': 21,
        'num_classes': 4
    }
    
    # 加载模型
    model = DLOmixAdapter(**config)
    model.load_model('/data0/wangb/cd/duibi0826/0826comparison/model_weights/dlomix_best.pth')
    model.to(device)
    model.eval()
    print("Model loaded successfully")
    
    # 测试集列表
    test_sets = [
        'DeepDetect_ecoli.csv',
        'DeepDetect_human.csv', 
        'DeepDetect_mouse.csv',
        'DeepDetect_yeast.csv',
        'DeepRescore2_HCC.csv',
        'DeepRescore2_label_free.csv',
        'DeepRescore2_UCEC.csv'
    ]
    
    test_data_dir = '/data0/wangb/wbscy/PhosSight/data/test'
    
    results = []
    
    for test_file in test_sets:
        print(f"\nTesting on: {test_file}")
        
        # 加载数据
        file_path = os.path.join(test_data_dir, test_file)
        df = pd.read_csv(file_path, header=None, names=['sequence', 'label'])
        
        sequences = df['sequence'].tolist()
        labels = df['label'].tolist()
        
        print(f"  Loaded {len(sequences)} samples")
        print(f"  Positive samples: {sum(labels)}, Negative samples: {len(labels) - sum(labels)}")
        
        # 创建数据集和dataloader
        dataset = PeptideDataset(sequences, labels)
        dataloader = DataLoader(dataset, batch_size=128, shuffle=False)
        
        # 评估
        all_predictions = []
        all_labels = []
        
        with torch.no_grad():
            for batch_sequences, batch_labels in tqdm(dataloader, desc=f"Evaluating {test_file}"):
                batch_sequences = batch_sequences.to(device)
                batch_labels = batch_labels.to(device)
                
                # 获取预测
                predictions = model.predict(batch_sequences)
                predictions = predictions.cpu().numpy().flatten()
                
                all_predictions.extend(predictions)
                all_labels.extend(batch_labels.cpu().numpy().flatten())
        
        # 计算AUC
        auc = roc_auc_score(all_labels, all_predictions)
        print(f"  AUC: {auc:.4f}")
        
        results.append({
            'test_set': test_file,
            'auc': auc,
            'num_samples': len(all_labels),
            'positive_samples': sum(all_labels)
        })
    
    # 打印汇总结果
    print("\n" + "="*60)
    print("SUMMARY RESULTS")
    print("="*60)
    print(f"{'Test Set':<25} {'AUC':<10} {'Samples':<10} {'Positive':<10}")
    print("-" * 60)
    
    for result in results:
        print(f"{result['test_set']:<25} {result['auc']:<10.4f} {result['num_samples']:<10} {result['positive_samples']:<10}")
    
    # 计算平均AUC
    avg_auc = np.mean([r['auc'] for r in results])
    print("-" * 60)
    print(f"{'Average AUC':<25} {avg_auc:<10.4f}")
    print("="*60)
    
    return results

if __name__ == "__main__":
    test_dlomix_on_datasets() 