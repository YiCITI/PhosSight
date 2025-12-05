#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试微调后的DeepDetect模型在DeepRescore2_UCEC数据集上的性能
"""

import torch
import pandas as pd
import numpy as np
import os
import sys
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score
from torch.utils.data import DataLoader

# 添加路径
sys.path.append('/data0/wangb/cd/duibi0826/model_adapters')
sys.path.append('/data0/wangb/cd/duibi')
from deepdetect_adapter import DeepDetectAdapter
from data_preprocessor import DataPreprocessor, PeptideDataset

def test_finetuned_model():
    """测试微调后的模型"""
    
    # 模型和数据路径
    model_path = '/data0/wangb/cd/duibi0826/0826comparison/finetuned_models/deepdetect_finetuned_DeepRescore2_UCEC_auc0.9483.pth'
    test_data_path = '/data0/wangb/wbscy/PhosSight/data/test/DeepRescore2_UCEC.csv'
    
    # 设置设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    # 加载模型
    print(f"Loading model from: {model_path}")
    model = DeepDetectAdapter()
    
    if os.path.exists(model_path):
        model.load_model(model_path)
        model.to(device)
        model.eval()
        print("✓ Model loaded successfully")
    else:
        print(f"✗ Model file not found: {model_path}")
        return
    
    # 加载测试数据
    print(f"Loading test data from: {test_data_path}")
    try:
        preprocessor = DataPreprocessor(test_data_path, model_type='deepdetect')
        preprocessor.load_data()
        sequences, labels = preprocessor.preprocess_sequences()
        
        print(f"Loaded {len(sequences)} samples")
        print(f"Label distribution: {pd.Series(labels).value_counts().to_dict()}")
        
        # 创建数据集和数据加载器
        dataset = PeptideDataset(sequences, labels, model_type='deepdetect')
        dataloader = DataLoader(dataset, batch_size=128, shuffle=False)
        
        # 运行推理
        predictions = []
        true_labels = []
        
        print("Running inference...")
        with torch.no_grad():
            for batch_idx, batch in enumerate(dataloader):
                inputs = batch[0].to(device)
                outputs = model.predict(inputs)
                predictions.extend(outputs.flatten().cpu().numpy().tolist())
                true_labels.extend(batch[1].numpy().tolist())
        
        # 计算指标
        auc = roc_auc_score(true_labels, predictions)
        pred_labels = [1 if p > 0.5 else 0 for p in predictions]
        accuracy = accuracy_score(true_labels, pred_labels)
        precision = precision_score(true_labels, pred_labels, zero_division=0)
        recall = recall_score(true_labels, pred_labels, zero_division=0)
        f1 = f1_score(true_labels, pred_labels, zero_division=0)
        
        print("\n" + "="*60)
        print("TEST RESULTS")
        print("="*60)
        print(f"Dataset: DeepRescore2_UCEC")
        print(f"Model: Finetuned DeepDetect")
        print(f"AUC: {auc:.4f}")
        print(f"Accuracy: {accuracy:.4f}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall: {recall:.4f}")
        print(f"F1-Score: {f1:.4f}")
        print("="*60)
        
        return {
            'dataset': 'DeepRescore2_UCEC',
            'model': 'Finetuned_DeepDetect',
            'auc': auc,
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        }
        
    except Exception as e:
        print(f"✗ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    results = test_finetuned_model()
    if results:
        print("\n✓ Testing completed successfully!")
    else:
        print("\n✗ Testing failed!") 