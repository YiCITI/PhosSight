#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
在七个测试集上微调DeepDetect_ori模型
参考evaluate_deepdetect_ori.py的代码
"""

import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
import os
import argparse
from sklearn.metrics import roc_auc_score
from torch.utils.data import DataLoader, Dataset
import sys
from tqdm import tqdm

# 添加DeepDetect模型路径
sys.path.append('/data0/wangb/cd/duibi0820')
from original_deepdetect_adapter import OriginalDeepDetectAdapter
from data_preprocessor import DataPreprocessor, PeptideDataset

class DeepDetectOriFinetuner:
    """DeepDetect_ori模型微调器"""
    
    def __init__(self, pretrained_model_path, device='cuda'):
        self.pretrained_model_path = pretrained_model_path
        self.device = torch.device(device if torch.cuda.is_available() else 'cpu')
        
        # 微调参数 - 参考PhosSight微调配置
        self.finetune_config = {
            'learning_rate': 1e-5,  # 与PhosSight保持一致
            'epochs': 20,  # 与PhosSight保持一致
            'batch_size': 16,  # 与PhosSight保持一致，更小的batch size
        }
        
        print(f"Using device: {self.device}")
        print(f"Fine-tune config: {self.finetune_config}")
    
    def load_pretrained_model(self):
        """加载预训练模型"""
        print(f"Loading pretrained model from {self.pretrained_model_path}")
        
        try:
            # 创建模型实例 - 使用OriginalDeepDetectAdapter
            model = OriginalDeepDetectAdapter()
            
            # 加载预训练权重
            if os.path.exists(self.pretrained_model_path):
                # 使用OriginalDeepDetectAdapter的load_model方法
                success = model.load_model(self.pretrained_model_path)
                if success:
                    model.to(self.device)
                    print("Pretrained weights loaded successfully")
                else:
                    print("Failed to load pretrained weights")
                    return None
            else:
                print(f"Pretrained model not found at {self.pretrained_model_path}")
                return None
            
            return model
        except Exception as e:
            print(f"Model loading failed: {e}")
            return None
    
    def prepare_data(self, test_file_path):
        """准备数据"""
        print(f"Loading data from {test_file_path}")
        
        try:
            # 使用DataPreprocessor加载数据 - 与成功评估代码保持一致
            preprocessor = DataPreprocessor(test_file_path, model_type='deepdetect')
            preprocessor.load_data()
            sequences, labels = preprocessor.preprocess_sequences()
            
            print(f"Loaded {len(sequences)} samples")
            print(f"Label distribution: {pd.Series(labels).value_counts().to_dict()}")
            
            # 创建数据集和数据加载器 - 与成功评估代码保持一致
            dataset = PeptideDataset(sequences, labels, model_type='deepdetect')
            dataloader = DataLoader(dataset, batch_size=self.finetune_config['batch_size'], shuffle=True)
            
            return dataloader, sequences, labels
        except Exception as e:
            print(f"Data preparation failed: {e}")
            return None, None, None
    
    def evaluate_model(self, model, sequences, labels):
        """评估模型性能"""
        model.eval()
        
        try:
            # 使用DataPreprocessor创建数据集 - 与成功评估代码保持一致
            dataset = PeptideDataset(sequences, labels, model_type='deepdetect')
            dataloader = DataLoader(dataset, batch_size=128, shuffle=False)
            
            all_predictions = []
            all_labels = []
            
            with torch.no_grad():
                for batch in dataloader:
                    inputs = batch[0].to(self.device)
                    predictions = model.predict(inputs)
                    all_predictions.extend(predictions.flatten().cpu().numpy().tolist())
                    all_labels.extend(batch[1].numpy().tolist())
            
            auc = roc_auc_score(all_labels, all_predictions)
            return auc
        except Exception as e:
            print(f"Evaluation failed: {e}")
            return 0.0
    
    def finetune_model(self, model, dataloader, sequences, labels, test_name):
        """微调模型"""
        print(f"Starting fine-tuning on {test_name}")
        
        # 记录初始性能
        initial_auc = self.evaluate_model(model, sequences, labels)
        print(f"Initial AUC: {initial_auc:.4f}")
        
        # 优化器和损失函数 - 参考PhosSight配置
        optimizer = optim.AdamW(
            model.parameters(), 
            lr=self.finetune_config['learning_rate'], 
            weight_decay=1e-4,  # 增加权重衰减
            betas=(0.9, 0.999)
        )
        criterion = nn.BCELoss()
        
        best_auc = initial_auc
        best_model_state = model.state_dict().copy()
        patience_counter = 0
        patience = 8  # 与PhosSight保持一致，减少早停耐心
        
        for epoch in range(self.finetune_config['epochs']):
            # 训练阶段
            model.train()
            train_loss = 0
            
            for batch in tqdm(dataloader, desc=f"Epoch {epoch+1}"):
                inputs = batch[0].to(self.device)
                targets = batch[1].to(self.device)
                
                optimizer.zero_grad()
                outputs = model(inputs)
                
                if outputs.dim() > 1:
                    outputs = outputs.squeeze()
                
                loss = criterion(outputs, targets)
                
                # 反向传播
                loss.backward()
                
                # 更严格的梯度裁剪
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5)
                
                optimizer.step()
                
                train_loss += loss.item()
            
            # 在每个epoch后评估完整测试集
            print(f"  Epoch {epoch+1}: Train Loss: {train_loss/len(dataloader):.4f}")
            
            current_auc = self.evaluate_model(model, sequences, labels)
            print(f"  Epoch {epoch+1}: Current AUC: {current_auc:.4f}")
            
            # 保存最佳模型
            if current_auc > best_auc:
                best_auc = current_auc
                best_model_state = model.state_dict().copy()
                patience_counter = 0
                print(f"  New best AUC: {best_auc:.4f}")
            else:
                patience_counter += 1
                print(f"  No improvement, patience: {patience_counter}/{patience}")
                
                # 早停检查
                if patience_counter >= patience:
                    print(f"  Early stopping triggered")
                    break
        
        # 加载最佳模型状态
        model.load_state_dict(best_model_state)
        print(f"  Final best AUC: {best_auc:.4f}")
        print(f"  AUC improvement: {best_auc - initial_auc:.4f}")
        
        return model, best_auc, initial_auc
    
    def finetune_on_seven_testsets(self):
        """在七个测试集上微调模型"""
        
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
        output_dir = '/data0/wangb/cd/duibi0826/0826comparison/finetuned_models'
        os.makedirs(output_dir, exist_ok=True)
        
        results = []
        
        for test_file in test_sets:
            print(f"\n{'='*60}")
            print(f"Processing: {test_file}")
            print(f"{'='*60}")
            
            # 加载预训练模型
            model = self.load_pretrained_model()
            if model is None:
                print(f"Failed to load model for {test_file}")
                continue
            
            # 准备数据
            test_file_path = os.path.join(test_data_dir, test_file)
            dataloader, sequences, labels = self.prepare_data(test_file_path)
            if dataloader is None:
                print(f"Failed to prepare data for {test_file}")
                continue
            
            # 微调模型
            finetuned_model, best_auc, initial_auc = self.finetune_model(model, dataloader, sequences, labels, test_file)
            
            # 保存微调后的模型（只有当性能超过原始性能时才保存）
            test_name = test_file.replace('.csv', '')
            if best_auc > initial_auc:
                model_save_path = os.path.join(output_dir, f'deepdetect_ori_finetuned_{test_name}_auc{best_auc:.4f}.pth')
                finetuned_model.save_model(model_save_path)
                print(f"  ✅ Performance improved, model saved to: {model_save_path}")
            else:
                model_save_path = os.path.join(output_dir, f'deepdetect_ori_finetuned_{test_name}_no_improvement.pth')
                finetuned_model.save_model(model_save_path)
                print(f"  ❌ No improvement, model saved to: {model_save_path}")
            
            results.append({
                'test_set': test_file,
                'initial_auc': initial_auc,
                'best_auc': best_auc,
                'improvement': best_auc - initial_auc,
                'model_path': model_save_path
            })
        
        # 打印汇总结果
        print(f"\n{'='*80}")
        print("FINETUNING SUMMARY")
        print(f"{'='*80}")
        print(f"{'Test Set':<25} {'Initial AUC':<12} {'Best AUC':<12} {'Improvement':<12} {'Model Path'}")
        print("-" * 80)
        
        for result in results:
            print(f"{result['test_set']:<25} {result['initial_auc']:<12.4f} {result['best_auc']:<12.4f} {result['improvement']:<12.4f} {result['model_path']}")
        
        # 计算平均指标
        avg_initial_auc = np.mean([r['initial_auc'] for r in results])
        avg_best_auc = np.mean([r['best_auc'] for r in results])
        avg_improvement = np.mean([r['improvement'] for r in results])
        
        print("-" * 80)
        print(f"{'Average':<25} {avg_initial_auc:<12.4f} {avg_best_auc:<12.4f} {avg_improvement:<12.4f}")
        print("=" * 80)
        
        return results

def main():
    parser = argparse.ArgumentParser(description='Finetune DeepDetect_ori model on seven test sets')
    parser.add_argument('--pretrained_model', type=str, 
                       default='/data0/wangb/cd/duibi0826/0826comparison/model_weights/deepdetect_ori.pth',
                       help='Path to pretrained DeepDetect_ori model')
    parser.add_argument('--device', type=str, default='cuda', help='Device to use (cuda/cpu)')
    
    args = parser.parse_args()
    
    # 创建微调器
    finetuner = DeepDetectOriFinetuner(args.pretrained_model, args.device)
    
    # 开始微调
    results = finetuner.finetune_on_seven_testsets()
    
    print("\nFinetuning completed successfully!")

if __name__ == "__main__":
    main() 