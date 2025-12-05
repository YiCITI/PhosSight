#!/usr/bin/env python3
"""
在七个测试集上微调DLOmix模型 - 正确版本
在完整测试集上微调，然后在完整测试集上评估AUC
"""

import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from tqdm import tqdm
import argparse

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

class DLOmixFinetuner:
    """DLOmix模型微调器 - 正确版本"""
    
    def __init__(self, pretrained_model_path, device='cuda'):
        self.device = torch.device(device if torch.cuda.is_available() else 'cpu')
        self.pretrained_model_path = pretrained_model_path
        
        # 模型配置
        self.config = {
            'num_units': 128,
            'alphabet_size': 21,
            'num_classes': 4
        }
        
        # 微调参数 - 参考PhosSight微调配置
        self.finetune_config = {
            'learning_rate': 1e-5,  # 与PhosSight保持一致
            'epochs': 20,  # 与PhosSight保持一致
            'batch_size': 16,  # 与PhosSight保持一致，更小的batch size
        }
        
        print(f"Using device: {self.device}")
    
    def load_pretrained_model(self):
        """加载预训练模型"""
        model = DLOmixAdapter(**self.config)
        model.load_model(self.pretrained_model_path)
        model.to(self.device)
        print(f"Loaded pretrained model from: {self.pretrained_model_path}")
        return model
    
    def prepare_data(self, test_file_path):
        """准备数据"""
        print(f"Loading data from: {test_file_path}")
        
        # 读取数据
        df = pd.read_csv(test_file_path, header=None, names=['sequence', 'label'])
        sequences = df['sequence'].tolist()
        labels = df['label'].tolist()
        
        print(f"  Total samples: {len(sequences)}")
        print(f"  Positive samples: {sum(labels)}, Negative samples: {len(labels) - sum(labels)}")
        
        # 创建数据集和数据加载器
        dataset = PeptideDataset(sequences, labels)
        dataloader = DataLoader(dataset, batch_size=self.finetune_config['batch_size'], shuffle=True)
        
        return dataloader, sequences, labels
    
    def evaluate_model(self, model, sequences, labels):
        """在完整测试集上评估模型"""
        model.eval()
        
        # 创建评估数据集
        eval_dataset = PeptideDataset(sequences, labels)
        eval_loader = DataLoader(eval_dataset, batch_size=128, shuffle=False)
        
        all_predictions = []
        all_labels = []
        
        with torch.no_grad():
            for batch_sequences, batch_labels in tqdm(eval_loader, desc="Evaluating"):
                batch_sequences = batch_sequences.to(self.device)
                batch_labels = batch_labels.to(self.device)
                
                # 获取预测
                predictions = model.predict(batch_sequences)
                predictions = predictions.cpu().numpy().flatten()
                
                all_predictions.extend(predictions)
                all_labels.extend(batch_labels.cpu().numpy().flatten())
        
        # 计算AUC
        auc = roc_auc_score(all_labels, all_predictions)
        return auc
    
    def finetune_model(self, model, dataloader, sequences, labels, test_name):
        """在完整测试集上微调模型"""
        print(f"\nStarting finetuning on {test_name}")
        
        # 首先评估微调前的性能
        print("Evaluating before finetuning...")
        initial_auc = self.evaluate_model(model, sequences, labels)
        print(f"  Initial AUC: {initial_auc:.4f}")
        
        # 优化器和损失函数 - 参考PhosSight配置
        optimizer = optim.AdamW(
            model.parameters(), 
            lr=self.finetune_config['learning_rate'], 
            weight_decay=1e-4,  # 增加权重衰减
            betas=(0.9, 0.999)
        )
        criterion = nn.CrossEntropyLoss()
        
        best_auc = initial_auc
        best_model_state = model.state_dict().copy()
        patience_counter = 0
        patience = 8  # 与PhosSight保持一致，减少早停耐心
        
        for epoch in range(self.finetune_config['epochs']):
            # 训练阶段
            model.train()
            train_loss = 0.0
            
            for batch_sequences, batch_labels in tqdm(dataloader, desc=f"Epoch {epoch+1}/{self.finetune_config['epochs']}"):
                batch_sequences = batch_sequences.to(self.device)
                batch_labels = batch_labels.to(self.device)
                
                # 前向传播
                optimizer.zero_grad()
                outputs = model(batch_sequences)
                
                # 计算损失（需要将二分类标签转换为4分类标签）
                # 根据您的模型配置：0->0 (Non-Flyer), 1->3 (Strong Flyer，代表磷酸化)
                # 因为您的模型训练时是将1映射到3的
                batch_labels_4class = torch.where(batch_labels == 1, torch.tensor(3, device=self.device), torch.tensor(0, device=self.device))
                loss = criterion(outputs, batch_labels_4class.long())
                
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
    
    def finetune_on_all_testsets(self):
        """在所有测试集上微调模型"""
        
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
            
            # 准备数据
            test_file_path = os.path.join(test_data_dir, test_file)
            dataloader, sequences, labels = self.prepare_data(test_file_path)
            
            # 微调模型
            finetuned_model, best_auc, initial_auc = self.finetune_model(model, dataloader, sequences, labels, test_file)
            
            # 保存微调后的模型（只有当性能超过原始性能时才保存）
            test_name = test_file.replace('.csv', '')
            if best_auc > initial_auc:
                model_save_path = os.path.join(output_dir, f'dlomix_finetuned_{test_name}_auc{best_auc:.4f}.pth')
                finetuned_model.save_model(model_save_path)
                print(f"  ✅ Performance improved, model saved to: {model_save_path}")
            else:
                model_save_path = os.path.join(output_dir, f'dlomix_finetuned_{test_name}_no_improvement.pth')
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
    parser = argparse.ArgumentParser(description='Finetune DLOmix model on seven test sets')
    parser.add_argument('--pretrained_model', type=str, 
                       default='/data0/wangb/cd/duibi0826/0826comparison/model_weights/dlomix_best.pth',
                       help='Path to pretrained DLOmix model')
    parser.add_argument('--device', type=str, default='cuda', help='Device to use (cuda/cpu)')
    
    args = parser.parse_args()
    
    # 创建微调器
    finetuner = DLOmixFinetuner(args.pretrained_model, args.device)
    
    # 开始微调
    results = finetuner.finetune_on_all_testsets()
    
    print("\nFinetuning completed successfully!")

if __name__ == "__main__":
    main() 