#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DLOmix模型训练脚本 - 专门用于balanced dataset训练
"""

import os
import sys
import argparse
import logging
import subprocess
from pathlib import Path
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from tqdm import tqdm

# 添加项目路径
sys.path.append('/data0/wangb/cd/duibi0826')
from model_adapters.dlomix_adapter import DLOmixAdapter

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PeptideDataset(Dataset):
    """肽段数据集"""
    
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

class DLOmixTrainer:
    """DLOmix模型训练器"""
    
    def __init__(self, config):
        self.config = config
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logger.info(f"Using device: {self.device}")
        
        # 创建模型
        self.model = DLOmixAdapter(
            num_units=config['num_units'],
            alphabet_size=config['alphabet_size'],
            num_classes=config['num_classes']
        ).to(self.device)
        
        # 损失函数和优化器
        self.criterion = nn.CrossEntropyLoss()
        self.optimizer = optim.Adam(self.model.parameters(), lr=config['learning_rate'], weight_decay=config['weight_decay'])
        
        # 学习率调度器
        self.scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer, mode='max', factor=0.7, patience=8, verbose=True, min_lr=1e-6
        )
        
        # 训练历史
        self.train_losses = []
        self.val_losses = []
        self.train_aucs = []
        self.val_aucs = []
        
        # 最佳模型
        self.best_val_auc = 0.0
        self.best_model_state = None
        
    def load_data(self):
        """加载数据"""
        logger.info("Loading balanced dataset...")
        
        # 读取数据
        df = pd.read_csv(self.config['data_path'], header=None, names=['sequence', 'label'])
        
        sequences = df['sequence'].tolist()
        labels = df['label'].tolist()
        
        # 数据统计
        total_samples = len(sequences)
        positive_samples = sum(labels)
        negative_samples = total_samples - positive_samples
        
        logger.info(f"Dataset statistics:")
        logger.info(f"  Total samples: {total_samples:,}")
        logger.info(f"  Positive samples: {positive_samples:,} ({positive_samples/total_samples*100:.1f}%)")
        logger.info(f"  Negative samples: {negative_samples:,} ({negative_samples/total_samples*100:.1f}%)")
        
        # 分割数据
        train_sequences, val_sequences, train_labels, val_labels = train_test_split(
            sequences, labels, test_size=0.2, random_state=42, stratify=labels
        )
        
        # 创建数据集
        train_dataset = PeptideDataset(train_sequences, train_labels)
        val_dataset = PeptideDataset(val_sequences, val_labels)
        
        # 创建数据加载器
        train_loader = DataLoader(train_dataset, batch_size=self.config['batch_size'], shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=self.config['batch_size'], shuffle=False)
        
        return train_loader, val_loader
    
    def train_epoch(self, train_loader):
        """训练一个epoch"""
        self.model.train()
        total_loss = 0
        all_predictions = []
        all_labels = []
        
        for batch_idx, (data, target) in enumerate(tqdm(train_loader, desc="Training")):
            data, target = data.to(self.device), target.to(self.device)
            
            self.optimizer.zero_grad()
            output = self.model(data)
            loss = self.criterion(output, target.long())
            loss.backward()
            self.optimizer.step()
            
            total_loss += loss.item()
            
            # 转换为二分类预测用于AUC计算
            binary_pred = self.model.predict(data)
            all_predictions.extend(binary_pred.detach().cpu().numpy().flatten())
            all_labels.extend(target.cpu().numpy().flatten())
        
        # 计算AUC
        auc = roc_auc_score(all_labels, all_predictions)
        
        return total_loss / len(train_loader), auc
    
    def validate_epoch(self, val_loader):
        """验证一个epoch"""
        self.model.eval()
        total_loss = 0
        all_predictions = []
        all_labels = []
        
        with torch.no_grad():
            for data, target in tqdm(val_loader, desc="Validation"):
                data, target = data.to(self.device), target.to(self.device)
                output = self.model(data)
                loss = self.criterion(output, target.long())
                
                total_loss += loss.item()
                
                # 转换为二分类预测
                binary_pred = self.model.predict(data)
                all_predictions.extend(binary_pred.cpu().numpy().flatten())
                all_labels.extend(target.cpu().numpy().flatten())
        
        # 计算指标
        auc = roc_auc_score(all_labels, all_predictions)
        predictions_binary = (np.array(all_predictions) > 0.5).astype(int)
        accuracy = accuracy_score(all_labels, predictions_binary)
        precision = precision_score(all_labels, predictions_binary)
        recall = recall_score(all_labels, predictions_binary)
        f1 = f1_score(all_labels, predictions_binary)
        
        return {
            'loss': total_loss / len(val_loader),
            'auc': auc,
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1': f1
        }
    
    def train(self, train_loader, val_loader):
        """训练模型"""
        logger.info("Starting DLOmix model training...")
        
        for epoch in range(self.config['epochs']):
            logger.info(f"Epoch {epoch+1}/{self.config['epochs']}")
            
            # 训练
            train_loss, train_auc = self.train_epoch(train_loader)
            
            # 验证
            val_metrics = self.validate_epoch(val_loader)
            
            # 记录历史
            self.train_losses.append(train_loss)
            self.val_losses.append(val_metrics['loss'])
            self.train_aucs.append(train_auc)
            self.val_aucs.append(val_metrics['auc'])
            
            # 学习率调度
            self.scheduler.step(val_metrics['auc'])
            
            # 打印结果
            logger.info(f"  Train Loss: {train_loss:.4f}, Train AUC: {train_auc:.4f}")
            logger.info(f"  Val Loss: {val_metrics['loss']:.4f}, Val AUC: {val_metrics['auc']:.4f}")
            logger.info(f"  Val Accuracy: {val_metrics['accuracy']:.4f}, Val F1: {val_metrics['f1']:.4f}")
            
            # 保存最佳模型
            if val_metrics['auc'] > self.best_val_auc:
                self.best_val_auc = val_metrics['auc']
                self.best_model_state = self.model.state_dict().copy()
                logger.info(f"  New best model! Val AUC: {self.best_val_auc:.4f}")
        
        # 保存最佳模型
        if self.best_model_state is not None:
            self.model.load_state_dict(self.best_model_state)
            torch.save(self.best_model_state, self.config['model_save_path'])
            logger.info(f"Best model saved to: {self.config['model_save_path']}")
        
        return self.best_val_auc
    
    def plot_training_history(self):
        """绘制训练历史"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # 损失曲线
        ax1.plot(self.train_losses, label='Train Loss')
        ax1.plot(self.val_losses, label='Val Loss')
        ax1.set_title('Training and Validation Loss')
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('Loss')
        ax1.legend()
        ax1.grid(True)
        
        # AUC曲线
        ax2.plot(self.train_aucs, label='Train AUC')
        ax2.plot(self.val_aucs, label='Val AUC')
        ax2.set_title('Training and Validation AUC')
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('AUC')
        ax2.legend()
        ax2.grid(True)
        
        # 学习率曲线
        ax3.plot([opt['lr'] for opt in self.optimizer.param_groups], label='Learning Rate')
        ax3.set_title('Learning Rate Schedule')
        ax3.set_xlabel('Epoch')
        ax3.set_ylabel('Learning Rate')
        ax3.legend()
        ax3.grid(True)
        
        # 最终验证指标
        final_val_metrics = {
            'AUC': self.val_aucs[-1],
            'Loss': self.val_losses[-1],
            'Best AUC': self.best_val_auc
        }
        
        ax4.bar(final_val_metrics.keys(), final_val_metrics.values())
        ax4.set_title('Final Validation Metrics')
        ax4.set_ylabel('Value')
        
        plt.tight_layout()
        plt.savefig(self.config['plot_save_path'], dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Training history plot saved to: {self.config['plot_save_path']}")

def main():
    parser = argparse.ArgumentParser(description='Train DLOmix model on balanced dataset')
    parser.add_argument('--data_path', type=str, 
                       default='/data0/wangb/wbscy/PhosSight/data/train/balanced_dataset_1.csv',
                       help='Path to balanced dataset')
    parser.add_argument('--output_dir', type=str, 
                       default='/data0/wangb/cd/duibi0826/0826comparison',
                       help='Output directory')
    
    args = parser.parse_args()
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "model_weights").mkdir(exist_ok=True)
    (output_dir / "plots").mkdir(exist_ok=True)
    
    # 配置
    config = {
        'data_path': args.data_path,
        'num_units': 128,
        'alphabet_size': 21,
        'num_classes': 4,
        'batch_size': 64,
        'learning_rate': 0.001,
        'weight_decay': 1e-5,
        'epochs': 50,
        'model_save_path': output_dir / "model_weights" / "dlomix_best.pth",
        'plot_save_path': output_dir / "plots" / "dlomix_training_history.png"
    }
    
    # 创建训练器
    trainer = DLOmixTrainer(config)
    
    # 加载数据
    train_loader, val_loader = trainer.load_data()
    
    # 训练模型
    best_auc = trainer.train(train_loader, val_loader)
    
    # 绘制训练历史
    trainer.plot_training_history()
    
    logger.info(f"DLOmix training completed! Best validation AUC: {best_auc:.4f}")

if __name__ == "__main__":
    main() 