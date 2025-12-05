#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试deepdetect_ori.pth模型在三个磷酸化测试集上的性能
"""

import torch
import pandas as pd
import numpy as np
import os
import sys
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score
from torch.utils.data import DataLoader

# 添加路径
sys.path.append('/data0/wangb/cd/duibi0820')
from original_deepdetect_adapter import OriginalDeepDetectAdapter
from data_preprocessor import DataPreprocessor, PeptideDataset

class DeepDetectOriTester:
    """DeepDetect_ori模型测试器"""
    
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model_path = '/data0/wangb/cd/duibi0826/0826comparison/model_weights/deepdetect_ori.pth'
        self.model = None
        
    def load_model(self):
        """加载Deepdetect_ori模型"""
        print("正在加载Deepdetect_ori模型...")
        try:
            self.model = OriginalDeepDetectAdapter()
            
            # 直接加载权重文件
            if os.path.exists(self.model_path):
                checkpoint = torch.load(self.model_path, map_location='cpu')
                
                # 检查checkpoint的格式并处理键名不匹配的问题
                if isinstance(checkpoint, dict):
                    if 'model_state_dict' in checkpoint:
                        state_dict = checkpoint['model_state_dict']
                    else:
                        state_dict = checkpoint
                else:
                    state_dict = checkpoint
                
                # 创建新的state_dict，添加"model."前缀
                new_state_dict = {}
                for key, value in state_dict.items():
                    if not key.startswith('model.'):
                        new_key = f'model.{key}'
                    else:
                        new_key = key
                    new_state_dict[new_key] = value
                
                # 加载修改后的state_dict
                self.model.load_state_dict(new_state_dict, strict=False)
                self.model.to(self.device)
                self.model.eval()
                print(f"✓ 模型加载成功 (设备: {self.device})")
                return True
            else:
                print(f"✗ 模型文件不存在: {self.model_path}")
                return False
        except Exception as e:
            print(f"✗ 模型加载失败: {e}")
            return False
    
    def evaluate_dataset(self, dataset_path, dataset_name):
        """评估单个数据集"""
        print(f"\n正在评估数据集: {dataset_name}")
        
        try:
            # 加载数据
            preprocessor = DataPreprocessor(dataset_path, model_type='deepdetect')
            preprocessor.load_data()  # 先加载数据
            sequences, labels = preprocessor.preprocess_sequences()
            
            # 创建数据集和数据加载器
            dataset = PeptideDataset(sequences, labels, model_type='deepdetect')
            dataloader = DataLoader(dataset, batch_size=64, shuffle=False)
            
            # 运行推理
            predictions = []
            true_labels = []
            
            with torch.no_grad():
                for batch_idx, batch in enumerate(dataloader):
                    inputs = batch[0].to(self.device)
                    outputs = self.model.predict(inputs)
                    predictions.extend(outputs.flatten().tolist())
                    true_labels.extend(batch[1].tolist())
            
            # 计算指标
            auc = roc_auc_score(true_labels, predictions)
            # 将预测值转换为二分类标签
            pred_labels = [1 if p > 0.5 else 0 for p in predictions]
            accuracy = accuracy_score(true_labels, pred_labels)
            precision = precision_score(true_labels, pred_labels, zero_division=0)
            recall = recall_score(true_labels, pred_labels, zero_division=0)
            f1 = f1_score(true_labels, pred_labels, zero_division=0)
            
            results = {
                'Dataset': dataset_name,
                'Dataset_Type': 'phosphorylation',
                'Model': 'Deepdetect_ori',
                'AUC': auc,
                'Accuracy': accuracy,
                'Precision': precision,
                'Recall': recall,
                'F1_Score': f1
            }
            
            print(f"✓ {dataset_name} 评估完成")
            print(f"  AUC: {auc:.4f}")
            print(f"  Accuracy: {accuracy:.4f}")
            print(f"  Precision: {precision:.4f}")
            print(f"  Recall: {recall:.4f}")
            print(f"  F1-Score: {f1:.4f}")
            
            return results
            
        except Exception as e:
            print(f"✗ 评估数据集 {dataset_name} 时出错: {e}")
            return None
    
    def evaluate_three_testsets(self):
        """评估三个磷酸化数据集"""
        datasets = [
            ('/data0/wangb/wbscy/PhosSight/data/test/DeepRescore2_HCC.csv', 'DeepRescore2_HCC'),
            ('/data0/wangb/wbscy/PhosSight/data/test/DeepRescore2_label_free.csv', 'DeepRescore2_label_free'),
            ('/data0/wangb/wbscy/PhosSight/data/test/DeepRescore2_UCEC.csv', 'DeepRescore2_UCEC')
        ]
        
        all_results = []
        
        for dataset_path, dataset_name in datasets:
            if os.path.exists(dataset_path):
                result = self.evaluate_dataset(dataset_path, dataset_name)
                if result:
                    all_results.append(result)
            else:
                print(f"数据集不存在: {dataset_path}")
        
        return all_results
    
    def save_results(self, results):
        """保存结果"""
        if results:
            df = pd.DataFrame(results)
            
            # 保存详细结果
            output_dir = "results/evaluation_results"
            os.makedirs(output_dir, exist_ok=True)
            
            # 保存到deepdetect_ori_phosphorylation_results.csv
            output_file = os.path.join(output_dir, "deepdetect_ori_phosphorylation_results.csv")
            df.to_csv(output_file, index=False)
            print(f"\n✓ 结果已保存到: {output_file}")
            
            # 计算平均性能
            avg_performance = df[['AUC', 'Accuracy', 'Precision', 'Recall', 'F1_Score']].mean()
            avg_df = pd.DataFrame([avg_performance])
            avg_df['Dataset_Type'] = 'phosphorylation'
            avg_df['Model'] = 'Deepdetect_ori'
            avg_df = avg_df[['Model', 'Dataset_Type', 'AUC', 'Accuracy', 'Precision', 'Recall', 'F1_Score']]
            
            # 保存平均性能
            avg_file = os.path.join(output_dir, "deepdetect_ori_phosphorylation_average_performance.csv")
            avg_df.to_csv(avg_file, index=False)
            print(f"✓ 平均性能已保存到: {avg_file}")
            
            return df
        else:
            print("没有结果可保存")
            return None

def main():
    tester = DeepDetectOriTester()
    
    if tester.load_model():
        results = tester.evaluate_three_testsets()
        tester.save_results(results)
    else:
        print("模型加载失败，无法进行评估")

if __name__ == "__main__":
    main() 