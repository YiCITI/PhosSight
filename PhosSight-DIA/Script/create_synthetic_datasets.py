#!/usr/bin/env python3
"""
为A549和syn数据集创建合成训练数据集
合并磷酸化肽段和非修饰肽段作为正样本，从full_dataset.csv中抽取等量负样本
"""

import os
import pandas as pd
import random
from typing import List, Tuple

def load_positive_samples(phospho_file: str, non_modified_file: str) -> List[Tuple[str, int]]:
    """加载正样本：磷酸化肽段 + 非修饰肽段"""
    positive_samples = []
    
    # 加载磷酸化肽段
    if os.path.exists(phospho_file):
        print(f"加载磷酸化肽段: {phospho_file}")
        with open(phospho_file, 'r', encoding='utf-8') as f:
            phospho_sequences = [line.strip() for line in f if line.strip()]
        print(f"  磷酸化肽段数量: {len(phospho_sequences)}")
        positive_samples.extend([(seq, 1) for seq in phospho_sequences])
    else:
        print(f"警告: 磷酸化肽段文件不存在: {phospho_file}")
    
    # 加载非修饰肽段
    if os.path.exists(non_modified_file):
        print(f"加载非修饰肽段: {non_modified_file}")
        with open(non_modified_file, 'r', encoding='utf-8') as f:
            non_modified_sequences = [line.strip() for line in f if line.strip()]
        print(f"  非修饰肽段数量: {len(non_modified_sequences)}")
        positive_samples.extend([(seq, 1) for seq in non_modified_sequences])
    else:
        print(f"警告: 非修饰肽段文件不存在: {non_modified_file}")
    
    print(f"正样本总数: {len(positive_samples)}")
    return positive_samples

def load_negative_samples(full_dataset_file: str, target_count: int) -> List[Tuple[str, int]]:
    """从full_dataset.csv中抽取负样本"""
    print(f"从 {full_dataset_file} 中抽取 {target_count} 个负样本...")
    
    if not os.path.exists(full_dataset_file):
        print(f"错误: 文件不存在 {full_dataset_file}")
        return []
    
    # 读取CSV文件
    df = pd.read_csv(full_dataset_file)
    print(f"  原始数据集大小: {len(df)}")
    
    # 检查列名
    print(f"  可用列: {list(df.columns)}")
    
    # 假设序列在第一列，标签在最后一列
    sequence_col = df.columns[0]
    label_col = df.columns[-1]
    
    print(f"  使用序列列: {sequence_col}")
    print(f"  使用标签列: {label_col}")
    
    # 过滤出负样本（标签为0）
    negative_df = df[df[label_col] == 0]
    print(f"  负样本数量: {len(negative_df)}")
    
    if len(negative_df) < target_count:
        print(f"  警告: 负样本数量不足，只有 {len(negative_df)} 个，需要 {target_count} 个")
        target_count = len(negative_df)
    
    # 随机抽取指定数量的负样本
    if target_count > 0:
        sampled_df = negative_df.sample(n=target_count, random_state=42)
        negative_samples = [(str(seq), 0) for seq in sampled_df[sequence_col]]
        print(f"  实际抽取负样本数量: {len(negative_samples)}")
    else:
        negative_samples = []
        print("  没有可用的负样本")
    
    return negative_samples

def create_synthetic_dataset(positive_samples: List[Tuple[str, int]], 
                           negative_samples: List[Tuple[str, int]], 
                           output_file: str) -> None:
    """创建合成数据集并保存"""
    print(f"\n创建合成数据集...")
    
    # 合并正负样本
    all_samples = positive_samples + negative_samples
    
    # 打乱顺序
    random.shuffle(all_samples)
    
    print(f"总样本数: {len(all_samples)}")
    print(f"正样本数: {len(positive_samples)}")
    print(f"负样本数: {len(negative_samples)}")
    
    # 统计标签分布
    label_counts = {}
    for _, label in all_samples:
        label_counts[label] = label_counts.get(label, 0) + 1
    
    print(f"标签分布: {label_counts}")
    
    # 保存到文件
    print(f"保存到: {output_file}")
    with open(output_file, 'w', encoding='utf-8') as f:
        for sequence, label in all_samples:
            f.write(f"{sequence}\t{label}\n")
    
    print("保存完成！")

def create_dataset_for_one(dataset_name: str, base_dir: str, full_dataset_file: str):
    """为单个数据集创建合成数据集"""
    print(f"\n{'='*60}")
    print(f"为 {dataset_name} 数据集创建合成数据集")
    print(f"{'='*60}")
    
    dataset_dir = os.path.join(base_dir, dataset_name, "for_finetuning_res")
    
    if not os.path.exists(dataset_dir):
        print(f"错误: 数据集目录不存在: {dataset_dir}")
        return
    
    phospho_file = os.path.join(dataset_dir, "phospho_peptides_for_finetuning.txt")
    non_modified_file = os.path.join(dataset_dir, "non_modified_peptides.txt")
    output_file = os.path.join(dataset_dir, "synthetic_dataset.txt")
    
    # 加载正样本
    positive_samples = load_positive_samples(phospho_file, non_modified_file)
    
    if not positive_samples:
        print(f"错误: {dataset_name} 数据集没有加载到正样本")
        return
    
    # 加载负样本（等量）
    negative_samples = load_negative_samples(full_dataset_file, len(positive_samples))
    
    if not negative_samples:
        print(f"错误: {dataset_name} 数据集没有加载到负样本")
        return
    
    # 创建合成数据集
    create_synthetic_dataset(positive_samples, negative_samples, output_file)
    
    print(f"\n✅ {dataset_name} 数据集合成完成！")
    print(f"文件保存在: {output_file}")

def main():
    """主函数"""
    base_dir = "/data0/wangb/wbscy/finetune0108"
    full_dataset_file = "/data0/wangb/wbscy/PhosSight/data/full_dataset.csv"
    
    # 如果full_dataset.csv不存在，尝试其他路径
    if not os.path.exists(full_dataset_file):
        # 尝试从czy目录中查找
        alternative_paths = [
            "/data0/wangb/wbscy/czy/for_finetuning_res/full_dataset.csv",
            "/data0/wangb/wbscy/PhosSight/data/full_dataset.csv"
        ]
        for alt_path in alternative_paths:
            if os.path.exists(alt_path):
                full_dataset_file = alt_path
                break
        else:
            print(f"警告: 找不到full_dataset.csv文件，将使用等量的随机负样本")
            full_dataset_file = None
    
    # 为A549数据集创建合成数据集
    create_dataset_for_one("A549", base_dir, full_dataset_file)
    
    # 为syn数据集创建合成数据集
    create_dataset_for_one("syn", base_dir, full_dataset_file)
    
    print(f"\n{'='*60}")
    print("所有合成数据集创建完成！")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()

