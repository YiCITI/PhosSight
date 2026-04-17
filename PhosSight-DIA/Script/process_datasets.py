#!/usr/bin/env python3
"""
处理A549和syn数据集，提取磷酸化肽段和非修饰肽段
"""

import os
import pandas as pd
import re
from typing import List, Set
import glob

def process_modified_sequence(modified_seq: str) -> str:
    """
    处理修饰序列：
    1. 去除(UniMod:4)固定修饰
    2. 将(UniMod:21)磷酸化修饰转换为小写sty
    3. 确保只有一个磷酸化位点
    """
    if pd.isna(modified_seq) or not isinstance(modified_seq, str):
        return None
    
    # 去除固定修饰(UniMod:4)
    seq = modified_seq.replace('(UniMod:4)', '')
    
    # 检查是否有磷酸化修饰(UniMod:21)
    phospho_matches = re.findall(r'([STY])\(UniMod:21\)', seq)
    
    # 确保只有一个磷酸化位点
    if len(phospho_matches) == 0:
        return None  # 没有磷酸化修饰
    elif len(phospho_matches) > 1:
        return None  # 多个磷酸化位点，跳过
    
    # 将磷酸化修饰转换为小写
    seq = re.sub(r'([STY])\(UniMod:21\)', lambda m: m.group(1).lower(), seq)
    
    return seq

def is_non_modified_sequence(modified_seq: str) -> bool:
    """
    判断序列是否为非修饰序列：
    1. 不包含(UniMod:4)固定修饰
    2. 不包含(UniMod:21)磷酸化修饰
    3. 序列不为空且有效
    """
    if pd.isna(modified_seq) or not isinstance(modified_seq, str):
        return False
    
    # 检查是否包含任何修饰
    has_unimod4 = '(UniMod:4)' in modified_seq
    has_unimod21 = '(UniMod:21)' in modified_seq
    
    # 只有不包含任何修饰的序列才是非修饰序列
    return not (has_unimod4 or has_unimod21)

def process_parquet_file_for_phospho(file_path: str) -> List[str]:
    """处理单个parquet文件，提取磷酸化肽段"""
    print(f"正在处理文件: {file_path}")
    
    try:
        df = pd.read_parquet(file_path)
        print(f"  文件行数: {len(df)}")
        
        if 'Modified.Sequence' not in df.columns:
            print(f"  警告: 文件中没有Modified.Sequence列")
            return []
        
        phospho_sequences = []
        total_sequences = 0
        valid_phospho = 0
        
        for idx, modified_seq in enumerate(df['Modified.Sequence']):
            total_sequences += 1
            
            if idx % 10000 == 0 and idx > 0:
                print(f"  已处理 {idx}/{len(df)} 序列...")
            
            processed_seq = process_modified_sequence(modified_seq)
            
            if processed_seq is not None:
                phospho_sequences.append(processed_seq)
                valid_phospho += 1
        
        print(f"  处理完成: 总序列数={total_sequences}, 有效单磷酸化序列={valid_phospho}")
        return phospho_sequences
        
    except Exception as e:
        print(f"  处理文件时出错: {str(e)}")
        return []

def process_parquet_file_for_non_modified(file_path: str) -> List[str]:
    """处理单个parquet文件，提取非修饰肽段"""
    print(f"正在处理文件: {file_path}")
    
    try:
        df = pd.read_parquet(file_path)
        print(f"  文件行数: {len(df)}")
        
        if 'Modified.Sequence' not in df.columns:
            print(f"  警告: 文件中没有Modified.Sequence列")
            return []
        
        non_modified_sequences = []
        total_sequences = 0
        valid_non_modified = 0
        
        for idx, modified_seq in enumerate(df['Modified.Sequence']):
            total_sequences += 1
            
            if idx % 10000 == 0 and idx > 0:
                print(f"  已处理 {idx}/{len(df)} 序列...")
            
            if is_non_modified_sequence(modified_seq):
                non_modified_sequences.append(modified_seq)
                valid_non_modified += 1
        
        print(f"  处理完成: 总序列数={total_sequences}, 非修饰序列={valid_non_modified}")
        return non_modified_sequences
        
    except Exception as e:
        print(f"  处理文件时出错: {str(e)}")
        return []

def process_dataset(base_dir: str, dataset_name: str):
    """处理单个数据集"""
    print(f"\n{'='*60}")
    print(f"处理数据集: {dataset_name}")
    print(f"基础目录: {base_dir}")
    print(f"{'='*60}")
    
    # 查找所有report.parquet文件
    parquet_files = glob.glob(os.path.join(base_dir, "*/report.parquet"))
    
    print(f"找到 {len(parquet_files)} 个report.parquet文件")
    
    # 提取磷酸化肽段
    all_phospho_sequences = []
    for file_path in sorted(parquet_files):
        folder_name = os.path.basename(os.path.dirname(file_path))
        print(f"\n处理文件夹: {folder_name} (磷酸化肽段)")
        sequences = process_parquet_file_for_phospho(file_path)
        all_phospho_sequences.extend(sequences)
        print(f"  当前总计: {len(all_phospho_sequences)} 条磷酸化肽段")
    
    # 提取非修饰肽段
    all_non_modified_sequences = []
    for file_path in sorted(parquet_files):
        folder_name = os.path.basename(os.path.dirname(file_path))
        print(f"\n处理文件夹: {folder_name} (非修饰肽段)")
        sequences = process_parquet_file_for_non_modified(file_path)
        all_non_modified_sequences.extend(sequences)
        print(f"  当前总计: {len(all_non_modified_sequences)} 条非修饰肽段")
    
    # 去重并保存
    unique_phospho = sorted(list(set(all_phospho_sequences)))
    unique_non_modified = sorted(list(set(all_non_modified_sequences)))
    
    # 创建输出目录（base_dir 已经是 for_finetuning_res 目录）
    output_dir = base_dir
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存磷酸化肽段
    phospho_file = os.path.join(output_dir, "phospho_peptides_for_finetuning.txt")
    with open(phospho_file, 'w', encoding='utf-8') as f:
        for seq in unique_phospho:
            f.write(f"{seq}\n")
    print(f"\n✅ 保存磷酸化肽段: {phospho_file} ({len(unique_phospho)} 条唯一序列)")
    
    # 保存非修饰肽段
    non_modified_file = os.path.join(output_dir, "non_modified_peptides.txt")
    with open(non_modified_file, 'w', encoding='utf-8') as f:
        for seq in unique_non_modified:
            f.write(f"{seq}\n")
    print(f"✅ 保存非修饰肽段: {non_modified_file} ({len(unique_non_modified)} 条唯一序列)")
    
    # 统计信息
    print(f"\n=== {dataset_name} 数据集统计 ===")
    print(f"磷酸化肽段: {len(unique_phospho)} 条唯一序列")
    print(f"非修饰肽段: {len(unique_non_modified)} 条唯一序列")
    
    # 磷酸化位点分布
    if unique_phospho:
        s_count = sum(1 for seq in unique_phospho if 's' in seq)
        t_count = sum(1 for seq in unique_phospho if 't' in seq)
        y_count = sum(1 for seq in unique_phospho if 'y' in seq)
        print(f"\n磷酸化位点分布:")
        print(f"  Serine (s): {s_count} ({s_count/len(unique_phospho)*100:.1f}%)")
        print(f"  Threonine (t): {t_count} ({t_count/len(unique_phospho)*100:.1f}%)")
        print(f"  Tyrosine (y): {y_count} ({y_count/len(unique_phospho)*100:.1f}%)")
    
    return phospho_file, non_modified_file

def main():
    """主函数"""
    base_dir = "/data0/wangb/wbscy/finetune0108"
    
    # 处理A549数据集
    a549_dir = os.path.join(base_dir, "A549", "for_finetuning_res")
    if os.path.exists(a549_dir):
        process_dataset(a549_dir, "A549")
    else:
        print(f"警告: A549目录不存在: {a549_dir}")
    
    # 处理syn数据集
    syn_dir = os.path.join(base_dir, "syn", "for_finetuning_res")
    if os.path.exists(syn_dir):
        process_dataset(syn_dir, "syn")
    else:
        print(f"警告: syn目录不存在: {syn_dir}")
    
    print(f"\n{'='*60}")
    print("所有数据集处理完成！")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()

