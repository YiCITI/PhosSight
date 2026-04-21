#!/usr/bin/env python3
"""
从 Yeast Phospho PeptideAtlas 构建评测数据集
- 正样本：GSB 中检测到的磷酸化肽段（磷酸化位点小写）
- 负样本：酵母蛋白组理论消化后，含 S/T/Y 的肽段但未在 GSB 中检测到

作者: PhosSight Team
日期: 2026-04-16
"""

import re
import pandas as pd
from pathlib import Path
import random

random.seed(42)

# 路径配置
DATA_DIR = Path("/data0/wangb/cd/duibi0826/0826comparison/data/yeast/public_dataset")
OUTPUT_DIR = DATA_DIR

# 磷酸化氨基酸
PHOSPHO_AA = {'S', 'T', 'Y'}


def trypsin_digest(protein_sequence, min_len=7, max_len=51, max_missed=2):
    """
    trypsin 消化蛋白质
    切割位点：K, R（但不在 P 后）
    """
    peptides = []
    
    # 找到所有切割位点
    cleavage_sites = [0]  # N端
    
    for i, aa in enumerate(protein_sequence):
        if aa in 'KR' and i + 1 < len(protein_sequence):
            next_aa = protein_sequence[i + 1]
            if next_aa != 'P':  # 不在 P 后切割
                cleavage_sites.append(i + 1)
    
    cleavage_sites.append(len(protein_sequence))  # C端
    
    # 生成肽段
    for i in range(len(cleavage_sites) - 1):
        start = cleavage_sites[i]
        for j in range(i + 1, min(i + 1 + max_missed + 1, len(cleavage_sites))):
            end = cleavage_sites[j]
            peptide = protein_sequence[start:end]
            if min_len <= len(peptide) <= max_len:
                peptides.append(peptide)
    
    return peptides


def parse_phosphopeptide(mod_seq):
    """
    从修饰序列中提取磷酸化肽段
    将 [Phospho] 前的氨基酸转为小写
    """
    if pd.isna(mod_seq) or not mod_seq:
        return None, 0
    
    mod_seq = str(mod_seq)
    
    # 找所有 [Phospho] 的位置
    phos_positions = []
    idx = 0
    while idx < len(mod_seq):
        pos = mod_seq.find('[Phospho]', idx)
        if pos == -1:
            break
        if pos > 0:
            aa = mod_seq[pos - 1]
            if aa in 'STY':
                phos_positions.append(pos - 1)
        idx = pos + len('[Phospho]')
    
    if not phos_positions:
        return None, 0
    
    # 移除所有 [xxx] 修饰标记
    clean_seq = re.sub(r'\[.*?\]', '', mod_seq)
    clean_seq = re.sub(r'-\d+$', '', clean_seq)
    clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', clean_seq)
    
    if not clean_seq:
        return None, 0
    
    # 将磷酸化位点转为小写
    result = list(clean_seq)
    clean_idx = 0
    idx = 0
    while idx < len(mod_seq):
        if mod_seq[idx] == '[':
            end = mod_seq.find(']', idx)
            if end != -1:
                if mod_seq[idx+1:end] == 'Phospho' and clean_idx > 0:
                    result[clean_idx - 1] = result[clean_idx - 1].lower()
                idx = end + 1
        else:
            clean_idx += 1
            idx += 1
    
    return ''.join(result), len(phos_positions)


def extract_positive_from_gsb(csv_path, min_len=7, max_len=51):
    """从 GSB.csv 提取磷酸化肽段（正样本）"""
    print(f"读取磷酸化位点数据: {csv_path}")
    df = pd.read_csv(csv_path, low_memory=False)
    print(f"总记录数: {len(df)}")
    
    # 只保留 Gold 和 Silver
    high_conf = df[df['PTM_FLR_category'].isin(['Gold', 'Silver'])].copy()
    print(f"高置信度 (Gold/Silver) 记录数: {len(high_conf)}")
    
    # 肽段列
    peptide_cols = [col for col in df.columns if '_peptide_mod_pos' in col]
    print(f"肽段来源列: {len(peptide_cols)}")
    
    # 提取磷酸化肽段
    phosphopeptides = {}
    for col in peptide_cols:
        for mod_seq in high_conf[col]:
            peptide, num_phos = parse_phosphopeptide(mod_seq)
            if peptide and min_len <= len(peptide) <= max_len and num_phos > 0:
                if peptide not in phosphopeptides:
                    phosphopeptides[peptide] = True
    
    print(f"提取的磷酸化肽段数 (去重): {len(phosphopeptides)}")
    return list(phosphopeptides.keys()), phosphopeptides


def generate_negative_peptides(fasta_path, positive_set):
    """
    从蛋白组生成负样本
    负样本：理论酶切后，含 S/T/Y 的肽段，但未在 GSB 中检测到
    负样本同样用小写字母标注磷酸化位点（理论磷酸化形式）
    """
    print(f"\n读取蛋白组: {fasta_path}")
    
    # 正样本集合（磷酸化形式，如 "ABCsDEF"）
    negative_peptides = {}
    
    with open(fasta_path, 'r') as f:
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if sequence:
                    # 消化当前蛋白
                    peptides = trypsin_digest(sequence)
                    for pep in peptides:
                        if 7 <= len(pep) <= 51:
                            # 检查是否含磷酸化位点（S/T/Y）
                            for i, aa in enumerate(pep):
                                if aa in PHOSPHO_AA:
                                    # 生成该位点磷酸化的版本（小写）
                                    phos_pep = pep[:i] + aa.lower() + pep[i+1:]
                                    # 排除已检测到的磷酸化肽段
                                    if phos_pep not in positive_set:
                                        negative_peptides[phos_pep] = True
                sequence = ''
            else:
                sequence += line.strip()
        
        # 处理最后一个蛋白
        if sequence:
            peptides = trypsin_digest(sequence)
            for pep in peptides:
                if 7 <= len(pep) <= 51:
                    for i, aa in enumerate(pep):
                        if aa in PHOSPHO_AA:
                            phos_pep = pep[:i] + aa.lower() + pep[i+1:]
                            if phos_pep not in positive_set:
                                negative_peptides[phos_pep] = True
    
    print(f"负样本候选数 (含 s/t/y 但未检测到): {len(negative_peptides)}")
    return list(negative_peptides.keys())
    return negative_peptides


def analyze_composition(peptides, name):
    """分析组成"""
    if not peptides:
        print(f"\n{name}: 无数据")
        return
    
    total = len(peptides)
    pS = sum(1 for p in peptides if 's' in p)
    pT = sum(1 for p in peptides if 't' in p)
    pY = sum(1 for p in peptides if 'y' in p)
    
    print(f"\n{name} 组成:")
    print(f"  总数: {total}")
    print(f"  含 pSer (s): {pS} ({pS/total*100:.1f}%)")
    print(f"  含 pThr (t): {pT} ({pT/total*100:.1f}%)")
    print(f"  含 pTyr (y): {pY} ({pY/total*100:.1f}%)")


def main():
    print("=" * 70)
    print("Yeast Phospho PeptideAtlas 评测数据集构建")
    print("(正样本: GSB 中检测到的磷酸化肽段)")
    print("(负样本: 理论酶切后含 s/t/y 但未在 GSB 中检测到)")
    print("=" * 70)
    
    gsb_path = DATA_DIR / "Yeast_phospho_GSB.csv"
    fasta_path = DATA_DIR / "Saccharomyces_cerevisiae.fasta"
    
    # 1. 提取正样本
    print("\n" + "=" * 50)
    print("步骤 1: 提取 GSB 中检测到的磷酸化肽段（正样本）")
    print("=" * 50)
    
    positive_peptides, positive_set = extract_positive_from_gsb(gsb_path)
    analyze_composition(positive_peptides, "正样本")
    
    print("\n正样本示例:")
    for pep in positive_peptides[:5]:
        print(f"  {pep}")
    
    # 2. 生成负样本
    print("\n" + "=" * 50)
    print("步骤 2: 生成负样本（理论酶切后含 S/T/Y 但未检测到）")
    print("=" * 50)
    
    negative_peptides = generate_negative_peptides(fasta_path, positive_set)
    analyze_composition(negative_peptides, "负样本")
    
    print("\n负样本示例:")
    for pep in negative_peptides[:5]:
        print(f"  {pep}")
    
    # 3. 创建平衡测试集
    print("\n" + "=" * 50)
    print("步骤 3: 创建平衡测试集")
    print("=" * 50)
    
    max_samples = min(len(positive_peptides), len(negative_peptides), 5000)
    print(f"将创建 {max_samples} 正样本 + {max_samples} 负样本")
    
    pos_samples = random.sample(positive_peptides, max_samples)
    neg_samples = random.sample(negative_peptides, max_samples)
    
    data = []
    for pep in pos_samples:
        data.append({'peptide': pep, 'label': 1})
    for pep in neg_samples:
        data.append({'peptide': pep, 'label': 0})
    
    test_df = pd.DataFrame(data)
    test_df = test_df.sample(frac=1, random_state=42).reset_index(drop=True)
    
    # 保存
    output_path = OUTPUT_DIR / "peptideatlas_yeast_test.csv"
    test_df.to_csv(output_path, index=False)
    
    print("\n" + "=" * 50)
    print("数据集统计")
    print("=" * 50)
    print(f"总样本数: {len(test_df)}")
    print(f"正样本: {sum(test_df['label'] == 1)}")
    print(f"负样本: {sum(test_df['label'] == 0)}")
    print(f"\n文件已保存: {output_path}")
    
    print("\n测试集样本示例:")
    print(test_df.head(10))
    
    # 保存完整正负样本
    pd.DataFrame({
        'peptide': positive_peptides,
        'label': 1,
        'source': 'PeptideAtlas_Yeast_GSB_Detected'
    }).to_csv(OUTPUT_DIR / "peptideatlas_yeast_positive.csv", index=False)
    
    pd.DataFrame({
        'peptide': negative_peptides,
        'label': 0,
        'source': 'Theoretical_Not_Detected_STY'
    }).to_csv(OUTPUT_DIR / "peptideatlas_yeast_negative.csv", index=False)
    
    print("\n" + "=" * 70)
    print("完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
