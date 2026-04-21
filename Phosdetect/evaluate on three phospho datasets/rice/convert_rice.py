#!/usr/bin/env python3
"""
从 Rice Phospho PeptideAtlas 构建评测数据集
- 正样本：GSB 中检测到的磷酸化肽段（磷酸化位点小写）
- 负样本：水稻蛋白组理论消化后，含 S/T/Y 的肽段但未在 GSB 中检测到

参考 Yeast 的处理方式
"""

import re
import os
import pandas as pd
from pathlib import Path
import random
import requests

random.seed(42)

# 路径配置
DATA_DIR = Path("/data0/wangb/cd/duibi0826/0826comparison/data/rice")
OUTPUT_DIR = DATA_DIR
DATA_DIR.mkdir(parents=True, exist_ok=True)

# PeptideAtlas 下载链接
GSB_URL = "https://peptideatlas.org/builds/rice/Rice_phospho_GSB_all_proteins.csv"
FASTA_URL = "https://peptideatlas.org/builds/rice/202204/Rice.fasta"

# 磷酸化氨基酸
PHOSPHO_AA = {'S', 'T', 'Y'}


def download_file(url, local_path):
    """下载文件"""
    if local_path.exists():
        print(f"  文件已存在: {local_path}")
        return True
    
    print(f"  下载: {url}")
    try:
        response = requests.get(url, timeout=300, stream=True)
        if response.status_code == 200:
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"  下载完成: {local_path}")
            return True
        else:
            print(f"  下载失败: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"  下载失败: {e}")
        return False


def trypsin_digest(protein_sequence, min_len=7, max_len=51, max_missed=2):
    """trypsin 消化蛋白质"""
    peptides = []
    cleavage_sites = [0]
    
    for i, aa in enumerate(protein_sequence):
        if aa in 'KR' and i + 1 < len(protein_sequence):
            next_aa = protein_sequence[i + 1]
            if next_aa != 'P':
                cleavage_sites.append(i + 1)
    
    cleavage_sites.append(len(protein_sequence))
    
    for i in range(len(cleavage_sites) - 1):
        start = cleavage_sites[i]
        for j in range(i + 1, min(i + 1 + max_missed + 1, len(cleavage_sites))):
            end = cleavage_sites[j]
            peptide = protein_sequence[start:end]
            if min_len <= len(peptide) <= max_len:
                peptides.append(peptide)
    
    return peptides


def parse_phosphopeptide(mod_seq):
    """从修饰序列中提取磷酸化肽段"""
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
    
    # 打印列名
    print(f"列名: {df.columns.tolist()[:10]}...")
    
    # 尝试找到正确的列
    peptide_cols = [col for col in df.columns if '_peptide_mod_pos' in col.lower() or 'peptide' in col.lower()]
    print(f"肽段相关列: {peptide_cols[:5]}")
    
    # 使用所有可能的肽段列
    if not peptide_cols:
        # 尝试查找包含修饰信息的列
        mod_cols = [col for col in df.columns if 'mod' in col.lower()]
        print(f"修饰相关列: {mod_cols[:5]}")
        peptide_cols = mod_cols
    
    # 提取磷酸化肽段
    phosphopeptides = {}
    for col in peptide_cols:
        for mod_seq in df[col]:
            peptide, num_phos = parse_phosphopeptide(mod_seq)
            if peptide and min_len <= len(peptide) <= max_len and num_phos > 0:
                if peptide not in phosphopeptides:
                    phosphopeptides[peptide] = True
    
    print(f"提取的磷酸化肽段数 (去重): {len(phosphopeptides)}")
    return list(phosphopeptides.keys()), phosphopeptides


def generate_negative_peptides(fasta_path, positive_set):
    """从蛋白组生成负样本"""
    print(f"\n读取蛋白组: {fasta_path}")
    
    negative_peptides = {}
    
    with open(fasta_path, 'r') as f:
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if sequence:
                    peptides = trypsin_digest(sequence)
                    for pep in peptides:
                        if 7 <= len(pep) <= 51:
                            for i, aa in enumerate(pep):
                                if aa in PHOSPHO_AA:
                                    phos_pep = pep[:i] + aa.lower() + pep[i+1:]
                                    if phos_pep not in positive_set:
                                        negative_peptides[phos_pep] = True
                sequence = ''
            else:
                sequence += line.strip()
        
        if sequence:
            peptides = trypsin_digest(sequence)
            for pep in peptides:
                if 7 <= len(pep) <= 51:
                    for i, aa in enumerate(pep):
                        if aa in PHOSPHO_AA:
                            phos_pep = pep[:i] + aa.lower() + pep[i+1:]
                            if phos_pep not in positive_set:
                                negative_peptides[phos_pep] = True
    
    print(f"负样本候选数: {len(negative_peptides)}")
    return list(negative_peptides.keys())


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
    print("Rice Phospho PeptideAtlas 评测数据集构建")
    print("=" * 70)
    
    gsb_path = DATA_DIR / "Rice_phospho_GSB.csv"
    fasta_path = DATA_DIR / "Rice.fasta"
    
    # 下载文件
    print("\n下载数据文件...")
    if not download_file(GSB_URL, gsb_path):
        print("下载 GSB 数据失败!")
        return
    
    if not download_file(FASTA_URL, fasta_path):
        print("下载 FASTA 文件失败!")
        return
    
    # 1. 提取正样本
    print("\n" + "=" * 50)
    print("步骤 1: 提取磷酸化肽段（正样本）")
    print("=" * 50)
    
    positive_peptides, positive_set = extract_positive_from_gsb(gsb_path)
    analyze_composition(positive_peptides, "正样本")
    
    print("\n正样本示例:")
    for pep in positive_peptides[:5]:
        print(f"  {pep}")
    
    if len(positive_peptides) < 50:
        print("\n正样本太少，尝试其他列...")
        # 尝试直接读取原始数据
        df = pd.read_csv(gsb_path, low_memory=False, nrows=100)
        print(f"前100行列名: {df.columns.tolist()}")
    
    # 2. 生成负样本
    if len(positive_peptides) >= 50:
        print("\n" + "=" * 50)
        print("步骤 2: 生成负样本")
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
        output_path = OUTPUT_DIR / "rice_test.csv"
        test_df.to_csv(output_path, index=False)
        
        print(f"\n数据集统计:")
        print(f"  总样本数: {len(test_df)}")
        print(f"  正样本: {sum(test_df['label'] == 1)}")
        print(f"  负样本: {sum(test_df['label'] == 0)}")
        print(f"  文件已保存: {output_path}")
        
        # 保存完整正负样本
        pd.DataFrame({
            'peptide': positive_peptides,
            'label': 1
        }).to_csv(OUTPUT_DIR / "rice_positive.csv", index=False)
        
        pd.DataFrame({
            'peptide': negative_peptides[:100000],
            'label': 0
        }).to_csv(OUTPUT_DIR / "rice_negative.csv", index=False)
    
    print("\n" + "=" * 70)
    print("完成!")
    print("=" * 70)


if __name__ == "__main__":
    main()
