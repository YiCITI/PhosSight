#!/usr/bin/env python3
"""Mouse 磷酸化数据处理，生成测试集"""
import re, pandas as pd, random
from pathlib import Path

random.seed(42)

DATA_DIR = Path("/data0/wangb/cd/duibi0826/0826comparison/data/mouse")
PHOSPHO_AA = {'S', 'T', 'Y'}

def trypsin_digest(protein_sequence, min_len=7, max_len=51, max_missed=2):
    """胰蛋白酶消化"""
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
    """解析磷酸化修饰序列"""
    if pd.isna(mod_seq) or not mod_seq:
        return None, 0
    mod_seq = str(mod_seq)
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
    clean_seq = re.sub(r'\[.*?\]', '', mod_seq)
    clean_seq = re.sub(r'-\d+$', '', clean_seq)
    clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', clean_seq)
    if not clean_seq:
        return None, 0
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

# 读取 GSB 数据
print("读取 Mouse_phospho_GSB.csv...")
df = pd.read_csv(DATA_DIR / "Mouse_phospho_GSB.csv", low_memory=False)
print(f"总记录数: {len(df)}")
print(f"列名: {df.columns.tolist()[:10]}")

# 提取磷酸化肽段
phosphopeptides = {}
peptide_cols = [col for col in df.columns if 'peptide' in col.lower()]
print(f"肽段列: {peptide_cols[:5]}")

for col in peptide_cols:
    for mod_seq in df[col]:
        peptide, num_phos = parse_phosphopeptide(mod_seq)
        if peptide and 7 <= len(peptide) <= 51 and num_phos > 0:
            phosphopeptides[peptide] = True

print(f"磷酸化肽段数: {len(phosphopeptides)}")

positive_peptides = list(phosphopeptides.keys())
positive_set = phosphopeptides

# 生成负样本
negative_peptides = {}
print("\n生成负样本...")
with open(DATA_DIR / "Mus_musculus.fasta", 'r') as f:
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

negative_peptides = list(negative_peptides.keys())
print(f"负样本数: {len(negative_peptides)}")

# 创建测试集
max_samples = min(len(positive_peptides), len(negative_peptides), 5000)
pos_samples = random.sample(positive_peptides, max_samples)
neg_samples = random.sample(negative_peptides, max_samples)

data = [(pep, 1) for pep in pos_samples] + [(pep, 0) for pep in neg_samples]
random.shuffle(data)

test_df = pd.DataFrame(data, columns=['peptide', 'label'])
test_df.to_csv(DATA_DIR / "mouse_test.csv", index=False)
print(f"\n测试集: {len(test_df)} (正: {max_samples}, 负: {max_samples})")

# 统计磷酸化位点类型
pSer = sum(1 for p in test_df[test_df['label']==1]['peptide'] if 's' in p.lower())
pThr = sum(1 for p in test_df[test_df['label']==1]['peptide'] if 't' in p.lower())
pTyr = sum(1 for p in test_df[test_df['label']==1]['peptide'] if 'y' in p.lower())
print(f"pSer: {pSer}, pThr: {pThr}, pTyr: {pTyr}")

pd.DataFrame({'peptide': positive_peptides, 'label': 1}).to_csv(DATA_DIR / "mouse_positive.csv", index=False)
pd.DataFrame({'peptide': negative_peptides[:100000], 'label': 0}).to_csv(DATA_DIR / "mouse_negative.csv", index=False)
print("完成!")
