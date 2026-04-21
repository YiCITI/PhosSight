#!/usr/bin/env python3
"""Mouse 数据下载和处理"""
import re, os, pandas as pd, random, requests
from pathlib import Path

random.seed(42)

DATA_DIR = Path("/data0/wangb/cd/duibi0826/0826comparison/data/mouse")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# 下载 GSB 数据
gsb_path = DATA_DIR / "Mouse_phospho_GSB.csv"
if not gsb_path.exists():
    print("下载 Mouse_phospho_GSB.csv...")
    url = "https://peptideatlas.org/builds/mouse/phospho/Mouse_phospho_GSB_all_proteins.csv"
    r = requests.get(url, timeout=600, stream=True)
    if r.status_code == 200:
        with open(gsb_path, 'wb') as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
        print(f"下载完成: {gsb_path}")
    else:
        print(f"下载失败: {r.status_code}")

# 下载 FASTA
fasta_path = DATA_DIR / "Mus_musculus.fasta"
if not fasta_path.exists():
    print("下载 Mus_musculus.fasta...")
    url = "https://peptideatlas.org/builds/mouse/phospho/202402/Mus_musculus.fasta"
    r = requests.get(url, timeout=600, stream=True)
    if r.status_code == 200:
        with open(fasta_path, 'wb') as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
        print(f"下载完成: {fasta_path}")
    else:
        print(f"下载失败: {r.status_code}")

# 检查文件
print(f"\n文件列表:")
for f in DATA_DIR.iterdir():
    size_mb = f.stat().st_size / 1024 / 1024
    print(f"  {f.name}: {size_mb:.1f} MB")

# 读取 GSB 数据
print(f"\n读取 GSB 数据...")
try:
    df = pd.read_csv(gsb_path, low_memory=False)
    print(f"列名: {df.columns.tolist()[:10]}")
    print(f"行数: {len(df)}")
except Exception as e:
    print(f"读取失败: {e}")
