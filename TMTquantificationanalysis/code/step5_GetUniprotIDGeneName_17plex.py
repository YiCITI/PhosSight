#!/usr/bin/env python3
"""
Step5: 获取UniprotID和GeneName（17 plex版本）
"""

import pandas as pd
import argparse
from pathlib import Path
import re
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def extract_uniprot_id(protein_str):
    """从蛋白质字符串中提取UniProt ID"""
    if pd.isna(protein_str):
        return None
    
    # 匹配 sp|XXXXX|YYYY 格式
    match = re.search(r'sp\|([A-Z0-9]+)\|', str(protein_str))
    if match:
        return match.group(1)
    
    # 匹配 tr|XXXXX|YYYY 格式
    match = re.search(r'tr\|([A-Z0-9]+)\|', str(protein_str))
    if match:
        return match.group(1)
    
    return None


def extract_gene_name(protein_str):
    """从蛋白质字符串中提取基因名"""
    if pd.isna(protein_str):
        return None
    
    # 匹配 sp|XXXXX|GENENAME_HUMAN 格式
    match = re.search(r'[sptr]+\|[A-Z0-9]+\|([A-Z0-9]+)_', str(protein_str))
    if match:
        return match.group(1)
    
    return None


def main():
    parser = argparse.ArgumentParser(description='获取UniprotID和GeneName（17 plex版本）')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', help='输出文件路径')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        logging.error(f"输入文件不存在: {input_path}")
        return
    
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.parent / f"{input_path.stem.replace('_Merged', '')}_UniprotID_GeneName.tsv"
    
    logging.info(f"读取输入文件: {input_path}")
    df = pd.read_csv(input_path, sep='\t')
    logging.info(f"原始行数: {len(df)}")
    
    # 从GeneID或Protein列提取UniprotID
    if 'GeneID' in df.columns:
        df['UniprotID'] = df['GeneID'].apply(extract_uniprot_id)
        df['GeneName'] = df['GeneID'].apply(extract_gene_name)
    elif 'Protein' in df.columns:
        df['UniprotID'] = df['Protein'].apply(extract_uniprot_id)
        df['GeneName'] = df['Protein'].apply(extract_gene_name)
    
    # 统计
    logging.info(f"成功提取UniprotID: {df['UniprotID'].notna().sum()}")
    logging.info(f"成功提取GeneName: {df['GeneName'].notna().sum()}")
    
    # 保存结果
    logging.info(f"保存到: {output_path}")
    df.to_csv(output_path, sep='\t', index=False)
    logging.info("完成！")


if __name__ == "__main__":
    main()












