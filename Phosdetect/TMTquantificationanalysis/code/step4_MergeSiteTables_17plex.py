#!/usr/bin/env python3
"""
Step4: 合并位点表（17 plex版本）
- 对于已经包含所有plex的表，主要是整理格式
"""

import pandas as pd
import argparse
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def main():
    parser = argparse.ArgumentParser(description='合并位点表（17 plex版本）')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', help='输出文件路径（默认: 在输入文件名后加_Merged）')
    parser.add_argument('-m', '--method', default='PhosSight', choices=['PhosSight', 'PhosphoRS', 'DeepRescore2'],
                       help='方法类型')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    if not input_path.exists():
        logging.error(f"输入文件不存在: {input_path}")
        return
    
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.parent / f"{input_path.stem}_Merged{input_path.suffix}"
    
    logging.info(f"读取输入文件: {input_path}")
    df = pd.read_csv(input_path, sep='\t')
    logging.info(f"原始行数: {len(df)}")
    
    # 对于17 plex版本，数据已经合并，只需要重命名
    logging.info(f"保存到: {output_path}")
    df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"完成！最终行数: {len(df)}")


if __name__ == "__main__":
    main()












