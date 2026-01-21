#!/usr/bin/env python3
"""
Step3: 调整位点表格式（17 plex版本）
- 规范化列名
- 确保输出格式一致
"""

import pandas as pd
import argparse
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def main():
    parser = argparse.ArgumentParser(description='调整位点表格式（17 plex版本）')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径')
    parser.add_argument('-o', '--output', help='输出文件路径（默认: 在输入文件名后加_Adjust）')
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
        output_path = input_path.parent / f"{input_path.stem}_Adjust{input_path.suffix}"
    
    logging.info(f"读取输入文件: {input_path}")
    df = pd.read_csv(input_path, sep='\t')
    logging.info(f"原始行数: {len(df)}")
    logging.info(f"原始列数: {len(df.columns)}")
    logging.info(f"原始列名: {list(df.columns)}")
    
    # 检查是否是长表格式（完整版sitequant输出）
    if 'label' in df.columns and 'quantity' in df.columns:
        logging.info("检测到长表格式，需要pivot成宽表...")
        
        # 组合Plex和label，生成新的列名
        if 'Plex' in df.columns:
            df['label_with_plex'] = df['Plex'] + '_' + df['label'].astype(str)
            logging.info(f"创建组合列 label_with_plex")
            logging.info(f"Label_with_plex示例: {df['label_with_plex'].unique()[:10]}")
        else:
            # 如果没有Plex列，尝试从basename推断
            if 'basename' in df.columns:
                # 从basename提取plex信息（如果basename包含TMT01等）
                df['label_with_plex'] = df['label'].astype(str)
                logging.warning("没有Plex列，尝试从basename推断...")
            else:
                df['label_with_plex'] = df['label'].astype(str)
        
        # 确定pivot的索引列（除了label、quantity、Plex、label_with_plex之外的所有列）
        id_vars = [col for col in df.columns if col not in ['label', 'quantity', 'Plex', 'label_with_plex']]
        
        logging.info(f"Pivot索引列数: {len(id_vars)}")
        logging.info(f"Label_with_plex唯一值数: {df['label_with_plex'].nunique()}")
        
        # Pivot: label_with_plex -> 列，quantity -> 值
        # 注意：id_vars 中可能包含 basename，这会导致同一个 Site 在不同 plex 中产生多行
        # 我们需要先按 Site 去重，然后再 pivot
        logging.info(f"Pivot前唯一Site数: {df['Site'].nunique() if 'Site' in df.columns else 'N/A'}")
        logging.info(f"Pivot前总行数: {len(df)}")
        
        # 关键修复：pivot 时只使用 Site 相关的列作为索引，不包括 basename
        # 因为 basename 在不同 plex 中可能不同，会导致同一个 Site 产生多行
        site_id_vars = ['Site', 'quality', 'AApos', 'AA', 'Modi', 'GeneID', 'Protein']
        site_id_vars = [col for col in site_id_vars if col in df.columns]
        
        logging.info(f"使用以下列作为 pivot 索引: {site_id_vars}")
        
        # 先保存 basename 的映射（每个 Site 取第一个 basename）
        basename_map = None
        if 'basename' in df.columns:
            basename_map = df.groupby('Site')['basename'].first().to_dict()
            logging.info(f"保存 basename 映射，共 {len(basename_map)} 个 Site")
        
        # Pivot: 使用 site_id_vars 作为索引，这样同一个 Site 只会产生一行
        df = df.pivot_table(
            index=site_id_vars,
            columns='label_with_plex',
            values='quantity',
            aggfunc='first'  # 如果有重复，取第一个
        ).reset_index()
        
        # 列名从MultiIndex变成普通列名
        df.columns.name = None
        
        # 恢复 basename 列
        if basename_map is not None:
            df['basename'] = df['Site'].map(basename_map)
            logging.info("已恢复 basename 列")
        
        logging.info(f"Pivot后行数: {len(df)}")
        logging.info(f"Pivot后唯一Site数: {df['Site'].nunique() if 'Site' in df.columns else 'N/A'}")
        logging.info(f"Pivot后列数: {len(df.columns)}")
        logging.info(f"Pivot后列名示例: {list(df.columns[:15])}")
        
        # 再次按 Site 去重（防止意外重复）
        if 'Site' in df.columns:
            before_dedup = len(df)
            unique_sites = df['Site'].nunique()
            if before_dedup > unique_sites:
                logging.warning(f"发现重复的 Site: {before_dedup} 行，{unique_sites} 个唯一位点")
                # 按 Site 分组，对信息列取第一个，对 TMT 列取第一个非空值
                info_cols_for_group = [col for col in ['Site', 'quality', 'AApos', 'AA', 'Modi', 'GeneID', 'Protein', 'basename'] if col in df.columns]
                # 识别所有 TMT 列（TMT01-TMT17）
                plex_list_for_group = [f'TMT{i:02d}' for i in range(1, 18)]  # TMT01-TMT17
                tmt_cols_in_df = [col for col in df.columns if any(col.startswith(f'{plex}_') for plex in plex_list_for_group)]
                
                # 对于信息列，取第一个值
                agg_dict = {col: 'first' for col in info_cols_for_group if col != 'Site'}
                # 对于 TMT 列，如果有非空值，取第一个非空值；否则取 0
                # 使用 'first' 聚合函数，pandas 会自动处理
                for col in tmt_cols_in_df:
                    agg_dict[col] = 'first'
                
                df = df.groupby('Site', as_index=False).agg(agg_dict)
                after_dedup = len(df)
                logging.info(f"按 Site 去重: {before_dedup} -> {after_dedup} 行")
            else:
                logging.info(f"没有重复的 Site，保持 {before_dedup} 行")
    
    # 列名应该已经是 "TMT01_126" 格式了（从label_with_plex来的）
    # 只需要确保格式正确，不需要重命名
    logging.info("列名应该已经是正确格式（TMT01_126等）")
    
    # 确保基本信息列存在（pivot后Plex列可能不存在，因为信息已经在列名里了）
    info_cols = ['Site', 'quality', 'AApos', 'AA', 'Modi', 'GeneID', 'Protein', 'basename']
    
    # 获取TMT通道列
    tmt_cols = []
    plex_list = [f'TMT{i:02d}' for i in range(1, 18)]  # TMT01-TMT17
    for col in df.columns:
        for plex in plex_list:
            if col.startswith(f'{plex}_') and col not in info_cols:
                tmt_cols.append(col)
                break
    
    logging.info(f"TMT通道列数: {len(tmt_cols)}")
    if len(tmt_cols) > 0:
        logging.info(f"TMT通道列示例: {tmt_cols[:10]}")
    
    # 重新排列列顺序
    existing_info_cols = [col for col in info_cols if col in df.columns]
    final_cols = existing_info_cols + sorted(tmt_cols)
    
    # 只保留存在的列
    final_cols = [col for col in final_cols if col in df.columns]
    df = df[final_cols]
    
    # 可选：只保留磷酸化位点（如果用户只关心磷酸化）
    # 注意：完整版sitequant会识别所有修饰类型，包括Carbamidomethyl和Oxidation
    # 如果需要只保留磷酸化位点，取消下面的注释
    # if 'Modi' in df.columns:
    #     before_filter = len(df)
    #     df = df[df['Modi'] == 'Phospho']
    #     after_filter = len(df)
    #     logging.info(f"过滤非磷酸化位点: {before_filter} -> {after_filter} 行")
    
    logging.info(f"保存到: {output_path}")
    df.to_csv(output_path, sep='\t', index=False)
    logging.info(f"完成！最终行数: {len(df)}, 列数: {len(df.columns)}")
    if 'Modi' in df.columns:
        logging.info(f"修饰类型统计: {df['Modi'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()

