#!/usr/bin/env python3
"""
Step2: 完整版sitequant位点定量（17 plex版本）
TMT01-TMT17共17个批次，每个批次单独调用原版sitequant，最后合并

核心逻辑：
1. 从 Intensity_17plex.txt 里按 Plex 拆分出各批次的 PSM
2. 为每个批次构造原版 sitequant 需要的 psm_* 宽表格式
3. 调用原版 sitequant (preprocess + runner + modisite_quant)
4. 合并17个批次的位点表，输出最终结果
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import warnings
from functools import partial

warnings.filterwarnings('ignore', category=FutureWarning)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 添加原版sitequant模块路径
sitequant_path = Path(__file__).parent.parent / "TMTQuantification" / "step2.sitequant" / "sitequant"
sys.path.insert(0, str(sitequant_path))

try:
    from runner import runner
    from modisite import modisite_quant
    SITEQUANT_AVAILABLE = True
except ImportError as e:
    logging.warning(f"无法导入sitequant模块: {e}")
    SITEQUANT_AVAILABLE = False


def load_fasta_simple(fasta_path: str) -> pd.DataFrame:
    """简单的FASTA解析"""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    fa = pd.DataFrame([
        {"header": h, "sequence": s, "geneid": h, "ref": h}
        for h, s in sequences.items()
    ])
    return fa


def maybe_split_on_proteinid(df):
    """处理多个蛋白质ID的情况"""
    SEP = ";"
    if "Proteins" not in df.columns:
        return df
    if not df["Proteins"].str.contains(SEP, na=False).any():
        if "GeneID" not in df.columns:
            df["GeneID"] = df["Proteins"]
        return df

    glstsplitter = (
        df["Proteins"]
        .str.split(SEP)
        .apply(pd.Series, 1)
        .stack()
        .to_frame(name="GeneID")
    )
    glstsplitter.index = glstsplitter.index.droplevel(-1)
    df = df.join(glstsplitter).reset_index(drop=True)
    df["GeneID"] = df.GeneID.fillna("-1")
    df.loc[df.GeneID == "", "GeneID"] = "-1"
    df["GeneID"] = df.GeneID.astype(str)
    return df


def remove_contaminants(df):
    """移除污染物"""
    if "GeneID" not in df.columns:
        return df
    contaminant_pattern = r"^(Contaminant_|CON__|CONT_|REV__|XXX_)"
    mask = ~df["GeneID"].str.contains(contaminant_pattern, regex=True, na=False)
    return df[mask]


def wide_to_long_table(df):
    """将宽表转换为长表（sitequant需要的格式）"""
    # 找出所有TMT通道列（去掉psm_前缀后的列名）
    tmt_cols = [col for col in df.columns if col in ['126', '127N', '127C', '128N', '128C', '129N', '129C', '130N', '130C', '131']]
    
    if not tmt_cols:
        logging.error("没有找到TMT通道列")
        return df
    
    id_vars = [col for col in df.columns if col not in tmt_cols]
    
    res = df.melt(
        id_vars=id_vars,
        value_vars=tmt_cols,
        value_name="SequenceArea",
        var_name="LabelFLAG"
    )
    return res


def prepare_plex_psm_table(df_plex, plex_name, ion_cols):
    """
    为单个plex准备原版sitequant需要的PSM表格式
    
    输入: df_plex 是已经筛选好的某个plex的PSM数据
    输出: 符合原版sitequant输入格式的DataFrame
    """
    # 通道名映射: Ion_126.128 -> psm_126, Ion_127.125 -> psm_127N, ...
    channel_map = {
        'Ion_126.128': '126',
        'Ion_127.125': '127N',
        'Ion_127.131': '127C',
        'Ion_128.128': '128N',
        'Ion_128.134': '128C',
        'Ion_129.131': '129N',
        'Ion_129.138': '129C',
        'Ion_130.135': '130N',
        'Ion_130.141': '130C',
        'Ion_131.138': '131'
    }

    def convert_modifications(raw_mod_str: str) -> str:
        """
        将 PhosSight / PhosphoRS 的修饰字符串转换为 modisite 期望的格式
        例如:
          原始: "3,Phospho[S];1,TMT6plex[AnyN-term];8,TMT6plex[K];"
          转换: "3S(Phospho)"

        只保留磷酸化等非TMT修饰, 丢弃 TMT6plex 等定量标签。
        """
        if pd.isna(raw_mod_str) or not str(raw_mod_str).strip():
            return ""

        parts = str(raw_mod_str).split(";")
        items = []
        for part in parts:
            part = part.strip()
            if not part:
                continue
            # 预期格式: "<pos>,<ModName>[<AA>]"
            # 例如: "3,Phospho[S]"
            try:
                pos_and_rest = part.split(",", 1)
                if len(pos_and_rest) != 2:
                    continue
                pos_str, rest = pos_and_rest
                pos_str = pos_str.strip()
                # 跳过TMT等标记
                if "TMT" in rest or "TMT6plex" in rest:
                    continue
                # 从 rest 中提取 ModName 和 氨基酸
                # rest 形如 "Phospho[S]" 或 "Oxidation[M]"
                mod_name = rest
                aa = ""
                if "[" in rest and "]" in rest:
                    mod_name = rest.split("[", 1)[0]
                    aa = rest.split("[", 1)[1].split("]", 1)[0]
                mod_name = mod_name.strip()
                if not mod_name or not pos_str.isdigit():
                    continue
                pos = int(pos_str)
                if aa:
                    items.append(f"{pos}{aa}({mod_name})")
                else:
                    items.append(f"{pos}({mod_name})")
            except Exception:
                continue

        return ", ".join(items)

    # 构建新的DataFrame
    result = pd.DataFrame()
    
    # 复制基本信息列
    result['Title'] = df_plex['Title']
    result['Charge'] = df_plex['Charge']
    result['Sequence'] = df_plex['Peptide']
    result['Proteins'] = df_plex['Proteins']
    
    # 处理修饰信息: 先取原始列, 再转换为 modisite 期望格式
    if 'Modifications_abbrev' in df_plex.columns:
        raw_mod = df_plex['Modifications_abbrev']
    elif 'Modification' in df_plex.columns:
        raw_mod = df_plex['Modification']
    else:
        raw_mod = ""

    result['Modifications'] = [convert_modifications(x) for x in raw_mod]

    # 处理SequenceModi（isoform 信息, 供 modisite 区分不同修饰序列）
    if 'IsoformSequence_PhosSight' in df_plex.columns:
        result['SequenceModi'] = df_plex['IsoformSequence_PhosSight']
    elif 'IsoformSequence_PhosphoRS' in df_plex.columns:
        result['SequenceModi'] = df_plex['IsoformSequence_PhosphoRS']
    elif 'IsoformSequence_DeepRescore2' in df_plex.columns:
        result['SequenceModi'] = df_plex['IsoformSequence_DeepRescore2']
    elif 'IsoformSequence' in df_plex.columns:
        result['SequenceModi'] = df_plex['IsoformSequence']
    elif 'Mod_Sequence_for_phosphoRS' in df_plex.columns:
        result['SequenceModi'] = df_plex['Mod_Sequence_for_phosphoRS']
    else:
        result['SequenceModi'] = df_plex['Peptide']
    
    # 添加PeakArea和ParentIonIntensity
    if f'{plex_name}_PeakArea' in df_plex.columns:
        result['PeakArea'] = df_plex[f'{plex_name}_PeakArea']
    else:
        result['PeakArea'] = 0
    
    if f'{plex_name}_ParentIonIntensity' in df_plex.columns:
        result['ParentIonIntensity'] = df_plex[f'{plex_name}_ParentIonIntensity']
    else:
        result['ParentIonIntensity'] = 0
    
    # 添加TMT通道列（使用原版sitequant需要的列名格式）
    for ion_col, short_name in channel_map.items():
        full_col = f'{plex_name}_{ion_col}'
        if full_col in df_plex.columns:
            result[short_name] = df_plex[full_col]
        else:
            result[short_name] = np.nan
    
    return result


def run_sitequant_for_plex(df_plex, plex_name, fasta_df, basename):
    """
    对单个plex运行完整版sitequant
    """
    logging.info(f"  运行sitequant for {plex_name}...")
    
    # 预处理
    df = df_plex.copy()
    df = maybe_split_on_proteinid(df)
    df = remove_contaminants(df)
    
    if "GeneID" not in df.columns:
        df["GeneID"] = df["Proteins"]
    
    df["GeneID"] = df["GeneID"].astype(str)
    df["Modifications"] = df["Modifications"].astype(str).fillna("")
    
    if "PSM_UseFLAG" not in df.columns:
        df["PSM_UseFLAG"] = 1
    if "IonScore" not in df.columns:
        df["IonScore"] = -1
    
    # 转换为长表
    df_long = wide_to_long_table(df)
    
    if len(df_long) == 0:
        logging.warning(f"  {plex_name}: 转换后数据为空")
        return None
    
    logging.info(f"    预处理后PSM数: {len(df)}, 长表行数: {len(df_long)}")
    
    # 运行sitequant
    unique_geneids = df_long.GeneID.unique()
    total_genes = len(unique_geneids)
    logging.info(f"    待处理蛋白质数: {total_genes}")
    
    ALL_RESULTS = []
    processed = 0
    failed_count = 0
    non_empty_df_count = 0
    
    for geneid in unique_geneids:
        try:
            res = runner(
                geneids=geneid,
                df=df_long,
                fasta=fasta_df,
                basename=f"{basename}_{plex_name}",
            )
            if res:
                # 统计每个geneid返回的非空DataFrame个数，便于调试
                for df_res in res:
                    if df_res is not None and len(df_res) > 0:
                        non_empty_df_count += 1
                ALL_RESULTS.extend(res)
            processed += 1
        except Exception as e:
            failed_count += 1
            processed += 1
    
    logging.info(
        f"    完成: 成功 {processed - failed_count}, 失败 {failed_count}, "
        f"runner返回结果对象数 {len(ALL_RESULTS)}, 其中非空DataFrame个数 {non_empty_df_count}"
    )
    
    if not ALL_RESULTS:
        logging.warning(f"    {plex_name}: ALL_RESULTS 为空")
        return None
    
    # 合并结果 - runner返回的是DataFrame列表
    # 过滤掉空的DataFrame，并增加更多调试信息
    valid_results = []
    for idx, df_res in enumerate(ALL_RESULTS):
        if df_res is None:
            continue
        rows = len(df_res)
        if rows > 0:
            valid_results.append(df_res)
        # 对前几个结果打印列信息，便于检查格式
        if idx < 3:
            logging.info(
                f"    调试: 结果[{idx}] 行数={rows}, 列数={df_res.shape[1] if df_res is not None else 'NA'}"
            )
            if rows > 0:
                logging.info(f"    调试: 结果[{idx}] 列名示例: {list(df_res.columns)[:10]}")

    logging.info(f"    有效结果数(非空DataFrame): {len(valid_results)}")
    
    if not valid_results:
        logging.warning(f"    {plex_name}: 没有有效的位点结果")
        return None
    
    result_df = pd.concat(valid_results, ignore_index=True)
    result_df['Plex'] = plex_name
    logging.info(f"    合并后位点数: {len(result_df)}")
    
    return result_df


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='完整版sitequant位点定量（17 plex版本）')
    parser.add_argument('-i', '--input', required=True, help='输入文件路径 (Intensity_17plex.txt)')
    parser.add_argument('-f', '--fasta', required=True, help='FASTA文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件路径')
    parser.add_argument('-m', '--method', default='PhosSight', choices=['PhosSight', 'PhosphoRS', 'DeepRescore2'],
                       help='方法类型')
    parser.add_argument('--plexes', default='1,2,6,7', help='要处理的plex范围，如 "1-17" 或 "1,2,6,7"')
    
    args = parser.parse_args()
    
    if not SITEQUANT_AVAILABLE:
        logging.error("sitequant模块不可用，无法运行完整版")
        return
    
    input_path = Path(args.input)
    fasta_path = Path(args.fasta)
    output_path = Path(args.output)
    
    if not input_path.exists():
        logging.error(f"输入文件不存在: {input_path}")
        return
    
    if not fasta_path.exists():
        logging.error(f"FASTA文件不存在: {fasta_path}")
        return
    
    # 解析plex范围
    if '-' in args.plexes:
        start, end = map(int, args.plexes.split('-'))
        plex_nums = list(range(start, end + 1))
    else:
        plex_nums = [int(x) for x in args.plexes.split(',')]
    
    plex_list = [f'TMT{num:02d}' for num in plex_nums]
    
    logging.info("=" * 60)
    logging.info(f"完整版sitequant位点定量（{args.method}）")
    logging.info("=" * 60)
    logging.info(f"输入文件: {input_path}")
    logging.info(f"FASTA文件: {fasta_path}")
    logging.info(f"输出文件: {output_path}")
    logging.info(f"处理批次: {', '.join(plex_list)}")
    logging.info("=" * 60)
    
    # TMT通道列名
    ion_cols = ['Ion_126.128', 'Ion_127.125', 'Ion_127.131', 'Ion_128.128', 'Ion_128.134',
                'Ion_129.131', 'Ion_129.138', 'Ion_130.135', 'Ion_130.141', 'Ion_131.138']
    
    # 读取数据
    logging.info("\n[1] 读取输入数据...")
    df = pd.read_csv(input_path, sep='\t')
    logging.info(f"  总PSM数: {len(df)}")
    
    # 加载FASTA
    logging.info("\n[2] 加载FASTA文件...")
    fasta_df = load_fasta_simple(str(fasta_path))
    logging.info(f"  FASTA序列数: {len(fasta_df)}")
    
    # 按plex处理
    logging.info("\n[3] 按批次运行sitequant...")
    
    all_site_tables = []
    basename = output_path.stem
    
    for plex in plex_list:
        logging.info(f"\n处理 {plex}...")
        
        # 筛选该plex的PSM
        # 优先检查该plex的PeakArea列是否有非NA数据
        plex_peakarea_col = f'{plex}_PeakArea'
        if plex_peakarea_col in df.columns:
            # 检查是否有非NA数据
            has_data = df[plex_peakarea_col].notna().any()
            if not has_data:
                logging.warning(f"  {plex}: PeakArea列全为NA，跳过")
                continue
            
            # 筛选有该plex数据的行
            df_plex = df[df[plex_peakarea_col].notna()].copy()
        elif 'Plex' in df.columns:
            # 使用Plex列筛选
            df_plex = df[df['Plex'] == plex].copy()
        else:
            # 根据Title前缀判断
            plex_num = int(plex[3:])
            prefix = f'{plex_num:02d}'
            df_plex = df[df['Title'].str.startswith(prefix)].copy()
        
        if len(df_plex) == 0:
            logging.warning(f"  {plex}: 没有找到PSM数据，跳过")
            continue
        
        logging.info(f"  {plex} PSM数: {len(df_plex)}")
        
        # 准备PSM表
        df_psm = prepare_plex_psm_table(df_plex, plex, ion_cols)
        
        # 运行sitequant
        site_table = run_sitequant_for_plex(df_psm, plex, fasta_df, basename)
        
        if site_table is not None and len(site_table) > 0:
            all_site_tables.append(site_table)
            logging.info(f"  {plex} 位点数: {len(site_table)}")
        else:
            logging.warning(f"  {plex}: sitequant没有产出位点")
    
    # 合并所有批次的位点表
    logging.info("\n[4] 合并所有批次的位点表...")
    
    if not all_site_tables:
        logging.error("没有任何批次产出位点，请检查输入数据格式")
        return
    
    final_df = pd.concat(all_site_tables, ignore_index=True)
    logging.info(f"  合并后总位点数: {len(final_df)}")
    
    # 保存结果
    logging.info(f"\n[5] 保存结果到: {output_path}")
    final_df.to_csv(output_path, sep='\t', index=False)
    
    logging.info("\n" + "=" * 60)
    logging.info("完成！")
    logging.info(f"  输出文件: {output_path}")
    logging.info(f"  总位点数: {len(final_df)}")
    logging.info("=" * 60)


if __name__ == "__main__":
    main()

