#!/usr/bin/env python3
"""
从FASTA文件中读取肽段序列，使用PhosSight模型进行可检测性打分
"""

import os
import sys
import argparse
from pathlib import Path
import time
from typing import List, Tuple, TYPE_CHECKING, Any

if TYPE_CHECKING:
    import torch  # noqa: F401
    import numpy as np  # noqa: F401
    from model import biGRU_Detect_Improved_V2  # noqa: F401

# 氨基酸到ID的映射
AA_TO_ID = {
    'Z': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
    'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
    'V': 18, 'W': 19, 'Y': 20, 's': 21, 't': 22, 'y': 23
}

def pad_and_encode(sequences: List[str]) -> Tuple["torch.Tensor", List[int]]:
    """将序列填充并编码为张量"""
    max_len = max(len(s) for s in sequences)
    encoded: List[List[int]] = []
    lengths: List[int] = []
    for seq in sequences:
        lengths.append(len(seq))
        padded = seq.ljust(max_len, 'Z')
        encoded.append([AA_TO_ID.get(ch, 0) for ch in padded])
    return torch.tensor(encoded, dtype=torch.long), lengths

def load_model(model_path: str, in_features: int, out_features: int, num_layers: int, dropout: float, device: "torch.device"):
    """加载训练好的模型"""
    model = biGRU_Detect_Improved_V2(
        in_features=in_features,
        out_features=out_features,
        num_layers=num_layers,
        dropout=dropout,
    )
    state = torch.load(model_path, map_location=device)
    model.load_state_dict(state)
    model.to(device)
    model.eval()
    return model

def infer(model: Any, sequences: List[str], batch_size: int, device: "torch.device") -> "np.ndarray":
    """使用模型进行推理"""
    probs: List[float] = []
    model.eval()
    total_batches = (len(sequences) + batch_size - 1) // batch_size
    
    print(f"开始推理，总共 {len(sequences)} 条序列，分为 {total_batches} 个批次")
    
    start_time = time.time()
    
    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch_idx = i // batch_size + 1
            batch = sequences[i:i + batch_size]
            
            # 显示进度和时间估算
            progress = (batch_idx / total_batches) * 100
            elapsed_time = time.time() - start_time
            
            if batch_idx > 1:
                avg_time_per_batch = elapsed_time / (batch_idx - 1)
                remaining_batches = total_batches - batch_idx + 1
                estimated_remaining_time = remaining_batches * avg_time_per_batch
                print(f"处理批次 {batch_idx}/{total_batches} ({progress:.1f}%) - 序列 {i+1}-{min(i+batch_size, len(sequences))}")
                print(f"  已用时: {elapsed_time:.1f}s, 预计剩余: {estimated_remaining_time:.1f}s")
            else:
                print(f"处理批次 {batch_idx}/{total_batches} ({progress:.1f}%) - 序列 {i+1}-{min(i+batch_size, len(sequences))}")
            
            inputs, _ = pad_and_encode(batch)
            inputs = inputs.to(device)
            outputs = model(inputs)
            if outputs.dim() > 1 and outputs.size(1) == 1:
                outputs = outputs.squeeze(1)
            probs.extend(outputs.detach().cpu().numpy().tolist())
            
            # 每10个批次显示一次详细进度
            if batch_idx % 10 == 0 or batch_idx == total_batches:
                avg_score = np.mean(probs[-len(batch):])
                print(f"  当前批次平均分数: {avg_score:.6f}")
    
    total_time = time.time() - start_time
    print(f"推理完成！总共处理了 {len(probs)} 条序列，用时 {total_time:.1f} 秒")
    print(f"平均处理速度: {len(probs)/total_time:.1f} 序列/秒")
    return np.array(probs)

def read_fasta_sequences(file_path: str) -> List[str]:
    """从FASTA文件中读取肽段序列"""
    sequences = []
    print(f"正在读取FASTA文件: {file_path}")
    
    try:
        record_count = 0
        valid_sequences = 0
        current_sequence = ""
        
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('>'):  # 序列标识符行
                    # 如果之前有序列，先保存
                    if current_sequence:
                        sequences.append(current_sequence)
                        valid_sequences += 1
                        current_sequence = ""
                    
                    record_count += 1
                    
                    # 显示读取进度
                    if record_count % 100000 == 0:
                        print(f"读取进度: {record_count} 条记录 - 有效序列: {valid_sequences}")
                
                else:  # 序列行
                    if line:  # 非空行
                        current_sequence += line
        
        # 处理最后一个序列
        if current_sequence:
            sequences.append(current_sequence)
            valid_sequences += 1
        
        print(f"FASTA文件读取完成！")
        print(f"总记录数: {record_count}")
        print(f"有效肽段序列: {len(sequences)}")
        return sequences
        
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        return []
    except Exception as e:
        print(f"读取FASTA文件时出错: {str(e)}")
        return []
    
def _build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="从FASTA文件中读取肽段序列，使用PhosSight模型进行可检测性打分",
    )
    parser.add_argument(
        "-w",
        "--weight",
        required=True,
        help="模型权重路径（.pth 文件）",
    )
    parser.add_argument(
        "-p",
        "--peptides",
        required=True,
        help="输入FASTA文件路径（肽段序列）",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="输出打分结果文件路径（tsv）",
    )
    return parser

def main(argv: List[str] | None = None):
    """主函数"""
    args = _build_argparser().parse_args(argv)

    # 延迟导入依赖，确保 `-h/--help` 在未装依赖时也可用
    try:
        import torch  # type: ignore
        import numpy as np  # type: ignore
        from model import biGRU_Detect_Improved_V2  # type: ignore
    except ModuleNotFoundError as e:
        missing = getattr(e, "name", None) or str(e)
        print(f"错误: 缺少依赖模块 `{missing}`。请先安装运行所需依赖（例如 torch、numpy）。")
        raise

    globals()["torch"] = torch
    globals()["np"] = np
    globals()["biGRU_Detect_Improved_V2"] = biGRU_Detect_Improved_V2

    # 配置参数（由命令行参数提供）
    model_path = args.weight
    assert os.path.exists(model_path), f"模型文件不存在: {model_path}"
    fasta_file = args.peptides
    output_file = args.output
    
    # 模型超参数（根据训练时使用的参数）
    in_features = 10
    out_features = 20
    num_layers = 2
    dropout = 0.3
    batch_size = 128
    
    print("=== PhosSight FASTA肽段可检测性打分 ===")
    print(f"模型路径: {model_path}")
    print(f"FASTA文件: {fasta_file}")
    print(f"输出文件: {output_file}")
    
    # 检查文件是否存在
    if not os.path.exists(model_path):
        print(f"错误: 模型文件不存在 {model_path}")
        return
    
    if not os.path.exists(fasta_file):
        print(f"错误: FASTA文件不存在 {fasta_file}")
        return
    
    # 读取肽段序列
    sequences = read_fasta_sequences(fasta_file)
    if not sequences:
        print("没有读取到有效的肽段序列")
        return
    
    # 验证序列字符
    valid_chars = set(AA_TO_ID.keys())
    invalid_sequences = []
    for i, seq in enumerate(sequences):
        invalid_chars = [ch for ch in seq if ch not in valid_chars]
        if invalid_chars:
            invalid_sequences.append((i, seq, invalid_chars))
    
    if invalid_sequences:
        print(f"警告: 发现 {len(invalid_sequences)} 条序列包含不支持字符")
        print("前5个示例:")
        for i, (idx, seq, chars) in enumerate(invalid_sequences[:5]):
            print(f"  序列 {idx}: {seq} (无效字符: {set(chars)})")
        print("这些序列将被跳过...")
    
    # 过滤掉包含无效字符的序列
    print("正在验证序列字符...")
    valid_sequences = []
    valid_indices = []
    
    for i, seq in enumerate(sequences):
        if i % 100000 == 0 and i > 0:
            print(f"验证进度: {i}/{len(sequences)} ({i/len(sequences)*100:.1f}%)")
        
        if all(ch in valid_chars for ch in seq):
            valid_sequences.append(seq)
            valid_indices.append(i)
    
    print(f"序列验证完成！")
    print(f"总序列数: {len(sequences)}")
    print(f"有效序列数量: {len(valid_sequences)}")
    print(f"无效序列数量: {len(sequences) - len(valid_sequences)}")
    
    if not valid_sequences:
        print("没有有效的序列可以处理")
        return
    
    # 设置设备
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"使用设备: {device}")
    
    # 加载模型
    print("正在加载模型...")
    try:
        model = load_model(
            model_path=model_path,
            in_features=in_features,
            out_features=out_features,
            num_layers=num_layers,
            dropout=dropout,
            device=device,
        )
        print("模型加载成功")
    except Exception as e:
        print(f"加载模型时出错: {str(e)}")
        return
    
    # 进行推理
    print("正在进行可检测性预测...")
    try:
        probs = infer(model, valid_sequences, batch_size, device)
        print(f"预测完成，共处理 {len(probs)} 条序列")
    except Exception as e:
        print(f"推理过程中出错: {str(e)}")
        return
    
    # 保存结果
    print("正在保存结果...")
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            # Output format must match `generate_pep_fasta/process_fasta/filter_pep_fasta.py`,
            # which parses each line as: peptide<comma-or-tab>score (no header required).
            for i, (seq, score) in enumerate(zip(valid_sequences, probs)):
                if i % 50000 == 0 and i > 0:
                    progress = (i / len(valid_sequences)) * 100
                    print(f"保存进度: {i}/{len(valid_sequences)} ({progress:.1f}%)")
                
                # Use comma-separated to be unambiguous.
                f.write(f"{seq},{score:.6f}\n")
        
        print(f"结果保存完成！")
        print(f"结果已保存到: {output_file}")
        
        # 打印统计信息
        print(f"\n=== 预测结果统计 ===")
        print(f"平均可检测性分数: {np.mean(probs):.6f}")
        print(f"最高分数: {np.max(probs):.6f}")
        print(f"最低分数: {np.min(probs):.6f}")
        print(f"标准差: {np.std(probs):.6f}")
        
        # 显示前10个结果
        print(f"\n前10个结果:")
        print("序列\t\t\t可检测性分数")
        print("-" * 40)
        for i in range(min(10, len(valid_sequences))):
            seq = valid_sequences[i]
            score = probs[i]
            print(f"{seq[:20]:<20}\t{score:.6f}")
        
    except Exception as e:
        print(f"保存结果时出错: {str(e)}")
        return
    
    print("\n处理完成！")

if __name__ == "__main__":
    main()
