import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import argparse
import os
import sys

# 导入必要的模块
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from read_fasta import read_fasta
from in_silico_digestion import digestion
from model import biGRU_Detect

# 设置设备
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# 编码函数
def encode_sequence(seq):
    dic = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
           'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
           'V': 18, 'W': 19, 'Y': 20, 'Z': 0}  # Z for padding
    return [dic.get(aa, 0) for aa in seq]

# 预测函数
def predict(model, sequences):
    model.eval()
    encoded_seqs = [encode_sequence(seq) for seq in sequences]
    X = torch.tensor(encoded_seqs, dtype=torch.long).to(device)
    dataset = TensorDataset(X)
    loader = DataLoader(dataset, batch_size=128, shuffle=False)
    predictions = []
    with torch.no_grad():
        for batch in loader:
            inputs = batch[0].to(device)
            outputs = model(inputs)
            predictions.extend(outputs.cpu().numpy().flatten().tolist())
    return predictions

def main():
    parser = argparse.ArgumentParser(description='Predict detectability of peptides from FASTA file.')
    parser.add_argument('--fasta', type=str, required=True, help='Path to the input FASTA file')
    parser.add_argument('--model', type=str, default='best_model.pth', help='Path to the trained model file')
    parser.add_argument('--output', type=str, default='peptide_predictions.txt', help='Path to the output file')
    parser.add_argument('--enzyme', type=str, default='Trypsin', help='Enzyme for digestion')
    parser.add_argument('--missed_cleavages', type=int, default=2, help='Number of allowed missed cleavages')
    parser.add_argument('--min_len', type=int, default=7, help='Minimum peptide length')
    parser.add_argument('--max_len', type=int, default=47, help='Maximum peptide length')
    args = parser.parse_args()

    # 加载模型
    model = biGRU_Detect(in_features=24, out_features=64, num_layers=2, dropout=0.2).to(device)
    state_dict = torch.load(args.model, map_location=device)
    model.load_state_dict(state_dict)

    # 读取 FASTA 文件
    sequences = read_fasta(args.fasta)
    print(f"Read {len(sequences)} protein sequences from {args.fasta}")

    # 模拟消化生成肽段
    all_peptides = []
    all_mers = []
    for seq_id, seq in sequences:
        if args.enzyme == 'Trypsin':
            sites = [-1] + [i for i, aa in enumerate(seq) if aa in 'KR' and i < len(seq) - 1 and seq[i + 1] != 'P'] + [len(seq) - 1]
        else:
            raise ValueError(f"Unsupported enzyme: {args.enzyme}")
        peptides, mers = digestion(seq, sites, 'C', args.missed_cleavages, args.min_len, args.max_len)
        all_peptides.extend([(seq_id, pep) for pep in peptides])
        all_mers.extend(mers)

    print(f"Generated {len(all_peptides)} peptides")

    # 预测可检测性概率
    probabilities = predict(model, all_mers)

    # 输出结果
    with open(args.output, 'w') as f:
        f.write("Sequence_ID\tPeptide\tDetectability_Probability\n")
        for (seq_id, peptide), prob in zip(all_peptides, probabilities):
            f.write(f"{seq_id}\t{peptide}\t{prob:.4f}\n")

    print(f"Predictions saved to {args.output}")

if __name__ == "__main__":
    main() 