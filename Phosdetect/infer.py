import os
import sys
import argparse
import csv
from typing import List, Tuple

import torch
import numpy as np

from model import biGRU_Detect_Improved_V2


AA_TO_ID = {
    'Z': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
    'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
    'V': 18, 'W': 19, 'Y': 20, 's': 21, 't': 22, 'y': 23
}


def pad_and_encode(sequences: List[str]) -> Tuple[torch.Tensor, List[int]]:
    max_len = max(len(s) for s in sequences)
    encoded: List[List[int]] = []
    lengths: List[int] = []
    for seq in sequences:
        lengths.append(len(seq))
        padded = seq.ljust(max_len, 'Z')
        encoded.append([AA_TO_ID.get(ch, 0) for ch in padded])
    return torch.tensor(encoded, dtype=torch.long), lengths


def load_model(model_path: str, in_features: int, out_features: int, num_layers: int, dropout: float, device: torch.device):
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


def infer(model, sequences: List[str], batch_size: int, device: torch.device) -> np.ndarray:
    probs: List[float] = []
    model.eval()
    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i + batch_size]
            inputs, _ = pad_and_encode(batch)
            inputs = inputs.to(device)
            outputs = model(inputs)
            if outputs.dim() > 1 and outputs.size(1) == 1:
                outputs = outputs.squeeze(1)
            probs.extend(outputs.detach().cpu().numpy().tolist())
    return np.array(probs)


def read_sequences_from_file(path: str) -> List[str]:
    sequences: List[str] = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            sequences.append(line)
    return sequences


def write_results(output_path: str, sequences: List[str], probs: np.ndarray) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["sequence", "probability"])
        for s, p in zip(sequences, probs):
            writer.writerow([s, float(p)])


def parse_args():
    parser = argparse.ArgumentParser(description="PhosSight V2 inference on custom sequences")
    parser.add_argument('--model_path', type=str, required=True, help='Path to checkpoint, e.g., models/xxxx/best_model.pth')
    parser.add_argument('--seq', type=str, nargs='*', default=None, help='Sequences provided directly on CLI')
    parser.add_argument('--seq_file', type=str, default=None, help='Path to a text file of sequences (one per line)')
    parser.add_argument('--batch_size', type=int, default=128, help='Batch size for inference')
    # Model hyperparameters used in training
    parser.add_argument('--in_features', type=int, default=10, help='Embedding size used in training')
    parser.add_argument('--out_features', type=int, default=20, help='GRU hidden size per direction used in training')
    parser.add_argument('--num_layers', type=int, default=2, help='Number of GRU layers used in training')
    parser.add_argument('--dropout', type=float, default=0.3, help='Dropout probability used in training')
    parser.add_argument('--output_csv', type=str, default=None, help='Optional: path to save results CSV')
    return parser.parse_args()


def main():
    args = parse_args()

    # Collect sequences
    sequences: List[str] = []
    if args.seq:
        sequences.extend(args.seq)
    if args.seq_file:
        if not os.path.exists(args.seq_file):
            print(f"Error: seq_file not found: {args.seq_file}")
            sys.exit(1)
        sequences.extend(read_sequences_from_file(args.seq_file))

    if not sequences:
        print("No sequences provided. Use --seq or --seq_file.")
        sys.exit(1)

    # Validate characters
    valid_chars = set(AA_TO_ID.keys())
    for s in sequences:
        invalid = [ch for ch in s if ch not in valid_chars]
        if invalid:
            print(f"Warning: sequence contains unsupported chars {set(invalid)}; they will be mapped to 'Z' if encountered.")

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    model = load_model(
        model_path=args.model_path,
        in_features=args.in_features,
        out_features=args.out_features,
        num_layers=args.num_layers,
        dropout=args.dropout,
        device=device,
    )

    probs = infer(model, sequences, args.batch_size, device)

    # Print to stdout
    for s, p in zip(sequences, probs):
        print(f"{s}\t{p:.6f}")

    # Optional CSV
    if args.output_csv:
        write_results(args.output_csv, sequences, probs)
        print(f"Saved results to {args.output_csv}")


if __name__ == '__main__':
    main()

