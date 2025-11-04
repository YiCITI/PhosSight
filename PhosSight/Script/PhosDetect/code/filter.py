import os
import torch
from torch import nn
import numpy as np
from model import biGRU_Detect as Detect
from torch.utils.data import DataLoader, Dataset
import random
import pandas as pd
from generate_samples import digest_fasta



def set_seed(seed):
    """Set random seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_device():
    """Get the device to be used for training."""
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")

def Coding(mers):
    dic = {'Z':0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
           'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
           'V': 18, 'W': 19, 'Y': 20, 's': 21, 't':22, 'y':23} # Z is used for padding, s, t, y is phosphorylation.
    
    coded_mer = [dic.get(aa, 0) for aa in mers]  # return list instead of array
    return coded_mer

def Decoding(coded_mers):
    dic = {0: 'Z', 1: 'A', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 9: 'K',
           10: 'L', 11: 'M', 12: 'N', 13: 'P', 14: 'Q', 15: 'R', 16: 'S', 17: 'T',
           18: 'V', 19: 'W', 20: 'Y', 21: 's', 22: 't', 23: 'y'} # Z is used for padding, s, t, y is phosphorylation.
    
    decoded_mer = ''.join([dic.get(code, 'Z') for code in coded_mers])
    decoded_mer = decoded_mer.rstrip('Z')  # remove padding 'Z' at the end
    return decoded_mer

class EncodedSequenceDataset(Dataset):
    def __init__(self, sequences):
        self.sequences = torch.tensor(sequences, dtype=torch.long)

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        return self.sequences[idx]
    
def prepare_test_data(digest_result, batch_size):
    # max_len = max(len(seq) for seq in digest_result)
    max_len = 53
    padded_sequences = [seq.ljust(max_len, 'Z') for seq in digest_result]  # padding with Z
    encoded_sequences = [Coding(seq) for seq in padded_sequences]
    test_dataset = EncodedSequenceDataset(encoded_sequences)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    return test_loader

def filter(model, test_loader):
    """Filter the peptites using the trained model. Keep the positive samples."""
    model.eval()
    positive_samples = []
    with torch.no_grad():
        for sequences in test_loader:
            sequences = sequences.to(get_device())
            outputs = model(sequences)
            preds = (outputs > 0.5).float()
            # keep the positive samples
            for seq, pred in zip(sequences, preds):
                if pred.item() == 1:
                    positive_samples.append(Decoding(seq.tolist()))
    
    print(f"Number of positive samples: {len(positive_samples)}")
    
    return positive_samples

def write_to_file(positive_samples, output_file):
    """Write the filtered peptides to a file."""
    with open(output_file, 'w') as f:
        for peptide in positive_samples:
            f.write(f"{peptide}\n")
    print(f"Filtered peptides written to {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Filter peptides using a trained model.")
    parser.add_argument('-m', '--model_path', type=str, default='../model/best_model.pth', help='Path to the trained model.')
    parser.add_argument('-f', '--fasta_path', type=str, default='/home/yue/wbscy/PhosSight-main/PhosSight-main/PhosSight/UP.fasta', help='Path to the FASTA file containing peptides.')
    parser.add_argument('-o', '--output_path', type=str, default='./output.txt', help='Path to the output file.')
    parser.add_argument('-b', '--batch_size', type=int, default=128, help='Batch size for testing.')


    set_seed(42)
    device = get_device()
    print("Using device:", device)
    args = parser.parse_args()
    model_path = args.model_path
    fasta_path = args.fasta_path
    output_path = args.output_path
    batch_size = args.batch_size
    model = Detect(in_features=10, out_features=20, num_layers=2, dropout=0.4)

    model.load_state_dict(torch.load(model_path, map_location=device))

    model.to(device)
    # Read the FASTA file and perform in silico digestion
    digest_result = digest_fasta(fasta_path)
    # digest_result is sperated by '\t', only keep the first column
    digest_result = [seq.split('\t')[0] for seq in digest_result if seq.strip()]  # Remove empty sequences
    digest_result = list(set(digest_result))  # Remove duplicates
    test_loader = prepare_test_data(digest_result, batch_size=batch_size)

    # Filter the peptides using the trained model
    positive_samples = filter(model, test_loader)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    # Write the filtered peptides to a file
    write_to_file(positive_samples, output_path)
