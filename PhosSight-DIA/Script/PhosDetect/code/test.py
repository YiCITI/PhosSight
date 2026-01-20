import time
import os
import torch
from torch import nn
import numpy as np
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, precision_score, recall_score
import random
from model import biGRU_Detect as Detect
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
import pandas as pd
import argparse
from datetime import datetime

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

class EncodedSequenceDataset(Dataset):
    def __init__(self, sequences, labels):
        self.sequences = torch.tensor(sequences, dtype=torch.long)
        self.labels = torch.tensor(labels, dtype=torch.float)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.sequences[idx], self.labels[idx]
    
def prepare_test_data(path, batch_size, max_len=51):
    data = pd.read_csv(path, header=None, sep=",")
    # the max len is set as the maximum length of the first column
    max_len = max(data[0].apply(len))
    data[0] = data[0].apply(lambda x: x.ljust(max_len, 'Z')) # padding with Z
    test_dataset = EncodedSequenceDataset(data[0].apply(lambda x: Coding(x)).tolist(), data[1].tolist())
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    return test_loader

def test(model, test_loader, device):
    model.eval()
    all_labels = []
    all_preds = []
    all_outputs = []
    total_loss = 0.0
    criterion = nn.BCELoss(reduction='sum')  # 注意这里用sum，为了累加loss

    with torch.no_grad():
        for sequences, labels in test_loader:
            sequences = sequences.to(device)
            labels = labels.to(device)

            outputs = model(sequences).squeeze()

            total_loss += criterion(outputs, labels.float()).item()  # 直接加总loss

            all_labels.append(labels.cpu())
            all_outputs.append(outputs.cpu())
            preds = (outputs > 0.5).float()
            all_preds.append(preds.cpu())

    # 整合所有batch的数据
    all_labels = torch.cat(all_labels).numpy()
    all_preds = torch.cat(all_preds).numpy()
    all_outputs = torch.cat(all_outputs).numpy()

    # 总loss除以样本数
    avg_loss = total_loss / len(all_labels)
    acc = accuracy_score(all_labels, all_preds)
    f1 = f1_score(all_labels, all_preds)
    precision = precision_score(all_labels, all_preds)
    recall = recall_score(all_labels, all_preds)

    try:
        auc = roc_auc_score(all_labels, all_outputs)
    except ValueError:
        auc = float('nan')  # 如果只有一个类别，AUC无法计算

    return avg_loss, acc, f1, precision, recall, auc

def save_results_to_file(output_path, model_path, data_path, test_loss, test_acc, test_f1, test_precision, test_recall, test_auc, device_info):
    """Save test results to a text file."""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("PhosSight Model Test Results\n")
        f.write("=" * 80 + "\n")
        f.write(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Device: {device_info}\n")
        f.write(f"Model Path: {model_path}\n")
        f.write(f"Data Path: {data_path}\n")
        f.write("-" * 80 + "\n")
        f.write("Test Metrics:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Loss: {test_loss:.6f}\n")
        f.write(f"Accuracy: {test_acc:.6f}\n")
        f.write(f"F1 Score: {test_f1:.6f}\n")
        f.write(f"Precision: {test_precision:.6f}\n")
        f.write(f"Recall: {test_recall:.6f}\n")
        f.write(f"AUC: {test_auc:.6f}\n")
        f.write("-" * 80 + "\n")
        f.write("Summary:\n")
        f.write(f"The model achieved {test_acc:.2%} accuracy and {test_auc:.4f} AUC on the test dataset.\n")
        f.write("=" * 80 + "\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Test PhosSight model and save results.")
    parser.add_argument('-m', '--model_path', type=str, default='../model/best_model.pth', 
                       help='Path to the trained model weights.')
    parser.add_argument('-p', '--path', type=str, default='../data/test/DeepDetect_human.csv', 
                       help='Path to the dataset CSV file.')
    parser.add_argument('-o', '--output', type=str, default='./log/log_phosight_test_balanced_dataset.txt',
                       help='Path to save the test results (txt file).')
    parser.add_argument('-b', '--batch_size', type=int, default=64,
                       help='Batch size for testing (default: 64).')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42).')
    return parser.parse_args()

def main():
    args = parse_args() 

    # Set random seed
    set_seed(args.seed)
    
    # Get device
    device = get_device()
    device_info = f"{device}"
    if torch.cuda.is_available():
        device_info += f" ({torch.cuda.get_device_name(0)})"
    
    print("Using device:", device_info)
    print(f"Model path: {args.model_path}")
    print(f"Data path: {args.path}")
    print(f"Output path: {args.output}")
    print(f"Batch size: {args.batch_size}")

    # Load model
    print("Loading model...")
    model = Detect(in_features=10, out_features=20, num_layers=2, dropout=0)
    model.load_state_dict(torch.load(args.model_path, map_location=device))
    model.to(device)
    print("Model loaded successfully!")

    # Prepare test data
    print("Preparing test data...")
    test_loader = prepare_test_data(args.path, batch_size=args.batch_size)
    print(f"Test data prepared: {len(test_loader.dataset)} samples")

    # Run test
    print("Running test...")
    start_time = time.time()
    test_loss, test_acc, test_f1, test_precision, test_recall, test_auc = test(model, test_loader, device)
    test_time = time.time() - start_time

    # Print results
    print("\n" + "="*60)
    print("TEST RESULTS")
    print("="*60)
    print(f"Test Loss: {test_loss:.6f}")
    print(f"Test Accuracy: {test_acc:.6f}")
    print(f"Test F1 Score: {test_f1:.6f}")
    print(f"Test Precision: {test_precision:.6f}")
    print(f"Test Recall: {test_recall:.6f}")
    print(f"Test AUC: {test_auc:.6f}")
    print(f"Test Time: {test_time:.2f} seconds")
    print("="*60)

    # Save results to file
    print(f"Saving results to: {args.output}")
    save_results_to_file(
        args.output, 
        args.model_path, 
        args.path, 
        test_loss, 
        test_acc, 
        test_f1, 
        test_precision, 
        test_recall, 
        test_auc, 
        device_info
    )
    print("Results saved successfully!")

if __name__ == "__main__":
    main()