import time
import os
import torch
import numpy as np
import random
from model import biGRU_Detect as Detect
from torch.utils.data import DataLoader, Dataset
import pandas as pd
import argparse
from datetime import datetime
from tqdm import tqdm
def set_seed(seed):
    """Set random seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def Coding(mers):
    dic = {'Z':0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
           'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
           'V': 18, 'W': 19, 'Y': 20, 's': 21, 't':22, 'y':23} # Z is used for padding, s, t, y is phosphorylation.
    
    coded_mer = [dic.get(aa, 0) for aa in mers]  # return list instead of array
    return coded_mer

class EncodedSequenceDataset(Dataset):
    def __init__(self, sequences):
        self.sequences = torch.tensor(sequences, dtype=torch.long)

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        return self.sequences[idx]
    
def prepare_prediction_data(path, batch_size, max_len=53):
    """Prepare data for prediction - only sequences, no labels."""
    data = pd.read_csv(path, sep='\t')
    # the max len is set as the maximum length of the first column
    peptides = data['peptide']
    # max_len = max(peptides.apply(len))
    peptides_padding = peptides.apply(lambda x: x.ljust(max_len, 'Z')) # padding with Z
    prediction_dataset = EncodedSequenceDataset(peptides_padding.apply(lambda x: Coding(x)).tolist())
    prediction_loader = DataLoader(prediction_dataset, batch_size=batch_size, shuffle=False)
    return prediction_loader, peptides.tolist(), peptides_padding.tolist()  # Return both loader and original sequences

def predict(model, prediction_loader, device):
    """Run prediction on the data."""
    model.eval()
    all_predictions = []
    
    with torch.no_grad():
        for sequences in tqdm(prediction_loader, desc='Predicting'):
            sequences = sequences.to(device)
            outputs = model(sequences).squeeze()
            
            # Ensure outputs is 1D
            if outputs.dim() == 0:
                outputs = outputs.unsqueeze(0)
            
            all_predictions.append(outputs.cpu())
    
    # Concatenate all predictions
    all_predictions = torch.cat(all_predictions).numpy()
    return all_predictions

def save_predictions_to_csv(output_path, sequences, predictions):
    """Save predictions to CSV file with two columns: sequence and prediction."""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create DataFrame with sequences and predictions
    results_df = pd.DataFrame({
        'peptide': sequences,
        'prediction': predictions
    })
    
    # Save to CSV
    results_df.to_csv(output_path, index=False, sep='\t')

def parse_args():
    parser = argparse.ArgumentParser(description="Run PhosSight model prediction on unlabeled data.")
    parser.add_argument('-m', '--model_path', type=str, required=True,
                       help='Path to the trained model weights.')
    parser.add_argument('-p', '--path', type=str, required=True,
                       help='Path to the dataset CSV file (single column with peptide sequences).')
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='Path to save the prediction results (CSV file).')
    parser.add_argument('-b', '--batch_size', type=int, default=64,
                       help='Batch size for prediction (default: 64).')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility (default: 42).')
    parser.add_argument('--device', type=str, default='cuda',
                       help='Device to use for prediction (default: cuda).')
    parser.add_argument('--max_len', type=int, default=53, help='Maximum sequence length for Transformer')
    return parser.parse_args()

def main():
    args = parse_args() 

    # Set random seed
    set_seed(args.seed)
    
    # Get device
    device = torch.device(args.device)
    print("Using device:", device)
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

    # Prepare prediction data
    print("Preparing prediction data...")
    prediction_loader, sequences, sequences_padding = prepare_prediction_data(args.path, batch_size=args.batch_size, max_len=args.max_len)
    print(f"Prediction data prepared: {len(sequences_padding)} samples")

    # Run prediction
    print("Running prediction...")
    predictions = predict(model, prediction_loader, device)
    # Save predictions to CSV
    print(f"Saving predictions to: {args.output}")
    save_predictions_to_csv(args.output, sequences, predictions)
    
if __name__ == "__main__":
    main() 