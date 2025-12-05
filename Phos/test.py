import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import os
import csv
import logging
from datetime import datetime
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, precision_score, recall_score
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import argparse
from model import biGRU_Detect, biGRU_Detect_Improved, biGRU_Detect_Improved_V2, Transformer_Detect

def setup_logging(log_level=logging.INFO, log_file=None):
    """Setup logging configuration."""
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = f"testing_log_{timestamp}.log"
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def get_device():
    """Get the device to be used for testing."""
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

def load_test_data(path, batch_size=32):
    """Load and prepare test data."""
    data = pd.read_csv(path, header=None, sep=",")
    # the max len is set as the maximum length of the first column
    max_len = max(data[0].apply(len))
    data[0] = data[0].apply(lambda x: x.ljust(max_len, 'Z')) # padding with Z
    data[0] = data[0].apply(lambda x: Coding(x)) # coding the sequence
    
    sequences = data[0].tolist()
    labels = data[1].tolist()

    sequences = np.array(sequences)
    labels = np.array(labels)

    test_dataset = EncodedSequenceDataset(sequences, labels)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    print(f"Test data loaded: {len(test_dataset)} samples")
    return test_loader

def create_model(model_type, **kwargs):
    """
    Create model based on model type.
    
    Args:
        model_type: 'bigru', 'bigru_improved', 'bigru_improved_v2' or 'transformer'
        **kwargs: Model parameters
        
    Returns:
        model: PyTorch model
    """
    if model_type.lower() == 'bigru':
        model = biGRU_Detect(
            in_features=kwargs.get('in_features', 10),
            out_features=kwargs.get('out_features', 20),
            num_layers=kwargs.get('num_layers', 2),
            dropout=kwargs.get('dropout', 0.3)
        )
    elif model_type.lower() == 'bigru_improved':
        model = biGRU_Detect_Improved(
            in_features=kwargs.get('in_features', 10),
            out_features=kwargs.get('out_features', 20),
            num_layers=kwargs.get('num_layers', 2),
            dropout=kwargs.get('dropout', 0.3)
        )
    elif model_type.lower() == 'bigru_improved_v2':
        model = biGRU_Detect_Improved_V2(
            in_features=kwargs.get('in_features', 10),
            out_features=kwargs.get('out_features', 20),
            num_layers=kwargs.get('num_layers', 2),
            dropout=kwargs.get('dropout', 0.3)
        )
    elif model_type.lower() == 'transformer':
        model = Transformer_Detect(
            vocab_size=kwargs.get('vocab_size', 25),
            hidden_size=kwargs.get('hidden_size', 256),
            num_layers=kwargs.get('num_layers', 6),
            num_heads=kwargs.get('num_heads', 8),
            dropout=kwargs.get('dropout', 0.3),
            max_length=kwargs.get('max_length', 512)
        )
    else:
        raise ValueError(f"Unknown model type: {model_type}. Supported types: 'bigru', 'bigru_improved', 'bigru_improved_v2', 'transformer'")
    
    return model

def Test(model, test_loader, weight_path, logger=None, verbose=False):
    """Test the model."""
    if logger is None:
        logger = logging.getLogger(__name__)
    
    device = get_device()
    test_loss = 0 
    test_correct = 0
    test_batch_num = len(test_loader)
    
    # Load model weights
    try:
        state_dict = torch.load(weight_path, map_location=device)
        model.load_state_dict(state_dict)
        logger.info(f"Model weights loaded from {weight_path}")
    except Exception as e:
        logger.error(f"Failed to load model weights from {weight_path}: {e}")
        return None, None, None, None, None, None
    
    model.to(device)
    model.eval()
    criterion = nn.BCELoss()
    labels_list = []
    outputs_list = []
    predictions_list = []
    
    logger.info("Starting model testing...")
    
    with torch.no_grad():
        # Create progress bar for testing
        test_pbar = tqdm(test_loader, desc="Testing", disable=not verbose)
        for batch_idx, (sequences, labels) in enumerate(test_pbar):
            sequences, labels = sequences.to(device), labels.float().to(device)
            outputs = model(sequences)
            if outputs.dim() > 1 and outputs.size(1) == 1:
                outputs = outputs.squeeze(1)
            loss = criterion(outputs, labels)

            test_loss += loss.item()
            test_correct += ((outputs > 0.5).float() == labels).sum().item()
            preds = (outputs > 0.5).float()
            labels_list.extend(labels.cpu().numpy())
            outputs_list.extend(outputs.cpu().numpy())
            predictions_list.extend(preds.cpu().numpy())
            
            # Update progress bar with current loss
            current_loss = loss.item()
            current_acc = ((outputs > 0.5).float() == labels).sum().item() / labels.size(0)
            test_pbar.set_postfix({
                'Loss': f'{current_loss:.4f}',
                'Acc': f'{current_acc:.4f}'
            })
    
    test_loss /= test_batch_num
    test_acc = test_correct / len(test_loader.dataset)
    test_auc = roc_auc_score(labels_list, outputs_list)
    test_f1 = f1_score(labels_list, predictions_list)
    test_precision = precision_score(labels_list, predictions_list)
    test_recall = recall_score(labels_list, predictions_list)

    logger.info(f"Test Results - Loss: {test_loss:.4f}, Acc: {test_acc:.4f}, F1: {test_f1:.4f}, Precision: {test_precision:.4f}, Recall: {test_recall:.4f}, AUC: {test_auc:.4f}")

    return test_loss, test_acc, test_f1, test_precision, test_recall, test_auc

def parse_args():
    parser = argparse.ArgumentParser(description="Test model on phosphorylation datasets.")
    parser.add_argument('-p', '--path', type=str, required=True, help='Path to the test dataset CSV file.')
    parser.add_argument('-m', '--model_path', type=str, required=True, help='Path to the trained model weights.')
    parser.add_argument('-b', '--batch_size', type=int, default=128, help='Batch size for testing.')
    parser.add_argument('--model_type', type=str, default='bigru_improved_v2', choices=['bigru', 'bigru_improved', 'bigru_improved_v2', 'transformer'], help='Model type')
    parser.add_argument('--log_level', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help='Logging level')
    parser.add_argument('--log_file', type=str, default=None, help='Log file path (optional)')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--in_features', type=int, default=10, help='Input features for biGRU')
    parser.add_argument('--out_features', type=int, default=20, help='Output features for biGRU')
    parser.add_argument('--num_layers', type=int, default=2, help='Number of layers for biGRU')
    parser.add_argument('--dropout', type=float, default=0.3, help='Dropout rate for the model')
    parser.add_argument('--vocab_size', type=int, default=25, help='Vocabulary size for Transformer')
    parser.add_argument('--hidden_size', type=int, default=256, help='Hidden size for Transformer')
    parser.add_argument('--num_heads', type=int, default=8, help='Number of heads for Transformer')
    parser.add_argument('--max_length', type=int, default=512, help='Maximum sequence length for Transformer')
    parser.add_argument('--output_file', type=str, default=None, help='Output file to save results (optional)')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    
    # Setup logging
    if args.log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.log_file = f"test_log_{timestamp}.log"
    
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    logger = setup_logging(log_level=getattr(logging, args.log_level), log_file=args.log_file)
    
    device = get_device()
    logger.info(f"Using device: {device}")

    # Load test data
    test_loader = load_test_data(args.path, args.batch_size)

    # Create model
    if args.model_type.lower() == 'bigru':
        model = create_model(
            model_type='bigru',
            in_features=args.in_features,
            out_features=args.out_features,
            num_layers=args.num_layers,
            dropout=args.dropout
        )
    elif args.model_type.lower() == 'bigru_improved':
        model = create_model(
            model_type='bigru_improved',
            in_features=args.in_features,
            out_features=args.out_features,
            num_layers=args.num_layers,
            dropout=args.dropout
        )
    elif args.model_type.lower() == 'bigru_improved_v2':
        model = create_model(
            model_type='bigru_improved_v2',
            in_features=args.in_features,
            out_features=args.out_features,
            num_layers=args.num_layers,
            dropout=args.dropout
        )
    elif args.model_type.lower() == 'transformer':
        model = create_model(
            model_type='transformer',
            vocab_size=args.vocab_size,
            hidden_size=args.hidden_size,
            num_layers=args.num_layers,
            num_heads=args.num_heads,
            dropout=args.dropout,
            max_length=args.max_length
        )
    
    logger.info(f"Created {args.model_type.upper()} model")

    # Test the model
    test_loss, test_acc, test_f1, test_precision, test_recall, test_auc = Test(
        model, test_loader, args.model_path, logger, args.verbose
    )

    if test_loss is not None:
        logger.info(f"Final Test Results - Loss: {test_loss:.4f}, Acc: {test_acc:.4f}, F1: {test_f1:.4f}, Precision: {test_precision:.4f}, Recall: {test_recall:.4f}, AUC: {test_auc:.4f}")
        
        # Save results to file if specified
        if args.output_file:
            results = {
                'dataset': os.path.basename(args.path),
                'model_type': args.model_type,
                'model_path': args.model_path,
                'test_loss': test_loss,
                'test_acc': test_acc,
                'test_f1': test_f1,
                'test_precision': test_precision,
                'test_recall': test_recall,
                'test_auc': test_auc,
                'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
            
            # Create results directory if it doesn't exist
            os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
            
            # Check if file exists to determine if we need to write header
            file_exists = os.path.exists(args.output_file)
            
            with open(args.output_file, mode='a', newline='') as f:
                fieldnames = ['dataset', 'model_type', 'model_path', 'test_loss', 'test_acc', 'test_f1', 'test_precision', 'test_recall', 'test_auc', 'timestamp']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                
                if not file_exists:
                    writer.writeheader()
                
                writer.writerow(results)
            
            logger.info(f"Results saved to {args.output_file}")
    else:
        logger.error("Testing failed due to model loading error") 