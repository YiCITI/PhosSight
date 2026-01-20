import time
import os
import torch
from torch import nn
import numpy as np
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, precision_score, recall_score
import random
from model import biGRU_Detect_Improved_V2 as Detect
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
import pandas as pd
import csv
import logging
from datetime import datetime
from tqdm import tqdm
import math

def augment_batch_sequences(sequences: torch.Tensor, probability: float, augmentation_type: str = 'mask'):
    """Apply lightweight token-level augmentation to a batch of integer-encoded sequences.

    - mask: randomly replace tokens (excluding padding 0) with padding token 0 with given probability
    - replace: randomly replace tokens (excluding padding 0) with a random amino-acid id in [1, 23]

    Args:
        sequences: LongTensor of shape (batch_size, seq_len)
        probability: probability to alter each token independently
        augmentation_type: 'mask' or 'replace'

    Returns:
        Augmented LongTensor of the same shape as input
    """
    if probability <= 0.0:
        return sequences
    if probability >= 1.0:
        probability = 1.0

    augmented = sequences.clone()
    non_pad_mask = augmented != 0
    rand_mask = (torch.rand_like(augmented, dtype=torch.float) < probability) & non_pad_mask

    if augmentation_type == 'mask':
        augmented[rand_mask] = 0
    elif augmentation_type == 'replace':
        random_tokens = torch.randint(low=1, high=24, size=augmented.shape, device=augmented.device)
        augmented[rand_mask] = random_tokens[rand_mask]
    else:
        return sequences

    return augmented

def setup_logging(log_level=logging.INFO):
    """Setup logging configuration."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

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
    def __init__(self, sequences, labels):
        self.sequences = torch.tensor(sequences, dtype=torch.long)
        self.labels = torch.tensor(labels, dtype=torch.float)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.sequences[idx], self.labels[idx]

class EarlyStopping:
    """Early stopping utility to prevent overfitting."""
    
    def __init__(self, patience=10, min_delta=0.0, restore_best_weights=True):
        self.patience = patience
        self.min_delta = min_delta
        self.restore_best_weights = restore_best_weights
        self.best_score = None
        self.counter = 0
        self.best_weights = None
        
    def __call__(self, val_score, model):
        if self.best_score is None:
            self.best_score = val_score
            self.best_weights = model.state_dict().copy()
        elif val_score > self.best_score + self.min_delta:
            self.best_score = val_score
            self.counter = 0
            self.best_weights = model.state_dict().copy()
        else:
            self.counter += 1
            
        if self.counter >= self.patience:
            if self.restore_best_weights:
                model.load_state_dict(self.best_weights)
            return True
        return False

# split the data into train, validation and test sets. Returns the data loaders for each set.
def split_data(path, batch_size=32, val_size=0.05, test_size=0.05, max_len=53):
    data = pd.read_csv(path, sep='\t')
    # the max len is set as the maximum length of the first column
    # max_len = max(data.iloc[:, 0].apply(len))
    data.iloc[:,0] = data.iloc[:,0].apply(lambda x: x[:max_len] if len(x) > max_len else x)
    data.iloc[:,0] = data.iloc[:,0].apply(lambda x: x.ljust(max_len, 'Z')) # padding with Z
    data.iloc[:,0] = data.iloc[:,0].apply(lambda x: Coding(x)) # coding the sequence
    
    # split the data into train, validation and test sets
    sequences = data.iloc[:,0].tolist()
    labels = data.iloc[:,1].tolist()

    sequences = np.array(sequences)
    labels = np.array(labels)

    X_train, X_test, y_train, y_test = train_test_split(
        sequences, labels, test_size=test_size, random_state=42, stratify=labels
    )
    X_train, X_val, y_train, y_val = train_test_split(
        X_train, y_train, test_size=val_size / (1 - test_size), random_state=42, stratify=y_train
    )
    print("train_data size:",len(X_train))
    print("val_data size:", len(X_val))
    print("test_data size:",len(X_test))
    for i in range(len(y_train)):
        if y_train[i]== None:
            print(i)
    train_dataset = EncodedSequenceDataset(X_train, y_train)
    val_dataset = EncodedSequenceDataset(X_val, y_val)
    test_dataset = EncodedSequenceDataset(X_test, y_test)
    # create the data loaders
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    return train_loader, val_loader, test_loader


def get_lr_scheduler(optimizer, scheduler_type, learning_rate, num_epochs, num_batches_per_epoch):
    """Get learning rate scheduler based on scheduler type."""
    if scheduler_type == 'cosine':
        # Cosine annealing with warmup
        warmup_steps = min(2000, num_batches_per_epoch * 3)
        total_steps = num_batches_per_epoch * num_epochs
        
        def lr_lambda(step):
            if step < warmup_steps:
                return float(step) / float(max(1, warmup_steps))
            else:
                progress = float(step - warmup_steps) / float(max(1, total_steps - warmup_steps))
                return max(0.0, 0.5 * (1.0 + math.cos(math.pi * progress)))
        
        scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)
        return scheduler, warmup_steps
    elif scheduler_type == 'plateau':
        # ReduceLROnPlateau
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=0.7, patience=8, min_lr=1e-6
        )
        return scheduler, 0
    elif scheduler_type == 'step':
        # StepLR
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.7)
        return scheduler, 0
    else:
        raise ValueError(f"Unknown scheduler type: {scheduler_type}")

def Train(args, model, train_loader, val_loader, learning_rate, num_epochs, model_path, device, patience=10, logger=None):
    """Train the model with unified strategy."""
    if logger is None:
        logger = logging.getLogger(__name__)
    
    model.to(device)
    
    logger.info(f"Training on device: {device}")
    logger.info(f"Model type: {type(model).__name__}")
    logger.info(f"Total parameters: {sum(p.numel() for p in model.parameters()):,}")
    logger.info(f"Trainable parameters: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")

    criterion = nn.BCELoss()
    
    # Unified optimizer settings
    if args.optimizer == 'adam':
        optimizer = torch.optim.Adam(
            params=model.parameters(), 
            lr=learning_rate, 
            weight_decay=args.weight_decay,
            betas=(args.beta1, args.beta2)
        )
    elif args.optimizer == 'adamw':
        optimizer = torch.optim.AdamW(
            params=model.parameters(), 
            lr=learning_rate, 
            weight_decay=args.weight_decay,
            betas=(args.beta1, args.beta2),
            eps=args.eps
        )
    elif args.optimizer == 'sgd':
        optimizer = torch.optim.SGD(
            params=model.parameters(),
            lr=learning_rate,
            momentum=args.momentum,
            weight_decay=args.weight_decay
        )
    else:
        raise ValueError(f"Unknown optimizer: {args.optimizer}")
    
    # Unified scheduler
    scheduler, warmup_steps = get_lr_scheduler(
        optimizer, args.scheduler, learning_rate, num_epochs, len(train_loader)
    )
    
    best_val_acc = 0.0
    best_model_path = model_path
    
    # Unified early stopping
    early_stopping = EarlyStopping(patience=patience, min_delta=args.min_delta)
    
    logger.info(f"Training started with learning_rate={learning_rate}, epochs={num_epochs}, patience={patience}")
    logger.info(f"Optimizer: {args.optimizer}, Scheduler: {args.scheduler}, Weight decay: {args.weight_decay}")

    log_data = []  
    log_header = ['epoch', 'train_loss', 'train_acc', 'val_loss', 'val_acc', 'lr']

    # Select best-score direction according to monitored metric
    monitor_metric = getattr(args, 'early_stopping_metric', 'acc')
    best_score = float('-inf') if monitor_metric != 'loss' else float('-inf')

    for epoch in range(num_epochs):
        model.train()
        start = time.perf_counter()
        train_loss = 0
        train_correct = 0
        train_batch_num = len(train_loader)
        
        logger.info('-'*30 + f' Epoch {epoch+1}/{num_epochs}, LR={optimizer.param_groups[0]["lr"]:.6f} ' + '-'*30)
        
        # Create progress bar for training
        train_pbar = tqdm(train_loader, desc=f"Training Step", disable=not args.verbose)
        for batch_idx, (sequences, labels) in enumerate(train_pbar):
            sequences, labels = sequences.to(device), labels.float().to(device)
            # On-the-fly sequence augmentation
            if getattr(args, 'seq_aug_prob', 0.0) > 0.0:
                sequences = augment_batch_sequences(sequences, args.seq_aug_prob, getattr(args, 'seq_aug_type', 'mask'))
            optimizer.zero_grad()
            outputs = model(sequences)
            if outputs.dim() > 1 and outputs.size(1) == 1:
                outputs = outputs.squeeze(1)
            # Apply label smoothing for BCE
            if getattr(args, 'label_smoothing', 0.0) > 0.0:
                smoothed_labels = labels * (1.0 - args.label_smoothing) + 0.5 * args.label_smoothing
            else:
                smoothed_labels = labels
            loss = criterion(outputs, smoothed_labels)
            loss.backward()
            
            # Gradient clipping
            if args.gradient_clip > 0:
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=args.gradient_clip)
            
            optimizer.step()

            train_loss += loss.item()
            train_correct += ((outputs > 0.5).float() == labels).sum().item()
            
            # Update progress bar with current loss
            current_loss = loss.item()
            current_acc = ((outputs > 0.5).float() == labels).sum().item() / labels.size(0)
            train_pbar.set_postfix({
                'Loss': f'{current_loss:.4f}',
                'Acc': f'{current_acc:.4f}'
            })

        train_loss /= train_batch_num
        train_acc = train_correct / len(train_loader.dataset)

        model.eval()
        val_loss = 0
        val_correct = 0
        val_batch_num = len(val_loader)
        with torch.no_grad():
            # Create progress bar for validation
            val_pbar = tqdm(val_loader, desc=f"Validation Step", disable=not args.verbose)
            for batch_idx, (sequences, labels) in enumerate(val_pbar):
                sequences, labels = sequences.to(device), labels.float().to(device)
                outputs = model(sequences)
                if outputs.dim() > 1 and outputs.size(1) == 1:
                    outputs = outputs.squeeze(1)
                loss = criterion(outputs, labels)

                val_loss += loss.item()
                val_correct += ((outputs > 0.5).float() == labels).sum().item()
                
                # Update progress bar with current loss
                current_loss = loss.item()
                current_acc = ((outputs > 0.5).float() == labels).sum().item() / labels.size(0)
                val_pbar.set_postfix({
                    'Loss': f'{current_loss:.4f}',
                    'Acc': f'{current_acc:.4f}'
                })

        val_loss /= val_batch_num
        val_acc = val_correct / len(val_loader.dataset)
        end_time = time.perf_counter()
        epoch_time = end_time - start
        
        # Update learning rate scheduler
        if args.scheduler == 'plateau':
            scheduler.step(val_acc)
        else:
            scheduler.step()
        
        current_lr = optimizer.param_groups[0]['lr']
        
        logger.info(f"Epoch [{epoch+1}/{num_epochs}], Train Loss: {train_loss:.4f}, Train Acc: {train_acc:.4f}, Val Loss: {val_loss:.4f}, Val Acc: {val_acc:.4f}, Time: {epoch_time:.2f}s, LR: {current_lr:.6f}")
        
        # save the log data for this epoch
        log_data.append([epoch + 1, train_loss, train_acc, val_loss, val_acc, current_lr])

        # Monitor metric: maximize val_acc or minimize val_loss via maximizing negative loss
        if monitor_metric == 'loss':
            current_score = -val_loss
        else:
            current_score = val_acc

        if current_score > best_score:
            best_score = current_score
            best_val_acc = val_acc
            os.makedirs(os.path.dirname(best_model_path), exist_ok=True)
            torch.save(model.state_dict(), best_model_path)
            logger.info(f"New best model saved (monitor={monitor_metric}) with Val Acc: {val_acc:.4f}, Val Loss: {val_loss:.4f}")

        # Check for early stopping
        if early_stopping(current_score, model):
            logger.info(f"Early stopping triggered at epoch {epoch + 1}")
            break

    logger.info(f"Training completed. Best model saved with val acc: {best_val_acc:.4f}")

    # write the log data to a CSV file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = os.path.join(os.path.dirname(args.model_path), f'training_log_{timestamp}.csv')
    with open(csv_path, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(log_header)
        writer.writerows(log_data)
    logger.info(f"Training log saved to {csv_path}")

    return None

def Test(model, test_loader, weight_path, device, logger=None):
    """Test the model."""
    if logger is None:
        logger = logging.getLogger(__name__)
    
    test_loss = 0 
    test_correct = 0
    test_batch_num = len(test_loader)
    state_dict = torch.load(weight_path, map_location=device)
    model.load_state_dict(state_dict)
    model.to(device)
    model.eval()
    criterion = nn.BCELoss()
    labels_list = []
    outputs_list = []
    predictions_list = []
    
    logger.info("Starting model testing...")
    
    with torch.no_grad():
        # Create progress bar for testing
        test_pbar = tqdm(test_loader, desc="Testing", disable=not args.verbose)
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

def create_model(**kwargs):
    """
    Create model based on model type.
    
    Args:
        **kwargs: Model parameters
        
    Returns:
        model: PyTorch model
    """
    model = Detect(
        in_features=kwargs.get('in_features', 10),
        out_features=kwargs.get('out_features', 20),
        num_layers=kwargs.get('num_layers', 2),
        dropout=kwargs.get('dropout', 0.3)
    )
    
    if 'pretrained_model_path' in kwargs and kwargs.get('pretrained_model_path') is not None:
        pretrained_model_path = kwargs.get('pretrained_model_path')
        state_dict = torch.load(pretrained_model_path, map_location=device)
        model.load_state_dict(state_dict)
    
    return model

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Input parameters.")
    parser.add_argument('-b', '--batch_size', type=int, default=256, help='Batch size for training and testing.')  # Reduced from 256
    parser.add_argument('-e', '--epochs', type=int, default=20, help='Number of epochs for training.')  # Reduced from 100
    parser.add_argument('-lr', '--learning_rate', type=float, default=1e-5, help='Learning rate for the optimizer.')  # Reduced from 0.001
    parser.add_argument('-d', '--dropout', type=float, default=0.3, help='Dropout rate for the model.')  # Increased from 0.3
    parser.add_argument('-p','--path', type=str, default='../data/train/balanced_dataset_1.csv', help='Path to the dataset CSV file.')
    parser.add_argument('-m', '--model_path', type=str, default='./best_model.pth', help='Path to save the best model.')
    parser.add_argument('-w', '--weight_path', type=str, default=None, help='Pretain weight path.')
    parser.add_argument('--patience', type=int, default=5, help='Patience for early stopping.')  # Reduced from 15
    parser.add_argument('--log_level', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help='Logging level')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--optimizer', type=str, default='adam', choices=['adam', 'adamw', 'sgd'], help='Optimizer type')
    parser.add_argument('--weight_decay', type=float, default=1e-5, help='Weight decay for optimizer')
    parser.add_argument('--beta1', type=float, default=0.9, help='Beta1 for Adam/AdamW')
    parser.add_argument('--beta2', type=float, default=0.999, help='Beta2 for Adam/AdamW')
    parser.add_argument('--eps', type=float, default=1e-8, help='Epsilon for Adam/AdamW')
    parser.add_argument('--momentum', type=float, default=0.9, help='Momentum for SGD')
    parser.add_argument('--scheduler', type=str, default='step', choices=['plateau', 'cosine', 'step'], help='Learning rate scheduler type')
    parser.add_argument('--gradient_clip', type=float, default=0.0, help='Gradient clipping norm (0 for no clipping)')
    parser.add_argument('--min_delta', type=float, default=0.001, help='Minimum delta for early stopping')
    parser.add_argument('--early_stopping_metric', type=str, default='acc', choices=['acc', 'loss'], help='Metric to monitor for early stopping and checkpointing')
    parser.add_argument('--label_smoothing', type=float, default=0.05, help='Label smoothing factor for BCE targets (0 disables)')
    parser.add_argument('--seq_aug_prob', type=float, default=0.05, help='Per-token probability for sequence augmentation during training (0 disables)')
    parser.add_argument('--seq_aug_type', type=str, default='mask', choices=['mask', 'replace'], help='Type of sequence augmentation applied during training')
    parser.add_argument('--in_features', type=int, default=10, help='Input features for biGRU')
    parser.add_argument('--out_features', type=int, default=20, help='Output features for biGRU')
    parser.add_argument('--num_layers', type=int, default=2, help='Number of layers for biGRU')
    parser.add_argument('--vocab_size', type=int, default=24, help='Vocabulary size for Transformer')
    parser.add_argument('--hidden_size', type=int, default=256, help='Hidden size for Transformer')
    parser.add_argument('--num_heads', type=int, default=8, help='Number of heads for Transformer')
    parser.add_argument('--max_length', type=int, default=512, help='Maximum sequence length for Transformer')
    parser.add_argument('--val_size', type=float, default=0.05, help='Validation set size')
    parser.add_argument('--test_size', type=float, default=0.05, help='Test set size')
    parser.add_argument('--sep', type=str, default='\t', help='Separator for the dataset CSV file')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    parser.add_argument('--max_len', type=int, default=53, help='Maximum sequence length for Transformer')
    parser.add_argument('--device', type=str, default='cuda', help='Device to use for training and testing')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    # Setup logging
    logger = setup_logging(log_level=getattr(logging, args.log_level))
    
    set_seed(args.seed)
    device = args.device
    logger.info(f"Using device: {device}")

    path = args.path
    batch_size = args.batch_size
    train_loader, val_loader, test_loader = split_data(path, batch_size=batch_size, val_size=args.val_size, test_size=args.test_size, max_len=args.max_len)

    dropout = args.dropout
    
    # Create model based on model type
    model = create_model(
        in_features=args.in_features,
        out_features=args.out_features,
        num_layers=args.num_layers,
        dropout=args.dropout,
        pretrained_model_path=args.weight_path
    )
    
    logger.info(f"Initialized model")

    learning_rate = args.learning_rate
    num_epochs = args.epochs
    patience = args.patience
    pretrained_model_path = args.weight_path
    model_path = args.model_path
    Train(args, model, train_loader, val_loader, learning_rate, num_epochs, model_path, device, patience, logger)
    model_path = args.model_path

    test_loss, test_acc, test_f1, test_precision, test_recall, test_auc = Test(model, test_loader, model_path, device, logger)

    logger.info(f"Final Test Results - Loss: {test_loss:.4f}, Acc: {test_acc:.4f}, F1: {test_f1:.4f}, Precision: {test_precision:.4f}, Recall: {test_recall:.4f}, AUC: {test_auc:.4f}")

