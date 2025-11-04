import time
import torch
import numpy as np
import random
import pandas as pd
import csv
from torch import nn
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from model import biGRU_Detect as Detect  

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def get_device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")

def Coding(mers):
    dic = {'Z': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
           'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
           'V': 18, 'W': 19, 'Y': 20, 's': 21, 't': 22, 'y': 23}
    return [dic.get(aa, 0) for aa in mers]

class EncodedSequenceDataset(Dataset):
    def __init__(self, sequences, labels):
        self.sequences = torch.tensor(sequences, dtype=torch.long)
        self.labels = torch.tensor(labels, dtype=torch.float)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.sequences[idx], self.labels[idx]

def split_data(path, batch_size=64, val_size=0.05, test_size=0.0):
    data = pd.read_csv(path, header=None, sep=",")
    max_len = max(data[0].apply(len))
    data[0] = data[0].apply(lambda x: x.ljust(max_len, 'Z'))
    data[0] = data[0].apply(Coding)

    sequences = np.array(data[0].tolist())
    labels = np.array(data[1].tolist())

    if test_size > 0:
        X_train, X_test, y_train, y_test = train_test_split(
            sequences, labels, test_size=test_size, random_state=42, stratify=labels
        )
    else:
        X_train, y_train = sequences, labels
        X_test, y_test = [], []

    X_train, X_val, y_train, y_val = train_test_split(
        X_train, y_train, test_size=val_size / (1 - test_size) if test_size < 1 else val_size,
        random_state=42, stratify=y_train
    )

    print("train_data:", len(X_train))
    print("val_data:", len(X_val))

    train_loader = DataLoader(EncodedSequenceDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(EncodedSequenceDataset(X_val, y_val), batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(EncodedSequenceDataset(X_test, y_test), batch_size=batch_size, shuffle=False) if test_size > 0 else None

    return train_loader, val_loader, test_loader

def Train_Finetune(model, train_loader, val_loader, learning_rate, num_epochs, finetune_path):
    device = get_device()
    model.to(device)

    criterion = nn.BCELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    best_val_acc = 0.0
    best_model_path = finetune_path

    for epoch in range(num_epochs):
        model.train()
        train_loss = 0
        train_correct = 0
        for sequences, labels in train_loader:
            sequences, labels = sequences.to(device), labels.float().to(device)
            outputs = model(sequences).squeeze(1)
            loss = criterion(outputs, labels)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
            train_correct += ((outputs > 0.5).float() == labels).sum().item()

        train_acc = train_correct / len(train_loader.dataset)

        model.eval()
        val_loss = 0
        val_correct = 0
        with torch.no_grad():
            for sequences, labels in val_loader:
                sequences, labels = sequences.to(device), labels.float().to(device)
                outputs = model(sequences).squeeze(1)
                loss = criterion(outputs, labels)
                val_loss += loss.item()
                val_correct += ((outputs > 0.5).float() == labels).sum().item()

        val_acc = val_correct / len(val_loader.dataset)
        print(f"Epoch [{epoch+1}/{num_epochs}], Train Acc: {train_acc:.4f}, Val Acc: {val_acc:.4f}")

        if val_acc > best_val_acc:
            best_val_acc = val_acc
            torch.save(model.state_dict(), best_model_path)

    print(f"Finetuned best model saved as {best_model_path}, Val Acc: {best_val_acc:.4f}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Input parameters.")
    parser.add_argument('-b', '--batch_size', type=int, default=512, help='Batch size for training and testing.')
    parser.add_argument('-e', '--epochs', type=int, default=5, help='Number of epochs for training.')
    parser.add_argument('-lr', '--learning_rate', type=float, default=0.0005, help='Learning rate for the optimizer.')
    parser.add_argument('-d', '--dropout', type=float, default=0.3, help='Dropout rate for the model.')
    parser.add_argument('-p','--path', type=str, default='../data/train/balanced_dataset_1.csv', help='Path to the dataset CSV file.')
    parser.add_argument('-m', '--model_path', type=str, default='../model/best_model.pth', help='Path to the model to be fine-tuned.')
    parser.add_argument('-t', '--finetune', type=str, default='./finetune.pth', help='Path to the finetuned model.')


    set_seed(42)
    print("Using device:", get_device())
    args = parser.parse_args()
    path = args.path
    batch_size = args.batch_size
    learning_rate = args.learning_rate
    num_epochs = args.epochs
    dropout = args.dropout

    model_path = args.model_path
    finetune_path = args.finetune


    train_loader, val_loader, _ = split_data(path, batch_size=batch_size, val_size=0.05, test_size=0.0)

    model = Detect(in_features=10, out_features=20, num_layers=2, dropout=dropout)
    model.load_state_dict(torch.load(model_path, map_location=get_device()))


    Train_Finetune(model, train_loader, val_loader, learning_rate, num_epochs, finetune_path)

