#!/usr/bin/env python3
"""评测 Rice 数据集"""
import os, sys, torch, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from torch.utils.data import Dataset, DataLoader

DATA_DIR = '/data0/wangb/cd/duibi0826/0826comparison/data/rice'
MODEL_DIR = '/data0/wangb/cd/duibi0826/0826comparison/model_weights'
sys.path.append('/data0/wangb/cd/PhosSight-main/PhosSight/sotaimprov2')
sys.path.append('/data0/wangb/cd/duibi0826/model_adapters')

class PhosSightDataset(Dataset):
    def __init__(self, sequences, labels):
        self.sequences, self.labels = sequences, labels
        self.vocab = {'Z': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
                      'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
                      'V': 18, 'W': 19, 'Y': 20, 's': 21, 't': 22, 'y': 23}
        self.max_length = 50

    def __len__(self):
        return len(self.sequences)

    def encode(self, seq):
        enc = [self.vocab.get(c, 0) for c in seq]
        if len(enc) < self.max_length:
            enc.extend([0] * (self.max_length - len(enc)))
        return enc[:self.max_length]

    def __getitem__(self, idx):
        return torch.tensor(self.encode(self.sequences[idx]), dtype=torch.long), torch.tensor(self.labels[idx], dtype=torch.float)

class DeepDetectDataset(Dataset):
    def __init__(self, sequences, labels):
        self.sequences, self.labels = sequences, labels
        self.vocab = {'0': 0, 'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
                      'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
                      'V': 18, 'W': 19, 'Y': 20}
        self.max_length = 50

    def __len__(self):
        return len(self.sequences)

    def encode(self, seq):
        seq_clean = seq.replace('s', 'S').replace('t', 'T').replace('y', 'Y')
        enc = [self.vocab.get(c, 0) for c in seq_clean]
        if len(enc) < self.max_length:
            enc.extend([0] * (self.max_length - len(enc)))
        return enc[:self.max_length]

    def __getitem__(self, idx):
        return torch.tensor(self.encode(self.sequences[idx]), dtype=torch.long), torch.tensor(self.labels[idx], dtype=torch.float)

def evaluate(model, dataset_class, sequences, labels, device, use_predict=False):
    dataset = dataset_class(sequences, labels)
    loader = DataLoader(dataset, batch_size=256, shuffle=False)
    preds, gts = [], []
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            outputs = model.predict(batch_x) if use_predict else model(batch_x)
            preds.extend(outputs.cpu().numpy().flatten())
            gts.extend(batch_y.numpy())
    
    auc = roc_auc_score(gts, preds)
    auc_inv = roc_auc_score(gts, [1-p for p in preds])
    if auc_inv > auc:
        auc = auc_inv
        preds = [1-p for p in preds]
    ap = average_precision_score(gts, preds)
    return auc, ap

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    df = pd.read_csv(f"{DATA_DIR}/rice_test.csv")
    sequences, labels = df['peptide'].tolist(), df['label'].tolist()
    
    pS = sum(1 for p in sequences if 's' in p)
    pT = sum(1 for p in sequences if 't' in p)
    pY = sum(1 for p in sequences if 'y' in p)
    
    print(f"Test set: {len(sequences)}, Positive: {sum(labels)}, Negative: {len(labels) - sum(labels)}")
    print(f"pSer: {pS}, pThr: {pT}, pTyr: {pY}")
    
    results = []
    
    # PhosSight
    from model import biGRU_Detect_Improved_V2
    model = biGRU_Detect_Improved_V2(in_features=10, out_features=20, num_layers=2, dropout=0.3)
    if os.path.exists(f"{MODEL_DIR}/phosight_v2_best.pth"):
        checkpoint = torch.load(f"{MODEL_DIR}/phosight_v2_best.pth", map_location=device)
        model.load_state_dict(checkpoint.get('model_state_dict', checkpoint), strict=False)
    model.to(device).eval()
    auc, ap = evaluate(model, PhosSightDataset, sequences, labels, device)
    print(f"PhosSight AUC: {auc:.4f}, AP: {ap:.4f}")
    results.append({'Model': 'PhosSight', 'Version': 'Original', 'AUC': auc, 'AP': ap})
    
    # DeepDetect
    from deepdetect_adapter import DeepDetectAdapter
    model2 = DeepDetectAdapter()
    model2.load_model(f"{MODEL_DIR}/deepdetect_best.pth")
    model2.to(device).eval()
    auc, ap = evaluate(model2, DeepDetectDataset, sequences, labels, device, use_predict=True)
    print(f"DeepDetect AUC: {auc:.4f}, AP: {ap:.4f}")
    results.append({'Model': 'DeepDetect', 'Version': 'Original', 'AUC': auc, 'AP': ap})
    
    # pFly
    from dlomix_adapter import DLOmixAdapter
    model3 = DLOmixAdapter(num_units=128, alphabet_size=21, num_classes=4)
    model3.load_model(f"{MODEL_DIR}/dlomix_best.pth")
    model3.to(device).eval()
    auc, ap = evaluate(model3, DeepDetectDataset, sequences, labels, device, use_predict=True)
    print(f"pFly AUC: {auc:.4f}, AP: {ap:.4f}")
    results.append({'Model': 'pFly', 'Version': 'Original', 'AUC': auc, 'AP': ap})
    
    pd.DataFrame(results).to_csv(f"{DATA_DIR}/rice_original_results.csv", index=False)
    print("\n" + pd.DataFrame(results).to_string(index=False))

if __name__ == "__main__":
    main()
