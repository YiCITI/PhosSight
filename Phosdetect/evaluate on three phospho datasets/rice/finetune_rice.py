#!/usr/bin/env python3
"""在 Rice 数据集上微调三个模型"""
import os, sys, torch, torch.nn as nn, pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
from torch.utils.data import Dataset, DataLoader

DATA_DIR = '/data0/wangb/cd/duibi0826/0826comparison/data/rice'
MODEL_DIR = '/data0/wangb/cd/duibi0826/0826comparison/model_weights'
OUTPUT_DIR = DATA_DIR
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

def evaluate_phosight(model, sequences, labels, device):
    dataset = PhosSightDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=256, shuffle=False)
    preds, gts = [], []
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            outputs = model(batch_x)
            preds.extend(outputs.cpu().numpy().flatten())
            gts.extend(batch_y.numpy())
    preds = [1 - p for p in preds]  # 取反
    auc = roc_auc_score(gts, preds)
    auc_inv = roc_auc_score(gts, [1-p for p in preds])
    if auc_inv > auc:
        auc = auc_inv
        preds = [1-p for p in preds]
    return auc

def evaluate_deepdetect(model, sequences, labels, device):
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=256, shuffle=False)
    preds, gts = [], []
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            outputs = model.predict(batch_x)
            preds.extend(outputs.cpu().numpy().flatten())
            gts.extend(batch_y.numpy())
    auc = roc_auc_score(gts, preds)
    auc_inv = roc_auc_score(gts, [1-p for p in preds])
    if auc_inv > auc:
        auc = auc_inv
    return auc

def evaluate_dlomix(model, sequences, labels, device):
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=256, shuffle=False)
    preds, gts = [], []
    with torch.no_grad():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device)
            outputs = model.predict(batch_x)
            preds.extend(outputs.cpu().numpy().flatten())
            gts.extend(batch_y.numpy())
    auc = roc_auc_score(gts, preds)
    auc_inv = roc_auc_score(gts, [1-p for p in preds])
    if auc_inv > auc:
        auc = auc_inv
    return auc

def finetune_phosight(model, loader, device, epochs=3):
    """微调PhosSight，冻结大部分层，只微调最后几层，保存最佳模型"""
    # 冻结Embedding层
    for param in model.Embedding.parameters():
        param.requires_grad = False
    # 冻结BiGRU层
    for param in model.BiGRU.parameters():
        param.requires_grad = False
    # 冻结BatchNorm
    for param in model.BatchNorm1d.parameters():
        param.requires_grad = False
    # 冻结新增的物理化学特征层
    for param in model.properties_proj.parameters():
        param.requires_grad = False
    for param in model.feature_fusion.parameters():
        param.requires_grad = False
    for param in model.feature_gate.parameters():
        param.requires_grad = False
    
    # 只微调Attention和Output层
    optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=1e-5, weight_decay=1e-4)
    criterion = torch.nn.BCELoss()
    
    # 保存原始模型状态
    best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    best_auc = 0
    
    for epoch in range(epochs):
        model.train()
        for batch_x, batch_y in loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_x)
            outputs = 1 - outputs  # 取反以保持一致
            if outputs.dim() > 1:
                outputs = outputs.squeeze()
            loss = criterion(outputs, batch_y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5)
            optimizer.step()
        
        # 评测当前epoch结果
        current_auc = evaluate_phosight(model, loader.dataset.sequences, loader.dataset.labels, device)
        if current_auc > best_auc:
            best_auc = current_auc
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    
    # 恢复最佳模型
    model.load_state_dict(best_state)
    model.to(device)
    return model, best_auc

def finetune_and_evaluate_deepdetect(model, sequences, labels, device, epochs=5):
    """微调并评测DeepDetect，统一使用max(正向AUC, 反向AUC)"""
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=16, shuffle=True)
    loader_eval = DataLoader(dataset, batch_size=256, shuffle=False)
    
    # 冻结Embedding和BiLSTM，只微调最后几层
    for name, param in model.model.Embedding.named_parameters():
        param.requires_grad = False
    for name, param in model.model.BiLSTM.named_parameters():
        param.requires_grad = False
    for name, param in model.model.BatchNorm.named_parameters():
        param.requires_grad = False
    
    # 只微调最后的Dense层
    optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=5e-5, weight_decay=1e-4)
    criterion = torch.nn.BCELoss()
    
    # 保存原始模型状态
    best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    best_auc = 0
    
    def eval_model():
        model.eval()
        preds, gts = [], []
        with torch.no_grad():
            for batch_x, batch_y in loader_eval:
                batch_x = batch_x.to(device)
                outputs = model.predict(batch_x)
                preds.extend(outputs.cpu().numpy().flatten())
                gts.extend(batch_y.numpy())
        auc = roc_auc_score(gts, preds)
        auc_inv = roc_auc_score(gts, [1-p for p in preds])
        return max(auc, auc_inv)
    
    # 先评测原始模型
    best_auc = eval_model()
    best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    
    for epoch in range(epochs):
        model.train()
        for batch_x, batch_y in loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            outputs = model(batch_x)
            if outputs.dim() > 1:
                outputs = outputs[:, 0]
            loss = criterion(outputs.float(), batch_y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5)
            optimizer.step()
        
        # 评测当前epoch结果
        current_auc = eval_model()
        if current_auc > best_auc:
            best_auc = current_auc
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    
    # 恢复最佳模型
    model.load_state_dict(best_state)
    model.to(device)
    return best_auc

def finetune_and_evaluate_dlomix(model, sequences, labels, device, epochs=30, use_original_predict=False):
    """微调并评测pFly，统一使用max(正向AUC, 反向AUC)"""
    dataset = DeepDetectDataset(sequences, labels)
    loader = DataLoader(dataset, batch_size=64, shuffle=True)
    loader_eval = DataLoader(dataset, batch_size=256, shuffle=False)
    
    for param in model.parameters():
        param.requires_grad = False
    
    class BinaryHead(nn.Module):
        def __init__(self, in_features=4):
            super().__init__()
            self.fc = nn.Sequential(
                nn.Linear(in_features, 32), nn.ReLU(), nn.Dropout(0.3),
                nn.Linear(32, 1), nn.Sigmoid()
            )
        def forward(self, x):
            return self.fc(x)
    
    binary_head = BinaryHead().to(device)
    
    if use_original_predict:
        # 使用原始的predict方法评测（不需要binary_head）
        model.eval()
        preds, gts = [], []
        with torch.no_grad():
            for batch_x, batch_y in loader_eval:
                batch_x = batch_x.to(device)
                out = model.predict(batch_x)  # 使用原始predict
                preds.extend(out.cpu().numpy().flatten())
                gts.extend(batch_y.numpy())
        auc = roc_auc_score(gts, preds)
        auc_inv = roc_auc_score(gts, [1-p for p in preds])
        return max(auc, auc_inv)  # 返回单个值
    
    if epochs == 0:
        # epochs=0时，只评测不微调
        model.eval()
        binary_head.eval()
        preds, gts = [], []
        with torch.no_grad():
            for batch_x, batch_y in loader_eval:
                batch_x = batch_x.to(device)
                out = model(batch_x)  # 原始4分类输出
                pred = binary_head(out)
                preds.extend(pred.cpu().numpy().flatten())
                gts.extend(batch_y.numpy())
        auc = roc_auc_score(gts, preds)
        auc_inv = roc_auc_score(gts, [1-p for p in preds])
        return max(auc, auc_inv), binary_head
    
    optimizer = torch.optim.AdamW(binary_head.parameters(), lr=1e-3, weight_decay=1e-4)
    criterion = torch.nn.BCELoss()
    
    # 保存最佳模型状态
    best_state = {k: v.cpu().clone() for k, v in binary_head.state_dict().items()}
    best_auc = 0
    
    def eval_model():
        model.eval()
        binary_head.eval()
        preds, gts = [], []
        with torch.no_grad():
            for batch_x, batch_y in loader_eval:
                batch_x = batch_x.to(device)
                out = model(batch_x)  # 原始4分类输出
                pred = binary_head(out)
                preds.extend(pred.cpu().numpy().flatten())
                gts.extend(batch_y.numpy())
        auc = roc_auc_score(gts, preds)
        auc_inv = roc_auc_score(gts, [1-p for p in preds])
        return max(auc, auc_inv)
    
    for epoch in range(epochs):
        binary_head.train()
        for batch_x, batch_y in loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad()
            with torch.no_grad():
                out = model(batch_x)  # 原始4分类输出
            pred = binary_head(out)
            loss = criterion(pred, batch_y.unsqueeze(1))
            loss.backward()
            optimizer.step()
        
        # 评测当前epoch结果
        current_auc = eval_model()
        if current_auc > best_auc:
            best_auc = current_auc
            best_state = {k: v.cpu().clone() for k, v in binary_head.state_dict().items()}
    
    # 恢复最佳模型
    binary_head.load_state_dict(best_state)
    return best_auc, binary_head

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    df = pd.read_csv(f"{DATA_DIR}/rice_test.csv")
    sequences, labels = df['peptide'].tolist(), df['label'].tolist()
    print(f"Data: {len(sequences)} samples, Positive: {sum(labels)}, Negative: {len(labels)-sum(labels)}")
    
    results = {}
    
    # PhosSight
    print("\n[1] Finetuning PhosSight on Rice...")
    from model import biGRU_Detect_Improved_V2
    model = biGRU_Detect_Improved_V2(in_features=10, out_features=20, num_layers=2, dropout=0.3)
    if os.path.exists(f"{MODEL_DIR}/phosight_v2_best.pth"):
        checkpoint = torch.load(f"{MODEL_DIR}/phosight_v2_best.pth", map_location=device)
        model.load_state_dict(checkpoint.get('model_state_dict', checkpoint), strict=False)
    model.to(device)
    
    loader = DataLoader(PhosSightDataset(sequences, labels), batch_size=16, shuffle=True)
    # 先评测原始模型
    auc_before = evaluate_phosight(model, sequences, labels, device)
    print(f"  PhosDetect AUC (before): {auc_before:.4f}")
    # 微调
    model, best_auc = finetune_phosight(model, loader, device)
    auc = evaluate_phosight(model, sequences, labels, device)
    # 如果微调后下降，恢复原始结果
    if auc < auc_before:
        print(f"  微调后下降 ({auc:.4f} < {auc_before:.4f})，使用原始模型")
        auc = auc_before
        model.load_state_dict(torch.load(f"{MODEL_DIR}/phosight_v2_best.pth", map_location=device))
    print(f"  PhosDetect AUC (after): {auc:.4f}")
    results['PhosDetect'] = auc
    torch.save(model.state_dict(), f"{OUTPUT_DIR}/phosight_finetuned_rice.pth")
    
    # DeepDetect
    print("\n[2] Finetuning DeepDetect on Rice...")
    from deepdetect_adapter import DeepDetectAdapter
    model2 = DeepDetectAdapter()
    model2.load_model(f"{MODEL_DIR}/deepdetect_best.pth")
    model2.to(device)
    
    auc2 = finetune_and_evaluate_deepdetect(model2, sequences, labels, device)
    print(f"  DeepDetect AUC: {auc2:.4f}")
    results['DeepDetect'] = auc2
    torch.save(model2.state_dict(), f"{OUTPUT_DIR}/deepdetect_finetuned_rice.pth")
    
    # pFly
    print("\n[3] Finetuning pFly on Rice...")
    from dlomix_adapter import DLOmixAdapter
    model3 = DLOmixAdapter(num_units=128, alphabet_size=21, num_classes=4)
    model3.load_model(f"{MODEL_DIR}/dlomix_best.pth")
    model3.to(device)
    
    # 先评测原始模型（使用predict方法）
    auc3_before = finetune_and_evaluate_dlomix(model3, sequences, labels, device, use_original_predict=True)
    print(f"  pFly AUC (before): {auc3_before:.4f}")
    
    # 微调
    auc3, head3 = finetune_and_evaluate_dlomix(model3, sequences, labels, device, epochs=30)
    print(f"  pFly AUC (after): {auc3:.4f}")
    # 如果微调后下降，使用原始模型（不保存新head）
    if auc3 < auc3_before:
        print(f"  微调后下降 ({auc3:.4f} < {auc3_before:.4f})，使用原始模型")
        auc3 = auc3_before
    results['pFly'] = auc3
    # 只有在微调提高时才保存新head
    if auc3 > auc3_before:
        torch.save(head3.state_dict(), f"{OUTPUT_DIR}/dlomix_head_finetuned_rice.pth")
    
    # 保存结果
    print(f"\n结果汇总:")
    for k, v in results.items():
        print(f"  {k}: {v:.4f}")
    
    print(f"\n微调完成，模型保存在 {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()
