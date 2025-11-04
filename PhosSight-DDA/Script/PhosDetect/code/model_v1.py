import torch.nn as nn
import torch
import re
from transformers import BertModel
import torch.nn.functional as F

def Embedding(in_features):
    return nn.Embedding(24, in_features)  # 20 aa + 3 phospho + 1 padding

def BatchNorm1d(out_features):
    return nn.BatchNorm1d(2 * out_features)  # 2 for BiGRU output size

def Dense(out_features):
    return nn.Linear(2 * out_features, 1)

def amino_acid_to_tokens(sequence, tokenizer, max_length=512):
    """
    Convert amino acid sequence to Transformer tokens.
    
    Args:
        sequence: Amino acid sequence string
        tokenizer: Transformer tokenizer
        max_length: Maximum sequence length
        
    Returns:
        dict: Tokenized inputs with input_ids and attention_mask
    """
    # Convert amino acids to space-separated tokens for Transformer
    # Transformer tokenizer works better with space-separated tokens
    spaced_sequence = ' '.join(list(sequence))
    
    # Tokenize the sequence
    tokenized = tokenizer(
        spaced_sequence,
        padding='max_length',
        truncation=True,
        max_length=max_length,
        return_tensors='pt'
    )
    
    return tokenized

class BiGRU(nn.Module):

    def __init__(self, in_features, out_features, num_layers, dropout):
        super().__init__()
        self.bigru = nn.GRU(
            in_features, out_features, num_layers=num_layers,
            dropout=dropout if num_layers > 1 else 0.0,
            bidirectional=True, batch_first=True
        )

    def forward(self, input):
        return self.bigru(input)

class Transformer_Detect(nn.Module):
    """Transformer-based model for phosphorylation site detection using amino acid sequences."""
    
    def __init__(self, 
                 vocab_size=24, 
                 hidden_size=256, 
                 num_layers=6, 
                 num_heads=8, 
                 dropout=0.3, 
                 max_length=512):
        super().__init__()
        
        # Custom Transformer-like architecture for amino acid sequences
        self.hidden_size = hidden_size
        self.max_length = max_length
        self.vocab_size = vocab_size
        
        # Embedding layer for amino acids (24 amino acids + 1 padding token)
        self.embedding = nn.Embedding(vocab_size, hidden_size)
        
        # Positional encoding
        self.positional_encoding = nn.Embedding(max_length, hidden_size)
        
        # Transformer encoder layers
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=hidden_size,
            nhead=num_heads,
            dim_feedforward=hidden_size * 4,
            dropout=dropout,
            batch_first=True
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # Classification head
        self.dropout = nn.Dropout(dropout)
        self.classifier = nn.Linear(hidden_size, 1)
        self.sigmoid = nn.Sigmoid()
        
        # Initialize parameters properly
        self._init_weights()
        
    def _init_weights(self):
        """Initialize model weights for better convergence."""
        # Initialize embedding layers with Xavier/Glorot initialization
        nn.init.xavier_uniform_(self.embedding.weight)
        nn.init.xavier_uniform_(self.positional_encoding.weight)
        
        # Initialize transformer encoder layers
        for module in self.transformer_encoder.modules():
            if isinstance(module, nn.Linear):
                # Initialize linear layers with Xavier/Glorot initialization
                nn.init.xavier_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.LayerNorm):
                # Initialize layer norm with standard initialization
                nn.init.ones_(module.weight)
                nn.init.zeros_(module.bias)
            elif isinstance(module, nn.MultiheadAttention):
                # Initialize attention weights
                nn.init.xavier_uniform_(module.in_proj_weight)
                nn.init.xavier_uniform_(module.out_proj.weight)
                if module.in_proj_bias is not None:
                    nn.init.zeros_(module.in_proj_bias)
                if module.out_proj.bias is not None:
                    nn.init.zeros_(module.out_proj.bias)
        
        # Initialize classifier with proper scaling
        nn.init.xavier_uniform_(self.classifier.weight)
        nn.init.zeros_(self.classifier.bias)
        
        # Scale the classifier output for better initial predictions
        # This helps with the sigmoid activation
        with torch.no_grad():
            self.classifier.weight.data *= 0.1
        
    def forward(self, input):
        """
        Forward pass through Transformer-like model.
        
        Args:
            input: Tensor of shape (batch_size, seq_len) with amino acid indices
            
        Returns:
            torch.Tensor: Binary classification output (0-1)
        """
        batch_size, seq_len = input.shape
        
        # Create positional indices
        positions = torch.arange(seq_len, device=input.device).unsqueeze(0).expand(batch_size, -1)
        
        # Embedding + positional encoding
        embedded = self.embedding(input) + self.positional_encoding(positions)
        
        # Scale embeddings by sqrt(hidden_size) for better training stability
        embedded = embedded * (self.hidden_size ** 0.5)
        
        # Create padding mask for transformer
        padding_mask = (input == 0)
        
        # Transformer encoder
        transformer_output = self.transformer_encoder(embedded, src_key_padding_mask=padding_mask)
        
        # Global average pooling over sequence dimension
        # Use mean pooling instead of masked pooling to avoid dimension issues
        pooled_output = transformer_output.mean(dim=1)
        
        # Classification
        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)
        output = self.sigmoid(logits)
        
        return output

class Bert_Detect(nn.Module):
    """BERT-based model for phosphorylation site detection with pre-training and fine-tuning support."""
    
    def __init__(self, 
                 model_name='google-bert/bert-base-cased',
                 cache_dir='/data0/wangb/wbscy/PhosSight/model',    
                 max_length=512,
                 dropout=0.3,
                 num_labels=1,
                 freeze_bert=False,
                 vocab_size=27,
                 hidden_size=256,
                 num_hidden_layers=6,
                 num_attention_heads=4,
                 intermediate_size=1024):  # New configurable transformer sizes
        super().__init__()
        
        self.max_length = max_length
        self.freeze_bert = freeze_bert
        self.vocab_size = vocab_size
        
        # Create new BERT model with custom vocab size instead of loading pre-trained
        from transformers import BertConfig, BertModel
        
        # Create BERT configuration
        config = BertConfig(
            vocab_size=vocab_size,
            hidden_size=hidden_size,
            num_hidden_layers=num_hidden_layers,
            num_attention_heads=num_attention_heads,
            intermediate_size=intermediate_size,
            hidden_dropout_prob=dropout,
            attention_probs_dropout_prob=dropout,
            max_position_embeddings=max_length,
            type_vocab_size=1,
            initializer_range=0.02,
            layer_norm_eps=1e-12,
            pad_token_id=0,
            cls_token_id=24,  # CLS token
            sep_token_id=25,  # SEP token
            mask_token_id=26,  # MASK token
        )
        
        # Create new BERT model
        self.bert = BertModel(config)
        
        # Freeze BERT layers if specified (for fine-tuning)
        if freeze_bert:
            for param in self.bert.parameters():
                param.requires_grad = False
        
        # Get BERT hidden size
        self.hidden_size = self.bert.config.hidden_size
        
        # Classification head for phosphorylation detection
        self.dropout = nn.Dropout(dropout)
        self.classifier = nn.Linear(self.hidden_size, num_labels)
        self.sigmoid = nn.Sigmoid()
        
        # MLM head (transform + decoder), with tied weights to embeddings
        self.mlm_transform = nn.Sequential(
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.GELU(),
            nn.LayerNorm(self.hidden_size, eps=1e-12)
        )
        self.mlm_decoder = nn.Linear(self.hidden_size, vocab_size, bias=False)
        self.mlm_bias = nn.Parameter(torch.zeros(vocab_size))
        # Tie decoder weights to word embeddings
        self.mlm_decoder.weight = self.bert.embeddings.word_embeddings.weight
        
        # Initialize classifier weights
        self._init_classifier_weights()
        
    def _init_classifier_weights(self):
        """Initialize classifier weights for better convergence."""
        nn.init.xavier_uniform_(self.classifier.weight)
        nn.init.zeros_(self.classifier.bias)
        
        # Scale the classifier output for better initial predictions
        with torch.no_grad():
            self.classifier.weight.data *= 0.1
    
    def forward(self, input_ids=None, attention_mask=None, token_type_ids=None, labels=None):
        """
        Forward pass through BERT model.
        
        Args:
            input_ids: Tokenized input sequences
            attention_mask: Attention mask for padding
            token_type_ids: Token type IDs (not used for amino acid sequences)
            labels: Ground truth labels (optional, for training)
            
        Returns:
            torch.Tensor: Binary classification output (0-1)
        """
        # BERT forward pass
        outputs = self.bert(
            input_ids=input_ids,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids
        )
        
        # Get pooled output (CLS token representation)
        pooled_output = outputs.pooler_output
        
        # Apply dropout
        pooled_output = self.dropout(pooled_output)
        
        # Classification
        logits = self.classifier(pooled_output)
        outputs = self.sigmoid(logits)
        
        return outputs
    
    def get_mlm_logits(self, hidden_states: torch.Tensor) -> torch.Tensor:
        """Compute MLM logits from hidden states using the persistent head."""
        transformed = self.mlm_transform(hidden_states)
        logits = self.mlm_decoder(transformed) + self.mlm_bias
        return logits
    
    def get_bert_embeddings(self, input_ids, attention_mask=None):
        """
        Get BERT embeddings for pre-training tasks.
        
        Args:
            input_ids: Tokenized input sequences
            attention_mask: Attention mask for padding
            
        Returns:
            torch.Tensor: BERT embeddings
        """
        with torch.no_grad():
            outputs = self.bert(
                input_ids=input_ids,
                attention_mask=attention_mask
            )
            return outputs.last_hidden_state

class biGRU_Detect(nn.Module):

    def __init__(self, in_features, out_features, num_layers, dropout):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.num_layers = num_layers
        self.Embedding = Embedding(in_features)
        self.BiGRU = BiGRU(in_features, out_features, num_layers, dropout)
        self.BatchNorm1d = BatchNorm1d(out_features)
        self.Dense = Dense(out_features)
        self.Flatten = nn.Flatten()
        self.Sigmoid = nn.Sigmoid()

    def forward(self, input):
        out = self.Embedding(input)
        out, h = self.BiGRU(out)  # out dim: (batch_size, seq_len, 2 * out_features), h dim: (2* num_layers, batch_size, out_features)
        h = h.view(self.num_layers, 2, input.size(0), -1)  # this code is to reshape h to (num_layers, 2, batch_size, out_features)
        h_last = h[-1]  # take the last layer's output, h_last dim: (2, batch_size, out_features)
        h_last = h_last.permute(1, 0, 2).contiguous().view(input.size(0), -1)  # h_last dim: (batch_size, 2 * out_features)
        out = self.BatchNorm1d(h_last) # out dim: (batch_size, 2 * out_features)
        out = self.Dense(out) # out dim: (batch_size, 1)
        out = self.Sigmoid(out)
        return out


class AminoAcidProperties:
    """氨基酸物理化学性质特征提取器"""
    
    def __init__(self):
        # 保守的物理化学性质特征（只选择最重要的几个）
        self.aa_properties = {
            'Z': {'hydrophobicity': 0.0, 'charge': 0.0, 'polarity': 0.0},  # padding
            'A': {'hydrophobicity': 1.8, 'charge': 0.0, 'polarity': 0.0},
            'C': {'hydrophobicity': 2.5, 'charge': 0.0, 'polarity': 0.0},
            'D': {'hydrophobicity': -3.5, 'charge': -1.0, 'polarity': 1.0},
            'E': {'hydrophobicity': -3.5, 'charge': -1.0, 'polarity': 1.0},
            'F': {'hydrophobicity': 2.8, 'charge': 0.0, 'polarity': 0.0},
            'G': {'hydrophobicity': -0.4, 'charge': 0.0, 'polarity': 0.0},
            'H': {'hydrophobicity': -3.2, 'charge': 0.1, 'polarity': 1.0},
            'I': {'hydrophobicity': 4.5, 'charge': 0.0, 'polarity': 0.0},
            'K': {'hydrophobicity': -3.9, 'charge': 1.0, 'polarity': 1.0},
            'L': {'hydrophobicity': 3.8, 'charge': 0.0, 'polarity': 0.0},
            'M': {'hydrophobicity': 1.9, 'charge': 0.0, 'polarity': 0.0},
            'N': {'hydrophobicity': -3.5, 'charge': 0.0, 'polarity': 1.0},
            'P': {'hydrophobicity': -1.6, 'charge': 0.0, 'polarity': 0.0},
            'Q': {'hydrophobicity': -3.5, 'charge': 0.0, 'polarity': 1.0},
            'R': {'hydrophobicity': -4.5, 'charge': 1.0, 'polarity': 1.0},
            'S': {'hydrophobicity': -0.8, 'charge': 0.0, 'polarity': 1.0},
            'T': {'hydrophobicity': -0.7, 'charge': 0.0, 'polarity': 1.0},
            'V': {'hydrophobicity': 4.2, 'charge': 0.0, 'polarity': 0.0},
            'W': {'hydrophobicity': -0.9, 'charge': 0.0, 'polarity': 0.0},
            'Y': {'hydrophobicity': -1.3, 'charge': 0.0, 'polarity': 1.0},
            's': {'hydrophobicity': -0.8, 'charge': -1.0, 'polarity': 1.0},  # phosphorylated S
            't': {'hydrophobicity': -0.7, 'charge': -1.0, 'polarity': 1.0},  # phosphorylated T
            'y': {'hydrophobicity': -1.3, 'charge': -1.0, 'polarity': 1.0}   # phosphorylated Y
        }
    
    def extract_properties(self, sequence):
        """提取序列的物理化学性质特征"""
        features = []
        for aa in sequence:
            if aa in self.aa_properties:
                props = self.aa_properties[aa]
                # 归一化特征值
                hydrophobicity = props['hydrophobicity'] / 5.0  # 归一化到[-1, 1]
                charge = props['charge'] / 1.0  # 已经是[-1, 1]
                polarity = props['polarity']  # 0或1
                features.append([hydrophobicity, charge, polarity])
            else:
                features.append([0.0, 0.0, 0.0])
        return torch.tensor(features, dtype=torch.float32)


class AttentionLayer(nn.Module):
    """注意力层"""
    
    def __init__(self, hidden_size):
        super().__init__()
        self.hidden_size = hidden_size
        self.attention = nn.Sequential(
            nn.Linear(hidden_size, hidden_size // 2),
            nn.Tanh(),
            nn.Linear(hidden_size // 2, 1)
        )
        
    def forward(self, hidden_states):
        attention_weights = self.attention(hidden_states)
        attention_weights = F.softmax(attention_weights, dim=1)
        attended_output = torch.sum(attention_weights * hidden_states, dim=1)
        return attended_output, attention_weights

class biGRU_Detect_Improved_V2(nn.Module):
    """改进的BiGRU模型V2 - 添加物理化学性质特征"""

    def __init__(self, in_features, out_features, num_layers, dropout):
        super().__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.num_layers = num_layers
        
        # 基础组件
        self.Embedding = Embedding(in_features)
        self.BiGRU = BiGRU(in_features, out_features, num_layers, dropout)
        self.Attention = AttentionLayer(2 * out_features)
        self.BatchNorm1d = BatchNorm1d(out_features)
        self.Dense = Dense(out_features)
        self.Sigmoid = nn.Sigmoid()
        self.Dropout = nn.Dropout(dropout)
        
        # 新增：物理化学性质特征处理
        self.aa_properties = AminoAcidProperties()
        self.properties_proj = nn.Linear(3, in_features)  # 3个性质特征投影到embedding维度
        self.feature_fusion = nn.Linear(in_features * 2, in_features)  # 特征融合
        
        # 新增：特征门控机制
        self.feature_gate = nn.Sequential(
            nn.Linear(in_features * 2, in_features),
            nn.Sigmoid()
        )
        
        # 新增：编码索引映射（优化性能）
        self.aa_properties_encoded = {
            0: {'hydrophobicity': 0.0, 'charge': 0.0, 'polarity': 0.0},  # Z
            1: {'hydrophobicity': 1.8, 'charge': 0.0, 'polarity': 0.0},  # A
            2: {'hydrophobicity': 2.5, 'charge': 0.0, 'polarity': 0.0},  # C
            3: {'hydrophobicity': -3.5, 'charge': -1.0, 'polarity': 1.0},  # D
            4: {'hydrophobicity': -3.5, 'charge': -1.0, 'polarity': 1.0},  # E
            5: {'hydrophobicity': 2.8, 'charge': 0.0, 'polarity': 0.0},  # F
            6: {'hydrophobicity': -0.4, 'charge': 0.0, 'polarity': 0.0},  # G
            7: {'hydrophobicity': -3.2, 'charge': 0.1, 'polarity': 1.0},  # H
            8: {'hydrophobicity': 4.5, 'charge': 0.0, 'polarity': 0.0},  # I
            9: {'hydrophobicity': -3.9, 'charge': 1.0, 'polarity': 1.0},  # K
            10: {'hydrophobicity': 3.8, 'charge': 0.0, 'polarity': 0.0},  # L
            11: {'hydrophobicity': 1.9, 'charge': 0.0, 'polarity': 0.0},  # M
            12: {'hydrophobicity': -3.5, 'charge': 0.0, 'polarity': 1.0},  # N
            13: {'hydrophobicity': -1.6, 'charge': 0.0, 'polarity': 0.0},  # P
            14: {'hydrophobicity': -3.5, 'charge': 0.0, 'polarity': 1.0},  # Q
            15: {'hydrophobicity': -4.5, 'charge': 1.0, 'polarity': 1.0},  # R
            16: {'hydrophobicity': -0.8, 'charge': 0.0, 'polarity': 1.0},  # S
            17: {'hydrophobicity': -0.7, 'charge': 0.0, 'polarity': 1.0},  # T
            18: {'hydrophobicity': 4.2, 'charge': 0.0, 'polarity': 0.0},  # V
            19: {'hydrophobicity': -0.9, 'charge': 0.0, 'polarity': 0.0},  # W
            20: {'hydrophobicity': -1.3, 'charge': 0.0, 'polarity': 1.0},  # Y
            21: {'hydrophobicity': -0.8, 'charge': -1.0, 'polarity': 1.0},  # s
            22: {'hydrophobicity': -0.7, 'charge': -1.0, 'polarity': 1.0},  # t
            23: {'hydrophobicity': -1.3, 'charge': -1.0, 'polarity': 1.0}   # y
        }

    def forward(self, input):
        # 1. 序列嵌入
        seq_embedding = self.Embedding(input)  # (batch_size, seq_len, in_features)
        
        # 2. 提取物理化学性质特征（优化版本）
        batch_size, seq_len = input.shape
        properties_features = []
        
        # 批量处理以提高效率
        for i in range(batch_size):
            # 直接使用编码索引获取性质，避免字符串转换
            batch_properties = []
            for j in range(seq_len):
                code = input[i, j].item()
                if code in self.aa_properties_encoded:
                    props = self.aa_properties_encoded[code]
                    batch_properties.append([props['hydrophobicity'], props['charge'], props['polarity']])
                else:
                    batch_properties.append([0.0, 0.0, 0.0])
            properties_features.append(torch.tensor(batch_properties, dtype=torch.float32))
        
        properties_tensor = torch.stack(properties_features).to(input.device)
        properties_proj = self.properties_proj(properties_tensor)  # (batch_size, seq_len, in_features)
        
        # 3. 特征融合
        combined_features = torch.cat([seq_embedding, properties_proj], dim=-1)  # (batch_size, seq_len, in_features*2)
        fused_features = self.feature_fusion(combined_features)  # (batch_size, seq_len, in_features)
        
        # 4. 特征门控
        gate_weights = self.feature_gate(combined_features)  # (batch_size, seq_len, in_features)
        gated_features = gate_weights * seq_embedding + (1 - gate_weights) * fused_features
        
        # 5. BiGRU处理
        out, h = self.BiGRU(gated_features)  # out: (batch_size, seq_len, 2*out_features)
        
        # 6. Dropout
        out = self.Dropout(out)
        
        # 7. 注意力机制
        attended_output, attention_weights = self.Attention(out)  # (batch_size, 2*out_features)
        
        # 8. Batch normalization
        out = self.BatchNorm1d(attended_output)  # (batch_size, 2*out_features)
        
        # 9. 最终分类
        out = self.Dense(out)  # (batch_size, 1)
        out = self.Sigmoid(out)
        
        return out
    
    def _decode_sequence(self, encoded_seq):
        """将数字编码转换回氨基酸序列"""
        aa_dict = {0: 'Z', 1: 'A', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 9: 'K',
                   10: 'L', 11: 'M', 12: 'N', 13: 'P', 14: 'Q', 15: 'R', 16: 'S', 17: 'T',
                   18: 'V', 19: 'W', 20: 'Y', 21: 's', 22: 't', 23: 'y'}
        
        sequence = ''
        for code in encoded_seq:
            if code.item() in aa_dict:
                sequence += aa_dict[code.item()]
            else:
                sequence += 'Z'
        return sequence

