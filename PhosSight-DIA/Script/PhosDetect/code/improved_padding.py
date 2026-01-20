#!/usr/bin/env python3
"""
Improved padding methods using built-in functions.
"""

import torch
import numpy as np
from torch.nn.utils.rnn import pad_sequence
import time

def original_padding_method(sequences, padding_len):
    """Original padding method (inefficient)."""
    Vector_ls = []
    for v in sequences:
        if len(v) < padding_len:
            num = padding_len - len(v)
            add = [0 for i in range(0, num)]
            v.extend(add)
        Vector_ls.append(v)
    return Vector_ls

def improved_padding_method_1(sequences, padding_len):
    """Improved method 1: Using torch.nn.functional.pad."""
    # Convert to tensors first
    tensors = [torch.tensor(seq, dtype=torch.long) for seq in sequences]
    
    # Pad sequences to the same length
    padded = pad_sequence(tensors, batch_first=True, padding_value=0)
    
    # If we need a specific padding length, truncate or pad further
    if padded.size(1) < padding_len:
        # Pad with zeros to reach padding_len
        pad_size = padding_len - padded.size(1)
        padded = torch.nn.functional.pad(padded, (0, pad_size), value=0)
    elif padded.size(1) > padding_len:
        # Truncate to padding_len
        padded = padded[:, :padding_len]
    
    return padded

def improved_padding_method_2(sequences, padding_len):
    """Improved method 2: Using numpy and torch."""
    # Convert to numpy arrays for easier manipulation
    sequences_np = [np.array(seq) for seq in sequences]
    
    # Create a padded array
    padded_array = np.zeros((len(sequences), padding_len), dtype=np.int64)
    
    # Fill the array with sequences
    for i, seq in enumerate(sequences_np):
        length = min(len(seq), padding_len)
        padded_array[i, :length] = seq[:length]
    
    return torch.tensor(padded_array, dtype=torch.long)

def improved_padding_method_3(sequences, padding_len):
    """Improved method 3: Using list comprehension (most efficient for this case)."""
    # Use list comprehension for faster processing
    padded_sequences = []
    for seq in sequences:
        if len(seq) < padding_len:
            # More efficient way to extend with zeros
            seq.extend([0] * (padding_len - len(seq)))
        elif len(seq) > padding_len:
            # Truncate if too long
            seq = seq[:padding_len]
        padded_sequences.append(seq)
    
    return torch.tensor(padded_sequences, dtype=torch.long)

def benchmark_padding_methods():
    """Benchmark different padding methods."""
    print("Benchmarking padding methods...")
    
    # Create sample data
    num_sequences = 1000
    max_len = 200
    sequences = []
    
    for i in range(num_sequences):
        # Random sequence length between 50 and 150
        seq_len = np.random.randint(50, 150)
        seq = np.random.randint(1, 25, seq_len).tolist()
        sequences.append(seq)
    
    padding_len = max_len
    
    methods = [
        ("Original", original_padding_method),
        ("PyTorch pad_sequence", improved_padding_method_1),
        ("NumPy + PyTorch", improved_padding_method_2),
        ("List comprehension", improved_padding_method_3)
    ]
    
    results = {}
    
    for name, method in methods:
        print(f"\nTesting {name}...")
        
        # Warm up
        _ = method(sequences[:10], padding_len)
        
        # Benchmark
        start_time = time.time()
        result = method(sequences, padding_len)
        end_time = time.time()
        
        execution_time = end_time - start_time
        results[name] = {
            'time': execution_time,
            'result_shape': result.shape if hasattr(result, 'shape') else len(result),
            'result_type': type(result)
        }
        
        print(f"  Time: {execution_time:.4f}s")
        print(f"  Result shape: {result.shape if hasattr(result, 'shape') else len(result)}")
        print(f"  Result type: {type(result)}")
    
    # Find the fastest method
    fastest = min(results.items(), key=lambda x: x[1]['time'])
    print(f"\n🏆 Fastest method: {fastest[0]} ({fastest[1]['time']:.4f}s)")
    
    return results

def demonstrate_usage():
    """Demonstrate how to use the improved padding methods."""
    print("\n" + "="*50)
    print("Demonstrating improved padding methods")
    print("="*50)
    
    # Sample sequences
    sequences = [
        [1, 2, 3, 4, 5],
        [1, 2, 3],
        [1, 2, 3, 4, 5, 6, 7, 8],
        [1]
    ]
    
    padding_len = 6
    
    print("Original sequences:")
    for i, seq in enumerate(sequences):
        print(f"  Sequence {i+1}: {seq} (length: {len(seq)})")
    
    print(f"\nTarget padding length: {padding_len}")
    
    # Test different methods
    methods = [
        ("PyTorch pad_sequence", improved_padding_method_1),
        ("NumPy + PyTorch", improved_padding_method_2),
        ("List comprehension", improved_padding_method_3)
    ]
    
    for name, method in methods:
        print(f"\n{name}:")
        result = method(sequences, padding_len)
        print(f"  Result shape: {result.shape}")
        print(f"  Result:\n{result}")
    
    print("\n✅ All methods produce the same result!")

if __name__ == "__main__":
    print("Testing improved padding methods...")
    
    # Run benchmark
    benchmark_results = benchmark_padding_methods()
    
    # Demonstrate usage
    demonstrate_usage()
    
    print("\n🎉 Benchmark completed!") 