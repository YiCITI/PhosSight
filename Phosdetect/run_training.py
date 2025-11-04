#!/usr/bin/env python3
"""
PhosSight V2 Training Script
训练改进的BiGRU模型V2（添加物理化学性质特征）
"""

import os
import sys
import subprocess
import argparse
from datetime import datetime

def run_training(args):
    """运行训练脚本"""
    
    # 创建日志和模型保存目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = f"logs/{timestamp}"
    model_dir = f"models/{timestamp}"
    
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(model_dir, exist_ok=True)
    
    # 构建训练命令
    cmd = [
        "python", "train.py",
        "-p", args.data_path,
        "-m", f"{model_dir}/best_model.pth",
        "-b", str(args.batch_size),
        "-e", str(args.epochs),
        "-lr", str(args.learning_rate),
        "-d", str(args.dropout),
        "--model_type", "bigru_improved_v2",
        "--log_file", f"{log_dir}/training.log",
        "--patience", str(args.patience),
        "--optimizer", args.optimizer,
        "--weight_decay", str(args.weight_decay),
        "--scheduler", args.scheduler,
        "--in_features", str(args.in_features),
        "--out_features", str(args.out_features),
        "--num_layers", str(args.num_layers)
    ]
    
    if args.verbose:
        cmd.append("--verbose")
    
    print("="*80)
    print("PhosSight V2 Training")
    print("="*80)
    print(f"Data path: {args.data_path}")
    print(f"Model type: bigru_improved_v2")
    print(f"Batch size: {args.batch_size}")
    print(f"Epochs: {args.epochs}")
    print(f"Learning rate: {args.learning_rate}")
    print(f"Dropout: {args.dropout}")
    print(f"Optimizer: {args.optimizer}")
    print(f"Scheduler: {args.scheduler}")
    print(f"Log directory: {log_dir}")
    print(f"Model directory: {model_dir}")
    print("="*80)
    
    # 运行训练
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"\n✅ Training completed successfully!")
        print(f"📁 Logs saved to: {log_dir}")
        print(f"📁 Model saved to: {model_dir}")
        return f"{model_dir}/best_model.pth"
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Training failed with error code: {e.returncode}")
        return None

def main():
    parser = argparse.ArgumentParser(description="PhosSight V2 Training Script")
    parser.add_argument("--data_path", type=str, default="../data/train/balanced_dataset_1.csv", 
                       help="Path to training data")
    parser.add_argument("--batch_size", type=int, default=128, help="Batch size")
    parser.add_argument("--epochs", type=int, default=100, help="Number of epochs")
    parser.add_argument("--learning_rate", type=float, default=0.0005, help="Learning rate")
    parser.add_argument("--dropout", type=float, default=0.3, help="Dropout rate")
    parser.add_argument("--patience", type=int, default=15, help="Early stopping patience")
    parser.add_argument("--optimizer", type=str, default="adam", choices=["adam", "adamw", "sgd"], 
                       help="Optimizer type")
    parser.add_argument("--weight_decay", type=float, default=1e-5, help="Weight decay")
    parser.add_argument("--scheduler", type=str, default="plateau", choices=["plateau", "cosine", "step"], 
                       help="Learning rate scheduler")
    parser.add_argument("--in_features", type=int, default=10, help="Input features")
    parser.add_argument("--out_features", type=int, default=20, help="Output features")
    parser.add_argument("--num_layers", type=int, default=2, help="Number of layers")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    # 检查数据文件是否存在
    if not os.path.exists(args.data_path):
        print(f"❌ Error: Data file not found: {args.data_path}")
        sys.exit(1)
    
    # 运行训练
    model_path = run_training(args)
    
    if model_path:
        print(f"\n🎯 Next steps:")
        print(f"1. Test the model: python run_testing.py --model_path {model_path}")
        print(f"2. Compare results: python compare_results.py")
        print(f"3. Visualize improvements: python plot_auc_comparison.py")

if __name__ == "__main__":
    main() 