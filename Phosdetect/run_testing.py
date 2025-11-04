#!/usr/bin/env python3


import os
import sys
import subprocess
import argparse
from datetime import datetime

def run_testing(model_path, test_datasets, args):
    """运行测试脚本"""
    
    # 创建结果保存目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = f"results/{timestamp}"
    log_dir = f"logs/{timestamp}"
    
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    results_file = f"{results_dir}/test_results.csv"
    
    print("="*80)
    print("PhosSight V2 Testing")
    print("="*80)
    print(f"Model path: {model_path}")
    print(f"Model type: bigru_improved_v2")
    print(f"Results directory: {results_dir}")
    print(f"Log directory: {log_dir}")
    print("="*80)
    
    all_results = []
    
    for dataset_name, dataset_path in test_datasets.items():
        print(f"\n🧪 Testing on {dataset_name}...")
        
        # 构建测试命令
        cmd = [
            "python", "test.py",
            "-p", dataset_path,
            "-m", model_path,
            "-b", str(args.batch_size),
            "--model_type", "bigru_improved_v2",
            "--log_file", f"{log_dir}/test_{dataset_name}.log",
            "--output_file", results_file,
            "--in_features", str(args.in_features),
            "--out_features", str(args.out_features),
            "--num_layers", str(args.num_layers),
            "--dropout", str(args.dropout)
        ]
        
        if args.verbose:
            cmd.append("--verbose")
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            print(f"✅ {dataset_name} testing completed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"❌ {dataset_name} testing failed with error code: {e.returncode}")
            continue
    
    print(f"\n📊 All testing completed!")
    print(f"📁 Results saved to: {results_file}")
    print(f"📁 Logs saved to: {log_dir}")
    
    return results_file

def main():
    parser = argparse.ArgumentParser(description="PhosSight V2 Testing Script")
    parser.add_argument("--model_path", type=str, required=True, 
                       help="Path to trained model weights")
    parser.add_argument("--batch_size", type=int, default=128, help="Batch size")
    parser.add_argument("--in_features", type=int, default=10, help="Input features")
    parser.add_argument("--out_features", type=int, default=20, help="Output features")
    parser.add_argument("--num_layers", type=int, default=2, help="Number of layers")
    parser.add_argument("--dropout", type=float, default=0.3, help="Dropout rate")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    # 检查模型文件是否存在
    if not os.path.exists(args.model_path):
        print(f"❌ Error: Model file not found: {args.model_path}")
        sys.exit(1)
    
    # 定义测试数据集
    test_datasets = {
        "DeepDetect_ecoli": "../data/test/DeepDetect_ecoli.csv",
        "DeepDetect_human": "../data/test/DeepDetect_human.csv",
        "DeepDetect_mouse": "../data/test/DeepDetect_mouse.csv",
        "DeepDetect_yeast": "../data/test/DeepDetect_yeast.csv",
        "DeepRescore2_HCC": "../data/test/DeepRescore2_HCC.csv",
        "DeepRescore2_label_free": "../data/test/DeepRescore2_label_free.csv",
        "DeepRescore2_UCEC": "../data/test/DeepRescore2_UCEC.csv"
    }
    
    # 检查数据集文件是否存在
    missing_datasets = []
    for name, path in test_datasets.items():
        if not os.path.exists(path):
            missing_datasets.append(name)
    
    if missing_datasets:
        print(f"⚠️  Warning: Some test datasets not found: {missing_datasets}")
        print("Available datasets will be tested.")
    
    # 运行测试
    results_file = run_testing(args.model_path, test_datasets, args)
    
    if results_file:
        print(f"\n🎯 Next steps:")
        print(f"1. Compare results: python compare_results.py")
        print(f"2. Visualize improvements: python plot_auc_comparison.py")

if __name__ == "__main__":
    main() 