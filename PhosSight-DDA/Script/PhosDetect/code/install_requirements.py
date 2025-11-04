#!/usr/bin/env python3
"""
Installation script for BERT model dependencies.

This script checks and installs the required packages for the BERT model training.
"""

import subprocess
import sys
import pkg_resources

def check_package(package_name):
    """Check if a package is installed."""
    try:
        pkg_resources.get_distribution(package_name)
        return True
    except pkg_resources.DistributionNotFound:
        return False

def install_package(package_name):
    """Install a package using pip."""
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package_name])
        return True
    except subprocess.CalledProcessError:
        return False

def check_torch_version():
    """Check PyTorch version and suggest upgrade if needed."""
    try:
        import torch
        version = torch.__version__
        print(f"PyTorch version: {version}")
        
        # Check if version is >= 2.1
        if version.startswith('2.0') or version.startswith('1.'):
            print("⚠️  Warning: PyTorch version < 2.1 detected")
            print("   Some transformers features may not work properly")
            print("   Consider upgrading PyTorch:")
            print("   pip install torch>=2.1.0")
            return False
        else:
            print("✓ PyTorch version is compatible")
            return True
    except ImportError:
        print("✗ PyTorch not found")
        return False

def main():
    """Main installation function."""
    print("BERT Model Dependencies Check")
    print("=" * 40)
    
    # Required packages
    required_packages = [
        "torch>=2.0.0",
        "transformers>=4.20.0",
        "tqdm",
        "pandas",
        "scikit-learn",
        "numpy"
    ]
    
    print("\nChecking required packages...")
    
    # Check PyTorch version first
    torch_ok = check_torch_version()
    
    # Check other packages
    missing_packages = []
    for package in required_packages:
        package_name = package.split('>=')[0].split('<=')[0]
        if check_package(package_name):
            print(f"✓ {package_name} is installed")
        else:
            print(f"✗ {package_name} is missing")
            missing_packages.append(package)
    
    # Install missing packages
    if missing_packages:
        print(f"\nInstalling {len(missing_packages)} missing packages...")
        for package in missing_packages:
            print(f"Installing {package}...")
            if install_package(package):
                print(f"✓ Successfully installed {package}")
            else:
                print(f"✗ Failed to install {package}")
    else:
        print("\n✓ All required packages are installed!")
    
    # Test imports
    print("\nTesting imports...")
    try:
        import torch
        print("✓ PyTorch imported successfully")
    except ImportError as e:
        print(f"✗ PyTorch import failed: {e}")
    
    try:
        from transformers import BertTokenizer, BertModel
        print("✓ Transformers imported successfully")
    except ImportError as e:
        print(f"✗ Transformers import failed: {e}")
    
    try:
        import tqdm
        print("✓ tqdm imported successfully")
    except ImportError as e:
        print(f"✗ tqdm import failed: {e}")
    
    try:
        import pandas
        print("✓ pandas imported successfully")
    except ImportError as e:
        print(f"✗ pandas import failed: {e}")
    
    try:
        import sklearn
        print("✓ scikit-learn imported successfully")
    except ImportError as e:
        print(f"✗ scikit-learn import failed: {e}")
    
    print("\n" + "=" * 40)
    print("Installation check completed!")
    print("=" * 40)
    
    if not torch_ok:
        print("\n⚠️  Note: Consider upgrading PyTorch for better compatibility")
        print("   pip install torch>=2.1.0")

if __name__ == "__main__":
    main() 