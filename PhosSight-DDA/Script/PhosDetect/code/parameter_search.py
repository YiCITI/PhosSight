#!/usr/bin/env python3
"""
Parameter search script for hyperparameter optimization.
"""

import subprocess
import itertools
import json
import os
from datetime import datetime

def run_training_experiment(params):
    """Run a single training experiment with given parameters."""
    cmd = [
        'python', 'train.py',
        '--model_type', params['model_type'],
        '--epochs', str(params['epochs']),
        '--batch_size', str(params['batch_size']),
        '--learning_rate', str(params['learning_rate']),
        '--dropout', str(params['dropout']),
        '--optimizer', params['optimizer'],
        '--scheduler', params['scheduler'],
        '--weight_decay', str(params['weight_decay']),
        '--patience', str(params['patience']),
        '--gradient_clip', str(params['gradient_clip']),
        '--log_file', f"./logs/experiment_{params['exp_id']}.log"
    ]
    
    # Add model-specific parameters
    if params['model_type'] == 'bigru':
        cmd.extend([
            '--in_features', str(params['in_features']),
            '--out_features', str(params['out_features']),
            '--num_layers', str(params['num_layers'])
        ])
    elif params['model_type'] == 'transformer':
        cmd.extend([
            '--vocab_size', str(params['vocab_size']),
            '--hidden_size', str(params['hidden_size']),
            '--num_layers', str(params['num_layers']),
            '--num_heads', str(params['num_heads']),
            '--max_length', str(params['max_length'])
        ])
    
    print(f"Running experiment {params['exp_id']}: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)  # 1 hour timeout
        return {
            'exp_id': params['exp_id'],
            'params': params,
            'success': result.returncode == 0,
            'stdout': result.stdout,
            'stderr': result.stderr
        }
    except subprocess.TimeoutExpired:
        return {
            'exp_id': params['exp_id'],
            'params': params,
            'success': False,
            'error': 'Timeout'
        }
    except Exception as e:
        return {
            'exp_id': params['exp_id'],
            'params': params,
            'success': False,
            'error': str(e)
        }

def generate_parameter_combinations():
    """Generate parameter combinations for grid search."""
    
    # Base parameters
    base_params = {
        'model_type': ['bigru', 'transformer'],
        'epochs': [50],
        'batch_size': [64, 128],
        'learning_rate': [0.001, 0.0005],
        'dropout': [0.3, 0.5],
        'optimizer': ['adam', 'adamw'],
        'scheduler': ['plateau', 'cosine'],
        'weight_decay': [1e-5, 1e-4],
        'patience': [10, 15],
        'gradient_clip': [0.0, 1.0]
    }
    
    # Model-specific parameters
    bigru_params = {
        'in_features': [10],
        'out_features': [16, 20],
        'num_layers': [1, 2]
    }
    
    transformer_params = {
        'vocab_size': [25],
        'hidden_size': [64, 128],
        'num_layers': [1, 2],
        'num_heads': [2, 4],
        'max_length': [128, 256]
    }
    
    experiments = []
    exp_id = 1
    
    # Generate all combinations
    for base_combo in itertools.product(*base_params.values()):
        base_dict = dict(zip(base_params.keys(), base_combo))
        
        if base_dict['model_type'] == 'bigru':
            for bigru_combo in itertools.product(*bigru_params.values()):
                bigru_dict = dict(zip(bigru_params.keys(), bigru_combo))
                experiment = {**base_dict, **bigru_dict, 'exp_id': exp_id}
                experiments.append(experiment)
                exp_id += 1
        else:  # transformer
            for transformer_combo in itertools.product(*transformer_params.values()):
                transformer_dict = dict(zip(transformer_params.keys(), transformer_combo))
                experiment = {**base_dict, **transformer_dict, 'exp_id': exp_id}
                experiments.append(experiment)
                exp_id += 1
    
    return experiments

def parse_results(log_file):
    """Parse training results from log file."""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract best validation accuracy
        lines = content.split('\n')
        best_val_acc = 0.0
        final_test_acc = 0.0
        final_test_auc = 0.0
        
        for line in lines:
            if 'New best model saved with val acc:' in line:
                try:
                    acc = float(line.split('val acc: ')[1].split()[0])
                    best_val_acc = max(best_val_acc, acc)
                except:
                    pass
            elif 'Final Test Results' in line:
                try:
                    parts = line.split('Acc: ')[1].split(',')
                    final_test_acc = float(parts[0])
                    final_test_auc = float(parts[-1].split('AUC: ')[1])
                except:
                    pass
        
        return {
            'best_val_acc': best_val_acc,
            'final_test_acc': final_test_acc,
            'final_test_auc': final_test_auc
        }
    except:
        return {
            'best_val_acc': 0.0,
            'final_test_acc': 0.0,
            'final_test_auc': 0.0
        }

def main():
    """Main function for parameter search."""
    print("Starting parameter search for model training...")
    
    # Create logs directory
    os.makedirs('./logs', exist_ok=True)
    
    # Generate parameter combinations
    experiments = generate_parameter_combinations()
    print(f"Generated {len(experiments)} experiments")
    
    # Run experiments
    results = []
    for i, params in enumerate(experiments):
        print(f"\n{'='*60}")
        print(f"Experiment {i+1}/{len(experiments)}")
        print(f"{'='*60}")
        
        result = run_training_experiment(params)
        
        if result['success']:
            # Parse results from log file
            log_file = f"./logs/experiment_{params['exp_id']}.log"
            metrics = parse_results(log_file)
            result['metrics'] = metrics
        else:
            result['metrics'] = {
                'best_val_acc': 0.0,
                'final_test_acc': 0.0,
                'final_test_auc': 0.0
            }
        
        results.append(result)
        
        # Save intermediate results
        with open('./logs/search_results.json', 'w') as f:
            json.dump(results, f, indent=2)
    
    # Find best results
    successful_results = [r for r in results if r['success']]
    if successful_results:
        best_by_val_acc = max(successful_results, key=lambda x: x['metrics']['best_val_acc'])
        best_by_test_acc = max(successful_results, key=lambda x: x['metrics']['final_test_acc'])
        best_by_auc = max(successful_results, key=lambda x: x['metrics']['final_test_auc'])
        
        print(f"\n{'='*60}")
        print("BEST RESULTS")
        print(f"{'='*60}")
        print(f"Best by validation accuracy: {best_by_val_acc['metrics']['best_val_acc']:.4f}")
        print(f"  Experiment ID: {best_by_val_acc['exp_id']}")
        print(f"  Parameters: {best_by_val_acc['params']}")
        
        print(f"\nBest by test accuracy: {best_by_test_acc['metrics']['final_test_acc']:.4f}")
        print(f"  Experiment ID: {best_by_test_acc['exp_id']}")
        print(f"  Parameters: {best_by_test_acc['params']}")
        
        print(f"\nBest by AUC: {best_by_auc['metrics']['final_test_auc']:.4f}")
        print(f"  Experiment ID: {best_by_auc['exp_id']}")
        print(f"  Parameters: {best_by_auc['params']}")
    
    # Save final results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    with open(f'./logs/final_results_{timestamp}.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nParameter search completed! Results saved to ./logs/final_results_{timestamp}.json")

if __name__ == "__main__":
    main() 