#!/usr/bin/env python3
"""
绘制 Mouse 数据集上三个模型的微调前后评测结果对比图
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'font.size': 7,
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.titlesize': 7,
    'axes.linewidth': 1,
    'grid.linewidth': 0.8,
    'lines.linewidth': 1,
    'patch.linewidth': 1,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'axes.unicode_minus': False
})

DATA_DIR = '/data0/wangb/cd/duibi0826/0826comparison/data/mouse'


def plot_comparison():
    # 原始评测结果
    original_results = {
        'PhosDetect': 0.9390,
        'DeepDetect': 0.8859,
        'pFly': 0.8722
    }
    
    # 微调后评测结果
    finetuned_results = {
        'PhosDetect': 0.9390,
        'DeepDetect': 0.9037,
        'pFly': 0.8722
    }
    
    models = ['PhosDetect', 'DeepDetect', 'pFly']
    original_scores = [original_results[m] for m in models]
    finetuned_scores = [finetuned_results[m] for m in models]
    
    colors_ori = ['#D4A5A5', '#B8D4E3', '#A8D4C8']
    colors_ft = ['#C59595', '#9CC5D8', '#90C5B0']
    
    fig, ax = plt.subplots(figsize=(4.5, 4))
    
    x = np.arange(len(models))
    width = 0.35
    
    for i, model in enumerate(models):
        ax.bar(x[i] - width/2, original_scores[i], width,
               color=colors_ori[i], edgecolor='white', linewidth=0.6, alpha=0.9)
        ax.bar(x[i] + width/2, finetuned_scores[i], width,
               color=colors_ft[i], edgecolor='white', linewidth=0.6, alpha=0.9)
    
    for i in range(len(models)):
        ax.annotate(f'{original_scores[i]:.3f}',
                   xy=(x[i] - width/2, original_scores[i]),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=6)
        ax.annotate(f'{finetuned_scores[i]:.3f}',
                   xy=(x[i] + width/2, finetuned_scores[i]),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=6, fontweight='bold')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#B8B8B8', edgecolor='white', label='Original'),
        Patch(facecolor='#888888', edgecolor='white', label='Finetuned')
    ]
    ax.legend(handles=legend_elements, ncol=1, frameon=False, loc='upper right', fontsize=7)
    
    ax.set_xlabel('Model', fontsize=8)
    ax.set_ylabel('AUC Score', fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=8)
    ax.set_ylim(0.5, 1.05)
    ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)
    
    plt.tight_layout()
    
    svg_path = os.path.join(DATA_DIR, 'mouse_comparison.svg')
    png_path = os.path.join(DATA_DIR, 'mouse_comparison.png')
    
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(png_path, format='png', bbox_inches='tight', dpi=300, facecolor='white', edgecolor='none')
    print(f"SVG saved: {svg_path}")
    print(f"PNG saved: {png_path}")
    
    plt.close()
    print("Done!")


if __name__ == "__main__":
    plot_comparison()
