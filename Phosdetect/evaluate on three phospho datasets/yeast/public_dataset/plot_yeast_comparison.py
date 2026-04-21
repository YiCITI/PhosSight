#!/usr/bin/env python3
"""
绘制 Yeast 数据集上三个模型的评测结果对比图
与 mouse 目录风格完全一致
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 标准样式 - 与mouse目录一致
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
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'axes.unicode_minus': False
})

DATA_DIR = '/data0/wangb/cd/duibi0826/0826comparison/data/yeast/public_dataset'


def plot_comparison():
    """绘制对比图 - 与mouse风格完全一致"""
    
    # 读取数据
    df = pd.read_csv(f"{DATA_DIR}/three_models_ft_v3_results.csv")
    print(f"读取数据:")
    print(df)
    
    models = df['Model'].tolist()
    original_scores = df['Original_AUC'].tolist()
    finetuned_scores = df['Finetuned_AUC'].tolist()
    
    # 配色方案 - 与mouse完全一致
    colors_ori = ['#D4A5A5', '#B8D4E3', '#A8D4C8']  # PhosDetect红色, DeepDetect蓝色, pFly绿色
    colors_ft = ['#C59595', '#9CC5D8', '#90C5B0']   # 对应的深色版本
    
    # 创建图形 - 与mouse风格完全一致 (figsize=(4.5, 4))
    fig, ax = plt.subplots(figsize=(4.5, 4))
    
    x = np.arange(len(models))
    width = 0.35
    
    # 绘制柱状图
    for i, model in enumerate(models):
        # Original 柱
        ax.bar(x[i] - width/2, original_scores[i], width, 
               color=colors_ori[i], 
               edgecolor='white', linewidth=0.6, alpha=0.9)
        
        # Finetuned 柱
        ax.bar(x[i] + width/2, finetuned_scores[i], width, 
               color=colors_ft[i], 
               edgecolor='white', linewidth=0.6, alpha=0.9)
    
    # 添加数值标签
    for i in range(len(models)):
        # Original 数值
        ax.annotate(f'{original_scores[i]:.3f}',
                   xy=(x[i] - width/2, original_scores[i]),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=6)
        
        # Finetuned 数值
        ax.annotate(f'{finetuned_scores[i]:.3f}',
                   xy=(x[i] + width/2, finetuned_scores[i]),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=6, fontweight='bold')
    
    # 添加图例 - 与mouse风格完全一致 (右上角)
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
    
    # 网格线
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # 边框样式 - 与mouse风格完全一致
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)
    
    plt.tight_layout()
    
    # 保存 - 文件名与mouse一致
    svg_path = os.path.join(DATA_DIR, 'three_models_ft_v3_comparison.svg')
    png_path = os.path.join(DATA_DIR, 'three_models_ft_v3_comparison.png')
    
    plt.savefig(svg_path, format='svg', bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(png_path, format='png', bbox_inches='tight', dpi=300, facecolor='white', edgecolor='none')
    print(f"\nSVG saved: {svg_path}")
    print(f"PNG saved: {png_path}")
    
    plt.close()
    print("Done!")


if __name__ == "__main__":
    plot_comparison()
