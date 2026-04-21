#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制磷酸化数据集的Original vs Finetuned AUC对比图
与plot_non_phosphorylation_auc_comparison.py保持一致的风格
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# 标准样式 - 与参考脚本完全一致
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'font.size': 9,
    'axes.titlesize': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 7,
    'axes.linewidth': 1,
    'grid.linewidth': 0.8,
    'lines.linewidth': 1,
    'patch.linewidth': 1,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'axes.unicode_minus': False,
    'svg.fonttype': 'none'
})

DATA_CSV = '/data0/wangb/cd/duibi0826/0826comparison/data/all_model_comparison_results.csv'
OUTDIR = '/data0/wangb/cd/duibi0826/0826comparison/data'


def load_data():
    """加载数据"""
    df = pd.read_csv(DATA_CSV)
    return df


def plot_phosphorylation_comparison(df, output_svg):
    """绘制磷酸化数据集的Original vs Finetuned对比图（左右拼接）"""
    datasets = ['Rice', 'Mouse', 'Yeast']
    models = ['PhosDetect', 'DeepDetect', 'pFly']
    model_labels = {'PhosDetect': 'PhosDetect', 'DeepDetect': 'DeepDetect', 'pFly': 'PFly'}
    
    # 配色与参考脚本一致
    original_colors = ['#D4A5A5', '#B8D4E3', '#A8D4C8']  # PhosDetect, DeepDetect, pFly
    finetuned_colors = ['#C59595', '#9CC5D8', '#90C5B0']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for ax, col_name, colors, ft_suffix in [
        (ax1, 'Original_AUC', original_colors, ''),
        (ax2, 'Finetuned_AUC', finetuned_colors, 'Fine-tuned ')
    ]:
        width = 0.2
        dataset_gap = 0.7
        x = np.arange(len(datasets)) * dataset_gap

        num_models = len(models)
        offsets = (np.arange(num_models) - (num_models - 1) / 2.0) * width

        # 绘制柱状图
        for i, (model, color) in enumerate(zip(models, colors)):
            values = []
            for ds in datasets:
                val = df[(df['Dataset'] == ds) & (df['Model'] == model)][col_name].values[0]
                values.append(val)
            
            bars = ax.bar(x + offsets[i], values, width, label=ft_suffix + model_labels[model], 
                         color=color, alpha=0.9, edgecolor='white', linewidth=0.6)

        ax.set_xlabel('Datasets')
        ax.set_ylabel('AUC Score')
        ax.set_xticks(x)
        ax.set_xticklabels(datasets)
        ax.set_ylim(0.82, 1.0)
        ax.set_yticks(np.arange(0.82, 1.01, 0.04))
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        # 图例放在顶部，字体放大到10
        ax.legend(ncol=3, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.02), fontsize=10)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.0)
        ax.spines['bottom'].set_linewidth(1.0)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_svg), exist_ok=True)
    plt.savefig(output_svg, bbox_inches='tight', facecolor='white', edgecolor='none', format='svg')
    plt.close(fig)


def main():
    """主函数"""
    print("开始创建磷酸化数据集AUC对比图...")
    
    df = load_data()
    
    output_svg = os.path.join(OUTDIR, 'phosphorylation_comparison_combined.svg')
    plot_phosphorylation_comparison(df, output_svg)
    
    print(f"✓ 已保存: {output_svg}")
    print("\n✅ 磷酸化数据集图表绘制完成!")


if __name__ == '__main__':
    main()
